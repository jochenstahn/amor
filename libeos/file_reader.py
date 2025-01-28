import logging
import os
import subprocess
import sys
from datetime import datetime, timezone
try:
    import zoneinfo
except ImportError:
    # for python versions < 3.9 try to use the backports version
    from backports import zoneinfo
from typing import List

import h5py
import numpy as np
from orsopy import fileio
from orsopy.fileio.model_language import SampleModel

from . import const
from .header import Header
from .instrument import Detector
from .options import ExperimentConfig, ReaderConfig

try:
    from . import nb_helpers
except Exception:
    nb_helpers = None

# Time zone used to interpret time strings
AMOR_LOCAL_TIMEZONE = zoneinfo.ZoneInfo(key='Europe/Zurich')

class AmorData:
    """read meta-data and event streams from .hdf file(s), apply filters and conversions"""
    chopperDetectorDistance: float
    chopperDistance: float
    chopperPhase: float
    chopperSpeed: float
    div: float
    data_file_numbers: List[int]
    delta_z: np.ndarray
    detZ_e: np.ndarray
    lamda_e: np.ndarray
    wallTime_e: np.ndarray
    kad: float
    kap: float
    lambdaMax: float
    lambda_e: np.ndarray
    #monitor: float
    mu: float
    nu: float
    tau: float
    tofCut: float
    start_date: str
    monitorType: str

    seriesStartTime = None

    #-------------------------------------------------------------------------------------------------
    def __init__(self, header: Header, reader_config: ReaderConfig, config: ExperimentConfig,
                 short_notation:str, norm=False):
        #self.startTime = reader_config.startTime
        self.header = header
        self.config = config
        self.reader_config = reader_config
        self.expand_file_list(short_notation)
        self.read_data(norm=norm)

    #-------------------------------------------------------------------------------------------------
    def read_data(self, norm=False):
        self.file_list = []
        for number in self.data_file_numbers:
            self.file_list.append(self.path_generator(number))
        ## read specific meta data and measurement from first file
        if norm:
            self.readHeaderInfo = False
        else:
            self.readHeaderInfo = True

        _detZ_e = []
        _lamda_e = []
        _wallTime_e = []
        #_monitor = 0
        _monitorPerPulse = []
        _pulseTimeS = []
        for file in self.file_list:
            self.read_individual_data(file, norm)
            _detZ_e = np.append(_detZ_e, self.detZ_e)
            _lamda_e = np.append(_lamda_e, self.lamda_e)
            _wallTime_e = np.append(_wallTime_e, self.wallTime_e)
            _monitorPerPulse = np.append(_monitorPerPulse, self.monitorPerPulse)
            _pulseTimeS = np.append(_pulseTimeS, self.pulseTimeS)
            #_monitor += self.monitor
        self.detZ_e = _detZ_e
        self.lamda_e = _lamda_e
        self.wallTime_e = _wallTime_e
        #self.monitor = _monitor
        self.monitorPerPulse = _monitorPerPulse   
        self.pulseTimeS    = _pulseTimeS

    #-------------------------------------------------------------------------------------------------
    #def path_generator(self, number):
    #    fileName = f'amor{self.reader_config.year}n{number:06d}.hdf'
    #    if os.path.exists(os.path.join(self.reader_config.dataPath,fileName)):
    #        path = self.reader_config.dataPath
    #    elif os.path.exists(fileName):
    #        path = '.'
    #    elif os.path.exists(os.path.join('.','raw', fileName)):
    #        path = os.path.join('.','raw')
    #    elif os.path.exists(os.path.join('..','raw', fileName)):
    #        path = os.path.join('..','raw')
    #    elif os.path.exists(f'/afs/psi.ch/project/sinqdata/{self.reader_config.year}/amor/{int(number/1000)}/{fileName}'):
    #        path = f'/afs/psi.ch/project/sinqdata/{self.reader_config.year}/amor/{int(number/1000)}'
    #    else:
    #        sys.exit(f'# ERROR: the file {fileName} is nowhere to be found!')
    #    return os.path.join(path, fileName)
    #-------------------------------------------------------------------------------------------------
    def path_generator(self, number):
        fileName = f'amor{self.reader_config.year}n{number:06d}.hdf'
        path = ''
        for rawd in self.reader_config.rawPath:
            if os.path.exists(os.path.join(rawd,fileName)):
                path = rawd
                break
        if not path:
            if os.path.exists(f'/afs/psi.ch/project/sinqdata/{self.reader_config.year}/amor/{int(number/1000)}/{fileName}'):
                path = f'/afs/psi.ch/project/sinqdata/{self.reader_config.year}/amor/{int(number/1000)}'
            else:
                sys.exit(f'# ERROR: the file {fileName} can not be found in {self.reader_config.rawPath}')
        return os.path.join(path, fileName)
    #-------------------------------------------------------------------------------------------------
    def expand_file_list(self, short_notation):
        """Evaluate string entry for file number lists"""
        #log().debug('Executing get_flist')
        file_list=[]
        for i in short_notation.split(','):
            if '-' in i:
                if ':' in i:
                    step = i.split(':', 1)[1]
                    file_list += range(int(i.split('-', 1)[0]), int((i.rsplit('-', 1)[1]).split(':', 1)[0])+1, int(step))
                else:
                    step = 1
                    file_list += range(int(i.split('-', 1)[0]), int(i.split('-', 1)[1])+1, int(step))
            else:
                file_list += [int(i)]
        self.data_file_numbers=sorted(file_list)
    #-------------------------------------------------------------------------------------------------
    def resolve_pixels(self):
        """determine spatial coordinats and angles from pixel number"""
        nPixel = Detector.nWires * Detector.nStripes * Detector.nBlades
        pixelID = np.arange(nPixel)
        (bladeNr, bPixel) = np.divmod(pixelID, Detector.nWires * Detector.nStripes)
        (bZi, detYi)      = np.divmod(bPixel, Detector.nStripes)                     # z index on blade, y index on detector
        detZi             = bladeNr * Detector.nWires + bZi                          # z index on detector
        detX              = bZi * Detector.dX                                        # x position in detector
        # detZ              = Detector.zero - bladeNr * Detector.bladeZ - bZi * Detector.dZ      # z position on detector
        bladeAngle        = np.rad2deg( 2. * np.arcsin(0.5*Detector.bladeZ / Detector.distance) )
        delta             = (Detector.nBlades/2. - bladeNr) * bladeAngle \
                            - np.rad2deg( np.arctan(bZi*Detector.dZ / ( Detector.distance + bZi * Detector.dX) ) )
        self.delta_z      = delta[detYi==1]
        return np.vstack((detYi.T, detZi.T, detX.T, delta.T)).T
    #-------------------------------------------------------------------------------------------------
    def read_individual_data(self, fileName, norm=False):
        self.hdf = h5py.File(fileName, 'r', swmr=True)

        if self.readHeaderInfo:
            self.read_header_info()

        logging.warning(f'    from file: {fileName}')
        self.read_individual_header()

        # add header content
        if self.readHeaderInfo:
            self.readHeaderInfo = False
            self.header.measurement_instrument_settings = fileio.InstrumentSettings(
                incident_angle = fileio.ValueRange(round(self.mu+self.kap+self.kad-0.5*self.div, 3),
                                                   round(self.mu+self.kap+self.kad+0.5*self.div, 3),
                                                   'deg'),
                wavelength = fileio.ValueRange(const.lamdaCut, self.config.lambdaRange[1], 'angstrom'),
                polarization = fileio.Polarization.unpolarized,
                )
            self.header.measurement_instrument_settings.mu = fileio.Value(round(self.mu, 3), 'deg', comment='sample angle to horizon')
            self.header.measurement_instrument_settings.nu = fileio.Value(round(self.nu, 3), 'deg', comment='detector angle to horizon')
            self.header.measurement_instrument_settings.div = fileio.Value(round(self.div, 3), 'deg', comment='incoming beam divergence')
            self.header.measurement_instrument_settings.kap = fileio.Value(round(self.kap, 3), 'deg', comment='incoming beam inclination')
            if abs(self.kad)>0.02:
                self.header.measurement_instrument_settings.kad = fileio.Value(round(self.kad, 3), 'deg', comment='incoming beam angular offset')
        if norm:
            self.header.measurement_additional_files.append(fileio.File(file=fileName.split('/')[-1], timestamp=self.fileDate))
        else:
            self.header.measurement_data_files.append(fileio.File(file=fileName.split('/')[-1], timestamp=self.fileDate))
        logging.info(f'      mu = {self.mu:6.3f}, nu = {self.nu:6.3f}, kap = {self.kap:6.3f}, kad = {self.kad:6.3f}')

        self.read_event_stream()
        totalNumber = np.shape(self.tof_e)[0]
        # check for empty event stream
        if totalNumber == 0:
             logging.error('empty event stream: can not determine end time')
             sys.exit()

        self.sort_pulses()

        self.associate_pulse_with_monitor()

        self.extract_walltime(norm)

        # following lines: debugging output to trace the time-offset of proton current and neutron pulses
        if self.config.monitorType == 'x':
            cpp, t_bins = np.histogram(self.wallTime_e, self.pulseTimeS)
            np.savetxt('tme.hst', np.vstack((self.pulseTimeS[:-1], cpp, self.monitorPerPulse[:-1])).T)

        #self.average_events_per_pulse() # for debugging only. VERY time consuming!!!

        self.monitor_threshold()

        self.filter_strange_times()

        self.merge_frames()

        self.filter_project_x()

        self.correct_for_chopper_opening()

        self.calculate_derived_properties()

        self.filter_qz_range(norm)

        logging.info(f'      number of events: total = {totalNumber:7d}, filtered = {np.shape(self.lamda_e)[0]:7d}')

    def sort_pulses(self):
        chopperPeriod = np.int64(2*self.tau*1e9)
        pulseTime = np.sort(self.dataPacketTime_p)
        pulseTime = pulseTime[np.abs(pulseTime[:]-np.roll(pulseTime, 1)[:])>5]

        pulseTime -= np.int64(self.seriesStartTime)
        self.stopTime = pulseTime[-1]
        pulseTime = pulseTime[pulseTime>=0]

        # fill in missing pulse times 
        # TODO: check for real end time
        try:
            # further files
            # TODO: use the first pulse of the respective measurement
            #nextPulseTime = startTime % np.int64(self.tau*2e9)
            #nextPulseTime = self.pulseTimeS[-1] + chopperPeriod
            nextPulseTime = pulseTime[0]
        except AttributeError:
            # first file
            nextPulseTime = pulseTime[0] % np.int64(self.tau*2e9)

        # calculate where time tiefference between pulses exceeds its time by more than 1/2
        # this yields the number of missing pulses
        pulseLengths = pulseTime[1:]-pulseTime[:-1]
        pulseExtra = (pulseLengths-np.int64(self.tau*1e9))//np.int64(self.tau*2e9)
        gap_indices = np.where(pulseExtra>0)[0]

        if len(gap_indices)==0:
            # no missing pulses, just use given array
            self.pulseTimeS = np.array(pulseTime, dtype=np.int64)
            return
        self.pulseTimeS = np.array(pulseTime[:gap_indices[0]+1], dtype=np.int64)
        last_index = gap_indices[0]
        for gapi in gap_indices[1:]:
            # insert missing pulses into each gap
            gap_pulses = pulseTime[last_index]+np.arange(1, pulseExtra[last_index]+1)*chopperPeriod
            self.pulseTimeS = np.append(self.pulseTimeS, gap_pulses)
            self.pulseTimeS = np.append(self.pulseTimeS, pulseTime[last_index+1:gapi+1])
            last_index = gapi
        if last_index<len(pulseTime):
            self.pulseTimeS = np.append(self.pulseTimeS, pulseTime[last_index:-1])

    def get_current_per_pulse(self, pulseTimeS, currentTimeS, currents):
        # add currents for early pulses and current time value after last pulse (j+1)
        currentTimeS = np.hstack([[0], currentTimeS, [pulseTimeS[-1]+1]])
        currents = np.hstack([[0], currents])
        pulseCurrentS = np.zeros(pulseTimeS.shape[0], dtype=float)
        j = 0
        for i, ti in enumerate(pulseTimeS):
            if ti >= currentTimeS[j+1]: 
                j += 1
            pulseCurrentS[i] = currents[j]
            #print(f' {i}  {pulseTimeS[i]}  {pulseCurrentS[i]}') 
        return pulseCurrentS

    def associate_pulse_with_monitor(self):
        if self.config.monitorType == 'p': # protonCharge
            self.currentTime -= np.int64(self.seriesStartTime)
            self.monitorPerPulse = self.get_current_per_pulse(self.pulseTimeS, self.currentTime, self.current) * 2*self.tau * 1e-3
            # filter low-current pulses
            self.monitorPerPulse = np.where(self.monitorPerPulse > 2*self.tau * self.config.lowCurrentThreshold * 1e-3, self.monitorPerPulse, 0)
        elif self.config.monitorType == 't': # countingTime
            self.monitorPerPulse = np.ones(np.shape(self.pulseTimeS)[0])*self.tau
        else: 
            self.monitorPerPulse = 1./np.shape(pulseTimeS)[1]

    def extract_walltime(self, norm):
        if nb_helpers:
            self.wallTime_e = nb_helpers.extract_walltime(self.tof_e, self.dataPacket_p, self.dataPacketTime_p)
        else:
            self.wallTime_e = np.empty(np.shape(self.tof_e)[0], dtype=np.int64)
            for i in range(len(self.dataPacket_p)-1):
                self.wallTime_e[self.dataPacket_p[i]:self.dataPacket_p[i+1]] = self.dataPacketTime_p[i]
            self.wallTime_e[self.dataPacket_p[-1]:] = self.dataPacketTime_p[-1]
        self.wallTime_e -= np.int64(self.seriesStartTime)
        logging.debug(f'      wall time from {self.wallTime_e[0]/1e9:6.1f} s to {self.wallTime_e[-1]/1e9:6.1f} s')

    def average_events_per_pulse(self):
        if self.config.monitorType == 'p':
            for i, time in enumerate(self.pulseTimeS):
                events = np.shape(self.wallTime_e[self.wallTime_e == time])[0]
                logging.info(f'pulse: {i:6.0f}, events: {events:6.0f}, monitor: {self.monitorPerPulse[i]:6.2f}')

    def monitor_threshold(self):
        if self.config.monitorType == 'p': # fix to check for file compatibility
            goodTimeS = self.pulseTimeS[self.monitorPerPulse!=0]
            filter_e = np.where(np.isin(self.wallTime_e, goodTimeS), True, False)
            self.tof_e = self.tof_e[filter_e]
            self.pixelID_e = self.pixelID_e[filter_e]
            self.wallTime_e = self.wallTime_e[filter_e]
            logging.info(f'      rejected {np.shape(self.monitorPerPulse)[0]-np.shape(goodTimeS)[0]} out of {np.shape(self.monitorPerPulse)[0]} pulses')
            logging.info(f'          with {np.shape(filter_e)[0]-np.shape(self.tof_e)[0]} events due to low beam current')
            logging.info(f'      average counts per pulse =  {np.shape(self.tof_e)[0] / np.shape(goodTimeS[goodTimeS!=0])[0]:7.1f}')

    def filter_qz_range(self, norm):
        if self.config.qzRange[1]<0.3 and not norm:
            self.mask_e = np.logical_and(self.mask_e,
                                         (self.config.qzRange[0]<=self.qz_e) & (self.qz_e<=self.config.qzRange[1]))
        self.detZ_e = self.detZ_e[self.mask_e]
        self.lamda_e = self.lamda_e[self.mask_e]
        self.wallTime_e = self.wallTime_e[self.mask_e]

    def calculate_derived_properties(self):
        self.lamdaMax = const.lamdaCut+1.e13*self.tau*const.hdm/(self.chopperDetectorDistance+124.)
        if nb_helpers:
            self.lamda_e, self.qz_e, self.mask_e = nb_helpers.calculate_derived_properties_focussing(
                    self.tof_e, self.detXdist_e, self.delta_e, self.mask_e,
                    self.config.lambdaRange[0], self.config.lambdaRange[1], self.nu, self.mu,
                    self.chopperDetectorDistance, const.hdm
                    )
            return
        # lambda
        self.lamda_e = (1.e13*const.hdm)*self.tof_e/(self.chopperDetectorDistance+self.detXdist_e)
        self.mask_e = np.logical_and(self.mask_e, (self.config.lambdaRange[0]<=self.lamda_e) & (
                    self.lamda_e<=self.config.lambdaRange[1]))
        # alpha_f
        # q_z
        if self.config.incidentAngle == 'alphaF':
            alphaF_e  = self.nu - self.mu + self.delta_e
            self.qz_e = 4*np.pi*(np.sin(np.deg2rad(alphaF_e))/self.lamda_e)
            # qx_e    = 0.
            self.header.measurement_scheme = 'angle- and energy-dispersive'
        elif self.config.incidentAngle == 'nu':
            alphaF_e  = (self.nu + self.delta_e + self.kap + self.kad) / 2.
            self.qz_e = 4*np.pi*(np.sin(np.deg2rad(alphaF_e))/self.lamda_e)
            # qx_e    = 0.
            self.header.measurement_scheme = 'energy-dispersive'
        else:
            alphaF_e  = self.nu - self.mu + self.delta_e
            alphaI    = self.kap + self.kad + self.mu
            self.qz_e = 2*np.pi * ((np.sin(np.deg2rad(alphaF_e)) + np.sin(np.deg2rad(alphaI)))/self.lamda_e)
            self.qx_e = 2*np.pi * ((np.cos(np.deg2rad(alphaF_e)) - np.cos(np.deg2rad(alphaI)))/self.lamda_e)
            self.header.measurement_scheme = 'energy-dispersive'

    def correct_for_chopper_opening(self):
        # correct tof for beam size effect at chopper:  t_cor = (delta / 180 deg) * tau
        if self.config.incidentAngle == 'alphaF':
            self.tof_e    -= ( self.delta_e / 180. ) * self.tau
        else:
            # TODO: check sign of correction
            self.tof_e    -= ( self.kad / 180. ) * self.tau

    def filter_project_x(self):
        pixelLookUp = self.resolve_pixels()
        if nb_helpers:
            (self.detZ_e, self.detXdist_e, self.delta_e, self.mask_e) = nb_helpers.filter_project_x(
                    pixelLookUp, self.pixelID_e.astype(np.int64), self.config.yRange[0], self.config.yRange[1]
                    )
        else:
            # resolve pixel ID into y and z indicees, x position and angle
            (detY_e, self.detZ_e, self.detXdist_e, self.delta_e) = pixelLookUp[np.int_(self.pixelID_e)-1, :].T
            # define mask and filter y range
            self.mask_e = (self.config.yRange[0]<=detY_e) & (detY_e<=self.config.yRange[1])

    def merge_frames(self):
        total_offset = self.tofCut+self.tau*self.config.chopperPhaseOffset/180.
        if nb_helpers:
            self.tof_e = nb_helpers.merge_frames(self.tof_e, self.tofCut, self.tau, total_offset)
        else:
            self.tof_e = np.remainder(self.tof_e-(self.tofCut-self.tau), self.tau)+total_offset  # tof shifted to 1 frame

    def filter_strange_times(self):
        # 'strange' tof times are those with t > 2 tau (originating from the efu)
        filter_e = (self.tof_e<=2*self.tau)
        self.tof_e = self.tof_e[filter_e]
        self.pixelID_e = self.pixelID_e[filter_e]
        self.wallTime_e = self.wallTime_e[filter_e]
        if np.shape(filter_e)[0]-np.shape(self.tof_e)[0]>0.5:
            logging.warning(f'        strange times: {np.shape(filter_e)[0]-np.shape(self.tof_e)[0]}')

    def read_event_stream(self):
        self.tof_e = np.array(self.hdf['/entry1/Amor/detector/data/event_time_offset'][:])/1.e9
        self.pixelID_e = np.array(self.hdf['/entry1/Amor/detector/data/event_id'][:], dtype=np.int64)
        self.dataPacket_p = np.array(self.hdf['/entry1/Amor/detector/data/event_index'][:], dtype=np.uint64)
        self.dataPacketTime_p = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][:], dtype=np.int64)
        if self.config.monitorType in ['auto', 'p']:
            try:
                self.currentTime = np.array(self.hdf['entry1/Amor/detector/proton_current/time'][:], dtype=np.int64)
                self.current = np.array(self.hdf['entry1/Amor/detector/proton_current/value'][:,0], dtype=float)
                if len(self.current)>4:
                    self.config.monitorType = 'p'
                else:
                    self.config.monitorType = 't'
            except(KeyError, IndexError):
                self.config.monitorType = 't'
        else:
            self.config.monitorType = 't'
        #TODO: protonMonitor

    def read_individual_header(self):
        self.chopperDistance = float(np.take(self.hdf['entry1/Amor/chopper/pair_separation'], 0))
        self.detectorDistance = float(np.take(self.hdf['entry1/Amor/detector/transformation/distance'], 0))
        self.chopperDetectorDistance = self.detectorDistance-float(np.take(self.hdf['entry1/Amor/chopper/distance'], 0))
        self.tofCut = const.lamdaCut*self.chopperDetectorDistance/const.hdm*1.e-13

        try:
            self.mu   = float(np.take(self.hdf['/entry1/Amor/master_parameters/mu/value'], 0))
            self.nu   = float(np.take(self.hdf['/entry1/Amor/master_parameters/nu/value'], 0))
            self.kap  = float(np.take(self.hdf['/entry1/Amor/master_parameters/kap/value'], 0))
            self.kad  = float(np.take(self.hdf['/entry1/Amor/master_parameters/kad/value'], 0))
            self.div  = float(np.take(self.hdf['/entry1/Amor/master_parameters/div/value'], 0))
            self.chopperSpeed = float(np.take(self.hdf['/entry1/Amor/chopper/rotation_speed/value'], 0))
            self.chopperPhase = float(np.take(self.hdf['/entry1/Amor/chopper/phase/value'], 0))
        except(KeyError, IndexError):
            logging.warning("     using parameters from nicos cache")
            year_date = str(self.start_date).replace('-', '/', 1)
            #cachePath = '/home/amor/nicosdata/amor/cache/'
            #cachePath = '/home/nicos/amorcache/'
            cachePath = '/home/amor/cache/'
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-mu/{year_date}')).split('\t')[-1]
            self.mu = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-nu/{year_date}')).split('\t')[-1]
            self.nu = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-kap/{year_date}')).split('\t')[-1]
            self.kap = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-kad/{year_date}')).split('\t')[-1]
            self.kad = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-div/{year_date}')).split('\t')[-1]
            self.div = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-ch1_speed/{year_date}')).split('\t')[-1]
            self.chopperSpeed = float(value)
            self.chopperPhase = self.config.chopperPhase
        self.tau     = 30. / self.chopperSpeed

        if self.config.muOffset:
            self.mu += self.config.muOffset
        if self.config.mu:
            self.mu = self.config.mu
        if self.config.nu:
            self.nu = self.config.nu

        # extract start time as unix time, adding UTC offset of 1h to time string
        dz = datetime.fromisoformat(self.hdf['/entry1/start_time'][0].decode('utf-8'))
        self.fileDate=dz.replace(tzinfo=AMOR_LOCAL_TIMEZONE)
        self.startTime = np.int64( (self.fileDate.timestamp() ) * 1e9 )
        if self.seriesStartTime is None:
            self.seriesStartTime = self.startTime 

    def read_header_info(self):
        # read general information and first data set
        logging.info(f'    meta data from: {self.file_list[0]}')
        self.hdf = h5py.File(self.file_list[0], 'r', swmr=True)
        title = self.hdf['entry1/title'][0].decode('utf-8')
        proposal_id = self.hdf['entry1/proposal_id'][0].decode('utf-8')
        user_name = self.hdf['entry1/user/name'][0].decode('utf-8')
        user_affiliation = 'unknown'
        user_email = self.hdf['entry1/user/email'][0].decode('utf-8')
        user_orcid = None
        sampleName = self.hdf['entry1/sample/name'][0].decode('utf-8')
        model = self.hdf['entry1/sample/model'][0].decode('utf-8')
        instrumentName = 'Amor'
        source = self.hdf['entry1/Amor/source/name'][0].decode('utf-8')
        sourceProbe = 'neutron'
        start_time = self.hdf['entry1/start_time'][0].decode('utf-8')
        self.start_date = start_time.split(' ')[0]
        if self.config.sampleModel:
            model = self.config.sampleModel
        # assembling orso header information
        self.header.owner = fileio.Person(
                name=user_name,
                affiliation=user_affiliation,
                contact=user_email,
                )
        if user_orcid:
            self.header.owner.orcid = user_orcid
        self.header.experiment = fileio.Experiment(
                title=title,
                instrument=instrumentName,
                start_date=self.start_date,
                probe=sourceProbe,
                facility=source,
                proposalID=proposal_id
                )
        self.header.sample = fileio.Sample(
                name=sampleName,
                model=SampleModel(stack=model),
                sample_parameters=None,
                )
        self.header.measurement_scheme = 'angle- and energy-dispersive'

