import logging
import os
import subprocess
import sys
import platform
from datetime import datetime, timezone
from dataclasses import dataclass
try:
    import zoneinfo
except ImportError:
    # for python versions < 3.9 try to use the backports version
    from backports import zoneinfo
from typing import List, Optional
from abc import ABC, abstractmethod

import h5py
import numpy as np
from orsopy import fileio
from orsopy.fileio.model_language import SampleModel

from . import const
from .header import Header
from .instrument import Detector
from .options import ExperimentConfig, IncidentAngle, MonitorType, ReaderConfig
from .helpers import merge_frames, extract_walltime, filter_project_x, calculate_derived_properties_focussing

# Time zone used to interpret time strings
AMOR_LOCAL_TIMEZONE = zoneinfo.ZoneInfo(key='Europe/Zurich')

if  platform.node().startswith('amor'):
    NICOS_CACHE_DIR = '/home/amor/nicosdata/amor/cache/'
    GREP = '/usr/bin/grep "%s"'
else:
    NICOS_CACHE_DIR = None

class AmorData:
    """read meta-data and event streams from .hdf file(s), apply filters and conversions"""
    chopperDetectorDistance: float
    chopperDistance: float
    chopperPhase: float
    chopperSpeed: float
    chopper1TriggerPhase: float
    chopper2TriggerPhase: float
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
                    file_list += range(int(i.split('-', 1)[0]),
                                       int((i.rsplit('-', 1)[1]).split(':', 1)[0])+1,
                                       int(step))
                else:
                    step = 1
                    file_list += range(int(i.split('-', 1)[0]),
                                       int(i.split('-', 1)[1])+1,
                                       int(step))
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
                #polarization = fileio.Polarization.unpolarized,
                polarization = fileio.Polarization(self.polarizationConfig)
                )
            self.header.measurement_instrument_settings.mu = fileio.Value(
                    round(self.mu, 3),
                    'deg',
                    comment='sample angle to horizon')
            self.header.measurement_instrument_settings.nu = fileio.Value(
                    round(self.nu, 3),
                    'deg',
                    comment='detector angle to horizon')
            self.header.measurement_instrument_settings.div = fileio.Value(
                    round(self.div, 3),
                    'deg',
                    comment='incoming beam divergence')
            self.header.measurement_instrument_settings.kap = fileio.Value(
                    round(self.kap, 3),
                    'deg',
                    comment='incoming beam inclination')
            if abs(self.kad)>0.02:
                self.header.measurement_instrument_settings.kad = fileio.Value(
                        round(self.kad, 3),
                        'deg',
                        comment='incoming beam angular offset')
        if norm:
            self.header.measurement_additional_files.append(fileio.File(
                file=fileName.split('/')[-1],
                timestamp=self.fileDate))
        else:
            self.header.measurement_data_files.append(fileio.File(
                file=fileName.split('/')[-1],
                timestamp=self.fileDate))

        logging.info(f'      mu = {self.mu:6.3f}, nu = {self.nu:6.3f}, kap = {self.kap:6.3f}, kad = {self.kad:6.3f}')

        self.read_event_stream()

        self.correct_for_chopper_phases()

        self.read_chopper_trigger_stream()

        self.extract_walltime(norm)

        self.read_proton_current_stream()

        self.associate_pulse_with_monitor()

        # following lines: debugging output to trace the time-offset of proton current and neutron pulses
        if self.config.monitorType == MonitorType.debug:
            cpp, t_bins = np.histogram(self.wallTime_e, self.pulseTimeS)
            np.savetxt('tme.hst', np.vstack((self.pulseTimeS[:-1], cpp, self.monitorPerPulse[:-1])).T)

        #self.average_events_per_pulse() # for debugging only. VERY time consuming!!!

        self.monitor_threshold()

        self.filter_strange_times()

        self.merge_time_frames()

        self.filter_project_x()

        self.correct_for_chopper_opening()

        self.calculate_derived_properties()

        self.filter_qz_range(norm)

        logging.info(f'      number of events: total = {self.totalNumber:7d}, filtered = {np.shape(self.lamda_e)[0]:7d}')

    def read_event_stream(self):
        self.tof_e = np.array(self.hdf['/entry1/Amor/detector/data/event_time_offset'][:])/1.e9
        self.pixelID_e = np.array(self.hdf['/entry1/Amor/detector/data/event_id'][:], dtype=np.int64)
        self.dataPacket_p = np.array(self.hdf['/entry1/Amor/detector/data/event_index'][:], dtype=np.uint64)
        self.dataPacketTime_p = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][:], dtype=np.int64)

    def correct_for_chopper_phases(self):
        #print(f'tof phase-offset: {self.ch1TriggerPhase - self.chopperPhase/2}')
        self.tof_e += self.tau * (self.ch1TriggerPhase - self.chopperPhase/2)/180

    def read_chopper_trigger_stream(self):
        self.chopper1TriggerTime = np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][:-2],
                                            dtype=np.int64)
        #self.chopper2TriggerTime = self.chopper1TriggerTime + np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time'][:-2], dtype=np.int64)
        #                           + np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_offset'][:], dtype=np.int64)
        if np.shape(self.chopper1TriggerTime)[0] > 2:
            self.startTime = self.chopper1TriggerTime[0]
            self.stopTime = self.chopper1TriggerTime[-1]
            self.pulseTimeS = self.chopper1TriggerTime
        else:
            logging.warn('     no chopper trigger data available, using event steram instead')
            self.startTime = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][0], dtype=np.int64)
            self.stopTime = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][-2], dtype=np.int64)
            self.pulseTimeS = np.arange(self.startTime, self.stopTime, self.tau*1e9) 
        if self.seriesStartTime is None:
            self.seriesStartTime = self.startTime
            logging.debug(f'      series start time (epoch): {self.seriesStartTime/1e9:13.2f} s')
        self.pulseTimeS -= self.seriesStartTime
        logging.debug(f'      epoch time from {self.startTime/1e9:13.2f} s to {self.stopTime/1e9:13.2f} s')
        logging.debug(f'      => counting time {self.stopTime/1e9-self.startTime/1e9:8.2f} s')

    def extract_walltime(self, norm):
        self.wallTime_e = extract_walltime(self.tof_e, self.dataPacket_p, self.dataPacketTime_p)
        self.wallTime_e -= np.int64(self.seriesStartTime)
        logging.debug(f'      wall time from {self.wallTime_e[0]/1e9:6.1f} s to {self.wallTime_e[-1]/1e9:6.1f} s')

    def read_proton_current_stream(self):
        self.currentTime = np.array(self.hdf['entry1/Amor/detector/proton_current/time'][:], dtype=np.int64)
        self.current = np.array(self.hdf['entry1/Amor/detector/proton_current/value'][:,0], dtype=float)
        if self.config.monitorType == MonitorType.auto:
            if self.current.sum() > 1:
                self.monitorType = MonitorType.proton_charge
                logging.warn('      monitor type set to "proton current"')
            else:
                self.monitorType = MonitorType.time
                logging.warn('      monitor type set to "time"')

    def associate_pulse_with_monitor(self):
        if self.config.monitorType == MonitorType.proton_charge or MonitorType.debug:
            self.currentTime -= np.int64(self.seriesStartTime)
            self.monitorPerPulse = self.get_current_per_pulse(self.pulseTimeS,
                                                              self.currentTime,
                                                              self.current)\
                                                              * 2*self.tau * 1e-3
            # filter low-current pulses
            self.monitorPerPulse = np.where(self.monitorPerPulse > 2*self.tau * self.config.lowCurrentThreshold * 1e-3,
                                            self.monitorPerPulse,
                                            0)
        elif self.config.monitorType == MonitorType.time:
            self.monitorPerPulse = np.ones(np.shape(self.pulseTimeS)[0])*2*self.tau
        else: # pulses
            self.monitorPerPulse = np.ones(np.shape(self.pulseTimeS)[0])

    def get_current_per_pulse(self, pulseTimeS, currentTimeS, currents):
        # add currents for early pulses and current time value after last pulse (j+1)
        currentTimeS = np.hstack([[0], currentTimeS, [pulseTimeS[-1]+1]])
        currents = np.hstack([[0], currents])
        pulseCurrentS = np.zeros(pulseTimeS.shape[0], dtype=float)
        j = 0
        for i, ti in enumerate(pulseTimeS):
            while ti >= currentTimeS[j+1]: 
                j += 1
            pulseCurrentS[i] = currents[j]
        return pulseCurrentS

    def average_events_per_pulse(self):
        if self.config.monitorType == MonitorType.proton_charge:
            for i, time in enumerate(self.pulseTimeS):
                events = np.shape(self.wallTime_e[self.wallTime_e == time])[0]
                logging.info(f'pulse: {i:6.0f}, events: {events:6.0f}, monitor: {self.monitorPerPulse[i]:6.2f}')

    def monitor_threshold(self):
        #if self.config.monitorType == MonitorType.proton_charge: # fix to check for file compatibility
        self.totalNumber = np.shape(self.tof_e[self.tof_e<=self.stopTime])[0]
        if True:
            goodTimeS = self.pulseTimeS[self.monitorPerPulse!=0]
            filter_e = np.where(np.isin(self.wallTime_e, goodTimeS), True, False)
            self.tof_e = self.tof_e[filter_e]
            self.pixelID_e = self.pixelID_e[filter_e]
            self.wallTime_e = self.wallTime_e[filter_e]
            logging.info(f'      low-beam (<{self.config.lowCurrentThreshold} mC) rejected pulses: {np.shape(self.monitorPerPulse)[0]-1-np.shape(goodTimeS)[0]} out of {np.shape(self.monitorPerPulse)[0]-1}')
            logging.info(f'          with {np.shape(filter_e)[0]-np.shape(self.tof_e)[0]} events')
            if np.shape(goodTimeS[goodTimeS!=0])[0]:
                logging.info(f'      average counts per pulse =  {np.shape(self.tof_e)[0] / np.shape(goodTimeS[goodTimeS!=0])[0]:7.1f}')
            else:
                logging.info(f'      average counts per pulse = undefined')

    def filter_strange_times(self):
        # 'strange' tof times are those with t > 2 tau (originating from the efu)
        filter_e = (self.tof_e<=2*self.tau)
        self.tof_e = self.tof_e[filter_e]
        self.pixelID_e = self.pixelID_e[filter_e]
        self.wallTime_e = self.wallTime_e[filter_e]
        if np.shape(filter_e)[0]-np.shape(self.tof_e)[0]>0.5:
            logging.warning(f'        strange times: {np.shape(filter_e)[0]-np.shape(self.tof_e)[0]}')

    # TODO: - handle each neutron pulse individually, - associate with correct monitor also for slow neutrons
    def merge_time_frames(self):
        total_offset = self.tofCut + self.tau * (self.ch1TriggerPhase + self.chopperPhase/2)/180
        self.tof_e = merge_frames(self.tof_e, self.tofCut, self.tau, total_offset)

    def filter_project_x(self):
        pixelLookUp = self.resolve_pixels()
        (self.detZ_e, self.detXdist_e, self.delta_e, self.mask_e) = filter_project_x(
                pixelLookUp, self.pixelID_e.astype(np.int64), self.config.yRange[0], self.config.yRange[1]
                )

    def correct_for_chopper_opening(self):
        # correct tof for beam size effect at chopper:  t_cor = (delta / 180 deg) * tau
        if self.config.incidentAngle == IncidentAngle.alphaF:
            self.tof_e    -= ( self.delta_e / 180. ) * self.tau
        else:
            # TODO: check sign of correction
            self.tof_e    -= ( self.kad / 180. ) * self.tau

    def calculate_derived_properties(self):
        self.lamdaMax = const.lamdaCut+1.e13*self.tau*const.hdm/(self.chopperDetectorDistance+124.)
        if False:
            self.lamda_e, self.qz_e, self.mask_e = calculate_derived_properties_focussing(
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
        if self.config.incidentAngle == IncidentAngle.alphaF:
            alphaF_e  = self.nu - self.mu + self.delta_e
            self.qz_e = 4*np.pi*(np.sin(np.deg2rad(alphaF_e))/self.lamda_e)
            # qx_e    = 0.
            self.header.measurement_scheme = 'angle- and energy-dispersive'
        elif self.config.incidentAngle == IncidentAngle.nu:
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

    def filter_qz_range(self, norm):
        if self.config.qzRange[1]<0.5 and not norm:
            self.mask_e = np.logical_and(self.mask_e,
                                         (self.config.qzRange[0]<=self.qz_e) & (self.qz_e<=self.config.qzRange[1]))
        self.detZ_e = self.detZ_e[self.mask_e]
        self.lamda_e = self.lamda_e[self.mask_e]
        self.wallTime_e = self.wallTime_e[self.mask_e]



    def read_individual_header(self):
        self.chopperDistance = float(np.take(self.hdf['entry1/Amor/chopper/pair_separation'], 0))
        self.detectorDistance = float(np.take(self.hdf['entry1/Amor/detector/transformation/distance'], 0))
        self.chopperDetectorDistance = self.detectorDistance-float(np.take(self.hdf['entry1/Amor/chopper/distance'], 0))
        self.tofCut = const.lamdaCut*self.chopperDetectorDistance/const.hdm*1.e-13

        #TODO: 'undefined' is not orso compatible - but should be.
        polarizationConfigs = ['undefined', 'unpolarized', 'po', 'mo', 'op', 'pp', 'mp', 'om', 'pm', 'mm']
        try:
            self.mu   = float(np.take(self.hdf['/entry1/Amor/instrument_control_parameters/mu'], 0))
            self.nu   = float(np.take(self.hdf['/entry1/Amor/instrument_control_parameters/nu'], 0))
            self.kap  = float(np.take(self.hdf['/entry1/Amor/instrument_control_parameters/kappa'], 0))
            self.kad  = float(np.take(self.hdf['/entry1/Amor/instrument_control_parameters/kappa_offset'], 0))
            #self.kap  = float(np.take(self.hdf['/entry1/Amor/instrument_control_parameters/kap'], 0))
            #self.kad  = float(np.take(self.hdf['/entry1/Amor/instrument_control_parameters/kad'], 0))
            self.div  = float(np.take(self.hdf['/entry1/Amor/instrument_control_parameters/div'], 0))
            self.ch1TriggerPhase = float(np.take(self.hdf['/entry1/Amor/chopper/ch1_trigger_phase'], 0))
            self.ch2TriggerPhase = float(np.take(self.hdf['/entry1/Amor/chopper/ch2_trigger_phase'], 0))
            try: 
                chopperTriggerTime = (float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][7])\
                                     - float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][0]))\
                                     / 7
                self.tau = int(1e-6*chopperTriggerTime/2+0.5)*(1e-3)
                self.chopperSpeed = 30/self.tau
                chopperTriggerTimeDiff =  float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_offset'][2])
                chopperTriggerPhase = 180e-9*chopperTriggerTimeDiff/self.tau
                #TODO: check the next line
                self.chopperPhase = chopperTriggerPhase + self.ch1TriggerPhase - self.ch2TriggerPhase
                #print(f'chopperTriggerPhase: {chopperTriggerPhase} + {self.ch1TriggerPhase} - {self.ch2TriggerPhase} chopper phase: {self.chopperPhase}')
            except(KeyError, IndexError):
                logging.debug('      chopper speed and phase taken from .hdf file')
                self.chopperSpeed = float(np.take(self.hdf['/entry1/Amor/chopper/rotation_speed'], 0))
                self.chopperPhase = float(np.take(self.hdf['/entry1/Amor/chopper/phase'], 0))
                self.tau = 30/self.chopperSpeed
            try:
                polarizationConfigLabel = int(self.hdf['/entry1/Amor/polarization/configuration/value'][0])
            except(KeyError, IndexError):
                polarizationConfigLabel = 0
        except(KeyError, IndexError):
            logging.warning("     using parameters from nicos cache")
            year_date = str(self.start_date).replace('-', '/', 1)
            # TODO: check new cache pathes
            cachePath = '/home/amor/nicosdata/amor/cache/'
            grp = '/usr/bin/grep "value"'
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-mu/{year_date}')).split('\t')[-1]
            self.mu = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-nu/{year_date}')).split('\t')[-1]
            self.nu = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-kappa/{year_date}')).split('\t')[-1]
            self.kap = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-kad/{year_date}')).split('\t')[-1]
            self.kad = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-div/{year_date}')).split('\t')[-1]
            self.div = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-ch1_speed/{year_date}')).split('\t')[-1]
            self.chopperSpeed = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-chopper_phase/{year_date}')).split('\t')[-1]
            self.chopperPhase = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-ch1_trigger_phase/{year_date}')).split('\t')[-1]
            self.ch1TriggerPhase = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-ch2_trigger_phase/{year_date}')).split('\t')[-1]
            self.ch2TriggerPhase = float(value)
            value = str(subprocess.getoutput(f'{grp} {cachePath}nicos-polarizer_config_label/{year_date}')).split('\t')[-1]
            self. polarizationConfigLabel = int(value)
            
            self.tau     = 30. / self.chopperSpeed

        self.polarizationConfig = polarizationConfigs[polarizationConfigLabel]
        logging.debug(f'      polarization configuration: {self.polarizationConfig} (index {polarizationConfigLabel})')

        logging.debug(f'        tau = {self.tau:5.3f} s')
        if self.config.muOffset:
            logging.debug(f'        set muOffset = {self.config.muOffset}')
            self.mu += self.config.muOffset
        if self.config.mu:
            logging.debug(f'        replaced mu = {self.mu} with {self.config.mu}')
            self.mu = self.config.mu
        if self.config.nu:
            logging.debug(f'        replaced nu = {self.nu} with {self.config.nu}')
            self.nu = self.config.nu
        if self.config.chopperPhaseOffset:
            logging.debug(f'        replaced ch1TriggerPhase = {self.ch1TriggerPhase} with {self.config.chopperPhaseOffset}')
            self.ch1TriggerPhase = self.config.chopperPhaseOffset

        # extract start time as unix time, adding UTC offset of 1h to time string
        dz = datetime.fromisoformat(self.hdf['/entry1/start_time'][0].decode('utf-8'))
        self.fileDate=dz.replace(tzinfo=AMOR_LOCAL_TIMEZONE)
        #self.startTime = np.int64( (self.fileDate.timestamp() ) * 1e9 )
        #if self.seriesStartTime is None:
        #    self.seriesStartTime = self.startTime 

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

