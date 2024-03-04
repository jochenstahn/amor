import logging
import os
import platform
import subprocess
import sys
from datetime import datetime
from dataclasses import dataclass
from typing import Optional, Tuple

import h5py
import numpy as np
from orsopy import fileio

from . import __version__, const
from .instrument import Detector

@dataclass
class DataReaderConfig:
    year: int
    dataPath: str
    sampleModel: str

    chopperPhase: float
    yRange: Tuple[float, float]
    lambdaRange: Tuple[float, float]
    qzRange: Tuple[float, float]

    chopperPhaseOffset: float = 0.0
    mu: Optional[float] = None
    nu: Optional[float] = None
    muOffset: Optional[float] = None
    offSpecular: bool = False

class Header:
    """orso compatible output file header content"""

    def __init__(self):
        self.owner                           = None
        self.experiment                      = None
        self.sample                          = None
        self.measurement_instrument_settings = None
        self.measurement_scheme              = None
        self.measurement_data_files          = []
        self.measurement_additional_files    = []

        self.reduction = fileio.Reduction(
            software    = fileio.Software('eos', version=__version__),
            call        = ' '.join(sys.argv),
            computer    = platform.node(),
            timestamp   = datetime.now(),
            creator     = None,
            corrections = ['histogramming in lambda and alpha_f',
                           'gravity'],
            )
    #-------------------------------------------------------------------------------------------------
    def data_source(self):
        return fileio.DataSource(
            self.owner,
            self.experiment,
            self.sample,
            fileio.Measurement(
                instrument_settings = self.measurement_instrument_settings,
                scheme              = self.measurement_scheme,
                data_files          = self.measurement_data_files,
                additional_files    = self.measurement_additional_files,
                ),
        )
    #-------------------------------------------------------------------------------------------------
    def columns(self):
        cols = [
            fileio.Column('Qz', '1/angstrom', 'normal momentum transfer'),
            fileio.Column('R', '', 'specular reflectivity'),
            fileio.ErrorColumn(error_of='R', error_type='uncertainty', distribution='gaussian', value_is='sigma'),
            fileio.ErrorColumn(error_of='Qz', error_type='resolution', distribution='gaussian', value_is='sigma'),
            ]
        return cols


class AmorData:
    """read meta-data and event streams from .hdf file(s), apply filters and conversions"""
    #-------------------------------------------------------------------------------------------------
    def __init__(self, startTime, header:Header, config: DataReaderConfig):
        self.startTime = startTime
        self.header = header
        self.config = config
    #-------------------------------------------------------------------------------------------------
    def read_data(self, short_notation, norm=False):
        self.data_file_numbers = self.expand_file_list(short_notation)
        #self.year = clas.year
        self.file_list = []
        for number in self.data_file_numbers:
            self.file_list.append(self.path_generator(number))
        ## read specific meta data and measurement from first file
        if norm:
            self.readHeaderInfo = False
        else:
            self.readHeaderInfo = True

        _detZ_e            = []
        _lamda_e           = []
        _wallTime_e        = []
        for file in self.file_list:
            self.read_individual_data(file, norm)
            _detZ_e        = np.append(_detZ_e,        self.detZ_e)
            _lamda_e       = np.append(_lamda_e,       self.lamda_e)
            _wallTime_e    = np.append(_wallTime_e,    self.wallTime_e)
        self.detZ_e        = _detZ_e
        self.lamda_e       = _lamda_e
        self.wallTime_e    = _wallTime_e

    #-------------------------------------------------------------------------------------------------
    def path_generator(self, number):
        fileName = f'amor{self.config.year}n{number:06d}.hdf'
        if   os.path.exists(f'{self.config.dataPath}/{fileName}'):
            path = self.config.dataPath
        elif os.path.exists(fileName):
            path = '.'
        elif os.path.exists(f'./raw/{fileName}'):
            path = './raw'
        elif os.path.exists(f'../raw/{fileName}'):
            path = '../raw'
        elif os.path.exists(f'/afs/psi.ch/project/sinqdata/{self.config.year}/amor/{int(number/1000)}/{fileName}'):
            path = '/afs/psi.ch/project/sinqdata/{self.config.year}/amor/{int(number/1000)}'
        else:
            sys.exit(f'# ERROR: the file {fileName} is nowhere to be found!')
        return f'{path}/{fileName}'
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
        return sorted(file_list)
    #-------------------------------------------------------------------------------------------------
    def resolve_pixels(self):
        """determine spatial coordinats and angles from pixel number"""
        det = Detector()
        nPixel = det.nWires * det.nStripes * det.nBlades
        pixelID = np.arange(nPixel)
        (bladeNr, bPixel) = np.divmod(pixelID, det.nWires * det.nStripes)
        (bZi, detYi)      = np.divmod(bPixel, det.nStripes)                     # z index on blade, y index on detector
        detZi             = bladeNr * det.nWires + bZi                          # z index on detector
        detX              = bZi * det.dX                                        # x position in detector
        # detZ              = det.zero - bladeNr * det.bladeZ - bZi * det.dZ      # z position on detector
        bladeAngle        = np.rad2deg( 2. * np.arcsin(0.5*det.bladeZ / det.distance) )
        delta             = (det.nBlades/2. - bladeNr) * bladeAngle \
                            - np.rad2deg( np.arctan(bZi*det.dZ / ( det.distance + bZi * det.dX) ) )
        self.delta_z      = delta[detYi==1]
        return np.vstack((detYi.T, detZi.T, detX.T, delta.T)).T
        #return matr
    #-------------------------------------------------------------------------------------------------
    def read_individual_data(self, fileName, norm=False):
        pixelLookUp = self.resolve_pixels()

        self.hdf = h5py.File(fileName, 'r', swmr=True)

        if self.readHeaderInfo:
            # read general information and first data set
            logging.info(f'#   meta data from: {self.file_list[0]}')
            self.hdf = h5py.File(self.file_list[0], 'r', swmr=True)

            title            = self.hdf['entry1/title'][0].decode('utf-8')
            proposal_id      = self.hdf['entry1/proposal_id'][0].decode('utf-8')
            user_name        = self.hdf['entry1/user/name'][0].decode('utf-8')
            user_affiliation = 'unknown'
            user_email       = self.hdf['entry1/user/email'][0].decode('utf-8')
            user_orcid       = None
            sampleName       = self.hdf['entry1/sample/name'][0].decode('utf-8')
            model            = self.hdf['entry1/sample/model'][0].decode('utf-8')
            instrumentName   = 'Amor'
            source           = self.hdf['entry1/Amor/source/name'][0].decode('utf-8')
            sourceProbe      = 'neutron'

            start_time       = self.hdf['entry1/start_time'][0].decode('utf-8')
            start_date       = start_time.split(' ')[0]

            if self.config.sampleModel:
                model = self.config.sampleModel
            else:
                model = None

            # assembling orso header information
            self.header.owner = fileio.Person(
                name = user_name,
                affiliation = user_affiliation,
                contact             = user_email,
                )
            if user_orcid:
                self.header.owner.orcid  = user_orcid
            self.header.experiment = fileio.Experiment(
                title               = title,
                instrument          = instrumentName,
                start_date          = start_date,
                probe               = sourceProbe,
                facility            = source,
                proposalID          = proposal_id
                )
            self.header.sample = fileio.Sample(
                name                = sampleName,
                model               = model,
                sample_parameters   = None,
                )
            self.header.measurement_scheme     = 'angle- and energy-dispersive'

        self.chopperDistance         = float(np.take(self.hdf['entry1/Amor/chopper/pair_separation'], 0))
        self.detectorDistance        = float(np.take(self.hdf['entry1/Amor/detector/transformation/distance'], 0))
        self.chopperDetectorDistance = self.detectorDistance - float(np.take(self.hdf['entry1/Amor/chopper/distance'], 0))
        self.tofCut                  = const.lamdaCut * self.chopperDetectorDistance / const.hdm * 1.e-13

        logging.info(f'#   data from file: {fileName}')
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
            cachePath = '/home/amor/nicosdata/amor/cache/'
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

        # TODO:  figure out real stop time....
        self.ctime=(self.hdf['/entry1/Amor/detector/data/event_time_zero'][-1]
            - self.hdf['/entry1/Amor/detector/data/event_time_zero'][0]) / 1.e9
        fileDate = datetime.fromisoformat( self.hdf['/entry1/start_time'][0].decode('utf-8') )

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
        if norm:
            self.header.measurement_additional_files.append(fileio.File(file=fileName.split('/')[-1], timestamp=fileDate))
        else:
            self.header.measurement_data_files.append(fileio.File(file=fileName.split('/')[-1], timestamp=fileDate))
        logging.info(f'#     mu = {self.mu:6.3f}, nu = {self.nu:6.3f}, kap = {self.kap:6.3f}, kad = {self.kap:6.3f}')

        # TODO: should extract monitor from counts or beam current times time
        self.monitor1 = self.ctime
        self.monitor2=self.monitor1

        # read data event streams
        tof_e        = np.array(self.hdf['/entry1/Amor/detector/data/event_time_offset'][:])/1.e9
        pixelID_e    = np.array(self.hdf['/entry1/Amor/detector/data/event_id'][:], dtype=int)
        totalNumber = np.shape(tof_e)[0]

        wallTime_e = np.empty(totalNumber)
        dataPacket_p = np.array(self.hdf['/entry1/Amor/detector/data/event_index'][:], dtype=np.uint64)
        dataPacketTime_p = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][:], dtype=np.uint64)/1e9
        #for i, index in enumerate(dataPacket_p):
        #    wallTime_e[index:] = dataPacketTime_p[i]
        for i in range(len(dataPacket_p)-1):
            wallTime_e[dataPacket_p[i]:dataPacket_p[i+1]] = dataPacketTime_p[i]
        wallTime_e[dataPacket_p[-1]:] = dataPacketTime_p[-1]
        if not self.startTime and not norm:
            self.startTime = wallTime_e[0]
        wallTime_e -= self.startTime
        logging.debug(f'wall time from {wallTime_e[0]} to {wallTime_e[-1]}')

        # filter 'strange' tof times > 2 tau
        if True:
            filter_e   = (tof_e <= 2*self.tau)
            tof_e      = tof_e[filter_e]
            pixelID_e  = pixelID_e[filter_e]
            wallTime_e = wallTime_e[filter_e]
            if np.shape(filter_e)[0]-np.shape(tof_e)[0] > 0.5 :
                logging.warning(f'#    strange times: {np.shape(filter_e)[0]-np.shape(tof_e)[0]}')
        tof_e        = np.remainder( tof_e - self.tofCut + self.tau, self.tau) + self.tofCut  # tof shifted to 1 frame
        tof_e        = tof_e + self.tau * self.config.chopperPhaseOffset / 180. # correction for time offset between chopper pulse and tof zero

        # resolve pixel ID into y and z indicees, x position and angle
        (detY_e, detZ_e, detXdist_e, delta_e) = pixelLookUp[np.int_(pixelID_e)-1,:].T

        # define mask and filter y range
        mask_e = (self.config.yRange[0] <= detY_e) & (detY_e <= self.config.yRange[1])

        # correct tof for beam size effect at chopper:  t_cor = (delta / 180 deg) * tau
        # TODO: check for correctness
        if not self.config.offSpecular:
            tof_e    -= ( delta_e / 180. ) * self.tau

        # lambda
        lamda_e = 1.e13 * tof_e * const.hdm / (self.chopperDetectorDistance + detXdist_e)
        self.lamdaMax = const.lamdaCut + 1.e13 * self.tau * const.hdm / (self.chopperDetectorDistance + 124.)
        mask_e = np.logical_and(mask_e, (self.config.lambdaRange[0] <= lamda_e) & (lamda_e <= self.config.lambdaRange[1]))

        # alpha_f
        alphaF_e = self.nu - self.mu + delta_e

        # q_z
        if self.config.offSpecular:
            alphaI = self.kap + self.kad + self.mu
            qz_e = 2 * np.pi * ( np.sin( np.deg2rad(alphaF_e) ) + np.sin( np.deg2rad( alphaI ) ) ) / lamda_e
            qx_e = 2 * np.pi * ( np.cos( np.deg2rad(alphaF_e) ) - np.cos( np.deg2rad( alphaI ) ) ) / lamda_e
            self.header.measurement_scheme     = 'energy-dispersive',
        else:
            qz_e = 4 * np.pi * np.sin( np.deg2rad(alphaF_e) ) / lamda_e
            # qx_e = 0.
            self.header.measurement_scheme     = 'angle- and energy-dispersive'

        # filter q_z range
        if self.config.qzRange[1] < 0.3 and not norm:
            mask_e = np.logical_and(mask_e, (self.config.qzRange[0] <= qz_e) & (qz_e <= self.config.qzRange[1]))

        self.detZ_e        = detZ_e[mask_e]
        self.lamda_e       = lamda_e[mask_e]
        self.wallTime_e    = wallTime_e[mask_e]

        logging.info(f'#     number of events: total = {totalNumber:7d}, filtered = {np.shape(self.lamda_e)[0]:7d}')

