"""
Reading of Amor NeXus data files to extract metadata and event stream.
"""
from typing import BinaryIO, List, Union

import h5py
import numpy as np
import platform
import logging
import subprocess
import sys
import os

from datetime import datetime

from orsopy import fileio
from orsopy.fileio.model_language import SampleModel

from . import const, event_handling as eh, event_analysis as ea
from .header import Header
from .helpers import extract_walltime
from .event_data_types import AmorGeometry, AmorTiming, AmorEventStream, PACKET_TYPE, EVENT_TYPE, PULSE_TYPE, PC_TYPE
from .options import ExperimentConfig, IncidentAngle, ReaderConfig

try:
    import zoneinfo
except ImportError:
    # for python versions < 3.9 try to use the backports version
    from backports import zoneinfo


# Time zone used to interpret time strings
AMOR_LOCAL_TIMEZONE = zoneinfo.ZoneInfo(key='Europe/Zurich')

if  platform.node().startswith('amor'):
    NICOS_CACHE_DIR = '/home/amor/nicosdata/amor/cache/'
    GREP = '/usr/bin/grep "%s"'
else:
    NICOS_CACHE_DIR = None

class PathResolver:
    def __init__(self, year, rawPath):
        self.year = year
        self.rawPath = rawPath

    def resolve(self, short_notation):
        return list(map(self.get_path, self.expand_file_list(short_notation)))

    def expand_file_list(self, short_notation):
        """Evaluate string entry for file number lists"""
        file_list = []
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
        return list(sorted(file_list))

    def get_path(self, number):
        fileName = f'amor{self.year}n{number:06d}.hdf'
        path = ''
        for rawd in self.rawPath:
            if os.path.exists(os.path.join(rawd, fileName)):
                path = rawd
                break
        if not path:
            if os.path.exists(
                    f'/afs/psi.ch/project/sinqdata/{self.year}/amor/{int(number/1000)}/{fileName}'):
                path = f'/afs/psi.ch/project/sinqdata/{self.year}/amor/{int(number/1000)}'
            else:
                sys.exit(f'# ERROR: the file {fileName} can not be found in {self.rawPath}')
        return os.path.join(path, fileName)

class AmorEventData:
    """
    Read one amor NeXus datafile and extract relevant header information.

    Implements EventDatasetProtocol
    """
    fileName: str
    first_index: int
    last_index: int = -1
    EOF: bool = False
    max_events: int
    owner: fileio.Person
    experiment: fileio.Experiment
    sample: fileio.Sample
    instrument_settings: fileio.InstrumentSettings
    geometry: AmorGeometry
    timing: AmorTiming
    data: AmorEventStream

    eventStartTime: np.int64

    def __init__(self, fileName:Union[str, h5py.File, BinaryIO], first_index:int=0, max_events:int=1_000_000):
        if type(fileName) is str:
            self.fileName = fileName
            self.hdf = h5py.File(fileName, 'r', swmr=True)
        elif type(fileName) is h5py.File:
            self.fileName = fileName.filename
            self.hdf = fileName
        else:
            self.fileName = repr(fileName)
            self.hdf = h5py.File(fileName, 'r')
        self.first_index = first_index
        self.max_events = max_events

        self.read_header_info()
        self.read_instrument_configuration()
        self.read_event_stream()

        # actions applied to any dataset
        self.correct_for_chopper_phases()
        self.read_chopper_trigger_stream()
        self.extract_walltime()
        self.read_proton_current_stream()

        if type(fileName) is str:
            # close the input file to free memory, only if the file was opened in this object
            self.hdf.close()
        del(self.hdf)

    def _replace_if_missing(self, key, nicos_key, dtype=float):
        try:
            return dtype(self.hdf[f'/entry1/Amor/{key}'][0])
        except(KeyError, IndexError):
            if NICOS_CACHE_DIR:
                try:
                    logging.warning(f"     using parameter {nicos_key} from nicos cache")
                    year_date = self.fileDate.strftime('%Y')
                    value = str(subprocess.getoutput(f'{GREP} {NICOS_CACHE_DIR}nicos-{nicos_key}/{year_date}')).split('\t')[-1]
                    return dtype(value)
                except Exception:
                    logging.error("Couldn't get value from nicos cache", exc_info=True)
                    return dtype(0)
            else:
                logging.warning(f"     parameter {key} not found, relpace by zero")
                return dtype(0)

    def read_header_info(self):
        # read general information and first data set
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
        # extract start time as unix time, adding UTC offset of 1h to time string
        start_date = datetime.fromisoformat(start_time)
        self.fileDate = start_date.replace(tzinfo=AMOR_LOCAL_TIMEZONE)

        self.owner = fileio.Person(
                name=user_name,
                affiliation=user_affiliation,
                contact=user_email,
                )
        if user_orcid:
            self.owner.orcid = user_orcid

        self.experiment = fileio.Experiment(
                title=title,
                instrument=instrumentName,
                start_date=start_date,
                probe=sourceProbe,
                facility=source,
                proposalID=proposal_id
                )
        self.sample = fileio.Sample(
                name=sampleName,
                model=SampleModel(stack=model),
                sample_parameters=None,
                )

    def read_instrument_configuration(self):
        chopperSeparation = float(np.take(self.hdf['entry1/Amor/chopper/pair_separation'], 0))
        detectorDistance = float(np.take(self.hdf['entry1/Amor/detector/transformation/distance'], 0))
        chopperDetectorDistance = detectorDistance-float(np.take(self.hdf['entry1/Amor/chopper/distance'], 0))

        polarizationConfigs = ['unpolarized', 'unpolarized', 'po', 'mo', 'op', 'pp', 'mp', 'om', 'pm', 'mm']

        mu = self._replace_if_missing('instrument_control_parameters/mu', 'mu', float)
        nu = self._replace_if_missing('instrument_control_parameters/nu', 'nu', float)
        kap = self._replace_if_missing('instrument_control_parameters/kappa', 'kappa', float)
        kad = self._replace_if_missing('instrument_control_parameters/kappa_offset', 'kad', float)
        div = self._replace_if_missing('instrument_control_parameters/div', 'div', float)
        ch1TriggerPhase = self._replace_if_missing('chopper/ch1_trigger_phase', 'ch1_trigger_phase', float)
        ch2TriggerPhase = self._replace_if_missing('chopper/ch2_trigger_phase', 'ch2_trigger_phase', float)
        try:
            chopperTriggerTime = (float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][7]) \
                                  -float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][0])) \
                                 /7
            chopperTriggerTimeDiff = float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_offset'][2])
        except (KeyError, IndexError):
            logging.debug('      chopper speed and phase taken from .hdf file')
            chopperSpeed = self._replace_if_missing('chopper/rotation_speed', 'chopper_phase', float)
            chopperPhase = self._replace_if_missing('chopper/phase', 'chopper_phase', float)
            tau = 30/chopperSpeed
        else:
            tau = int(1e-6*chopperTriggerTime/2+0.5)*(1e-3)
            chopperTriggerPhase = 180e-9*chopperTriggerTimeDiff/tau
            chopperSpeed = 30/tau
            chopperPhase = chopperTriggerPhase+ch1TriggerPhase-ch2TriggerPhase

        self.geometry = AmorGeometry(mu, nu, kap, kad, div,
                                     chopperSeparation, detectorDistance, chopperDetectorDistance)
        self.timing = AmorTiming(ch1TriggerPhase, ch2TriggerPhase, chopperSpeed, chopperPhase, tau)

        polarizationConfigLabel = self._replace_if_missing('polarization/configuration/value', 'polarizer_config_label', int)
        polarizationConfig = fileio.Polarization(polarizationConfigs[polarizationConfigLabel])
        logging.debug(f'      polarization configuration: {polarizationConfig} (index {polarizationConfigLabel})')


        self.instrument_settings = fileio.InstrumentSettings(
            incident_angle = fileio.ValueRange(round(mu+kap+kad-0.5*div, 3),
                                               round(mu+kap+kad+0.5*div, 3),
                                               'deg'),
            wavelength = fileio.ValueRange(const.lamdaCut, const.lamdaMax, 'angstrom'),
            #polarization = fileio.Polarization.unpolarized,
            polarization = fileio.Polarization(polarizationConfig)
            )
        self.instrument_settings.mu = fileio.Value(
                round(mu, 3),
                'deg',
                comment='sample angle to horizon')
        self.instrument_settings.nu = fileio.Value(
                round(nu, 3),
                'deg',
                comment='detector angle to horizon')
        self.instrument_settings.div = fileio.Value(
                round(div, 3),
                'deg',
                comment='incoming beam divergence')
        self.instrument_settings.kap = fileio.Value(
                round(kap, 3),
                'deg',
                comment='incoming beam inclination')
        if abs(kad)>0.02:
            self.instrument_settings.kad = fileio.Value(
                    round(kad, 3),
                    'deg',
                    comment='incoming beam angular offset')


    def update_header(self, header:Header):
        """
        Add dataset information into an existing header.
        """
        header.owner = self.owner
        header.experiment = self.experiment
        header.sample = self.sample
        header.measurement_instrument_settings = self.instrument_settings

    def read_event_stream(self):
        """
        Read the actual event data from file. If file is too large, find event index from packets
        that allow splitting of file smaller than self.max_events.
        """
        packets = np.recarray(self.hdf['/entry1/Amor/detector/data/event_index'].shape, dtype=PACKET_TYPE)
        packets.start_index = self.hdf['/entry1/Amor/detector/data/event_index'][:]
        packets.Time = self.hdf['/entry1/Amor/detector/data/event_time_zero'][:]
        try:
            # packet index that matches first event index
            start_packet = int(np.where(packets.start_index==self.first_index)[0][0])
        except IndexError:
            raise IndexError(f'No event packet found starting at event #{self.first_index}')
        packets = packets[start_packet:]

        nevts = self.hdf['/entry1/Amor/detector/data/event_time_offset'].shape[0]
        if (nevts-self.first_index)>self.max_events:
            end_packet = np.where(packets.start_index<=(self.first_index+self.max_events))[0][-1]
            self.last_index = packets.start_index[end_packet]-1
            packets = packets[:end_packet]
        else:
            self.last_index = nevts-1
            self.EOF = True
        nevts = self.last_index+1-self.first_index

        # adapte packet to event index relation
        packets.start_index -= self.first_index

        events = np.recarray(nevts, dtype=EVENT_TYPE)
        events.tof = np.array(self.hdf['/entry1/Amor/detector/data/event_time_offset'][self.first_index:self.last_index+1])/1.e9
        events.pixelID = self.hdf['/entry1/Amor/detector/data/event_id'][self.first_index:self.last_index+1]
        self.data = AmorEventStream(events, packets)

    def correct_for_chopper_phases(self):
        self.data.events.tof += self.timing.tau*(self.timing.ch1TriggerPhase-self.timing.chopperPhase/2)/180

    def read_chopper_trigger_stream(self):
        chopper1TriggerTime = np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][:-2],
                                            dtype=np.int64)
        #self.chopper2TriggerTime = self.chopper1TriggerTime + np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time'][:-2], dtype=np.int64)
        #                           + np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_offset'][:], dtype=np.int64)
        if np.shape(chopper1TriggerTime)[0] > 2:
            startTime = chopper1TriggerTime[0]
            pulseTimeS = chopper1TriggerTime
        else:
            logging.warn('     no chopper trigger data available, using event steram instead')
            startTime = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][0], dtype=np.int64)
            stopTime = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][-2], dtype=np.int64)
            pulseTimeS = np.arange(startTime, stopTime, self.timing.tau*1e9, dtype=np.int64)
        pulses = np.recarray(pulseTimeS.shape, dtype=PULSE_TYPE)
        pulses.time = pulseTimeS
        # apply filter in case the events were filtered
        if self.first_index>0 or not self.EOF:
            pulses = pulses[(pulses.time>=self.data.packets.Time[0])&(pulses.time<=self.data.packets.Time[-1])]
        self.data.pulses = pulses
        self.eventStartTime = startTime

    def extract_walltime(self):
        # TODO: fix numba type definition after refactor
        self.data.events.wallTime = extract_walltime(self.data.events.tof,
                                                        self.data.packets.start_index.astype(np.uint64),
                                                        self.data.packets.Time)

    def read_proton_current_stream(self):
        proton_current = np.recarray(self.hdf['entry1/Amor/detector/proton_current/time'].shape, dtype=PC_TYPE)
        proton_current.time = self.hdf['entry1/Amor/detector/proton_current/time'][:]
        proton_current.current = self.hdf['entry1/Amor/detector/proton_current/value'][:,0]
        if self.first_index>0 or not self.EOF:
            proton_current = proton_current[(proton_current.time>=self.data.packets.Time[0])&
                                            (proton_current.time<=self.data.packets.Time[-1])]
        self.data.proton_current = proton_current

    def info(self):
        output = ""
        for key in ['owner', 'experiment', 'sample', 'instrument_settings']:
            value = repr(getattr(self, key)).replace("\n","\n      ")
            output += f'\n{key}={value},'
        output += '\n'
        return output

    def append(self, other):
        """
        Append event streams from another file to this one. Adjusts the event indices in the
        packets to stay valid.
        """
        new_events = np.concatenate([self.data.events, other.data.events]).view(np.recarray)
        new_pulses = np.concatenate([self.data.pulses, other.data.pulses]).view(np.recarray)
        new_proton_current = np.concatenate([self.data.proton_current, other.data.proton_current]).view(np.recarray)
        new_packets = np.concatenate([self.data.packets, other.data.packets]).view(np.recarray)
        new_packets.start_index[self.data.packets.shape[0]:] += self.data.events.shape[0]
        self.data = AmorEventStream(new_events, new_packets, new_pulses, new_proton_current)
        # Indicate that this is amodified dataset, basically counts number of appends as negative indices
        self.last_index = min(self.last_index-1, -1)


    def __repr__(self):
        output = (f"AmorEventData({self.fileName!r}) # {self.data.events.shape[0]} events, "
                  f"{self.data.pulses.shape[0]} pulses")

        return output

class AmorData:
    """re-implement old AmorData class functionality until refactoring is complete"""
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
    # monitor: float
    mu: float
    nu: float
    tau: float
    tofCut: float
    start_date: str
    monitorType: str

    seriesStartTime = None

    # -------------------------------------------------------------------------------------------------
    def __init__(self, header: Header, reader_config: ReaderConfig, config: ExperimentConfig,
                 short_notation: str, norm=False):
        # self.startTime = reader_config.startTime
        self.header = header
        self.config = config
        self.reader_config = reader_config
        resolver = PathResolver(reader_config.year, reader_config.rawPath)
        self.file_list = resolver.resolve(short_notation)
        self.norm = norm
        self.prepare_actions()
        self.process()
        self.assign()

    def prepare_actions(self):
        # setup all actions performed in origin AmorData, time correction requires first dataset start time
        self.event_actions = eh.AssociatePulseWithMonitor(self.config.monitorType, self.config.lowCurrentThreshold)
        self.event_actions |= eh.FilterStrangeTimes()
        self.event_actions |= eh.MergeFrames()
        self.event_actions |= ea.AnalyzePixelIDs(self.config.yRange)
        self.event_actions |= ea.TofTimeCorrection(self.config.incidentAngle==IncidentAngle.alphaF)
        self.event_actions |= ea.WavelengthAndQ(self.config.lambdaRange, self.config.incidentAngle)
        self.event_actions |= ea.FilterQzRange(self.config.qzRange)
        self.event_actions |= ea.ApplyMask()

    def process(self):
        self.dataset = AmorEventData(self.file_list[0])
        time_correction = eh.CorrectSeriesTime(self.dataset.eventStartTime)
        time_correction(self.dataset)
        self.event_actions(self.dataset)
        if not self.norm:
            self.dataset.update_header(self.header)
            time_correction.update_header(self.header)
            self.event_actions.update_header(self.header)
        for fi in self.file_list[1:]:
            di = AmorEventData(fi)
            time_correction(di)
            self.event_actions(di)
            self.dataset.append(di)

        for fileName in self.file_list:
            if self.norm:
                self.header.measurement_additional_files.append(fileio.File(
                        file=fileName.split('/')[-1],
                        timestamp=self.dataset.fileDate))
            else:
                self.header.measurement_data_files.append(fileio.File(
                        file=fileName.split('/')[-1],
                        timestamp=self.dataset.fileDate))

    def assign(self):
        # assigne old class parameters from dataset object.
        ds = self.dataset
        ev = ds.data.events
        self.detZ_e = ev.detZ
        self.lamda_e = ev.lamda
        self.wallTime_e = ev.wallTime
        self.qz_e = ev.qz
        self.qx_e = ev.qx
        self.pulseTimeS = ds.data.pulses.time
        self.monitorPerPulse = ds.data.pulses.monitor

        for key, value in ds.geometry.__dict__.items():
            setattr(self, key, value)
        for key, value in ds.timing.__dict__.items():
            setattr(self, key, value)
        self.startTime = ds.eventStartTime
