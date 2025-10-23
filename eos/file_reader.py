"""
Reading of Amor NeXus data files to extract metadata and event stream.
"""
from typing import BinaryIO, List, Union

import h5py
import numpy as np
import platform
import logging
import subprocess

from datetime import datetime

from orsopy import fileio
from orsopy.fileio.model_language import SampleModel

from . import const
from .header import Header
from .event_data_types import AmorGeometry, AmorTiming, AmorEventStream, PACKET_TYPE, EVENT_TYPE, PULSE_TYPE, PC_TYPE

try:
    import zoneinfo
except ImportError:
    # for python versions < 3.9 try to use the backports version
    from backports import zoneinfo


# Time zone used to interpret time strings
AMOR_LOCAL_TIMEZONE = zoneinfo.ZoneInfo(key='Europe/Zurich')

if  platform.node().startswith('amor'):
    NICOS_CACHE_DIR = '/home/amor/nicosdata/amor/cache/'
    GREP = '/usr/bin/grep "value"'
else:
    NICOS_CACHE_DIR = None

class AmorHeader:
    """
    Collects header information from Amor NeXus fiel without reading event data.
    """

    def __init__(self, fileName:Union[str, h5py.File, BinaryIO]):
        if type(fileName) is str:
            logging.warning(f'    {fileName.split("/")[-1]}')
            self.hdf = h5py.File(fileName, 'r', swmr=True)
        elif type(fileName) is h5py.File:
            self.hdf = fileName
        else:
            self.hdf = h5py.File(fileName, 'r')

        self.read_header_info()
        self.read_instrument_configuration()

        if type(fileName) is str:
            # close the input file to free memory, only if the file was opened in this object
            self.hdf.close()
        del(self.hdf)

    def _replace_if_missing(self, key, nicos_key, dtype=float, suffix=''):
        try:
            return dtype(self.hdf[f'/entry1/Amor/{key}'][0])
        except(KeyError, IndexError):
            if NICOS_CACHE_DIR:
                try:
                    logging.info(f"     using parameter {nicos_key} from nicos cache")
                    year_date = self.fileDate.strftime('%Y')
                    call = f'{GREP} {NICOS_CACHE_DIR}nicos-{nicos_key}/{year_date}{suffix}'
                    value = str(subprocess.getoutput(call)).split('\t')[-1]
                    return dtype(value)
                except Exception:
                    logging.error(f"Couldn't get value from nicos cache {nicos_key}, {call}")
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
        if 'stack:' in model:
            import yaml
            model = yaml.safe_load(model)
        else:
            model = dict(stack=model)
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
                model=SampleModel.from_dict(model),
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
        kad = self._replace_if_missing('instrument_control_parameters/kappa_offset', 'kappa_offset', float)
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

        polarizationConfigLabel = self._replace_if_missing('polarization/configuration/average_value', 'polarization_config_label', int, suffix='/*')
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
        logging.info(f'    meta data from: {self.file_list[0]}')
        header.owner = self.owner
        header.experiment = self.experiment
        header.sample = self.sample
        header.measurement_instrument_settings = self.instrument_settings


class AmorEventData(AmorHeader):
    """
    Read one amor NeXus datafile and extract relevant header information.

    Implements EventDatasetProtocol
    """
    file_list: List[str]
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

    def __init__(self, fileName:Union[str, h5py.File, BinaryIO], first_index:int=0, max_events:int=100_000_000):
        if type(fileName) is str:
            logging.warning(f'    {fileName.split("/")[-1]}')
            self.file_list = [fileName]
            hdf = h5py.File(fileName, 'r', swmr=True)
        elif type(fileName) is h5py.File:
            self.file_list = [fileName.filename]
            hdf = fileName
        else:
            self.file_list = [repr(fileName)]
            hdf = h5py.File(fileName, 'r')
        self.first_index = first_index
        self.max_events = max_events

        super().__init__(hdf)
        self.hdf = hdf
        self.read_event_stream()

        if type(fileName) is str:
            # close the input file to free memory, only if the file was opened in this object
            self.hdf.close()
        del(self.hdf)


    def read_event_stream(self):
        """
        Read the actual event data from file. If file is too large, find event index from packets
        that allow splitting of file smaller than self.max_events.
        """
        packets = np.recarray(self.hdf['/entry1/Amor/detector/data/event_index'].shape, dtype=PACKET_TYPE)
        packets.start_index = self.hdf['/entry1/Amor/detector/data/event_index'][:]
        packets.time = self.hdf['/entry1/Amor/detector/data/event_time_zero'][:]
        try:
            # packet index that matches first event index
            start_packet = int(np.where(packets.start_index==self.first_index)[0][0])
        except IndexError:
            raise EOFError(f'No event packet found starting at event #{self.first_index}, '
                           f'number of events is {self.hdf["/entry1/Amor/detector/data/event_time_offset"].shape[0]}')
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
        events.mask = 0

        pulses = self.read_chopper_trigger_stream(packets)
        current = self.read_proton_current_stream(packets)
        self.data = AmorEventStream(events, packets, pulses, current)

        if self.first_index>0 and not self.EOF:
            # label the file name if not all events were used
            self.file_list[0] += f'[{self.first_index}:{self.last_index}]'

    def read_chopper_trigger_stream(self, packets):
        chopper1TriggerTime = np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][:-2], dtype=np.int64)
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
        pulses.monitor = 1. # default is monitor pulses as it requires no calculation
        # apply filter in case the events were filtered
        if self.first_index>0 or not self.EOF:
            pulses = pulses[(pulses.time>=packets.time[0])&(pulses.time<=packets.time[-1])]
        self.eventStartTime = startTime
        return pulses

    def read_proton_current_stream(self, packets):
        proton_current = np.recarray(self.hdf['entry1/Amor/detector/proton_current/time'].shape, dtype=PC_TYPE)
        proton_current.time = self.hdf['entry1/Amor/detector/proton_current/time'][:]
        proton_current.current = self.hdf['entry1/Amor/detector/proton_current/value'][:,0]
        if self.first_index>0 or not self.EOF:
            proton_current = proton_current[(proton_current.time>=packets.time[0])&
                                            (proton_current.time<=packets.time[-1])]
        return proton_current

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
        self.file_list += other.file_list

    def __repr__(self):
        output = (f"AmorEventData({self.file_list!r}) # {self.data.events.shape[0]} events, "
                  f"{self.data.pulses.shape[0]} pulses")

        return output

    def get_timeslice(self, start, end)->'AmorEventData':
        # return a new dataset with just events that occured in given time slice
        if not 'wallTime' in self.data.events.dtype.names:
            raise ValueError("This dataset is missing a wallTime that is required for time slicing")
        # convert from seconds to epoch integer values
        start , end = start*1e9, end*1e9
        event_filter = self.data.events.wallTime>=start
        event_filter &= self.data.events.wallTime<end
        pulse_filter = self.data.pulses.time>=start
        pulse_filter &= self.data.pulses.time<end
        output = super().__new__(AmorEventData)
        for key, value in self.__dict__.items():
            if key == 'data':
                continue
            else:
                setattr(output, key, value)
        # TODO: this is not strictly correct, as the packet/event relationship is lost
        output.data = AmorEventStream(self.data.events[event_filter], self.data.packets,
                                      self.data.pulses[pulse_filter], self.data.proton_current)
        return output
