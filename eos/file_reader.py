"""
Reading of Amor NeXus data files to extract metadata and event stream.
"""
from typing import BinaryIO, List, Union

import h5py
import numpy as np
import logging
import subprocess

from datetime import datetime

from orsopy import fileio
from orsopy.fileio.model_language import SampleModel

from . import const
from .header import Header
from .event_data_types import AmorGeometry, AmorTiming, AmorEventStream, LOG_TYPE, PACKET_TYPE, EVENT_TYPE, PULSE_TYPE, \
    PC_TYPE

try:
    import zoneinfo
except ImportError:
    # for python versions < 3.9 try to use the backports version
    from backports import zoneinfo


# Time zone used to interpret time strings
AMOR_LOCAL_TIMEZONE = zoneinfo.ZoneInfo(key='Europe/Zurich')

class AmorHeader:
    """
    Collects header information from Amor NeXus fiel without reading event data.
    """
    # mapping of names to (hdf_path, dtype, nicos_key[, suffix])
    hdf_paths = dict(
            title=('entry1/title', str),
            proposal_id=('entry1/proposal_id', str),
            user_name=('entry1/user/name', str),
            user_email=('entry1/user/email', str),
            sample_name=('entry1/sample/name', str),
            source_name=('entry1/Amor/source/name', str),
            sample_model=('entry1/sample/model', str),
            start_time=('entry1/start_time', str),
            start_time_fallback=('entry1/Amor/instrument_control_parameters/start_time', str),

            chopper_separation=('entry1/Amor/chopper/pair_separation', float),
            detector_distance=('entry1/Amor/detector/transformation/distance', float),
            chopper_distance=('entry1/Amor/chopper/distance', float),
            sample_temperature=('entry1/sample/temperature', float, 'ignore'),
            sample_magnetic_field=('entry1/sample/magnetic_field', float, 'ignore'),

            mu=('entry1/Amor/instrument_control_parameters/mu', float, 'mu'),
            nu=('entry1/Amor/instrument_control_parameters/nu', float, 'nu'),
            kap=('entry1/Amor/instrument_control_parameters/kappa', float, 'kappa'),
            kad=('entry1/Amor/instrument_control_parameters/kappa_offset', float, 'kappa_offset'),
            div=('entry1/Amor/instrument_control_parameters/div', float, 'div'),
            ch1_trigger_phase=('entry1/Amor/chopper/ch1_trigger_phase', float, 'ch1_trigger_phase'),
            ch2_trigger_phase=('entry1/Amor/chopper/ch2_trigger_phase', float, 'ch2_trigger_phase'),
            chopper_speed=('entry1/Amor/chopper/rotation_speed', float, 'chopper_phase'),
            chopper_phase=('entry1/Amor/chopper/phase', float, 'chopper_phase'),
            polarization_config_label=('entry1/Amor/polarization/configuration', int, 'polarization_config_label', '/*'),
            )

    def __init__(self, fileName:Union[str, h5py.File, BinaryIO]):
        if type(fileName) is str:
            logging.warning(f'    {fileName.split("/")[-1]}')
            self.hdf = h5py.File(fileName, 'r', swmr=True)
        elif type(fileName) is h5py.File:
            self.hdf = fileName
        else:
            self.hdf = h5py.File(fileName, 'r')

        self._log_keys = []

        self.read_header_info()
        self.read_instrument_configuration()

        if type(fileName) is str:
            # close the input file to free memory, only if the file was opened in this object
            self.hdf.close()
        del(self.hdf)

    def _replace_if_missing(self, key, nicos_key, dtype=float, suffix=''):
        from .nicos_interface import lookup_nicos_value
        year = self.fileDate.strftime('%Y')
        return lookup_nicos_value(key, nicos_key, dtype, suffix, year)

    def rv(self, key):
        """
        Generic read value methos based on key in hdf_paths dictionary.
        """
        hdf_path, dtype, *nicos = self.hdf_paths[key]
        try:
            hdfgrp = self.hdf[hdf_path]
            if hdfgrp.attrs.get('NX_class', None) == 'NXlog':
                self._log_keys.append(key)
                # use the last value that was recoreded before the count started
                time_column = hdfgrp['time'][:]
                try:
                    start_index = np.where(time_column<=self._start_time_ns)[0][0]
                except IndexError:
                    start_index = 0
                if hdfgrp['value'].ndim==1:
                    return dtype(hdfgrp['value'][start_index])
                else:
                    return dtype(hdfgrp['value'][start_index, 0])
            elif dtype is str:
                return self.read_string(hdf_path)
            else:
                if len(hdfgrp.shape)==1:
                    return dtype(hdfgrp[0])
                else:
                    return dtype(hdfgrp[()])
        except (KeyError, IndexError):
            if nicos:
                nicos_key = nicos[0]
                suffix = nicos[1] if len(nicos)>1 else ''
                return self._replace_if_missing(key, nicos_key, dtype, suffix)
            else:
                raise


    def read_string(self, path):
        if not np.shape(self.hdf[path]):
            return self.hdf[path][()].decode('utf-8')
        else:
            # format until 2025
            return self.hdf[path][0].decode('utf-8')

    def read_header_info(self):
        self._start_time_ns = np.uint64(0)
        try:
            start_time = self.rv('start_time')
        except KeyError:
            start_time = self.rv('start_time_fallback')

        # extract start time as unix time, adding UTC offset of 1h to time string
        start_date = datetime.fromisoformat(start_time)
        self.fileDate = start_date.replace(tzinfo=AMOR_LOCAL_TIMEZONE)
        self._start_time_ns = np.uint64(self.fileDate.timestamp()*1e9)

        # read general information and first data set
        title = self.rv('title')
        proposal_id = self.rv('proposal_id')
        user_name = self.rv('user_name')
        user_affiliation = 'unknown'
        user_email = self.rv('user_email')
        user_orcid = None
        sampleName = self.rv('sample_name')
        instrumentName = 'Amor'
        source = self.rv('source_name')
        sourceProbe = 'neutron'
        model = self.rv('sample_model')
        if 'stack' in model:
            import yaml
            model = yaml.safe_load(model)
        else:
            model = dict(stack=model)

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
        if model['stack'] == '':
            om = None
        else:
            om = SampleModel.from_dict(model)
        self.sample = fileio.Sample(
                name=sampleName,
                model=om,
                sample_parameters={},
                )
        # while event times are not evaluated, use average_value reported in file for SEE
        if self.hdf['entry1/sample'].get('temperature', None) is not None:
            sample_temperature = self.rv('sample_temperature')
            self.sample.sample_parameters['temperature'] = fileio.Value(sample_temperature, unit='K')
        if self.hdf['entry1/sample'].get('magnetic_field', None) is not None:
            sample_magnetic_field = self.rv('sample_magnetic_field')
            self.sample.sample_parameters['magnetic_field'] = fileio.Value(sample_magnetic_field, unit='T')

    def read_instrument_configuration(self):
        chopperSeparation = self.rv('chopper_separation')
        detectorDistance = self.rv('detector_distance')
        chopperDistance = self.rv('chopper_distance')
        chopperDetectorDistance = detectorDistance - chopperDistance

        polarizationConfigs = ['unpolarized', 'unpolarized', 'po', 'mo', 'op', 'pp', 'mp', 'om', 'pm', 'mm']

        mu = self.rv('mu')
        nu = self.rv('nu')
        kap = self.rv('kap')
        kad = self.rv('kad')
        div = self.rv('div')
        ch1TriggerPhase = self.rv('ch1_trigger_phase')
        ch2TriggerPhase = self.rv('ch2_trigger_phase')
        try:
            chopperTriggerTime = (float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][7]) \
                                  -float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][0])) \
                                 /7
            chopperTriggerTimeDiff = float(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_offset'][2])
        except (KeyError, IndexError):
            logging.debug('      chopper speed and phase taken from .hdf file')
            chopperSpeed = self.rv('chopper_speed')
            chopperPhase = self.rv('chopper_phase')
            tau = 30/chopperSpeed
        else:
            tau = int(1e-6*chopperTriggerTime/2+0.5)*(1e-3)
            chopperTriggerPhase = 180e-9*chopperTriggerTimeDiff/tau
            chopperSpeed = 30/tau
            chopperPhase = chopperTriggerPhase+ch1TriggerPhase-ch2TriggerPhase

        self.geometry = AmorGeometry(mu, nu, kap, kad, div,
                                     chopperSeparation, detectorDistance, chopperDetectorDistance)
        self.timing = AmorTiming(ch1TriggerPhase, ch2TriggerPhase, chopperSpeed, chopperPhase, tau)

        polarizationConfigLabel = self.rv('polarization_config_label')
        polarizationConfig = fileio.Polarization(polarizationConfigs[polarizationConfigLabel])
        logging.debug(f'      polarization configuration: {polarizationConfig} (index {polarizationConfigLabel})')


        self.instrument_settings = fileio.InstrumentSettings(
            incident_angle = fileio.ValueRange(round(mu+kap+kad-0.5*div, 3),
                                               round(mu+kap+kad+0.5*div, 3),
                                               'deg'),
            wavelength = fileio.ValueRange(const.lamdaCut, const.lamdaMax, 'angstrom'),
            polarization = fileio.Polarization(polarizationConfig)
            )
        self.instrument_settings.qz = fileio.ValueRange(round(4*np.pi*np.sin(np.deg2rad(mu+kap+kad-0.5*div))/const.lamdaMax, 3),
                                                        round(4*np.pi*np.sin(np.deg2rad(mu+kap+kad+0.5*div))/const.lamdaCut, 3),
                                                        '1/angstrom')
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
        try:
            self.read_event_stream()
        except EOFError:
            self.hdf.close()
            del(self.hdf)
            raise
        self.read_log_streams()

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
        if packets.shape[0]==0:
            raise EOFError(f'No more packets left after start_packet filter')

        nevts = self.hdf['/entry1/Amor/detector/data/event_time_offset'].shape[0]
        if (nevts-self.first_index)>self.max_events:
            end_packet = np.where(packets.start_index<=(self.first_index+self.max_events))[0][-1]
            end_packet = max(1, end_packet)
            if len(packets)==1:
                self.last_index = nevts-1
            else:
                self.last_index = packets.start_index[end_packet]-1
            packets = packets[:end_packet]
        else:
            self.last_index = nevts-1
            self.EOF = True

        if packets.shape[0]==0:
            raise EOFError(f'No more packets left after end_packet filter, first_index={self.first_index}, '
                           f'max_events={self.max_events}, nevts={nevts}')

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

    def read_log_streams(self):
        """
        Read the individual NXlog datasets that can later be used for filtering etc.
        """
        for key in self._log_keys:
            hdf_path, dtype, *_ = self.hdf_paths[key]
            hdfgroup = self.hdf[hdf_path]
            shape = hdfgroup['time'].shape
            data = np.recarray(shape, dtype=LOG_TYPE)
            data.time = hdfgroup['time'][:]
            if len(hdfgroup['value'].shape)==1:
                data.value = hdfgroup['value'][:]
            else:
                data.value = hdfgroup['value'][:, 0]
            self.data.device_logs[key] = data

    def read_chopper_trigger_stream(self, packets):
        chopper1TriggerTime = np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_zero'][:-2], dtype=np.int64)
        #self.chopper2TriggerTime = self.chopper1TriggerTime + np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time'][:-2], dtype=np.int64)
        #                           + np.array(self.hdf['entry1/Amor/chopper/ch2_trigger/event_time_offset'][:], dtype=np.int64)
        if np.shape(chopper1TriggerTime)[0] > 2:
            startTime = chopper1TriggerTime[0]
            pulseTimeS = chopper1TriggerTime
        else:
            logging.critical('     No chopper trigger data available, using event steram instead, pulse filtering will fail!')
            startTime = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][0], dtype=np.int64)
            stopTime = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][-2], dtype=np.int64)
            pulseTimeS = np.arange(startTime, stopTime, self.timing.tau*1e9, dtype=np.int64)
        pulses = np.recarray(pulseTimeS.shape, dtype=PULSE_TYPE)
        pulses.time = pulseTimeS
        pulses.monitor = 1. # default is monitor pulses as it requires no calculation
        # apply filter in case the events were filtered
        if (self.first_index>0 or not self.EOF):
            pulses = pulses[(pulses.time>=packets.time[0])&(pulses.time<=packets.time[-1])]
        self.eventStartTime = startTime
        return pulses

    def read_proton_current_stream(self, packets):
        proton_current = np.recarray(self.hdf['entry1/Amor/detector/proton_current/time'].shape, dtype=PC_TYPE)
        proton_current.time = self.hdf['entry1/Amor/detector/proton_current/time'][:]
        if self.hdf['entry1/Amor/detector/proton_current/value'].ndim==1:
            proton_current.current = self.hdf['entry1/Amor/detector/proton_current/value'][:]
        else:
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
