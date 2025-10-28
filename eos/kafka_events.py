"""
Collect AMOR detector events send via Kafka.
"""

import logging
import numpy as np
from threading import Thread, Event
from time import time

from .event_data_types import AmorGeometry, AmorTiming, AmorEventStream, PACKET_TYPE, EVENT_TYPE, PULSE_TYPE, PC_TYPE

from uuid import uuid4
from streaming_data_types.eventdata_ev44 import EventData
from streaming_data_types.logdata_f144 import ExtractedLogData
from streaming_data_types import deserialise_f144, deserialise_ev44
from confluent_kafka import Consumer

from .header import Header


try:
    from streaming_data_types.utils import get_schema
except ImportError:
    from streaming_data_types.utils import _get_schema as get_schema


KAFKA_BROKER = 'linkafka01.psi.ch:9092'
AMOR_EVENTS = 'amor_ev44'
AMOR_NICOS = 'AMOR_nicosForwarder'

class KafkaFrozenData:
    """
    Represents event stream data from Kafka at a given time.
    Will be returned by KafkaEventData to be use in conjunction
    with data processing and projections.

    Implements EventDatasetProtocol
    """
    geometry: AmorGeometry
    timing: AmorTiming
    data: AmorEventStream

    def __init__(self, geometry, timing, data, monitor=1.):
        self.geometry = geometry
        self.timing = timing
        self.data = data
        self.monitor = monitor

    def append(self, other):
        raise NotImplementedError("can't append live datastream to other event data")

    def update_header(self, header:Header):
        # maybe makes sense later, but for now just used for live vizualization
        ...

class KafkaEventData(Thread):
    """
    Read Nicos information and events from Kafka. Creates a background
    thread that listens to Kafka events and converts them to eos compatible information.
    """
    geometry: AmorGeometry
    timing: AmorTiming
    events: np.recarray

    def __init__(self):
        self.stop_event = Event()
        self.stop_counting = Event()
        self.new_events = Event()
        self.last_read = 0
        self.last_read_time = 0.
        self.start_time = time()
        self.consumer = Consumer(
                {'bootstrap.servers': 'linkafka01.psi.ch:9092',
                 'group.id': uuid4()})
        self.consumer.subscribe([AMOR_EVENTS, AMOR_NICOS])
        self.geometry = AmorGeometry(1.0, 2.0, 0., 0., 1.5, 10.0, 4.0, 10.0)
        self.timing = AmorTiming(0., 0., 500., 0., 30./500.)
        # create empty dataset
        self.events = np.recarray(0, dtype=EVENT_TYPE)
        super().__init__()

    def run(self):
        while not self.stop_event.is_set():
            messages = self.consumer.consume(10, timeout=1)
            for message in messages:
                self.process_message(message)

    def process_message(self, message):
        if message.error():
            logging.info(f"  received Kafka message with error: {message.error()}")
            return
        schema = get_schema(message.value())
        if message.topic()==AMOR_EVENTS and schema=='ev44':
            events:EventData = deserialise_ev44(message.value())
            self.add_events(events)
            self.new_events.set()
            logging.debug(f'  new events {events}')
        elif message.topic()==AMOR_NICOS and schema=='f144':
            nicos_data:ExtractedLogData = deserialise_f144(message.value())
            if nicos_data.source_name in self.nicos_mapping.keys():
                logging.debug(f'  {nicos_data.source_name} = {nicos_data.value}')
                self.update_instrument(nicos_data)

    def add_events(self, events:EventData):
        """
        Add new events to the Dataset. The object keeps raw events
        and only copies the latest set to the self.data object,
        this allows to run the event processing to be performed on a "clean"
        evnet stream each time.
        """
        if self.stop_counting.is_set():
            return
        prev_size = self.events.shape[0]
        new_events = events.pixel_id.shape[0]
        self.events.resize(prev_size+new_events, refcheck=False)
        self.events.pixelID[prev_size:] = events.pixel_id
        self.events.mask[prev_size:] = 0
        self.events.tof[prev_size:] = events.time_of_flight/1.e9

    nicos_mapping = {
        'mu': ('geometry', 'mu'),
        'nu':  ('geometry', 'nu'),
        'kappa':  ('geometry', 'kap'),
        'kappa_offset':  ('geometry', 'kad'),
        'ch1_trigger_phase':  ('timing', 'ch1TriggerPhase'),
        'ch2_trigger_phase':  ('timing', 'ch2TriggerPhase'),
        'ch2_speed':  ('timing', 'chopperSpeed'),
        'chopper_phase':  ('timing', 'chopperPhase'),
        }

    def update_instrument(self, nicos_data:ExtractedLogData):
        if nicos_data.source_name in self.nicos_mapping:
            attr, subattr = self.nicos_mapping[nicos_data.source_name]
            setattr(getattr(self, attr), subattr, nicos_data.value)
            if nicos_data.source_name=='ch2_speed':
                self.timing.tau = 30./self.timing.chopperSpeed

    def monitor(self):
        return time()-self.start_time

    def restart(self):
        # empty event buffer
        self.events = np.recarray(0, dtype=EVENT_TYPE)
        self.stop_counting.clear()
        self.last_read = 0
        self.start_time = time()
        self.new_events.clear()

    def get_events(self, total_counts=False):
        packets = np.recarray(0, dtype=PACKET_TYPE)
        pulses = np.recarray(0, dtype=PULSE_TYPE)
        pc = np.recarray(0, dtype=PC_TYPE)
        if total_counts:
            last_read = 0
        else:
            last_read = self.last_read
        if last_read>=self.events.shape[0]:
            raise EOFError("No new events arrived")
        data = AmorEventStream(self.events[last_read:].copy(), packets, pulses, pc)
        self.last_read = self.events.shape[0]
        self.new_events.clear()
        t_now = time()
        monitor = t_now-self.last_read_time
        self.last_read_time = t_now
        return KafkaFrozenData(self.geometry, self.timing, data, monitor=monitor)
