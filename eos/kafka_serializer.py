"""
Allows to send eos projections to Kafka using ESS histogram serialization.

For histogram_h01 the message is build using:

hist = {
    "source": "some_source",
    "timestamp": 123456,
    "current_shape": [2, 5],
    "dim_metadata": [
        {
            "length": 2,
            "unit": "a",
            "label": "x",
            "bin_boundaries": np.array([10, 11, 12]),
        },
        {
            "length": 5,
            "unit": "b",
            "label": "y",
            "bin_boundaries": np.array([0, 1, 2, 3, 4, 5]),
        },
    ],
    "last_metadata_timestamp": 123456,
    "data": np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]]),
    "errors": np.array([[5, 4, 3, 2, 1], [10, 9, 8, 7, 6]]),
    "info": "info_string",
}
"""
import logging
from typing import List, Tuple, Union
from threading import Thread, Event

import numpy as np
import json
from time import time
from dataclasses import dataclass, asdict
from streaming_data_types import histogram_hs01
from confluent_kafka import Producer, Consumer

from uuid import uuid4

from .projection import TofZProjection, YZProjection

KAFKA_BROKER = 'linkafka01.psi.ch:9092'
KAFKA_TOPICS = {
    'histogram': 'AMOR_histograms',
    'response': 'AMOR_histResponse',
    'command':  'AMOR_histCommands'
    }

def ktime():
    return int(time()*1_000)

@dataclass
class DimMetadata:
    length: int
    unit: str
    label: str
    bin_boundaries: np.ndarray

@dataclass
class HistogramMessage:
    source: str
    timestamp: int
    current_shape: Tuple[int, int]
    dim_metadata: Tuple[DimMetadata, DimMetadata]
    last_metadata_timestamp: int
    data: np.ndarray
    errors: np.ndarray
    info: str

    def serialize(self):
        return histogram_hs01.serialise_hs01(asdict(self))

@dataclass
class CommandMessage:
    msg_id: str

    cmd=None

    @classmethod
    def get_message(cls, data):
        """
        Uses the sub-class cmd attribute to select which message to retugn
        """
        msg = dict([(ci.cmd, ci) for ci in cls.__subclasses__()])
        return msg[data['cmd']](**data)


@dataclass
class Stop(CommandMessage):
    hist_id: str
    id: str
    cmd:str = 'stop'

@dataclass
class HistogramConfig:
    id: str
    type: str
    data_brokers: List[str]
    topic: str
    data_topics: List[str]
    tof_range: Tuple[float, float]
    det_range: Tuple[int, int]
    num_bins: int
    width: int
    height: int
    left_edges: list
    source: str

@dataclass
class ConfigureHistogram(CommandMessage):
    histograms: List[HistogramConfig]
    start: int
    cmd:str = 'config'

    def __post_init__(self):
        self.histograms = [HistogramConfig(**cfg) for cfg in self.histograms]


class ESSSerializer:

    def __init__(self):
        self.producer = Producer({
            'bootstrap.servers': KAFKA_BROKER,
            'message.max.bytes': 4_000_000,
            })
        self.consumer = Consumer({
            'bootstrap.servers': KAFKA_BROKER,
            "group.id": uuid4(),
            "default.topic.config": {"auto.offset.reset": "latest"},
            })
        self._active_histogram_yz = None
        self._active_histogram_tofz = None
        self.new_count_started = Event()
        self.count_stopped = Event()

        self.consumer.subscribe([KAFKA_TOPICS['command']])

    def process_message(self, message):
        if message.error():
            logging.error("Command Consumer Error: %s", message.error())
        else:
            command = json.loads(message.value().decode())
            try:
                command = CommandMessage.get_message(command)
            except Exception:
                logging.error(f'Could not interpret message: \n{command}', exc_info=True)
                return
            logging.info(command)
            resp = json.dumps({
                "msg_id":   getattr(command, "id", None) or command.msg_id,
                "response": "ACK",
                "message":  ""
                })
            self.producer.produce(
                    topic=KAFKA_TOPICS['response'],
                    value=resp
                    )
            self.producer.flush()
            if isinstance(command, Stop):
                if command.hist_id == self._active_histogram_yz:
                    self.count_stopped.set()
                else:
                    return
            elif isinstance(command, ConfigureHistogram):
                for hist in command.histograms:
                    if hist.topic == KAFKA_TOPICS['histogram']+'_YZ':
                        self._active_histogram_yz = hist.id
                        logging.debug(f"   histogram data_topic: {hist.data_topics}")
                        self._start = command.start
                        self.count_stopped.clear()
                        self.new_count_started.set()
                    if hist.topic == KAFKA_TOPICS['histogram']+'_TofZ':
                        self._active_histogram_tofz = hist.id

    def receive(self, timeout=5):
        rec = self.consumer.poll(timeout)
        if rec is not None:
            self.process_message(rec)
            return True
        else:
            return False

    def receive_loop(self):
        while not self._stop_receiving.is_set():
            try:
                self.receive()
            except Exception:
                logging.error("Exception while receiving", exc_info=True)

    def start_command_thread(self):
        self._stop_receiving = Event()
        self._command_thread = Thread(target=self.receive_loop)
        self._command_thread.start()

    def end_command_thread(self, event=None):
        self._stop_receiving.set()
        self._command_thread.join()

    def acked(self, err, msg):
        # We need to have callback to produce-method to catch server errors
        if err is not None:
            logging.warning("Failed to deliver message: %s: %s" % (str(msg), str(err)))
        else:
            logging.debug("Message produced: %s" % (str(msg)))

    def send(self, proj: Union[YZProjection, TofZProjection], final=False):
        if final:
            state = 'FINISHED'
        else:
            state = 'COUNTING'
        if isinstance(proj, YZProjection):
            if self._active_histogram_yz is None:
                return
            suffix = 'YZ'
            message = HistogramMessage(
                source='amor-eos',
                timestamp=ktime(),
                current_shape=(proj.y.shape[0]-1, proj.z.shape[0]-1),
                dim_metadata=(
                    DimMetadata(
                            length=proj.y.shape[0]-1,
                            unit="pixel",
                            label="Y",
                            bin_boundaries=proj.y,
                            ),
                    DimMetadata(
                            length=proj.z.shape[0]-1,
                            unit="pixel",
                            label="Z",
                            bin_boundaries=proj.z,
                            )
                ),
                last_metadata_timestamp=0,
                data=proj.data.cts,
                errors=np.sqrt(proj.data.cts),
                info=json.dumps({
                    "start": self._start,
                    "state": state,
                    "num events": proj.data.cts.sum()
                })
                )
            logging.info(f"   {state}: Sending {proj.data.cts.sum()} events to Nicos")
        elif isinstance(proj, TofZProjection):
            if self._active_histogram_tofz is None:
                return
            suffix = 'TofZ'
            message = HistogramMessage(
                source='amor-eos',
                timestamp=ktime(),
                current_shape=(proj.tof.shape[0]-1, proj.z.shape[0]-1),
                dim_metadata=(
                    DimMetadata(
                            length=proj.tof.shape[0]-1,
                            unit="ms",
                            label="ToF",
                            bin_boundaries=proj.tof,
                            ),
                    DimMetadata(
                            length=proj.z.shape[0]-1,
                            unit="pixel",
                            label="Z",
                            bin_boundaries=proj.z,
                            ),
                ),
                last_metadata_timestamp=0,
                data=proj.data.cts,
                errors=np.sqrt(proj.data.cts),
                info=json.dumps({
                    "start": self._start,
                    "state": state,
                    "num events": proj.data.I.sum()
                })
                )
        else:
            raise NotImplementedError(f"Histogram for {proj.__class__.__name__} not implemented")

        self.producer.produce(value=message.serialize(),
                              topic=KAFKA_TOPICS['histogram']+'_'+suffix,
                              callback=self.acked)
        self.producer.flush()
