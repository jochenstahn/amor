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
from typing import Tuple, Union

import numpy as np
import json
from time import time
from dataclasses import dataclass, asdict
from streaming_data_types import histogram_hs00
from confluent_kafka import Producer, Consumer, TopicPartition

from uuid import uuid4

from .projection import LZProjection, YZProjection

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
        return histogram_hs00.serialise_hs00(asdict(self))

class ESSSerializer:

    def __init__(self):
        self.producer = Producer({
            'bootstrap.servers': 'linkafka01.psi.ch:9092',
            })
        self.consumer = Consumer({
            'bootstrap.servers': 'linkafka01.psi.ch:9092',
            "group.id": uuid4(),
            "default.topic.config": {"auto.offset.reset": "latest"},
            })

        #tp = [TopicPartition( "SANSLLB_histCommands",0)]
        #self.consumer.assign(tp)
        self.consumer.subscribe(["SANSLLB_histCommands"])

    def process_message(self, message):
        if message.error():
            print("Command Consumer Error: %s", message.error())
        else:
            command = json.loads(message.value().decode())
            print(command)
            resp = json.dumps({
                "msg_id":   command.get("id") or command.get("msg_id"),
                "response": "ACK",
                "message":  ""
                })
            self.producer.produce(
                    topic="SANSLLB_histResponse",
                    value=resp
                    )
            self.producer.flush()

    def receive(self, timeout=5):
        rec = self.consumer.poll(5)
        if rec is not None:
            self.process_message(rec)
            return True
        else:
            return False


    def acked(self, err, msg):
        # We need to have callback to produce-method to catch server errors
        if err is not None:
            print("Failed to deliver message: %s: %s" % (str(msg), str(err)))
        else:
            print("Message produced: %s" % (str(msg)))

    def send(self, proj: Union[YZProjection, LZProjection]):
        if isinstance(proj, YZProjection):
            message = HistogramMessage(
                source='just-bin-it',
                timestamp=int(time()),
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
                data=proj.data.I,
                errors=proj.data.err,
                info=json.dumps({
                    "start": int(time()),
                    "state": 'COUNTING',
                    "num events": proj.data.cts.sum()
                })
                )
        elif isinstance(proj, LZProjection):
            message = HistogramMessage(
                source='just-bin-it',
                timestamp=int(time()),
                current_shape=(proj.lamda.shape[0]-1, proj.alphaF.shape[0]-1),
                dim_metadata=(
                    DimMetadata(
                            length=proj.lamda.shape[0]-1,
                            unit="Angstrom",
                            label="Lambda",
                            bin_boundaries=proj.lamda,
                            ),
                    DimMetadata(
                            length=proj.alphaF.shape[0]-1,
                            unit="degrees",
                            label="Theta",
                            bin_boundaries=proj.alphaF,
                            )
                ),
                last_metadata_timestamp=0,
                data=proj.data.ref,
                errors=proj.data.err,
                info=json.dumps({
                    "start": int(time()),
                    "state": 'COUNTING',
                    "num events": proj.data.I.sum()
                })
                )
        else:
            raise NotImplementedError(f"Histogram for {proj.__class__.__name__} not implemented")

        self.producer.produce(value=message.serialize(),
                              topic='SANSLLB_histograms',
                              callback=self.acked)
        self.producer.flush()
