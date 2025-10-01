"""
Specify the data type and protocol used for event datasets.
"""
from typing import Optional, Protocol
from dataclasses import dataclass
from .header import Header
from abc import ABC, abstractmethod
import  numpy as np
import logging

@dataclass
class AmorGeometry:
    mu:float
    nu:float
    kap:float
    kad:float
    div:float

    chopperSeparation: float
    detectorDistance: float
    chopperDetectorDistance: float

    delta_z: Optional[float] = None

@dataclass
class AmorTiming:
    ch1TriggerPhase: float
    ch2TriggerPhase: float
    chopperSpeed: float
    chopperPhase: float
    tau: float

# Structured datatypes used for event streams
EVENT_TYPE = np.dtype([('tof', np.float64),('pixelID', np.uint32), ('wallTime', np.int64)])
PACKET_TYPE = np.dtype([('start_index', np.uint32), ('Time', np.int64)])
PULSE_TYPE = np.dtype([('time', np.int64), ('monitor', np.float32)])
PC_TYPE = np.dtype([('current', np.float32), ('time', np.int64)])
# analyzed event sreams with extra attributes
ANA_EVENT_TYPE = np.dtype([('tof', np.float64),('pixelID', np.uint32), ('wallTime', np.int64),
                           ('detZ', np.float64), ('detXdist', np.float64), ('delta', np.float64),
                           ('mask', bool)])

FINAL_EVENT_TYPE = np.dtype([('tof', np.float64),('pixelID', np.uint32), ('wallTime', np.int64),
                           ('detZ', np.float64), ('detXdist', np.float64), ('delta', np.float64),
                           ('mask', bool),
                           ('lamda', np.float64), ('qz', np.float64), ('qx', np.float64), ])

@dataclass
class AmorEventStream:
    events: np.recarray # EVENT_TYPE
    packets: np.recarray # PACKET_TYPE
    pulses: Optional[np.recarray] = None # PULSE_TYPE
    proton_current: Optional[np.recarray] = None # PC_TYPE

class EventDatasetProtocol(Protocol):
    """
    Minimal attributes a dataset needs to provide to work with EventDataAction
    """
    geometry: AmorGeometry
    timing: AmorTiming
    data: AmorEventStream

class EventDataAction(ABC):
    """
    Abstract base class used for actions applied to an EventDatasetProtocol based objects.
    Each action can optionally modify the header information.
    """

    def __call__(self, dataset: EventDatasetProtocol)->None:
        logging.debug(f"        Enter action {self.__class__.__name__} on {dataset!r}")
        self.perform_action(dataset)

    @abstractmethod
    def perform_action(self, dataset: EventDatasetProtocol)->None: ...

    def update_header(self, header:Header)->None:
        if hasattr(self, 'action_name'):
            header.reduction.corrections.append(getattr(self, 'action_name'))

