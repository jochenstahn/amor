"""
Specify the data type and protocol used for event datasets.
"""
from typing import List, Optional, Protocol, Tuple
from dataclasses import dataclass
from .header import Header
from abc import ABC, abstractmethod
from hashlib import sha256
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

@dataclass
class AmorTiming:
    ch1TriggerPhase: float
    ch2TriggerPhase: float
    chopperSpeed: float
    chopperPhase: float
    tau: float

# Structured datatypes used for event streams
EVENT_TYPE = np.dtype([('tof', np.float64), ('pixelID', np.uint32), ('mask', np.int32)])
PACKET_TYPE = np.dtype([('start_index', np.uint32), ('time', np.int64)])
PULSE_TYPE = np.dtype([('time', np.int64), ('monitor', np.float32)])
PC_TYPE = np.dtype([('current', np.float32), ('time', np.int64)])

# define the bitmask for individual event filters
EVENT_BITMASKS = {
    'MonitorThreshold': 1,
    'StrangeTimes': 2,
    'yRange': 4,
    'LamdaRange': 8,
    'qRange': 16,
    }

def append_fields(input: np.recarray, new_fields: List[Tuple[str, np.dtype]]):
    # add one ore more fields to a recarray, numpy functions seems to fail
    flds = [(name, dtypei[0]) for name, dtypei in input.dtype.fields.items()]
    flds += new_fields
    output = np.recarray(len(input), dtype=flds)

    for field in input.dtype.fields.keys():
        output[field] = input[field]
    return output

@dataclass
class AmorEventStream:
    events: np.recarray # EVENT_TYPE
    packets: np.recarray # PACKET_TYPE
    pulses: np.recarray  # PULSE_TYPE
    proton_current: np.recarray # PC_TYPE

class EventDatasetProtocol(Protocol):
    """
    Minimal attributes a dataset needs to provide to work with EventDataAction
    """
    geometry: AmorGeometry
    timing: AmorTiming
    data: AmorEventStream

    def append(self, other):
        # Should define a way to add events from other to own
        ...

    def update_header(self, header:Header):
        # update a header with the information read from file
        ...

class EventDataAction(ABC):
    """
    Abstract base class used for actions applied to an EventDatasetProtocol based objects.
    Each action can optionally modify the header information.
    Actions can be combined using the pipe operator | (OR).
    """

    def __call__(self, dataset: EventDatasetProtocol)->None:
        logging.debug(f"        Enter action {self.__class__.__name__} on {dataset!r}")
        self.perform_action(dataset)

    @abstractmethod
    def perform_action(self, dataset: EventDatasetProtocol)->None: ...

    def update_header(self, header:Header)->None:
        if hasattr(self, 'action_name'):
            header.reduction.corrections.append(getattr(self, 'action_name'))

    def __or__(self, other:'EventDataAction')->'CombinedAction':
        return CombinedAction([self, other])

    def __repr__(self):
        output = self.__class__.__name__+'('
        for key,value in self.__dict__.items():
            output += f'{key}={value}, '
        return output.rstrip(', ')+')'

    def action_hash(self)->bytes:
        # generate a unique hash that encodes this action with its configuration parameters
        mh = sha256()
        mh.update(self.__class__.__name__.encode())
        for key,value in sorted(self.__dict__.items()):
            mh.update(repr(value).encode())
        return mh.hexdigest()

class CombinedAction(EventDataAction):
    """
    Used to perform multiple actions in one call. Stores a sequence of actions
    that are then performed individually one after the other.
    """
    def __init__(self, actions: List[EventDataAction]) -> None:
        self._actions = actions

    def perform_action(self, dataset: EventDatasetProtocol)->None:
        for action in self._actions:
            action(dataset)

    def update_header(self, header:Header)->None:
        for action in self._actions:
            action.update_header(header)

    def __or__(self, other:'EventDataAction')->'CombinedAction':
        return CombinedAction(self._actions+[other])

    def __repr__(self):
        output = repr(self._actions[0])
        for ai in self._actions[1:]:
            output += ' | '+repr(ai)
        return output

    def action_hash(self)->bytes:
        mh = sha256()
        for action in self._actions:
            mh.update(action.action_hash().encode())
        return mh.hexdigest()
