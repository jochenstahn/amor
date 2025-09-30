"""
Calculations performed on AmorEventData.
"""
import logging
import numpy as np

from abc import ABC, abstractmethod

from .header import Header
from .options import MonitorType
from .file_reader import AmorEventData

class EventDataAction(ABC):
    """
    Abstract base class used for actions applied to an AmorEventData object.
    Each action can optionally modify the header information.
    """

    @abstractmethod
    def __call__(self, dataset: AmorEventData)->None: ...

    def update_header(self, header:Header)->None:
        if hasattr(self, 'action_name'):
            header.reduction.corrections.append(getattr(self, 'action_name'))

class CorrectSeriesTime(EventDataAction):
    def __init__(self, seriesStartTime):
        self.seriesStartTime = np.int64(seriesStartTime)

    def __call__(self, dataset: AmorEventData)->None:
        dataset.data.pulses.time -= self.seriesStartTime
        dataset.data.events.wallTime -= self.seriesStartTime
        dataset.data.proton_current.time -= self.seriesStartTime
        start, stop = dataset.data.proton_current.time[0], dataset.data.proton_current.time[-1]
        logging.debug(f'      wall time from {start:6.1f} s to {stop/1e9:6.1f} s, '
                      f'series time = {self.seriesStartTime/1e9:6.1f}')

class AssociatePulseWithMonitor(EventDataAction):
    def __init__(self, monitorType:MonitorType, lowCurrentThreshold:float):
        self.monitorType = monitorType
        self.lowCurrentThreshold = lowCurrentThreshold

    def __call__(self, dataset: AmorEventData)->None:
        logging.debug(f'      using monitor type {self.monitorType}')
        if self.monitorType in [MonitorType.proton_charge or MonitorType.debug]:
            monitorPerPulse = self.get_current_per_pulse(dataset.data.pulses.time,
                                                              dataset.data.proton_current.time,
                                                              dataset.data.proton_current.current)\
                                                              * 2*dataset.timing.tau * 1e-3
            # filter low-current pulses
            dataset.data.pulses.monitor = np.where(
                    monitorPerPulse > 2*dataset.timing.tau * self.lowCurrentThreshold * 1e-3,
                    monitorPerPulse, 0)
        elif self.monitorType==MonitorType.time:
            dataset.data.pulses.monitor  = 2*dataset.timing.tau
        else:  # pulses
            dataset.data.pulses.monitor  = 1

    @staticmethod
    def get_current_per_pulse(pulseTimeS, currentTimeS, currents):
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

class FilterStrangeTimes(EventDataAction):
    def __call__(self, dataset: AmorEventData)->None:
        filter_e = (dataset.data.events.tof<=2*dataset.timing.tau)
        dataset.data.events = dataset.data.events[filter_e]
        if filter_e.any():
            logging.warning(f'        strange times: {filter_e.sum()}')
