"""
Calculations performed on AmorEventData.
"""
import logging
import numpy as np

from .options import MonitorType
from .event_data_types import EventDatasetProtocol, EventDataAction
from .helpers import merge_frames


class CorrectSeriesTime(EventDataAction):
    def __init__(self, seriesStartTime):
        self.seriesStartTime = np.int64(seriesStartTime)

    def perform_action(self, dataset: EventDatasetProtocol)->None:
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

    def perform_action(self, dataset: EventDatasetProtocol)->None:
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

        if self.monitorType == MonitorType.debug:
            cpp, t_bins = np.histogram(dataset.data.events.wallTime, dataset.data.pulses.time)
            np.savetxt('tme.hst', np.vstack((dataset.data.pulses.time[:-1], cpp, dataset.data.pulses.monitor[:-1])).T)

        if self.monitorType in [MonitorType.proton_charge or MonitorType.debug]:
            goodTimeS = dataset.data.pulses.time[dataset.data.pulses.monitor!=0]
            filter_e = np.isin(dataset.data.events.wallTime, goodTimeS)
            dataset.data.events = dataset.data.events[filter_e]
            logging.info(f'      low-beam (<{self.lowCurrentThreshold} mC) rejected pulses: '
                         f'{dataset.data.pulses.monitor.shape[0]-goodTimeS.shape[0]} '
                         f'out of {dataset.data.pulses.monitor.shape[0]}')
            logging.info(f'          with {filter_e.shape[0]-dataset.data.events.shape[0]} events')
            if goodTimeS.shape[0]:
                logging.info(f'      average counts per pulse =  {dataset.data.events.shape[0] / goodTimeS.shape[0]:7.1f}')
            else:
                logging.info(f'      average counts per pulse = undefined')

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
    def perform_action(self, dataset: EventDatasetProtocol)->None:
        filter_e = (dataset.data.events.tof<=2*dataset.timing.tau)
        dataset.data.events = dataset.data.events[filter_e]
        if not filter_e.all():
            logging.warning(f'        strange times: {np.logical_not(filter_e).sum()}')

class MergeFrames(EventDataAction):
    def __init__(self, tofCut:float):
        self.tofCut = tofCut

    def perform_action(self, dataset: EventDatasetProtocol)->None:
        total_offset = (self.tofCut +
                        dataset.timing.tau * (dataset.timing.ch1TriggerPhase + dataset.timing.chopperPhase/2)/180)
        dataset.data.events.tof = merge_frames(dataset.data.events.tof, self.tofCut, dataset.timing.tau, total_offset)
