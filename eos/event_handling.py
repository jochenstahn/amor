"""
Calculations performed on AmorEventData.
This module contains actions that do not need the numba base helper functions. Other actions are in event_analysis
"""
import logging
import os
import numpy as np

from .header import Header
from .options import ExperimentConfig, MonitorType
from .event_data_types import EventDatasetProtocol, EventDataAction, EVENT_BITMASKS

class ApplyPhaseOffset(EventDataAction):
    def __init__(self, chopperPhaseOffset: float):
        self.chopperPhaseOffset=chopperPhaseOffset

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        logging.debug(
            f'        replaced ch1TriggerPhase = {dataset.timing.ch1TriggerPhase} '
            f'with {self.chopperPhaseOffset}')
        dataset.timing.ch1TriggerPhase = self.chopperPhaseOffset

class ApplyParameterOverwrites(EventDataAction):
    def __init__(self, config: ExperimentConfig):
        self.config=config

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        if self.config.muOffset:
            logging.debug(f'        set muOffset = {self.config.muOffset}')
            dataset.geometry.mu += self.config.muOffset
        if self.config.mu:
            logging.debug(f'        replaced mu = {dataset.geometry.mu} with {self.config.mu}')
            dataset.geometry.mu = self.config.mu
        if self.config.nu:
            logging.debug(f'        replaced nu = {dataset.geometry.nu} with {self.config.nu}')
            dataset.geometry.nu = self.config.nu
        logging.info(f'      mu = {dataset.geometry.mu:6.3f}, '
                     f'nu = {dataset.geometry.nu:6.3f}, '
                     f'kap = {dataset.geometry.kap:6.3f}, '
                     f'kad = {dataset.geometry.kad:6.3f}')

    def update_header(self, header:Header) ->None:
        if self.config.sampleModel:
            import yaml
            from orsopy.fileio.model_language import SampleModel
            if self.config.sampleModel.endswith('.yml') or self.config.sampleModel.endswith('.yaml'):
                if os.path.isfile(self.config.sampleModel):
                    with open(self.config.sampleModel, 'r') as model_yml:
                        model = yaml.safe_load(model_yml)
                else:
                    logging.warning(f'  ! the file {self.config.sampleModel}.yml does not exist. Ignored!')
                    return
            else:
                model = dict(stack=self.config.sampleModel)
            logging.debug(f'        set sample.model = {self.config.sampleModel}')
            header.sample.model = SampleModel.from_dict(model)


class CorrectChopperPhase(EventDataAction):
    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        dataset.data.events.tof += dataset.timing.tau*(dataset.timing.ch1TriggerPhase-dataset.timing.chopperPhase/2)/180


class CorrectSeriesTime(EventDataAction):
    def __init__(self, seriesStartTime):
        self.seriesStartTime = np.int64(seriesStartTime)

    def perform_action(self, dataset: EventDatasetProtocol)->None:
        if not 'wallTime' in dataset.data.events.dtype.names:
            raise ValueError("CorrectTimeSeries requires walltTime to be extracted, please run ExtractWalltime first")
        dataset.data.pulses.time -= self.seriesStartTime
        dataset.data.events.wallTime -= self.seriesStartTime
        dataset.data.proton_current.time -= self.seriesStartTime
        start, stop = dataset.data.events.wallTime[0], dataset.data.events.wallTime[-1]
        logging.debug(f'      wall time from {start/1e9:6.1f} s to {stop/1e9:6.1f} s, '
                      f'series time = {self.seriesStartTime/1e9:6.1f}')


class AssociatePulseWithMonitor(EventDataAction):
    def __init__(self, monitorType:MonitorType):
        self.monitorType = monitorType

    def perform_action(self, dataset: EventDatasetProtocol)->None:
        logging.debug(f'      using monitor type {self.monitorType}')
        if self.monitorType in [MonitorType.proton_charge or MonitorType.debug]:
            dataset.data.pulses.monitor = self.get_current_per_pulse(dataset.data.pulses.time,
                                                              dataset.data.proton_current.time,
                                                              dataset.data.proton_current.current)\
                                                              * 2*dataset.timing.tau * 1e-3
        elif self.monitorType==MonitorType.time:
            dataset.data.pulses.monitor  = 2*dataset.timing.tau
        else:  # pulses
            dataset.data.pulses.monitor  = 1

        if self.monitorType == MonitorType.debug:
            if not 'wallTime' in dataset.data.events.dtype.names:
                raise ValueError(
                    "AssociatePulseWithMonitor requires walltTime for debugging, please run ExtractWalltime first")
            cpp, t_bins = np.histogram(dataset.data.events.wallTime, dataset.data.pulses.time)
            np.savetxt('tme.hst', np.vstack((dataset.data.pulses.time[:-1], cpp, dataset.data.pulses.monitor[:-1])).T)

    @staticmethod
    def get_current_per_pulse(pulseTimeS, currentTimeS, currents):
        # add currents for early pulses and current time value after last pulse (j+1)
        currentTimeS = np.hstack([[0], currentTimeS, [pulseTimeS[-1]+1]])
        currents = np.hstack([[0], currents])
        pulseCurrentS = np.zeros(pulseTimeS.shape[0], dtype=float)
        j = 0
        for i, ti in enumerate(pulseTimeS):
            # find the last current item that was before this pulse
            while ti >= currentTimeS[j+1]:
                j += 1
            pulseCurrentS[i] = currents[j]
        return pulseCurrentS

class FilterMonitorThreshold(EventDataAction):
    def __init__(self, lowCurrentThreshold:float):
        self.lowCurrentThreshold = lowCurrentThreshold

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        if not 'wallTime' in dataset.data.events.dtype.names:
            raise ValueError(
                    "FilterMonitorThreshold requires walltTime to be extracted, please run ExtractWalltime first")
        low_current_filter = dataset.data.pulses.monitor>2*dataset.timing.tau*self.lowCurrentThreshold*1e-3
        dataset.data.pulses.monitor[np.logical_not(low_current_filter)] = 0.
        goodTimeS = dataset.data.pulses.time[low_current_filter]
        filter_e = np.logical_not(np.isin(dataset.data.events.wallTime, goodTimeS))

        dataset.data.events.mask += EVENT_BITMASKS['MonitorThreshold']*filter_e
        logging.info(f'      low-beam (<{self.lowCurrentThreshold} mC) rejected pulses: '
                     f'{dataset.data.pulses.monitor.shape[0]-goodTimeS.shape[0]} '
                     f'out of {dataset.data.pulses.monitor.shape[0]}')
        logging.info(f'          with {filter_e.sum()} events')
        if goodTimeS.shape[0]:
            logging.info(f'      average counts per pulse =  {dataset.data.events.shape[0]/goodTimeS.shape[0]:7.1f}')
        else:
            logging.info(f'      average counts per pulse = undefined')

class FilterStrangeTimes(EventDataAction):
    def perform_action(self, dataset: EventDatasetProtocol)->None:
        filter_e = np.logical_not(dataset.data.events.tof<=2*dataset.timing.tau)
        dataset.data.events.mask += EVENT_BITMASKS['StrangeTimes']*filter_e
        if filter_e.any():
            logging.warning(f'        strange times: {filter_e.sum()}')


class TofTimeCorrection(EventDataAction):
    def __init__(self, correct_chopper_opening: bool = True):
        self.correct_chopper_opening = correct_chopper_opening

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        if not 'delta' in dataset.data.events.dtype.names:
            raise ValueError(
                    "TofTimeCorrection requires delta to be extracted, please run AnalyzePixelIDs first")
        d = dataset.data
        if self.correct_chopper_opening:
            d.events.tof -= ( d.events.delta / 180. ) * dataset.timing.tau
        else:
            d.events.tof -= ( dataset.geometry.kad / 180. ) * dataset.timing.tau

class ApplyMask(EventDataAction):
    def __init__(self, bitmask_filter=None):
        self.bitmask_filter = bitmask_filter

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        # TODO: why is this action time consuming?
        d = dataset.data
        pre_filter = d.events.shape[0]
        if logging.getLogger().level == logging.DEBUG:
            # only run this calculation if debug level is actually active
            filtered_by_mask = {}
            for key, value in EVENT_BITMASKS.items():
                filtered_by_mask[key] = ((d.events.mask & value)!=0).sum()
            logging.debug(f"        Removed by filters: {filtered_by_mask}")
        if self.bitmask_filter is None:
            d.events = d.events[d.events.mask==0]
        else:
            # remove the provided bitmask_filter bits from the events
            # this means that all bits that are set in bitmask_filter will NOT be used to filter events
            fltr = (d.events.mask & (~self.bitmask_filter)) == 0
            d.events = d.events[fltr]
        post_filter = d.events.shape[0]
        logging.info(f'      number of events: total = {pre_filter:7d}, filtered = {post_filter:7d}')
