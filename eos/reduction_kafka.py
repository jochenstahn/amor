"""
Events 2 histogram, quick reduction of single file to display during experiment.
Can be used as a live preview with automatic update when files are modified.
"""

import logging
import os

from time import sleep
from .file_reader import AmorEventData, AmorHeader
from .header import Header
from .options import E2HConfig
from . import event_handling as eh, event_analysis as ea
from .path_handling import PathResolver
from .projection import  TofZProjection,  YZProjection
from .kafka_serializer import ESSSerializer


class KafkaReduction:
    config: E2HConfig
    header: Header
    event_actions: eh.EventDataAction

    _last_mtime = 0.
    proj_yz: YZProjection
    proj_tofz = TofZProjection

    def __init__(self, config: E2HConfig):
        self.config = config

        self.header = Header()

        self.prepare_actions()

    def prepare_actions(self):
        """
        Does not do any actual reduction.
        """
        self.path_resolver = PathResolver(self.config.reader.year, self.config.reader.rawPath)
        self.current_file = self.path_resolver.resolve('0')[0]

        # Actions on datasets not used for normalization
        self.event_actions = eh.ApplyPhaseOffset(self.config.experiment.chopperPhaseOffset)
        self.event_actions |= eh.CorrectChopperPhase()
        self.event_actions |= ea.MergeFrames()
        self.event_actions |= eh.ApplyMask()

    def reduce(self):
        last_file_header = AmorHeader(self.current_file)

        self.proj_yz = YZProjection()
        self.proj_tofz = TofZProjection(last_file_header.timing.tau, foldback=True)

        self.read_data()
        self.add_data()

        self.serializer = ESSSerializer()
        self.serializer.start_command_thread()

        self.loop()


    def read_data(self):
        self.dataset = AmorEventData(self.current_file, max_events=self.config.reduction.max_events)
        self.event_actions(self.dataset)

    def add_data(self):
        self.monitor = self.dataset.data.pulses.monitor.sum()
        self.proj_yz.project(self.dataset, monitor=self.monitor)
        self.proj_tofz.project(self.dataset, monitor=self.monitor)

    def replace_dataset(self, latest):
        new_file = self.path_resolver.resolve('0')[0]
        if not os.path.exists(new_file):
            return
        try:
            # check that events exist in the new file
            AmorEventData(new_file, 0, max_events=1000)
        except Exception:
            logging.debug("Problem when trying to load new dataset", exc_info=True)
            return

        logging.warning(f"Preceding to next file {latest}")
        self.current_file = new_file
        self.read_data()
        self.add_data()

    def loop(self):
        while True:
            try:
                self.update()
                sleep(1.0)
            except KeyboardInterrupt:
                self.serializer.end_command_thread()
                return

    def update(self):
        logging.debug("    check for update")
        if self.config.reduction.fileIdentifier=='0':
            # if latest file was choosen, check if new one available and switch to it
            current = int(os.path.basename(self.current_file)[9:15])
            latest = self.path_resolver.search_latest(0)
            if latest>current:
                self.replace_dataset(latest)
                return
        # if all events were read last time, only load more if file was modified
        if self.dataset.EOF and os.path.getmtime(self.current_file)<=self._last_mtime:
            return

        self._last_mtime = os.path.getmtime(self.current_file)
        try:
            update_data = AmorEventData(self.current_file, self.dataset.last_index+1,
                                        max_events=self.config.reduction.max_events)
        except EOFError:
            return
        logging.info("    updating with new data")

        self.event_actions(update_data)
        self.dataset=update_data
        self.monitor = self.dataset.data.pulses.monitor.sum()
        self.proj_yz.project(update_data, self.monitor)
        self.proj_tofz.project(update_data, self.monitor)

        self.serializer.send(self.proj_yz)
        self.serializer.send(self.proj_tofz)
