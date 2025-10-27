"""
Events 2 histogram, quick reduction of single file to display during experiment.
Can be used as a live preview with automatic update when files are modified.
"""

import logging
import os

from time import sleep
from .kafka_events import KafkaEventData
from .header import Header
from .options import E2HConfig
from . import event_handling as eh, event_analysis as ea
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
        self.event_data = KafkaEventData()
        self.event_data.start()

        self.prepare_actions()

    def prepare_actions(self):
        """
        Does not do any actual reduction.
        """
        # Actions on datasets not used for normalization
        self.event_actions = eh.ApplyPhaseOffset(self.config.experiment.chopperPhaseOffset)
        self.event_actions |= eh.CorrectChopperPhase()
        self.event_actions |= ea.MergeFrames()
        self.event_actions |= eh.ApplyMask()

    def reduce(self):
        self.create_projections()
        self.read_data()
        self.add_data()

        self.serializer = ESSSerializer()
        self.serializer.start_command_thread()

        self.loop()

    def create_projections(self):
        self.proj_yz = YZProjection()
        self.proj_tofz = TofZProjection(self.event_data.timing.tau, foldback=True, combine=2)

    def read_data(self):
        # make sure the first events have arrived before starting analysis
        self.event_data.new_events.wait()
        self.dataset = self.event_data.get_events()
        self.event_actions(self.dataset)


    def add_data(self):
        self.monitor = self.dataset.monitor
        self.proj_yz.project(self.dataset, monitor=self.monitor)
        self.proj_tofz.project(self.dataset, monitor=self.monitor)

    def loop(self):
        self.wait_for = self.serializer.new_count_started
        while True:
            try:
                self.update()
                self.wait_for.wait(1.0)
            except KeyboardInterrupt:
                self.event_data.stop_event.set()
                self.event_data.join()
                self.serializer.end_command_thread()
                return

    def update(self):
        if self.serializer.new_count_started.is_set():
            logging.warning('Start new count, clearing event data')
            self.wait_for = self.serializer.count_stopped
            self.event_data.restart()
            self.serializer.new_count_started.clear()
            self.create_projections()
            return
        elif self.serializer.count_stopped.is_set() and not self.event_data.stop_counting.is_set():
            return self.finish_count()
        try:
            update_data = self.event_data.get_events()
        except EOFError:
            return
        logging.info("    updating with new data")

        self.event_actions(update_data)
        self.dataset=update_data
        self.monitor = self.dataset.monitor
        self.proj_yz.project(update_data, self.monitor)
        self.proj_tofz.project(update_data, self.monitor)

        self.serializer.send(self.proj_yz)
        self.serializer.send(self.proj_tofz)

    def finish_count(self):
        logging.debug("    stop event set, hold event collection and send final results")
        self.wait_for = self.serializer.new_count_started
        self.event_data.stop_counting.set()

        try:
            update_data = self.event_data.get_events()
        except EOFError:
            pass
        else:
            self.event_actions(update_data)
            self.dataset = update_data
            self.monitor = self.dataset.monitor
            self.proj_yz.project(update_data, self.monitor)
            self.proj_tofz.project(update_data, self.monitor)

        logging.warning(f'  stop counting, total events {int(self.proj_tofz.data.cts.sum())}')

        self.serializer.send(self.proj_yz, final=True)
        self.serializer.send(self.proj_tofz, final=True)
