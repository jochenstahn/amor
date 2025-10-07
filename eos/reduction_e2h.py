"""
Events 2 histogram, quick reduction of single file to display during experiment.
Can be used as a live preview with automatic update when files are modified.
"""

import logging

from .file_reader import AmorEventData
from .header import Header
from .instrument import LZGrid
from .normalization import LZNormalisation
from .options import E2HConfig, IncidentAngle, MonitorType, NormalisationMethod, MONITOR_UNITS, ExperimentConfig
from . import event_handling as eh, event_analysis as ea

class E2HReduction:
    config: E2HConfig
    header: Header
    event_actions: eh.EventDataAction

    def __init__(self, config: E2HConfig):
        self.config = config

        self.header = Header()

        self.prepare_actions()

    def prepare_actions(self):
        """
        Does not do any actual reduction.
        """
        # Actions on datasets not used for normalization
        self.event_actions = eh.ApplyPhaseOffset(self.config.experiment.chopperPhaseOffset)
        self.event_actions |= eh.CorrectChopperPhase()
        self.event_actions |= ea.ExtractWalltime()
        self.event_actions |= eh.AssociatePulseWithMonitor(self.config.experiment.monitorType)
        if self.config.experiment.monitorType in [MonitorType.proton_charge or MonitorType.debug]:
            # the filtering only makes sense if using actual monitor data, not time
            self.event_actions |= eh.FilterMonitorThreshold(self.config.experiment.lowCurrentThreshold)
        self.event_actions |= eh.FilterStrangeTimes()
        self.event_actions |= ea.MergeFrames()
        self.event_actions |= ea.AnalyzePixelIDs(self.config.experiment.yRange)
        self.event_actions |= eh.TofTimeCorrection(self.config.experiment.incidentAngle==IncidentAngle.alphaF)
        self.event_actions |= ea.CalculateWavelength(self.config.experiment.lambdaRange)
        self.event_actions |= eh.ApplyMask()

        self.grid = LZGrid([0.01], [0.0, 0.25])

    def reduce(self):
        self.norm = LZNormalisation.unity(self.grid)

