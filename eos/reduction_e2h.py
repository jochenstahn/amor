"""
Events 2 histogram, quick reduction of single file to display during experiment.
Can be used as a live preview with automatic update when files are modified.
"""

import logging
import os
import matplotlib.pyplot as plt
import numpy as np

from orsopy import fileio
from datetime import datetime

from .file_reader import AmorEventData, AmorHeader
from .header import Header
from .instrument import LZGrid
from .normalization import LZNormalisation
from .options import E2HConfig, E2HPlotArguments, IncidentAngle, MonitorType, E2HPlotSelection
from . import event_handling as eh, event_analysis as ea
from .path_handling import PathResolver
from .projection import LZProjection, ProjectionInterface


class E2HReduction:
    config: E2HConfig
    header: Header
    event_actions: eh.EventDataAction

    _last_mtime = 0.

    def __init__(self, config: E2HConfig):
        self.config = config

        self.header = Header()

        self.prepare_actions()

    def prepare_actions(self):
        """
        Does not do any actual reduction.
        """
        self.path_resolver = PathResolver(self.config.reader.year, self.config.reader.rawPath)
        self.file_list = self.path_resolver.resolve(self.config.reduction.fileIdentifier)
        self.file_index = 0
        self.plot_kwds = {}
        self.fig = plt.figure()

        if self.config.reduction.update:
            # live update implies plotting
            self.config.reduction.show_plot = True

        # Actions on datasets not used for normalization
        self.event_actions = eh.ApplyPhaseOffset(self.config.experiment.chopperPhaseOffset)
        if not self.config.reduction.fast:
            self.event_actions |= eh.CorrectChopperPhase()
            self.event_actions |= ea.ExtractWalltime()
        else:
            logging.info('    Fast reduction always uses time normalization')
            self.config.experiment.monitorType = MonitorType.time
        self.event_actions |= eh.AssociatePulseWithMonitor(self.config.experiment.monitorType)
        if self.config.experiment.monitorType in [MonitorType.proton_charge, MonitorType.debug]:
            # the filtering only makes sense if using actual monitor data, not time
            self.event_actions |= eh.FilterMonitorThreshold(self.config.experiment.lowCurrentThreshold)
        if not self.config.reduction.fast:
            self.event_actions |= eh.FilterStrangeTimes()
        # select needed actions in depenence of plots
        if self.config.reduction.plot in [E2HPlotSelection.All, E2HPlotSelection.LT, E2HPlotSelection.Q,
                                          E2HPlotSelection.L]:
            self.event_actions |= ea.MergeFrames()
            self.event_actions |= ea.AnalyzePixelIDs(self.config.experiment.yRange)
            self.event_actions |= eh.TofTimeCorrection(self.config.experiment.incidentAngle==IncidentAngle.alphaF)
            self.event_actions |= ea.CalculateWavelength(self.config.experiment.lambdaRange)
        self.event_actions |= eh.ApplyMask()

        # plot dependant options
        if self.config.reduction.plot in [E2HPlotSelection.All, E2HPlotSelection.LT, E2HPlotSelection.Q]:
            self.grid = LZGrid(0.01, [0.0, 0.25])

        if self.config.reduction.plot==E2HPlotSelection.LT:
            self.plot_kwds['colorbar'] = True
            self.plot_kwds['cmap'] = str(self.config.reduction.plot_colormap)

    def reduce(self):
        self.norm = LZNormalisation.unity(self.grid)

        self.prepare_graphs()

        while self.file_index < len(self.file_list):
            self.read_data()
            self.add_data()

        if self.config.reduction.plotArgs==E2HPlotArguments.OutputFile:
            self.create_file_output()
        if self.config.reduction.plotArgs!=E2HPlotArguments.OutputFile or self.config.reduction.show_plot:
            self.create_graph()

        if self.config.reduction.plotArgs==E2HPlotArguments.Default:
            plt.savefig(f'e2h_{self.config.reduction.plot}.png', dpi=300)
        if self.config.reduction.update:
            self.timer = self.fig.canvas.new_timer(1000)
            self.timer.add_callback(self.update)
            self.timer.start()
        if self.config.reduction.show_plot:
            plt.show()

    def prepare_graphs(self):
        last_file_header = AmorHeader(self.file_list[-1])
        tthh  = last_file_header.geometry.nu - last_file_header.geometry.mu


        if self.config.reduction.plot==E2HPlotSelection.LT:
            self.projection = LZProjection(tthh, self.grid)
            if not self.config.reduction.fast:
                self.projection.correct_gravity(last_file_header.geometry.detectorDistance)

    def read_data(self):
        fileName = self.file_list[self.file_index]
        self.file_index += 1
        self.dataset = AmorEventData(fileName)
        self.event_actions(self.dataset)
        self.dataset.update_header(self.header)

        self.header.measurement_data_files.append(fileio.File(file=fileName.split('/')[-1],
                                                              timestamp=self.dataset.fileDate))

    def add_data(self):
        self.monitor = self.dataset.data.pulses.monitor.sum()
        self.projection.project(self.dataset, monitor=self.monitor)

    def create_file_output(self):
        ...

    def create_title(self):
        output = "Events to Histogram - "
        output += ",".join(["#"+os.path.basename(fi)[9:15].lstrip('0') for fi in self.file_list])
        output += f"\n$\\mu$={self.dataset.geometry.mu:.2f}"
        output += f" $\\nu$={self.dataset.geometry.nu:.2f}"
        if self.config.reduction.update:
            output += f"\n at "+datetime.now().strftime("%m/%d/%Y %I:%M:%S")
        return output

    def create_graph(self):
        plt.title(self.create_title())
        self.projection.plot(**self.plot_kwds)
        plt.tight_layout()
        if self.config.reduction.plot==E2HPlotSelection.LT:
            plt.connect('button_press_event', self.draw_qline)

    def draw_qline(self, event):
        if event.button is plt.MouseButton.LEFT and self.fig.canvas.manager.toolbar.mode=='':
            slope = event.ydata/event.xdata
            xmax = 12.5
            plt.plot([0, xmax], [0, slope*xmax], '-', color='grey')
            plt.text(event.xdata, event.ydata, f'q={np.deg2rad(slope)*4.*np.pi:.3f}', backgroundcolor='white')
            plt.draw()
        if event.button is plt.MouseButton.RIGHT and self.fig.canvas.manager.toolbar.mode=='':
            for art in list(plt.gca().lines)+list(plt.gca().texts):
                art.remove()
            plt.draw()

    def replace_dataset(self, latest):
        self.file_list = self.path_resolver.resolve(f'{latest}')
        self.file_index = 0
        self.read_data()
        self.projection.clear()
        self.add_data()
        self.fig.clear()
        self.create_graph()

    def update(self):
        if self.config.reduction.fileIdentifier=='0':
            # if latest file was choosen, check if new one available and switch to it
            current = int(os.path.basename(self.file_list[-1])[9:15])
            latest = self.path_resolver.search_latest(0)
            if latest>current:
                self.replace_dataset(latest)
                return
        if os.path.getmtime(self.file_list[-1])<=self._last_mtime:
            return

        self._last_mtime = os.path.getmtime(self.file_list[-1])
        try:
            update_data = AmorEventData(self.file_list[-1], self.dataset.last_index+1)
        except EOFError:
            return


        self.event_actions(update_data)
        self.dataset=update_data
        self.monitor = self.dataset.data.pulses.monitor.sum()
        self.projection.project(update_data, self.monitor)

        self.projection.update_plot()
        plt.title(self.create_title())
        plt.draw()