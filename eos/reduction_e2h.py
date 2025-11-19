"""
Events 2 histogram, quick reduction of single file to display during experiment.
Can be used as a live preview with automatic update when files are modified.
"""

import logging
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from orsopy import fileio
from datetime import datetime

from .file_reader import AmorEventData, AmorHeader
from .header import Header
from .instrument import LZGrid
from .normalization import LZNormalisation
from .options import E2HConfig, E2HPlotArguments, IncidentAngle, MonitorType, E2HPlotSelection
from . import event_handling as eh
from .path_handling import PathResolver
from .projection import CombinedProjection, LProjection, LZProjection, ProjectionInterface, ReflectivityProjector, \
    TofProjection,  TofZProjection,  TProjection, YTProjection, YZProjection

NEEDS_LAMDA = (E2HPlotSelection.All, E2HPlotSelection.LT, E2HPlotSelection.Q, E2HPlotSelection.L)

class E2HReduction:
    config: E2HConfig
    header: Header
    event_actions: eh.EventDataAction

    _last_mtime = 0.
    projection: ProjectionInterface

    def __init__(self, config: E2HConfig):
        self.config = config

        self.header = Header()

        self.fig = plt.figure()
        self.register_colormap()
        self.prepare_actions()

    def prepare_actions(self):
        """
        Does not do any actual reduction.
        """
        self.path_resolver = PathResolver(self.config.reader.year, self.config.reader.rawPath)
        self.file_list = self.path_resolver.resolve(self.config.reduction.fileIdentifier)
        self.file_index = 0
        self.plot_kwds = {}
        plt.rcParams.update({'font.size': self.config.reduction.fontsize})

        self.overwrite = eh.ApplyParameterOverwrites(self.config.experiment) # some actions use instrument parameters, change before that
        if self.config.reduction.update:
            # live update implies plotting
            self.config.reduction.show_plot = True

        if self.config.reduction.plot==E2HPlotSelection.Raw:
            # Raw implies fast caculations
            self.config.reduction.fast = True
        if not self.config.experiment.is_default('lambdaRange'):
            # filtering wavelength requires frame analysis
            self.config.reduction.fast = False

        if not self.config.reduction.fast or self.config.reduction.plot in NEEDS_LAMDA:
            from . import event_analysis as ea

        # Actions on datasets not used for normalization
        self.event_actions = eh.ApplyPhaseOffset(self.config.experiment.chopperPhaseOffset)
        if not self.config.reduction.fast:
            self.event_actions |= self.overwrite
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
            if self.config.reduction.plot in [E2HPlotSelection.YT, E2HPlotSelection.YZ]:
                # perform time fold-back and apply yRange filter if not fast mode
                self.event_actions |= ea.MergeFrames()
                self.event_actions |= ea.AnalyzePixelIDs(self.config.experiment.yRange)
            if self.config.reduction.plot==E2HPlotSelection.YT:
                # perform corrections for tof if not fast mode
                self.event_actions |= eh.TofTimeCorrection(self.config.experiment.incidentAngle==IncidentAngle.alphaF)
        # select needed actions in depenence of plots
        if self.config.reduction.plot in NEEDS_LAMDA or not self.config.experiment.is_default('lambdaRange'):
            self.event_actions |= ea.MergeFrames(lamdaCut=self.config.experiment.lambdaRange[0])
            self.event_actions |= ea.AnalyzePixelIDs(self.config.experiment.yRange)
            self.event_actions |= eh.TofTimeCorrection(self.config.experiment.incidentAngle==IncidentAngle.alphaF)
            self.event_actions |= ea.CalculateWavelength(self.config.experiment.lambdaRange)
        self.event_actions |= eh.ApplyMask()

        # plot dependant options
        if self.config.reduction.plot in [E2HPlotSelection.All, E2HPlotSelection.LT, E2HPlotSelection.Q]:
            self.grid = LZGrid(0.05, [0.0, 0.25], lambda_overwrite=self.config.experiment.lambdaRange)
            self.grid.dldl = 0.01

        if self.config.reduction.plot in [E2HPlotSelection.All, E2HPlotSelection.Raw,
                                          E2HPlotSelection.LT, E2HPlotSelection.YT,
                                          E2HPlotSelection.YZ, E2HPlotSelection.TZ]:
            self.plot_kwds['colorbar'] = True
            self.plot_kwds['cmap'] = str(self.config.reduction.plot_colormap)
            if self.config.reduction.plotArgs==E2HPlotArguments.Linear:
                self.plot_kwds['norm'] = None

    def reduce(self):
        if self.config.reduction.plot in [E2HPlotSelection.All, E2HPlotSelection.LT, E2HPlotSelection.Q]:
            if self.config.reduction.normalizationModel:
                self.norm = LZNormalisation.model(self.grid)
            else:
                self.norm = LZNormalisation.unity(self.grid)

        self.prepare_graphs()

        while self.file_index < len(self.file_list):
            self.read_data()
            self.add_data()

        if self.config.reduction.plotArgs==E2HPlotArguments.OutputFile:
            self.create_file_output()
        if self.config.reduction.plotArgs!=E2HPlotArguments.OutputFile or self.config.reduction.show_plot:
            self.create_graph()

        if self.config.reduction.plotArgs==E2HPlotArguments.Default and not self.config.reduction.update:
            # safe to image file if not auto-updating graph
            plt.savefig(f'e2h_{self.config.reduction.plot}.png', dpi=300)
        if self.config.reduction.kafka:
            from .kafka_serializer import ESSSerializer
            self.serializer = ESSSerializer()
            self.fig.canvas.mpl_connect('close_event', self.serializer.end_command_thread)
            self.serializer.start_command_thread()
            self.serializer.send(self.projection)
        if self.config.reduction.update:
            self.timer = self.fig.canvas.new_timer(1000)
            self.timer.add_callback(self.update)
            self.timer.start()
        if self.config.reduction.show_plot:
            plt.show()


    def register_colormap(self):
        cmap = plt.colormaps['turbo'](np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = ListedColormap(cmap, name='jochen_deluxe', N=cmap.shape[0])
        #cmap.set_bad((1.,1.,0.9))
        plt.colormaps.register(cmap)

    def prepare_graphs(self):
        last_file_header = AmorHeader(self.file_list[-1])
        self.overwrite.perform_action(last_file_header)
        tthh  = last_file_header.geometry.nu - last_file_header.geometry.mu

        if not self.config.reduction.is_default('thetaRangeR'):
            # adjust range based on detector center
            thetaRange = [ti+tthh for ti in self.config.reduction.thetaRangeR]
        else:
            thetaRange = [tthh - last_file_header.geometry.div/2, tthh + last_file_header.geometry.div/2]

        if self.config.reduction.plot==E2HPlotSelection.LT:
            self.projection = LZProjection(tthh, self.grid)
            if not self.config.reduction.fast:
                self.projection.correct_gravity(last_file_header.geometry.detectorDistance)
            self.projection.apply_lamda_mask(self.config.experiment.lambdaRange)
            self.projection.apply_theta_mask(thetaRange)
            for thi in self.config.reduction.thetaFilters:
                self.projection.apply_theta_filter((thi[0]+tthh, thi[1]+tthh))
            self.projection.apply_norm_mask(self.norm)

        if self.config.reduction.plot==E2HPlotSelection.Q:
            plz = LZProjection(tthh, self.grid)
            if not self.config.reduction.fast:
                plz.correct_gravity(last_file_header.geometry.detectorDistance)
            plz.calculate_q()
            plz.apply_lamda_mask(self.config.experiment.lambdaRange)
            plz.apply_theta_mask(thetaRange)
            for thi in self.config.reduction.thetaFilters:
                self.projection.apply_theta_filter((thi[0]+tthh, thi[1]+tthh))
            plz.apply_norm_mask(self.norm)
            self.projection = ReflectivityProjector(plz, self.norm)

        if self.config.reduction.plot==E2HPlotSelection.YZ:
            self.projection = YZProjection()

        if self.config.reduction.plot==E2HPlotSelection.YT:
            self.projection = YTProjection(tthh)

        if self.config.reduction.plot==E2HPlotSelection.T:
            self.projection = TProjection(tthh)

        if self.config.reduction.plot==E2HPlotSelection.L:
            self.projection = LProjection()

        if self.config.reduction.plot==E2HPlotSelection.TZ:
            self.projection = TofZProjection(last_file_header.timing.tau, foldback=not self.config.reduction.fast)

        if self.config.reduction.plot==E2HPlotSelection.ToF:
            self.projection = TofProjection(last_file_header.timing.tau, foldback=not self.config.reduction.fast)

        if self.config.reduction.plot==E2HPlotSelection.All:
            plz = LZProjection(tthh, self.grid)
            if not self.config.reduction.fast:
                plz.correct_gravity(last_file_header.geometry.detectorDistance)
            plz.calculate_q()
            plz.apply_lamda_mask(self.config.experiment.lambdaRange)
            plz.apply_theta_mask(thetaRange)
            for thi in self.config.reduction.thetaFilters:
                plz.apply_theta_filter((thi[0]+tthh, thi[1]+tthh))
            plz.apply_norm_mask(self.norm)
            pr = ReflectivityProjector(plz, self.norm)
            pyz = YZProjection()
            self.projection = CombinedProjection(3, 2, [plz, pyz, pr],
                                                 [(0, 2, 0, 1), (0, 2, 1, 2), (2, 3, 0, 2)])

        if self.config.reduction.plot==E2HPlotSelection.Raw:
            del(self.plot_kwds['colorbar'])
            # A combined graph that does not require longer calculations
            plyt = YTProjection(tthh)
            pltofz = TofZProjection(last_file_header.timing.tau, foldback=not self.config.reduction.fast)
            pltof =  TofProjection(last_file_header.timing.tau, foldback=not self.config.reduction.fast)
            plt = TProjection(tthh)

            self.projection = CombinedProjection(3, 3, [plyt, pltofz, plt, pltof],
                                                 [(0,2, 0, 1), (0, 2, 1, 3), (2,3, 0,1),(2,3,1,3)])

    def read_data(self):
        fileName = self.file_list[self.file_index]
        self.dataset = AmorEventData(fileName, max_events=self.config.reduction.max_events)
        if self.dataset.EOF or fileName==self.file_list[-1]:
            self.file_index += 1
        self.event_actions(self.dataset)
        self.dataset.update_header(self.header)

        self.header.measurement_data_files.append(fileio.File(file=fileName.split('/')[-1],
                                                              timestamp=self.dataset.fileDate))

    def add_data(self):
        self.monitor = self.dataset.data.pulses.monitor.sum()
        self.projection.project(self.dataset, monitor=self.monitor)
        if self.config.reduction.plot==E2HPlotSelection.LT:
            self.projection.normalize_over_illuminated(self.norm)

    def create_file_output(self):
        raise NotImplementedError("Export to text output not yet implemented")

    def create_title(self):
        output = "Events to Histogram - "
        output += ",".join(["#"+os.path.basename(fi)[9:15].lstrip('0') for fi in self.file_list])
        output += f" ($\\mu$={self.dataset.geometry.mu:.2f} ;"
        output += f" $\\nu$={self.dataset.geometry.nu:.2f})"
        if self.config.reduction.update:
            output += f"\n at "+datetime.now().strftime("%m/%d/%Y %H:%M:%S")
        return output

    def create_graph(self):
        plt.suptitle(self.create_title())
        self.projection.plot(**self.plot_kwds)
        plt.tight_layout(pad=0.5)

    def replace_dataset(self, latest):
        new_files = self.path_resolver.resolve(f'{latest}')
        if not os.path.exists(new_files[-1]):
            return
        try:
            # check that events exist in the new file
            AmorEventData(new_files[-1], 0, max_events=1000)
        except Exception:
            logging.debug("Problem when trying to load new dataset", exc_info=True)
            return

        logging.warning(f"Preceding to next file {latest}")
        self.file_list = new_files
        self.file_index = 0
        self.prepare_actions()
        self.prepare_graphs()
        self.read_data()
        self.add_data()
        self.fig.clear()
        self.create_graph()
        plt.draw()

    def update(self):
        logging.debug("    check for update")
        if self.config.reduction.fileIdentifier=='0':
            # if latest file was choosen, check if new one available and switch to it
            current = int(os.path.basename(self.file_list[-1])[9:15])
            latest = self.path_resolver.search_latest(0)
            if latest>current:
                self.replace_dataset(latest)
                return
        # if all events were read last time, only load more if file was modified
        if self.dataset.EOF and os.path.getmtime(self.file_list[-1])<=self._last_mtime:
            return

        self._last_mtime = os.path.getmtime(self.file_list[-1])
        try:
            update_data = AmorEventData(self.file_list[-1], self.dataset.last_index+1,
                                        max_events=self.config.reduction.max_events)
        except EOFError:
            return
        logging.info("    updating with new data")

        self.event_actions(update_data)
        self.dataset=update_data
        self.monitor = self.dataset.data.pulses.monitor.sum()
        self.projection.project(update_data, self.monitor)

        self.projection.update_plot()
        plt.suptitle(self.create_title())
        plt.draw()

        if self.config.reduction.kafka:
            self.serializer.send(self.projection)
