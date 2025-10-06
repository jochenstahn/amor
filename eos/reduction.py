import logging
import os
import sys

import numpy as np
from orsopy import fileio

from .file_reader import AmorEventData
from .header import Header
from .path_handling import PathResolver
from .options import EOSConfig, IncidentAngle, MonitorType, NormalisationMethod
from .instrument import Detector, LZGrid
from .normalisation import LZNormalisation
from . import event_handling as eh, event_analysis as ea
from .projection import LZProjection


MONITOR_UNITS = {
    MonitorType.neutron_monitor: 'cnts',
    MonitorType.proton_charge: 'mC',
    MonitorType.time: 's',
    MonitorType.auto: 'various',
    MonitorType.debug: 'mC',
    }

class AmorReduction:
    def __init__(self, config: EOSConfig):
        self.experiment_config = config.experiment
        self.reader_config = config.reader
        self.reduction_config = config.reduction
        self.output_config = config.output

        self.header = Header()
        self.header.reduction.call = config.call_string()

        self.prepare_actions()

    def prepare_actions(self):
        """
        Does not do any actual reduction.
        """
        self.path_resolver = PathResolver(self.reader_config.year, self.reader_config.rawPath)

        # setup all actions performed on event datasets before projection on the grid
        # The order of these corrections matter as some rely on parameters modified before
        if self.reduction_config.normalisationFileIdentifier:
            # explicit steps performed on AmorEventDataset for normalization files
            self.normevent_actions = eh.ApplyPhaseOffset(self.experiment_config.chopperPhaseOffset)
            self.normevent_actions |= eh.CorrectChopperPhase()
            if self.experiment_config.monitorType in [MonitorType.proton_charge, MonitorType.debug]:
                self.normevent_actions |= ea.ExtractWalltime()
            self.normevent_actions |= eh.AssociatePulseWithMonitor(self.experiment_config.monitorType,
                                                                   self.experiment_config.lowCurrentThreshold)
            self.normevent_actions |= eh.FilterStrangeTimes()
            self.normevent_actions |= ea.MergeFrames()
            self.normevent_actions |= ea.AnalyzePixelIDs(self.experiment_config.yRange)
            self.normevent_actions |= eh.TofTimeCorrection(self.experiment_config.incidentAngle==IncidentAngle.alphaF)
            self.normevent_actions |= ea.CalculateWavelength(self.experiment_config.lambdaRange)
            self.normevent_actions |= eh.ApplyMask()
        # Actions on datasets not used for normalization
        self.dataevent_actions = eh.ApplyPhaseOffset(self.experiment_config.chopperPhaseOffset)
        self.dataevent_actions |= eh.ApplyParameterOverwrites(self.experiment_config) # some actions use instrument parameters, change before that
        self.dataevent_actions |= eh.CorrectChopperPhase()
        self.dataevent_actions |= ea.ExtractWalltime()
        self.dataevent_time_correction = eh.CorrectSeriesTime(0) # will be set from first dataset
        self.dataevent_actions |= self.dataevent_time_correction
        self.dataevent_actions |= eh.AssociatePulseWithMonitor(self.experiment_config.monitorType,
                                                               self.experiment_config.lowCurrentThreshold)
        self.dataevent_actions |= eh.FilterStrangeTimes()
        self.dataevent_actions |= ea.MergeFrames()
        self.dataevent_actions |= ea.AnalyzePixelIDs(self.experiment_config.yRange)
        self.dataevent_actions |= eh.TofTimeCorrection(self.experiment_config.incidentAngle==IncidentAngle.alphaF)
        self.dataevent_actions |= ea.CalculateWavelength(self.experiment_config.lambdaRange)
        self.dataevent_actions |= ea.CalculateQ(self.experiment_config.incidentAngle)
        self.dataevent_actions |= ea.FilterQzRange(self.reduction_config.qzRange)
        self.dataevent_actions |= eh.ApplyMask()

        self.grid = LZGrid(self.reduction_config.qResolution, self.reduction_config.qzRange)

    def reduce(self):
        if not os.path.exists(f'{self.output_config.outputPath}'):
            logging.debug(f'Creating destination path {self.output_config.outputPath}')
            os.system(f'mkdir {self.output_config.outputPath}')

        # load or create normalisation matrix
        if self.reduction_config.normalisationFileIdentifier:
            self.create_normalisation_map(self.reduction_config.normalisationFileIdentifier[0])
        else:
            self.norm = LZNormalisation.unity(self.grid)

        # load R(q_z) curve to be subtracted:
        if self.reduction_config.subtract:
            self.sq_q, self.sR_q, self.sdR_q, self.sFileName = self.loadRqz(self.reduction_config.subtract)
            logging.warning(f'loaded background file: {self.sFileName}')
            self.header.reduction.corrections.append(f'background from \'{self.sFileName}\' subtracted')
            self.subtract = True
        else:
            self.subtract = False

        # load measurement data and do the reduction
        self.datasetsRqz = []
        self.datasetsRlt = []
        for i, short_notation in enumerate(self.reduction_config.fileIdentifier):
            self.read_file_block(i, short_notation)

        # output
        logging.warning('output:')

        if 'Rqz.ort' in self.output_config.outputFormats:
            self.save_Rqz()

        if 'Rlt.ort' in self.output_config.outputFormats:
            self.save_Rtl()

    def read_file_block(self, i, short_notation):
        logging.warning('reading input:')
        file_list = self.path_resolver.resolve(short_notation)

        self.header.measurement_data_files = []

        self.dataset = AmorEventData(file_list[0])
        if self.experiment_config.monitorType==MonitorType.auto:
            if self.dataset.data.proton_current.current.sum()>1:
                self.experiment_config.monitorType = MonitorType.proton_charge
                logging.debug('      monitor type set to "proton current"')
            else:
                self.experiment_config.monitorType = MonitorType.time
                logging.debug('      monitor type set to "time"')
            # update actions to sue selected monitor
            self.prepare_actions()

        self.dataevent_time_correction.seriesStartTime = self.dataset.eventStartTime
        self.dataevent_actions(self.dataset)
        self.dataset.update_header(self.header)
        self.dataevent_actions.update_header(self.header)
        for fi in file_list[1:]:
            di = AmorEventData(fi)
            self.dataevent_actions(di)
            self.dataset.append(di)

        for fileName in file_list:
            self.header.measurement_data_files.append(fileio.File( file=fileName.split('/')[-1],
                                                                   timestamp=self.dataset.fileDate))


        if self.reduction_config.timeSlize:
            self.read_timeslices(i)
        else:
            self.read_unsliced(i)

    def read_unsliced(self, i):
        self.monitor = np.sum(self.dataset.data.pulses.monitor)
        logging.warning(f'    monitor = {self.monitor:8.2f} {MONITOR_UNITS[self.experiment_config.monitorType]}')

        proj:LZProjection = self.project_on_lz()
        try:
            scale = self.reduction_config.scale[i]
        except IndexError:
            scale = self.reduction_config.scale[-1]
        proj.scale(scale)

        if 'Rqz.ort' in self.output_config.outputFormats:
            headerRqz = self.header.orso_header()
            headerRqz.data_set = f'Nr {i} : mu = {self.dataset.geometry.mu:6.3f} deg'

            # projection on qz-grid
            result = proj.project_on_qz()

            if self.reduction_config.autoscale:
                if i==0:
                    result.autoscale(self.reduction_config.autoscale)
                else:
                    result.stitch(self.last_result)

            if self.subtract:
                if len(result.Q)==len(self.sq_q):
                    result.subtract(self.sR_q, self.sdR_q)
                else:
                    logging.warning(
                            f'backgroung file {self.sFileName} not compatible with q_z scale ({len(self.sq_q)} vs. {len(result.Q)})')

            orso_data = fileio.OrsoDataset(headerRqz, result.data)
            self.last_result = result
            self.datasetsRqz.append(orso_data)
        if 'Rlt.ort' in self.output_config.outputFormats:
            columns = [
                fileio.Column('Qz', '1/angstrom', 'normal momentum transfer'),
                fileio.Column('R', '', 'specular reflectivity'),
                fileio.ErrorColumn(error_of='R', error_type='uncertainty', value_is='sigma'),
                fileio.ErrorColumn(error_of='Qz', error_type='resolution', value_is='sigma'),
                fileio.Column('lambda', 'angstrom', 'wavelength'),
                fileio.Column('alpha_f', 'deg', 'final angle'),
                fileio.Column('l', '', 'index of lambda-bin'),
                fileio.Column('t', '', 'index of theta bin'),
                fileio.Column('intensity', '', 'filtered neutron events per pixel'),
                fileio.Column('norm', '', 'normalisation matrix'),
                fileio.Column('mask', '', 'pixels used for calculating R(q_z)'),
                fileio.Column('Qx', '1/angstrom', 'parallel momentum transfer'),
                ]

            ts, zs = proj.data.shape
            lindex_lz = np.tile(np.arange(1, ts+1), (zs, 1)).T
            tindex_lz = np.tile(np.arange(1, zs+1), (ts, 1))

            j = 0
            for item in zip(
                    proj.data.qz.T,
                    proj.data.ref.T,
                    proj.data.err.T,
                    proj.data.res.T,
                    proj.lamda.T,
                    proj.alphaF.T,
                    lindex_lz.T,
                    tindex_lz.T,
                    proj.data.I.T,
                    proj.data.norm.T,
                    proj.data.mask.T,
                    proj.data.qx.T,
                    ):
                data = np.array(list(item)).T
                headerRlt = self.header.orso_header(columns=columns)
                headerRlt.data_set = f'dataset_{i}_{j+1} : alpha_f = {proj.alphaF[0, j]:6.3f} deg'
                orso_data = fileio.OrsoDataset(headerRlt, data)
                self.datasetsRlt.append(orso_data)
                j += 1

    def read_timeslices(self, i):
        wallTime_e = np.float64(self.dataset.data.events.wallTime)/1e9
        pulseTimeS = np.float64(self.dataset.data.pulses.time)/1e9
        interval = self.reduction_config.timeSlize[0]
        try:
            start = self.reduction_config.timeSlize[1]
        except IndexError:
            start = 0
        try:
            stop = self.reduction_config.timeSlize[2]
        except IndexError:
            stop = wallTime_e[-1]
        # make overwriting log lines possible by removing newline at the end
        #logging.StreamHandler.terminator = "\r"
        logging.warning(f'    time slizing')
        logging.info('      slize  time  monitor')
        for ti, time in enumerate(np.arange(start, stop, interval)):
            slice = self.dataset.get_timeslice(time, time+interval)
            self.monitor = np.sum(slice.data.pulses.monitor)
            logging.info(f'      {ti:<4d}  {time:6.0f}  {self.monitor:7.2f} {MONITOR_UNITS[self.experiment_config.monitorType]}')

            proj: LZProjection = self.project_on_lz(slice)
            try:
                scale = self.reduction_config.scale[i]
            except IndexError:
                scale = self.reduction_config.scale[-1]
            proj.scale(scale)

            # projection on qz-grid
            result = proj.project_on_qz()

            if self.reduction_config.autoscale:
                if i==0:
                    result.autoscale(self.reduction_config.autoscale)
                else:
                    result.stitch(self.last_result)

            if self.subtract:
                if len(result.Q)==len(self.sq_q):
                    result.subtract(self.sR_q, self.sdR_q)
                else:
                    logging.warning(
                            f'backgroung file {self.sFileName} not compatible with q_z scale ({len(self.sq_q)} vs. {len(result.Q)})')

            headerRqz = self.header.orso_header(
                    extra_columns=[fileio.Column('time', 's', 'time relative to start of measurement series')])
            headerRqz.data_set = f'{i}_{ti}: time = {time:8.1f} s  to {time+interval:8.1f} s'
            orso_data = fileio.OrsoDataset(headerRqz, result.data_for_time(time))
            self.datasetsRqz.append(orso_data)
        self.last_result = result
        # reset normal logging behavior
        #logging.StreamHandler.terminator = "\n"
        logging.info(f'      done  {min(time+interval, pulseTimeS[-1]):5.0f}')

    def save_Rqz(self):
        fname = os.path.join(self.output_config.outputPath, f'{self.output_config.outputName}.Rqz.ort')
        logging.warning(f'    {fname}')
        theSecondLine = f' {self.header.experiment.title} | {self.header.experiment.start_date} | sample {self.header.sample.name} | R(q_z)'
        fileio.save_orso(self.datasetsRqz, fname, data_separator='\n', comment=theSecondLine)

    def save_Rtl(self):
        fname = os.path.join(self.output_config.outputPath, f'{self.output_config.outputName}.Rlt.ort')
        logging.warning(f'    {fname}')
        theSecondLine = f' {self.header.experiment.title} | {self.header.experiment.start_date} | sample {self.header.sample.name} | R(lambda, theta)'
        fileio.save_orso(self.datasetsRlt, fname, data_separator='\n', comment=theSecondLine)

    def loadRqz(self, name):
        fname = os.path.join(self.output_config.outputPath, name)
        if os.path.exists(fname):
            fileName = fname
        elif os.path.exists(f'{fname}.Rqz.ort'):
            fileName = f'{fname}.Rqz.ort'
        else:
            sys.exit(f'### the background file \'{fname}\' does not exist! => stopping')

        q_q, Sq_q, dS_q = np.loadtxt(fileName, usecols=(0, 1, 2), comments='#', unpack=True)

        return q_q, Sq_q, dS_q, fileName

    def create_normalisation_map(self, short_notation):
        outputPath = self.output_config.outputPath
        normalisation_list = self.path_resolver.expand_file_list(short_notation)
        name = '_'.join(map(str, normalisation_list))
        n_path = os.path.join(outputPath, f'{name}.norm')

        self.norm = None
        if os.path.exists(n_path):
            logging.debug(f'trying to load matrix from file {n_path}')
            try:
                self.norm = LZNormalisation.from_file(n_path, self.normevent_actions.action_hash())
            except (ValueError, EOFError):
                self.norm =None
            else:
                logging.warning(f'normalisation matrix: found and using {n_path}')
        if self.norm is None:
            # in case file does not exist or the action hash doesn't match, create new normalization
            logging.warning(f'normalisation matrix: using the files {normalisation_list}')
            normalization_files = list(map(self.path_resolver.get_path, normalisation_list))
            reference = AmorEventData(normalization_files[0])
            self.normevent_actions(reference)
            for nfi in normalization_files[1:]:
                toadd = AmorEventData(nfi)
                self.normevent_actions(toadd)
                reference.append(toadd)
            self.norm = LZNormalisation(reference, self.reduction_config.normalisationMethod, self.grid)
            if reference.data.events.shape[0] > 1e6:
                self.norm.safe(n_path, self.normevent_actions.action_hash())

        self.header.measurement_additional_files = self.norm.file_list
        self.header.reduction.corrections.append('normalisation with \'additional files\'')

    def project_on_lz(self, dataset=None):
        if dataset is None:
            dataset=self.dataset
        proj = LZProjection.from_dataset(dataset, self.grid,
                                         has_offspecular=(self.experiment_config.incidentAngle==IncidentAngle.alphaF))

        if self.reduction_config.is_default('thetaRangeR'):
            t0 = dataset.geometry.nu - dataset.geometry.mu
            # adjust range based on detector center
            thetaRange = [ti+t0 for ti in self.reduction_config.thetaRangeR]
            proj.apply_theta_mask(thetaRange)
        elif self.reduction_config.is_default('thetaRange'):
            proj.apply_theta_mask(self.reduction_config.thetaRange)
        else:
            # TODO: review what this is for an check if this should be removed or updated in case of thetaRangeR
            self.reduction_config.thetaRange = [dataset.geometry.nu - dataset.geometry.mu - dataset.geometry.div/2,
                                              dataset.geometry.nu - dataset.geometry.mu + dataset.geometry.div/2]
            proj.apply_theta_mask(self.reduction_config.thetaRange)

        if self.experiment_config.is_default('lambdaRange'):
            proj.apply_lamda_mask(self.experiment_config.lambdaRange)

        proj.apply_norm_mask(self.norm)
        proj.calculate_q()

        proj.project(dataset, self.monitor)

        if self.reduction_config.normalisationMethod == NormalisationMethod.over_illuminated:
            logging.debug('      assuming an overilluminated sample and correcting for the angle of incidence')
            proj.normalize_over_illuminated(self.norm)
        elif self.reduction_config.normalisationMethod==NormalisationMethod.under_illuminated:
            logging.debug('      assuming an underilluminated sample and ignoring the angle of incidence')
            proj.normalize_no_footprint(self.norm)
        elif self.reduction_config.normalisationMethod==NormalisationMethod.direct_beam:
            logging.debug('      assuming direct beam for normalisation and ignoring the angle of incidence')
            proj.normalize_no_footprint(self.norm)
        else:
            logging.error('unknown normalisation method! Use [u]nder, [o]ver or [d]irect illumination')
            proj.normalize_no_footprint(self.norm)
        if self.monitor<=1e-6:
            logging.info('                               low monitor -> nan output')
            proj.data.ref *= np.nan

        return proj
