import logging
import os
import sys

import numpy as np
from orsopy import fileio

from .command_line import expand_file_list
from .file_reader import AmorData
from .header import Header
from .options import EOSConfig
from .instrument import Grid

class AmorReduction:
    def __init__(self, config: EOSConfig):
        self.experiment_config = config.experiment
        self.reader_config = config.reader
        self.reduction_config = config.reduction
        self.output_config = config.output
        self.grid = Grid(config.reduction.qResolution, config.experiment.qzRange)

        self.header = Header()
        self.header.reduction.call = config.call_string()

        self.monitorUnit = {'n': 'cnts', 'p': 'mC', 't': 's'}

    def reduce(self):
        if not os.path.exists(f'{self.output_config.outputPath}'):
            logging.debug(f'Creating destination path {self.output_config.outputPath}')
            os.system(f'mkdir {self.output_config.outputPath}')

        # load or create normalisation matrix
        if self.reduction_config.normalisationFileIdentifier:
            self.create_normalisation_map(self.reduction_config.normalisationFileIdentifier[0])
        else:
            self.norm_lz = self.grid.lz()
            self.normAngle = 1.
            self.normMonitor = 1.

            logging.warning('normalisation matrix: none requested')

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
        self.header.measurement_data_files = []
        self.file_reader = AmorData(header=self.header,
                                    reader_config=self.reader_config,
                                    config=self.experiment_config,
                                    short_notation=short_notation)
        if self.reduction_config.timeSlize:
            self.read_timeslices(i)
        else:
            self.read_unsliced(i)

    def read_unsliced(self, i):
        lamda_e = self.file_reader.lamda_e
        detZ_e  = self.file_reader.detZ_e
        self.monitor = np.sum(self.file_reader.monitorPerPulse)
        logging.warning(f'    monitor = {self.monitor:8.2f} {self.monitorUnit[self.experiment_config.monitorType]}')
        qz_lz, qx_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, self.mask_lz = self.project_on_lz(
                self.file_reader, self.norm_lz, self.normAngle, lamda_e, detZ_e)
        #if self.monitor>1 :
        #    ref_lz /= self.monitor
        #    err_lz /= self.monitor
        try:
            ref_lz *= self.reduction_config.scale[i]
            err_lz *= self.reduction_config.scale[i]
        except IndexError:
            ref_lz *= self.reduction_config.scale[-1]
            err_lz *= self.reduction_config.scale[-1]
        if 'Rqz.ort' in self.output_config.outputFormats:
            headerRqz = self.header.orso_header()
            headerRqz.data_set = f'Nr {i} : mu = {self.file_reader.mu:6.3f} deg'

            if qz_lz[0,int(np.shape(qz_lz)[1]/2)]  < 0:
                # assuming a 'measurement from below' when center of detector at negative qz
                qz_lz *= -1

            # projection on qz-grid
            q_q, R_q, dR_q, dq_q = self.project_on_qz(qz_lz, ref_lz, err_lz, res_lz, self.norm_lz, self.mask_lz)

            # The filtering is now done by restricting the qz-grid
            #filter_q = np.where((self.experiment_config.qzRange[0]>q_q) & (q_q>self.experiment_config.qzRange[1]),
            #                    False, True)
            #q_q = q_q[filter_q]
            #R_q = R_q[filter_q]
            #dR_q = dR_q[filter_q]
            #dq_q = dq_q[filter_q]

            if self.reduction_config.autoscale:
                if i==0:
                    R_q, dR_q = self.autoscale(q_q, R_q, dR_q)
                else:
                    pRq_z = self.datasetsRqz[i-1].data[:, 1]
                    pdRq_z = self.datasetsRqz[i-1].data[:, 2]
                    R_q, dR_q = self.autoscale(q_q, R_q, dR_q, pRq_z, pdRq_z)

            if self.subtract:
                if len(q_q)==len(self.sq_q):
                    R_q -= self.sR_q
                    dR_q = np.sqrt(dR_q**2+self.sdR_q**2)
                else:
                    logging.warning(
                            f'backgroung file {self.sFileName} not compatible with q_z scale ({len(self.sq_q)} vs. {len(q_q)})')

            data = np.array([q_q, R_q, dR_q, dq_q]).T
            orso_data = fileio.OrsoDataset(headerRqz, data)
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
            # data_source = file_reader.data_source

            ts, zs = ref_lz.shape
            lindex_lz = np.tile(np.arange(1, ts+1), (zs, 1)).T
            tindex_lz = np.tile(np.arange(1, zs+1), (ts, 1))

            j = 0
            for item in zip(
                    qz_lz.T,
                    ref_lz.T,
                    err_lz.T,
                    res_lz.T,
                    lamda_lz.T,
                    theta_lz.T,
                    lindex_lz.T,
                    tindex_lz.T,
                    int_lz.T,
                    self.norm_lz.T,
                    np.where(self.mask_lz, 1, 0).T,
                    qx_lz.T,
                    ):
                data = np.array(list(item)).T
                headerRlt = self.header.orso_header(columns=columns)
                headerRlt.data_set = f'dataset_{i}_{j+1} : alpha_f = {theta_lz[0, j]:6.3f} deg'
                orso_data = fileio.OrsoDataset(headerRlt, data)
                self.datasetsRlt.append(orso_data)
                j += 1

    def read_timeslices(self, i):
        wallTime_e = np.float64(self.file_reader.wallTime_e)/1e9
        pulseTimeS = np.float64(self.file_reader.pulseTimeS)/1e9
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

            filter_e = np.where((time<wallTime_e) & (wallTime_e<time+interval), True, False)
            lamda_e = self.file_reader.lamda_e[filter_e]
            detZ_e = self.file_reader.detZ_e[filter_e]
            filter_m = np.where((time<pulseTimeS) & (pulseTimeS<time+interval), True, False)
            self.monitor = np.sum(self.file_reader.monitorPerPulse[filter_m])
            logging.info(f'      {ti:<4d}  {time:6.0f}  {self.monitor:7.2f} {self.monitorUnit[self.experiment_config.monitorType]}')

            qz_lz, qx_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz = self.project_on_lz(
                    self.file_reader, self.norm_lz, self.normAngle, lamda_e, detZ_e)
            try:
                ref_lz *= self.reduction_config.scale[i]
                err_lz *= self.reduction_config.scale[i]
            except IndexError:
                ref_lz *= self.reduction_config.scale[-1]
                err_lz *= self.reduction_config.scale[-1]
            q_q, R_q, dR_q, dq_q = self.project_on_qz(qz_lz, ref_lz, err_lz, res_lz, self.norm_lz, mask_lz)

            filter_q = np.where((self.experiment_config.qzRange[0]<q_q) & (q_q<self.experiment_config.qzRange[1]),
                                True, False)
            q_q = q_q[filter_q]
            R_q = R_q[filter_q]
            dR_q = dR_q[filter_q]
            dq_q = dq_q[filter_q]

            if self.reduction_config.autoscale:
                R_q, dR_q = self.autoscale(q_q, R_q, dR_q)

            if self.subtract:
                if len(q_q)==len(self.sq_q):
                    R_q -= self.sR_q
                    dR_q = np.sqrt(dR_q**2+self.sdR_q**2)
                else:
                    self.subtract = False
                    logging.warning(
                            f'background file {self.sFileName} not compatible with q_z scale ({len(self.sq_q)} vs. {len(q_q)})')

            tme_q = np.ones(np.shape(q_q))*time
            data = np.array([q_q, R_q, dR_q, dq_q, tme_q]).T
            headerRqz = self.header.orso_header(
                    extra_columns=[fileio.Column('time', 's', 'time relative to start of measurement series')])
            headerRqz.data_set = f'{i}_{ti}: time = {time:8.1f} s  to {time+interval:8.1f} s'
            orso_data = fileio.OrsoDataset(headerRqz, data)
            self.datasetsRqz.append(orso_data)
        # reset normal logging behavior
        #logging.StreamHandler.terminator = "\n"
        logging.info(f'      done  {time+interval:5.0f}')

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

    def autoscale(self, q_q, R_q, dR_q, pR_q=[], pdR_q=[]):
        autoscale = self.reduction_config.autoscale
        if len(pR_q) == 0:
            filter_q  = np.where((autoscale[0]<=q_q)&(q_q<=autoscale[1]), True, False)
            filter_q  = np.where(dR_q>0, filter_q, False)
            if len(filter_q[filter_q]) > 0:
                scale = np.sum(R_q[filter_q]**2/dR_q[filter_q]) / np.sum(R_q[filter_q]/dR_q[filter_q])
            else:
                logging.warning('      automatic scaling not possible')
                scale = 1.
        else:
            filter_q  = np.where(np.isnan(pR_q*R_q), False, True)
            filter_q  = np.where(R_q>0, filter_q, False)
            filter_q  = np.where(pR_q>0, filter_q, False)
            if len(filter_q[filter_q]) > 0:
                scale = np.sum(R_q[filter_q]**3 * pR_q[filter_q] / (dR_q[filter_q]**2 * pdR_q[filter_q]**2)) \
                      / np.sum(R_q[filter_q]**2 * pR_q[filter_q]**2 / (dR_q[filter_q]**2  * pdR_q[filter_q]**2))
            else:
                logging.warning('      automatic scaling not possible')
                scale = 1.
        R_q  /= scale
        dR_q /= scale
        logging.info(f'      scaling factor = {1/scale}')

        return R_q, dR_q

    def project_on_qz(self, q_lz, R_lz, dR_lz, dq_lz, norm_lz, mask_lz):
        q_q       = self.grid.q()
        mask_lzf  = mask_lz.flatten()
        q_lzf     = q_lz.flatten()[mask_lzf]
        R_lzf     = R_lz.flatten()[mask_lzf]
        dR_lzf    = dR_lz.flatten()[mask_lzf]
        dq_lzf    = dq_lz.flatten()[mask_lzf]
        norm_lzf  = norm_lz.flatten()[mask_lzf]

        weights_lzf = norm_lzf
        #weights_lzf = np.sqrt(norm_lzf)
        #weights_lzf = 1 / dR_lzf

        N_q       = np.histogram(q_lzf, bins = q_q, weights = weights_lzf )[0]
        N_q       = np.where(N_q > 0, N_q, np.nan)

        R_q       = np.histogram(q_lzf, bins = q_q, weights = weights_lzf * R_lzf )[0]
        R_q       = R_q / N_q

        dR_q      = np.histogram(q_lzf, bins = q_q, weights = (weights_lzf * dR_lzf)**2 )[0]
        dR_q      = np.sqrt( dR_q ) / N_q

        # TODO: different error propagations for dR and dq!
        # this is what should work:
        #dq_q      = np.histogram(q_lzf, bins = q_q, weights = (weights_lzf * dq_lzf)**2 )[0]
        #dq_q      = np.sqrt( dq_q ) / N_q 
        # and this actually works:
        N_q       = np.histogram(q_lzf, bins = q_q, weights = weights_lzf**2 )[0]
        N_q       = np.where(N_q > 0, N_q, np.nan)
        dq_q      = np.histogram(q_lzf, bins = q_q, weights = (weights_lzf * dq_lzf)**2 )[0]
        dq_q      = np.sqrt( dq_q / N_q )

        q_q       = 0.5 * (q_q + np.roll(q_q, 1))

        return q_q[1:], R_q, dR_q, dq_q

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
        normalisation_list = expand_file_list(short_notation)
        name = str(normalisation_list[0])
        for i in range(1, len(normalisation_list), 1):
            name = f'{name}_{normalisation_list[i]}'
        n_path = os.path.join(outputPath, f'{name}.norm')
        if os.path.exists(n_path):
            logging.warning(f'normalisation matrix: found and using {n_path}')
            with open(n_path, 'rb') as fh:
                self.normFileList = np.load(fh, allow_pickle=True)
                self.normAngle    = np.load(fh, allow_pickle=True)
                self.norm_lz      = np.load(fh, allow_pickle=True)
                self.normMonitor  = np.load(fh, allow_pickle=True)
            for i, entry in enumerate(self.normFileList):
                 self.normFileList[i] = entry.split('/')[-1]
            self.header.measurement_additional_files = self.normFileList
        else:
            logging.warning(f'normalisation matrix: using the files {normalisation_list}')
            fromHDF = AmorData(header=self.header,
                               reader_config=self.reader_config,
                               config=self.experiment_config,
                               short_notation=short_notation, norm=True)
            self.normAngle     = fromHDF.nu - fromHDF.mu
            lamda_e          = fromHDF.lamda_e
            detZ_e           = fromHDF.detZ_e
            self.normMonitor = np.sum(fromHDF.monitorPerPulse)
            norm_lz, bins_l, bins_z = np.histogram2d(lamda_e, detZ_e, bins = (self.grid.lamda(), self.grid.z()))
            norm_lz = np.where(norm_lz>2, norm_lz, np.nan)
            if self.reduction_config.normalisationMethod == 'd':
                # direct reference => invert map vertically
                self.norm_lz = np.flip(norm_lz, 1)
            else:
                # correct for reference sm reflectivity
                lamda_l  = self.grid.lamda()
                theta_z  = self.normAngle + fromHDF.delta_z
                lamda_lz = (self.grid.lz().T*lamda_l[:-1]).T
                theta_lz = self.grid.lz()*theta_z
                qz_lz    = 4.0*np.pi * np.sin(np.deg2rad(theta_lz)) / lamda_lz
                # TODO: introduce variable for `m` and propably for the slope
                Rsm_lz   = np.ones(np.shape(qz_lz))
                Rsm_lz   = np.where(qz_lz>0.0217, 1-(qz_lz-0.0217)*(0.0625/0.0217), Rsm_lz)
                Rsm_lz   = np.where(qz_lz>0.0217*5, np.nan, Rsm_lz)
                self.norm_lz  = norm_lz / Rsm_lz

            if len(lamda_e) > 1e6:
                with open(n_path, 'wb') as fh:
                    np.save(fh, np.array(fromHDF.file_list), allow_pickle=False)
                    np.save(fh, np.array(self.normAngle), allow_pickle=False)
                    np.save(fh, self.norm_lz, allow_pickle=False)
                    np.save(fh, self.normMonitor, allow_pickle=False)
            self.normFileList = fromHDF.file_list
        self.header.reduction.corrections.append('normalisation with \'additional files\'')

    def project_on_lz(self, fromHDF, norm_lz, normAngle, lamda_e, detZ_e):
        # projection on lambda-z-grid
        lamda_l  = self.grid.lamda()
        alphaF_z  = fromHDF.nu - fromHDF.mu + fromHDF.delta_z
        # TODO: implement various methods to obtain alpha_i.
        #if self.experiment_config.incidentAngle == 'alphaF':
        #  # for specular reflectometry with a highly divergent beam
        #  alphaF_z  = fromHDF.nu - fromHDF.mu + fromHDF.delta_z
        #elif self.experiment_config.incidentAngle == 'nu':
        #  # for specular reflectometry, using kappa nad nu but ignoring mu
        #  alphaF_z  = (fromHDF.nu + fromHDF.delta_z + fromHDF.kap + fromHDF.kad) / 2.
        #else:
        #  # using kappa, for a collimated incoming beam
        #  pass
        lamda_lz  = (self.grid.lz().T*lamda_l[:-1]).T
        alphaF_lz = self.grid.lz()*alphaF_z

        mask_lz   = np.where(np.isnan(norm_lz), False, True)
        mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(alphaF_lz)>5e-3, True, False))
        if self.reduction_config.thetaRangeR[1]<12:
          t0 = fromHDF.nu - fromHDF.mu
          mask_lz   = np.logical_and(mask_lz, np.where(alphaF_lz-t0 >= self.reduction_config.thetaRangeR[0], True, False))
          mask_lz   = np.logical_and(mask_lz, np.where(alphaF_lz-t0 <= self.reduction_config.thetaRangeR[1], True, False))
        elif self.reduction_config.thetaRange[1]<12:
          mask_lz   = np.logical_and(mask_lz, np.where(alphaF_lz >= self.reduction_config.thetaRange[0], True, False))
          mask_lz   = np.logical_and(mask_lz, np.where(alphaF_lz <= self.reduction_config.thetaRange[1], True, False))
        else:
          self.reduction_config.thetaRange = [fromHDF.nu - fromHDF.mu - fromHDF.div/2, 
                                              fromHDF.nu - fromHDF.mu + fromHDF.div/2]
          mask_lz   = np.logical_and(mask_lz, np.where(alphaF_lz >= self.reduction_config.thetaRange[0], True, False))
          mask_lz   = np.logical_and(mask_lz, np.where(alphaF_lz <= self.reduction_config.thetaRange[1], True, False))
        if self.experiment_config.lambdaRange[1]<15:
          mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz >= self.experiment_config.lambdaRange[0], True, False))
          mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz <= self.experiment_config.lambdaRange[1], True, False))

        #           gravity correction
        #alphaF_lz += np.rad2deg( np.arctan( 3.07e-10 * (fromHDF.detectorDistance + detXdist_e) * lamda_lz**2 ) )
        alphaF_lz += np.rad2deg( np.arctan( 3.07e-10 * fromHDF.detectorDistance * lamda_lz**2 ) )

        if self.experiment_config.incidentAngle == 'alphaF':
          #alphaI_lz = alphaF_lz
          qz_lz     = 4.0*np.pi * np.sin(np.deg2rad(alphaF_lz)) / lamda_lz
          qx_lz     = self.grid.lz() * 0.
        else:
          alphaI_lz = self.grid.lz()*(fromHDF.mu + fromHDF.kap + fromHDF.kad)
          qz_lz     = 2.0*np.pi * (np.sin(np.deg2rad(alphaF_lz)) + np.sin(np.deg2rad(alphaI_lz))) / lamda_lz
          qx_lz     = 2.0*np.pi * (np.cos(np.deg2rad(alphaF_lz)) - np.cos(np.deg2rad(alphaI_lz))) / lamda_lz

        int_lz, bins_l, bins_z  = np.histogram2d(lamda_e, detZ_e, bins = (lamda_l, self.grid.z()))
        #           cut normalisation sample horizon
        int_lz    = np.where(mask_lz, int_lz, np.nan)
        thetaF_lz = np.where(mask_lz, alphaF_lz, np.nan)

        if self.reduction_config.normalisationMethod == 'o':
            logging.debug('      assuming an overilluminated sample and correcting for the angle of incidence')
            thetaN_z  = fromHDF.delta_z + normAngle
            thetaN_lz = np.ones(np.shape(norm_lz))*thetaN_z
            thetaN_lz = np.where(np.absolute(thetaN_lz)>5e-3, thetaN_lz, np.nan)
            mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(thetaN_lz)>5e-3, True, False))
            ref_lz    = (int_lz * np.absolute(thetaN_lz)) / (norm_lz * np.absolute(thetaF_lz))
        elif self.reduction_config.normalisationMethod == 'u':
            logging.debug('      assuming an underilluminated sample and ignoring the angle of incidence')
            ref_lz    = (int_lz / norm_lz)
        elif self.reduction_config.normalisationMethod == 'd':
            logging.debug('      assuming direct beam for normalisation and ignoring the angle of incidence')
            ref_lz    = (int_lz / norm_lz)
        else:
            logging.error('unknown normalisation method! Use [u]nder, [o]ver or [d]irect illumination')
            ref_lz    = (int_lz / norm_lz)
        if self.monitor > 1e-6 :
            ref_lz   *= self.normMonitor / self.monitor
        else:
            logging.info('       too small monitor value for normalisation -> ignoring monitors')
        err_lz    = ref_lz * np.sqrt( 1/(int_lz+.1) + 1/norm_lz ) 

        # TODO: allow for non-ideal Delta lambda / lambda (rather than 2.2%)
        res_lz    = np.ones((np.shape(lamda_l[:-1])[0], np.shape(alphaF_z)[0])) * 0.022**2
        res_lz    = res_lz + (0.008/alphaF_lz)**2
        res_lz    = qz_lz * np.sqrt(res_lz)

        return qz_lz, qx_lz, ref_lz, err_lz, res_lz, lamda_lz, alphaF_lz, int_lz, mask_lz


    @staticmethod
    def histogram2d_lz(lamda_e, detZ_e, bins):
        """
        Perform binning operation equivalent to numpy bin2d for the sepcial case
        of the second dimension using integer positions (pre-defined pixels).
        Based on the devide_bin algorithm below.
        """
        dimension = bins[1].shape[0]-1
        if not (np.array(bins[1])==np.arange(dimension+1)).all():
            raise ValueError("histogram2d_lz requires second bin dimension to be contigous integer range")
        binning = AmorReduction.devide_bin(lamda_e, detZ_e.astype(np.int64), bins[0], dimension)
        return np.array(binning), bins[0], bins[1]

    @staticmethod
    def devide_bin(lambda_e, position_e, lamda_edges, dimension):
        '''
        Use a divide and conquer strategy to bin the data. For the actual binning the
        numpy bincount function is used, as it is much faster than histogram for
        counting of integer values.

        :param lambda_e: Array of wavelength for each event
        :param position_e: Array of positional indices for each event
        :param lamda_edges: The edges of bins to be used for the histogram
        :param dimension: position number of buckets in output arrray

        :return: 2D list of dimensions (lambda, x) of counts
        '''
        if len(lambda_e)==0:
            # no more events in range, return empty bins
            return [np.zeros(dimension, dtype=np.int64).tolist()]*(len(lamda_edges)-1)
        if len(lamda_edges)==2:
            # deepest recursion reached, all items should be within the two ToF edges
            return [np.bincount(position_e, minlength=dimension).tolist()]
        # split all events into two time of flight regions
        split_idx = len(lamda_edges)//2
        left_region = lambda_e<lamda_edges[split_idx]
        left_list = AmorReduction.devide_bin(lambda_e[left_region], position_e[left_region],
                                             lamda_edges[:split_idx+1], dimension)
        right_region = np.logical_not(left_region)
        right_list = AmorReduction.devide_bin(lambda_e[right_region], position_e[right_region],
                                              lamda_edges[split_idx:], dimension)
        return left_list+right_list
