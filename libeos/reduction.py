import logging
import os
import sys

import numpy as np
from orsopy import fileio

from .command_line import expand_file_list
from .dataset import AmorData, Header
from .options import DataReaderConfig, ReductionConfig, OutputConfig
from .instrument import Grid

class AmorReduction:
    def __init__(self,
                  data_reader_config: DataReaderConfig,
                  reduction_config: ReductionConfig,
                  output_config: OutputConfig):
        self.data_reader_config = data_reader_config
        self.reduction_config = reduction_config
        self.output_config = output_config

    def reduce(self):
        self.grid = Grid(self.reduction_config.qResolution)
        self.header = Header()
        self.startTime = 0
        if not os.path.exists(f'{self.data_reader_config.dataPath}'):
            os.system(f'mkdir {self.data_reader_config.dataPath}')
        fromHDF = AmorData(self.startTime, header=self.header, config=self.data_reader_config)
        logging.warning('\n######## eos - data reduction for Amor ########')

        # load or create normalisation matrix
        if self.reduction_config.normalisationFileIdentifier:
            normalise = True
            norm_lz, normAngle, normFileList = self.normalisation_map(self.reduction_config.normalisationFileIdentifier[0])
            self.header.reduction.corrections.append('normalisation with \'additional files\'')
        else:
            normalise = False
            norm_lz = self.grid.lz()
            normAngle = 1.

            logging.warning('normalisation matrix: none requested')

        # load R(q_z) curve to be subtracted:
        if self.reduction_config.subtract:
            sq_q, sR_q, sdR_q, sFileName = self.loadRqz(self.reduction_config.subtract)
            subtract = True
            logging.warning(f'loaded background file: {sFileName}')
            self.header.reduction.corrections.append(f'background from \'{sFileName}\' subtracted')
        else:
            subtract = False

        # load measurement data and do the reduction
        datasetsRqz = []
        datasetsRlt = []
        for i, short_notation in enumerate(self.reduction_config.fileIdentifier):
            logging.warning('reading input:')
            self.header.measurement_data_files = []
            fromHDF.read_data(short_notation)

            if self.reduction_config.timeSlize:
                wallTime_e = fromHDF.wallTime_e
                columns = self.header.columns()+[fileio.Column('time', 's', 'time relative to start of measurement series')]
                headerRqz = fileio.Orso(self.header.data_source(), self.header.reduction, columns)

                interval = self.reduction_config.timeSlize[0]
                try:
                    start = self.reduction_config.timeSlize[1]
                except:
                    start = 0
                try:
                    stop = self.reduction_config.timeSlize[2]
                except:
                    stop = wallTime_e[-1]
                # make overwriting log lines possible by removing newline at the end
                logging.StreamHandler.terminator = "\r"
                for i, time in enumerate(np.arange(start, stop, interval)):
                    logging.warning(f' time slize {i:4d}')

                    filter_e = np.where((time<wallTime_e) & (wallTime_e<time+interval), True, False)
                    lamda_e = fromHDF.lamda_e[filter_e]
                    detZ_e = fromHDF.detZ_e[filter_e]

                    qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz = self.project_on_lz(
                            fromHDF, norm_lz, normAngle, lamda_e, detZ_e)
                    q_q, R_q, dR_q, dq_q = self.project_on_qz(qz_lz, ref_lz, err_lz, res_lz, norm_lz, mask_lz)

                    filter_q = np.where((self.data_reader_config.qzRange[0]<q_q) & (q_q<self.data_reader_config.qzRange[1]),
                                        True, False)
                    q_q = q_q[filter_q]
                    R_q = R_q[filter_q]
                    dR_q = dR_q[filter_q]
                    dq_q = dq_q[filter_q]

                    if self.reduction_config.autoscale:
                        R_q, dR_q = self.autoscale(q_q, R_q, dR_q)

                    if subtract:
                        if len(q_q)==len(sq_q):
                            R_q -= sR_q
                            dR_q = np.sqrt(dR_q**2+sdR_q**2)
                        else:
                            subtract = False
                            logging.warning(
                                f'background file {sFileName} not compatible with q_z scale ({len(sq_q)} vs. {len(q_q)})')

                    tme_q = np.ones(np.shape(q_q))*time
                    data = np.array([q_q, R_q, dR_q, dq_q, tme_q]).T
                    headerRqz.data_set = f'{i}: time = {time:8.1f} s  to {time+interval:8.1f} s'
                    orso_data = fileio.OrsoDataset(headerRqz, data)
                    # make a copy of the header for the next iteration
                    headerRqz = fileio.Orso.from_dict(headerRqz.to_dict())
                    datasetsRqz.append(orso_data)
                # reset normal logging behavior
                logging.StreamHandler.terminator = "\n"
                logging.warning(f' time slize {i:4d}, done')
            else:
                lamda_e = fromHDF.lamda_e
                detZ_e = fromHDF.detZ_e

                qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz = self.project_on_lz(
                        fromHDF, norm_lz, normAngle, lamda_e, detZ_e)
                try:
                    ref_lz *= self.reduction_config.scale[i]
                    err_lz *= self.reduction_config.scale[i]
                except IndexError:
                    ref_lz *= self.reduction_config.scale[-1]
                    err_lz *= self.reduction_config.scale[-1]

                if 'Rqz.ort' in self.output_config.outputFormats:
                    headerRqz = fileio.Orso(self.header.data_source(), self.header.reduction, self.header.columns())
                    headerRqz.data_set = f'Nr {i} : mu = {fromHDF.mu:6.3f} deg'

                    # projection on q-grid
                    q_q, R_q, dR_q, dq_q = self.project_on_qz(qz_lz, ref_lz, err_lz, res_lz, norm_lz, mask_lz, self.grid)

                    filter_q = np.where((self.reduction_config.qzRange[0]<q_q) & (q_q<self.reduction_config.qzRange[1]), True, False)
                    q_q = q_q[filter_q]
                    R_q = R_q[filter_q]
                    dR_q = dR_q[filter_q]
                    dq_q = dq_q[filter_q]

                    if self.reduction_config.autoscale:
                        if i==0:
                            R_q, dR_q = self.autoscale(q_q, R_q, dR_q)
                        else:
                            pRq_z = datasetsRqz[i-1].data[:, 1]
                            pdRq_z = datasetsRqz[i-1].data[:, 2]
                            R_q, dR_q = self.autoscale(q_q, R_q, dR_q, pRq_z, pdRq_z)

                    if subtract:
                        if len(q_q)==len(sq_q):
                            R_q -= sR_q
                            dR_q = np.sqrt(dR_q**2+sdR_q**2)
                        else:
                            logging.warning(
                                f'backgroung file {sFileName} not compatible with q_z scale ({len(sq_q)} vs. {len(q_q)})')

                    data = np.array([q_q, R_q, dR_q, dq_q]).T
                    orso_data = fileio.OrsoDataset(headerRqz, data)
                    headerRqz = fileio.Orso(**headerRqz.to_dict())
                    datasetsRqz.append(orso_data)

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
                        ]
                    # data_source = fromHDF.data_source
                    headerRlt = fileio.Orso(self.header.data_source, self.header.reduction, columns)

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
                            norm_lz.T,
                            np.where(mask_lz, 1, 0).T
                            ):
                        data = np.array(list(item)).T
                        headerRlt = fileio.Orso(**headerRlt.to_dict())
                        headerRlt.data_set = f'dataset_{i}_{j+1} : alpha_f = {theta_lz[0, j]:6.3f} deg'
                        orso_data = fileio.OrsoDataset(headerRlt, data)
                        datasetsRlt.append(orso_data)
                        j += 1

        # output
        logging.warning('output:')

        if 'Rqz.ort' in self.output_config.outputFormats:
            logging.warning(f'  {self.data_reader_config.dataPath}/{self.output_config.outputName}.Rqz.ort')
            theSecondLine = f' {self.header.experiment.title} | {self.header.experiment.start_date} | sample {self.header.sample.name} | R(q_z)'
            fileio.save_orso(datasetsRqz, f'{self.data_reader_config.dataPath}/{self.output_config.outputName}.Rqz.ort', data_separator='\n',
                             comment=theSecondLine)

        if 'Rlt.ort' in self.output_config.outputFormats:
            logging.warning(f'  {self.data_reader_config.dataPath}/{self.output_config.outputName}.Rlt.ort')
            theSecondLine = f' {self.header.experiment.title} | {self.header.experiment.start_date} | sample {self.header.sample.name} | R(lambda, theta)'
            fileio.save_orso(datasetsRlt, f'{self.data_reader_config.dataPath}/{self.output_config.outputName}.Rlt.ort', data_separator='\n',
                             comment=theSecondLine)


    def autoscale(self, q_q, R_q, dR_q, pR_q=[], pdR_q=[]):
        autoscale = self.reduction_config.autoscale
        if len(pR_q) == 0:
            filter_q  = np.where((autoscale[0]<=q_q)&(q_q<=autoscale[1]), True, False)
            filter_q  = np.where(dR_q>0, filter_q, False)
            if len(filter_q[filter_q]) > 0:
                scale = np.sum(R_q[filter_q]**2/dR_q[filter_q]) / np.sum(R_q[filter_q]/dR_q[filter_q])
            else:
                logging.warning(f'#     automatic scaling not possible')
                scale = 1.
        else:
            filter_q  = np.where(np.isnan(pR_q*R_q), False, True)
            filter_q  = np.where(R_q>0, filter_q, False)
            filter_q  = np.where(pR_q>0, filter_q, False)
            if len(filter_q[filter_q]) > 0:
                scale = np.sum(R_q[filter_q]**3 * pR_q[filter_q] / (dR_q[filter_q]**2 * pdR_q[filter_q]**2)) \
                      / np.sum(R_q[filter_q]**2 * pR_q[filter_q]**2 / (dR_q[filter_q]**2  * pdR_q[filter_q]**2))
            else:
                logging.warning(f'#     automatic scaling not possible')
                scale = 1.
        R_q  /= scale
        dR_q /= scale
        logging.debug(f'#     scaling factor = {scale}')

        return R_q, dR_q

    def project_on_qz(self, q_lz, R_lz, dR_lz, dq_lz, norm_lz, mask_lz):
        q_q       = self.grid.q()
        mask_lzf  = mask_lz.flatten()
        q_lzf     = q_lz.flatten()[mask_lzf]
        R_lzf     = R_lz.flatten()[mask_lzf]
        dR_lzf    = dR_lz.flatten()[mask_lzf]
        dq_lzf    = dq_lz.flatten()[mask_lzf]
        norm_lzf  = norm_lz.flatten()[mask_lzf]

        N_q       = np.histogram(q_lzf, bins = q_q, weights = norm_lzf )[0]
        N_q       = np.where(N_q > 0, N_q, np.nan)

        R_q       = np.histogram(q_lzf, bins = q_q, weights = norm_lzf * R_lzf )[0]
        R_q       = R_q / N_q

        dR_q      = np.histogram(q_lzf, bins = q_q, weights = (norm_lzf * dR_lzf)**2 )[0]
        dR_q      = np.sqrt( dR_q ) / N_q

        # TODO: different error propagations for dR and dq!
        N_q       = np.histogram(q_lzf, bins = q_q, weights = norm_lzf**2 )[0]
        N_q       = np.where(N_q > 0, N_q, np.nan)
        dq_q      = np.histogram(q_lzf, bins = q_q, weights = (norm_lzf * dq_lzf)**2 )[0]
        dq_q      = np.sqrt( dq_q / N_q )

        q_q       = 0.5 * (q_q + np.roll(q_q, 1))

        return q_q[1:], R_q, dR_q, dq_q

    def loadRqz(self, name):
        if os.path.exists(f'{self.data_reader_config.dataPath}/{name}'):
            fileName = f'{self.data_reader_config.dataPath}/{name}'
        elif os.path.exists(f'{self.data_reader_config.dataPath}/{name}.Rqz.ort'):
            fileName = f'{self.data_reader_config.dataPath}/{name}.Rqz.ort'
        else:
            sys.exit(f'### the background file \'{self.data_reader_config.dataPath}/{name}\' does not exist! => stopping')

        q_q, Sq_q, dS_q = np.loadtxt(fileName, usecols=(0, 1, 2), comments='#', unpack=True)

        return q_q, Sq_q, dS_q, fileName

    def normalisation_map(self, short_notation):
        dataPath = self.data_reader_config.dataPath
        fromHDF = AmorData(self.startTime, self.header, self.data_reader_config)
        normalisation_list = expand_file_list(short_notation)
        name = str(normalisation_list[0])
        for i in range(1, len(normalisation_list), 1):
            name = f'{name}_{normalisation_list[i]}'
        if os.path.exists(f'{dataPath}/{name}.norm'):
            logging.info(f'# normalisation matrix: found and using {dataPath}/{name}.norm')
            norm_lz = np.loadtxt(f'{dataPath}/{name}.norm')
            fh = open(f'{dataPath}/{name}.norm', 'r')
            fh.readline()
            normFileList = fh.readline().split('[')[1].split(']')[0].replace('\'', '').split(', ')
            normAngle = float(fh.readline().split('= ')[1])
            fh.close()
            for i, entry in enumerate(normFileList):
                 normFileList[i] = entry.split('/')[-1]
            self.header.measurement_additional_files = normFileList
        else:
            logging.info(f'# normalisation matrix: using the files {normalisation_list}')
            fromHDF.read_data(short_notation, norm=True)
            normAngle     = fromHDF.nu - fromHDF.mu
            lamda_e       = fromHDF.lamda_e
            detZ_e        = fromHDF.detZ_e
            norm_lz, bins_l, bins_z = np.histogram2d(lamda_e, detZ_e, bins = (self.grid.lamda(), self.grid.z()))
            norm_lz = np.where(norm_lz>0, norm_lz, np.nan)
            # correct for the SM reflectivity
            lamda_l  = self.grid.lamda()
            theta_z  = normAngle + fromHDF.delta_z
            lamda_lz = (self.grid.lz().T*lamda_l[:-1]).T
            theta_lz = self.grid.lz()*theta_z
            qz_lz    = 4.0*np.pi * np.sin(np.deg2rad(theta_lz)) / lamda_lz
            Rsm_lz   = np.ones(np.shape(qz_lz))
            Rsm_lz   = np.where(qz_lz>0.0217, 1-(qz_lz-0.0217)*(0.0625/0.0217), Rsm_lz)
            Rsm_lz   = np.where(qz_lz>0.0217*5, np.nan, Rsm_lz)
            norm_lz  = norm_lz / Rsm_lz
            if len(lamda_e) > 1e6:
                head = ('normalisation matrix based on the measurements\n'
                       f'{fromHDF.file_list}\n'
                       f'nu - mu = {normAngle}\n'
                       f'shape= {np.shape(norm_lz)} (lambda, z)\n'
                       f'measured at mu = {fromHDF.mu:6.3f} deg\n'
                       f'N(l_lambda, z) = theta(z) / sum_i=-1..1 I(l_lambda+i, z)')
                head = head.replace('../', '')
                head = head.replace('./', '')
                head = head.replace('raw/', '')
                np.savetxt(f'{dataPath}/{name}.norm', norm_lz, header = head)
            normFileList = fromHDF.file_list
        return norm_lz, normAngle, normFileList


    def project_on_lz(self, fromHDF, norm_lz, normAngle, lamda_e, detZ_e):
        # projection on lambda-z-grid
        lamda_l  = self.grid.lamda()
        theta_z  = fromHDF.nu - fromHDF.mu + fromHDF.delta_z
        lamda_lz = (self.grid.lz().T*lamda_l[:-1]).T
        theta_lz = self.grid.lz()*theta_z

        thetaN_z  = fromHDF.delta_z + normAngle
        thetaN_lz = np.ones(np.shape(norm_lz))*thetaN_z
        thetaN_lz = np.where(np.absolute(thetaN_lz)>5e-3, thetaN_lz, np.nan)

        mask_lz   = np.where(np.isnan(norm_lz), False, True)
        mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(thetaN_lz)>5e-3, True, False))
        mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(theta_lz)>5e-3, True, False))
        if self.reduction_config.thetaRange[1]<12:
          mask_lz   = np.logical_and(mask_lz, np.where(theta_lz >= self.reduction_config.thetaRange[0], True, False))
          mask_lz   = np.logical_and(mask_lz, np.where(theta_lz <= self.reduction_config.thetaRange[1], True, False))
        if self.reduction_config.thetaRangeR[1]<12:
          t0 = fromHDF.nu - fromHDF.mu
          mask_lz   = np.logical_and(mask_lz, np.where(theta_lz-t0 >= self.reduction_config.thetaRangeR[0], True, False))
          mask_lz   = np.logical_and(mask_lz, np.where(theta_lz-t0 <= self.reduction_config.thetaRangeR[1], True, False))
        if self.reduction_config.lambdaRange[1]<15:
          mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz >= self.reduction_config.lambdaRange[0], True, False))
          mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz <= self.reduction_config.lambdaRange[1], True, False))

        #           gravity correction
        #theta_lz += np.rad2deg( np.arctan( 3.07e-10 * (fromHDF.detectorDistance + detXdist_e) * lamda_lz**2 ) )
        theta_lz += np.rad2deg( np.arctan( 3.07e-10 * fromHDF.detectorDistance * lamda_lz**2 ) )

        z_z       = enumerate(theta_z)
        qz_lz     = 4.0*np.pi * np.sin(np.deg2rad(theta_lz)) / lamda_lz
        int_lz, bins_l, bins_z  = np.histogram2d(lamda_e, detZ_e, bins = (lamda_l, self.grid.z()))
        #           cut normalisation sample horizon
        int_lz    = np.where(mask_lz, int_lz, np.nan)
        thetaF_lz  = np.where(mask_lz, theta_lz, np.nan)

        ref_lz    = (int_lz * np.absolute(thetaN_lz)) / (norm_lz * np.absolute(thetaF_lz))
        err_lz    = ref_lz * np.sqrt( 1/(int_lz+.1) + 1/norm_lz )

        res_lz    = np.ones((np.shape(lamda_l[:-1])[0], np.shape(theta_z)[0])) * 0.022**2
        res_lz    = res_lz + (0.008/theta_lz)**2
        res_lz    = qz_lz * np.sqrt(res_lz)

        return qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz

