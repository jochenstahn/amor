#!/usr/bin/env python
"""
eos reduces measurements performed on Amor@SINQ, PSI

Author: Jochen Stahn

conventions (not strictly followed, yet):
- array names end with the suffix '_x[y]' with the meaning
    _e  = events
    _tof
    _l  = lambda
    _t  = theta
    _z  = detector z
    _lz = (lambda, detector z)
    _q  = q_z
- to come
"""


import os
import sys
import logging
import logging.handlers
import numpy as np
from orsopy import fileio

import libeos.reduction
from libeos.command_line import commandLineArgs, output_format_list
from libeos.dataset import AmorData, DataReaderConfig, Header
from libeos.instrument import Grid
from libeos.logconfig import setup_logging, update_loglevel
from libeos.reduction import autoscale, normalisation_map, project_on_lz, project_on_qz


#=====================================================================================================
# TODO:
# - calculate resolution using the chopperPhase
# - deal with background correction
# - format of 'call' + add '-Y' if not supplied
#=====================================================================================================
def loadRqz(name):

    if os.path.exists(f'{clas.dataPath}/{name}'):
        fileName = f'{clas.dataPath}/{name}'
    elif os.path.exists(f'{clas.dataPath}/{name}.Rqz.ort'):
        fileName = f'{clas.dataPath}/{name}.Rqz.ort'
    else:
        sys.exit(f'### the background file \'{clas.dataPath}/{name}\' does not exist! => stopping')

    q_q, Sq_q, dS_q = np.loadtxt(fileName, usecols=(0, 1, 2), comments='#', unpack=True)

    return q_q, Sq_q, dS_q, fileName

def main():
    setup_logging()

    global startTime, grid, clas, header
    clas   = commandLineArgs()
    update_loglevel(clas.verbose, clas.debug)

    grid   = Grid(clas.qResolution)
    header = Header()
    startTime = 0
    if not os.path.exists(f'{clas.dataPath}'):
        os.system(f'mkdir {clas.dataPath}')
    fromHDF = AmorData(startTime, header=header, config=DataReaderConfig(
                            year=clas.year,
                            dataPath=clas.dataPath,
                            sampleModel=clas.sampleModel,
                            chopperPhase=clas.chopperPhase,
                            chopperPhaseOffset=clas.chopperPhaseOffset,
                            yRange=clas.yRange,
                            lambdaRange=clas.lambdaRange,
                            qzRange=clas.qzRange,
                            offSpecular=clas.offSpecular,
                            mu=clas.mu,
                            nu=clas.nu,
                            muOffset=clas.muOffset
                            ))
    logging.warning('\n######## eos - data reduction for Amor ########')
    
    # load or create normalisation matrix
    if clas.normalisationFileIdentifier:
        normalise = True
        norm_lz, normAngle, normFileList = normalisation_map(clas.normalisationFileIdentifier[0],
                                                             header, grid, clas.dataPath)
        header.reduction.corrections.append('normalisation with \'additional files\'')
    else:
        normalise = False
        norm_lz   = grid.lz()
        normAngle = 1.

        logging.warning('normalisation matrix: none requested')

    # load R(q_z) curve to be subtracted:
    if clas.subtract:
        sq_q, sR_q, sdR_q, sFileName = loadRqz(clas.subtract)
        subtract = True
        logging.warning(f'loaded background file: {sFileName}')
        header.reduction.corrections.append(f'background from \'{sFileName}\' subtracted')
    else:
        subtract = False
        
    # load measurement data and do the reduction
    datasetsRqz = []
    datasetsRlt = []
    for i, short_notation in enumerate(clas.fileIdentifier):
        logging.warning('reading input:')
        header.measurement_data_files = []
        fromHDF.read_data(short_notation)

        if clas.timeSlize:
            wallTime_e = fromHDF.wallTime_e
            columns = header.columns() + [fileio.Column('time', 's', 'time relative to start of measurement series')]
            headerRqz = fileio.Orso(header.data_source(), header.reduction, columns)

            interval = clas.timeSlize[0]
            try:
                start = clas.timeSlize[1]
            except: 
                start = 0 
            try:
                stop  = clas.timeSlize[2]
            except:
                stop  = wallTime_e[-1]
            # make overwriting log lines possible by removing newline at the end
            logging.StreamHandler.terminator="\r"
            for i, time in enumerate(np.arange(start, stop, interval)):
                logging.warning(f' time slize {i:4d}')

                filter_e = np.where((time < wallTime_e) & (wallTime_e < time+interval), True, False)
                lamda_e  = fromHDF.lamda_e[filter_e] 
                detZ_e   = fromHDF.detZ_e[filter_e]
            
                qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz = project_on_lz(
                        fromHDF, norm_lz, normAngle, lamda_e, detZ_e,
                        grid, clas.thetaRange, clas.thetaRangeR, clas.lambdaRange)
                q_q, R_q, dR_q, dq_q = project_on_qz(qz_lz, ref_lz, err_lz, res_lz, norm_lz, mask_lz, grid)

                filter_q = np.where((clas.qzRange[0] < q_q) & (q_q < clas.qzRange[1]), True, False)
                q_q = q_q[filter_q]
                R_q = R_q[filter_q]
                dR_q = dR_q[filter_q]
                dq_q = dq_q[filter_q]

                if clas.autoscale:
                    R_q, dR_q = autoscale(q_q, R_q, dR_q, clas.autoscale)

                if subtract:
                    if len(q_q) == len(sq_q):
                        R_q  -= sR_q
                        dR_q = np.sqrt( dR_q**2 + sdR_q**2 )
                    else:
                        subtract = False
                        logging.warning(f'background file {sFileName} not compatible with q_z scale ({len(sq_q)} vs. {len(q_q)})')

                tme_q              = np.ones(np.shape(q_q))*time
                data               = np.array([q_q, R_q, dR_q, dq_q, tme_q]).T
                headerRqz.data_set = f'{i}: time = {time:8.1f} s  to {time+interval:8.1f} s'
                orso_data          = fileio.OrsoDataset(headerRqz, data)
                # make a copy of the header for the next iteration
                headerRqz          = fileio.Orso.from_dict(headerRqz.to_dict())
                datasetsRqz.append(orso_data)
            # reset normal logging behavior
            logging.StreamHandler.terminator="\n"
            logging.warning(f' time slize {i:4d}, done')
        else:
            lamda_e  = fromHDF.lamda_e
            detZ_e   = fromHDF.detZ_e

            qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz = project_on_lz(
                    fromHDF, norm_lz, normAngle, lamda_e, detZ_e,
                    grid, clas.thetaRange, clas.thetaRangeR, clas.lambdaRange
                    )
            try:
                ref_lz *= clas.scale[i]
                err_lz *= clas.scale[i]
            except:
                ref_lz *= clas.scale[-1]
                err_lz *= clas.scale[-1]

            if 'Rqz.ort' in output_format_list(clas.outputFormat):
                headerRqz = fileio.Orso(header.data_source(), header.reduction, header.columns())
                headerRqz.data_set = f'Nr {i} : mu = {fromHDF.mu:6.3f} deg'

                # projection on q-grid 
                q_q, R_q, dR_q, dq_q = project_on_qz(qz_lz, ref_lz, err_lz, res_lz, norm_lz, mask_lz, grid)

                filter_q = np.where((clas.qzRange[0] < q_q) & (q_q < clas.qzRange[1]), True, False)
                q_q = q_q[filter_q]
                R_q = R_q[filter_q]
                dR_q = dR_q[filter_q]
                dq_q = dq_q[filter_q]

                if libeos.reduction.autoscale:
                    if i == 0:
                        R_q, dR_q = autoscale(q_q, R_q, dR_q, clas.autoscale)
                    else: 
                        pRq_z     = datasetsRqz[i-1].data[:,1]
                        pdRq_z    = datasetsRqz[i-1].data[:,2]
                        R_q, dR_q = autoscale(q_q, R_q, dR_q, clas.autoscale, pRq_z, pdRq_z)

                if subtract:
                    if len(q_q) == len(sq_q):
                        R_q  -= sR_q
                        dR_q = np.sqrt( dR_q**2 + sdR_q**2 )
                    else:
                        logging.warning(f'backgroung file {sFileName} not compatible with q_z scale ({len(sq_q)} vs. {len(q_q)})')

                data = np.array([q_q, R_q, dR_q, dq_q]).T
                orso_data       = fileio.OrsoDataset(headerRqz, data)
                headerRqz       = fileio.Orso(**headerRqz.to_dict())
                datasetsRqz.append(orso_data)
   
            if 'Rlt.ort' in output_format_list(clas.outputFormat):
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
                #data_source = fromHDF.data_source
                headerRlt = fileio.Orso(header.data_source, header.reduction, columns)
        
                ts, zs=ref_lz.shape
                lindex_lz=np.tile(np.arange(1, ts+1), (zs, 1)).T
                tindex_lz=np.tile(np.arange(1, zs+1), (ts, 1))
  
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
                    data               = np.array(list(item)).T
                    headerRlt          = fileio.Orso(**headerRlt.to_dict())
                    headerRlt.data_set = f'dataset_{i}_{j+1} : alpha_f = {theta_lz[0,j]:6.3f} deg'
                    orso_data          = fileio.OrsoDataset(headerRlt, data)
                    datasetsRlt.append(orso_data)
                    j += 1

    # output
    logging.warning('output:')

    if 'Rqz.ort' in output_format_list(clas.outputFormat):
        logging.warning(f'  {clas.dataPath}/{clas.outputName}.Rqz.ort')
        theSecondLine = f' {header.experiment.title} | {header.experiment.start_date} | sample {header.sample.name} | R(q_z)'
        fileio.save_orso(datasetsRqz, f'{clas.dataPath}/{clas.outputName}.Rqz.ort', data_separator='\n', comment=theSecondLine)

    if 'Rlt.ort' in output_format_list(clas.outputFormat):
        logging.warning(f'  {clas.dataPath}/{clas.outputName}.Rlt.ort')
        theSecondLine = f' {header.experiment.title} | {header.experiment.start_date} | sample {header.sample.name} | R(lambda, theta)'
        fileio.save_orso(datasetsRlt, f'{clas.dataPath}/{clas.outputName}.Rlt.ort', data_separator='\n', comment=theSecondLine)

if __name__ == '__main__':
    main()
