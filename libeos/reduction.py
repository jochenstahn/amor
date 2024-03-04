import logging
import os

import numpy as np

from libeos.command_line import expand_file_list
from libeos.dataset import AmorData


def normalisation_map(short_notation, header, grid, dataPath):
    fromHDF = AmorData()
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
        header.measurement_additional_files = normFileList
    else:
        logging.info(f'# normalisation matrix: using the files {normalisation_list}')
        fromHDF.read_data(short_notation, norm=True)
        normAngle     = fromHDF.nu - fromHDF.mu
        lamda_e       = fromHDF.lamda_e
        detZ_e        = fromHDF.detZ_e
        norm_lz, bins_l, bins_z = np.histogram2d(lamda_e, detZ_e, bins = (grid.lamda(), grid.z()))
        norm_lz = np.where(norm_lz>0, norm_lz, np.nan)
        # correct for the SM reflectivity
        lamda_l  = grid.lamda()
        theta_z  = normAngle + fromHDF.delta_z
        lamda_lz = (grid.lz().T*lamda_l[:-1]).T
        theta_lz = grid.lz()*theta_z
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


def project_on_lz(fromHDF, norm_lz, normAngle, lamda_e, detZ_e, grid, thetaRange, thetaRangeR, lambdaRange):
    # projection on lambda-z-grid
    lamda_l  = grid.lamda()
    theta_z  = fromHDF.nu - fromHDF.mu + fromHDF.delta_z
    lamda_lz = (grid.lz().T*lamda_l[:-1]).T
    theta_lz = grid.lz()*theta_z

    thetaN_z  = fromHDF.delta_z + normAngle
    thetaN_lz = np.ones(np.shape(norm_lz))*thetaN_z
    thetaN_lz = np.where(np.absolute(thetaN_lz)>5e-3, thetaN_lz, np.nan)

    mask_lz   = np.where(np.isnan(norm_lz), False, True)
    mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(thetaN_lz)>5e-3, True, False))
    mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(theta_lz)>5e-3, True, False))
    if thetaRange[1]<12:
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz >= thetaRange[0], True, False))
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz <= thetaRange[1], True, False))
    if thetaRangeR[1]<12:
      t0 = fromHDF.nu - fromHDF.mu
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz-t0 >= thetaRangeR[0], True, False))
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz-t0 <= thetaRangeR[1], True, False))
    if lambdaRange[1]<15:
      mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz >= lambdaRange[0], True, False))
      mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz <= lambdaRange[1], True, False))

    #           gravity correction
    #theta_lz += np.rad2deg( np.arctan( 3.07e-10 * (fromHDF.detectorDistance + detXdist_e) * lamda_lz**2 ) )
    theta_lz += np.rad2deg( np.arctan( 3.07e-10 * fromHDF.detectorDistance * lamda_lz**2 ) )

    z_z       = enumerate(theta_z)
    qz_lz     = 4.0*np.pi * np.sin(np.deg2rad(theta_lz)) / lamda_lz
    int_lz, bins_l, bins_z  = np.histogram2d(lamda_e, detZ_e, bins = (lamda_l, grid.z()))
    #           cut normalisation sample horizon
    int_lz    = np.where(mask_lz, int_lz, np.nan)
    thetaF_lz  = np.where(mask_lz, theta_lz, np.nan)

    ref_lz    = (int_lz * np.absolute(thetaN_lz)) / (norm_lz * np.absolute(thetaF_lz))
    err_lz    = ref_lz * np.sqrt( 1/(int_lz+.1) + 1/norm_lz )

    res_lz    = np.ones((np.shape(lamda_l[:-1])[0], np.shape(theta_z)[0])) * 0.022**2
    res_lz    = res_lz + (0.008/theta_lz)**2
    res_lz    = qz_lz * np.sqrt(res_lz)

    return qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz


def project_on_qz(q_lz, R_lz, dR_lz, dq_lz, norm_lz, mask_lz, grid):
    q_q       = grid.q()
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


def autoscale(q_q, R_q, dR_q, autoscale, pR_q=[], pdR_q=[]):
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
