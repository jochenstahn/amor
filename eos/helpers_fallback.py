"""
Equivalent function as in helpers_numba.py but using just numpy functionality.
"""

import numpy as np

def merge_frames(tof_e, tofCut, tau, total_offset):
    # tof shifted to 1 frame
    return np.remainder(tof_e-(tofCut-tau), tau)+total_offset

def extract_walltime(tof_e, dataPacket_p, dataPacketTime_p):
    output = np.empty(np.shape(tof_e)[0], dtype=np.int64)
    for i in range(len(dataPacket_p)-1):
        output[dataPacket_p[i]:dataPacket_p[i+1]] = dataPacketTime_p[i]
    output[dataPacket_p[-1]:] = dataPacketTime_p[-1]
    return output

def filter_project_x(pixelLookUp, pixelID_e, ymin, ymax):
    (detY_e, detZ_e, detXdist_e, delta_e) = pixelLookUp[np.int_(pixelID_e)-1, :].T
    # define mask and filter y range
    mask_e = (ymin<=detY_e) & (detY_e<=ymax)
    return (detZ_e, detXdist_e, delta_e, mask_e)

def calculate_derived_properties_focussing(tof_e, detXdist_e, delta_e, mask_e,
                                           lmin, lmax, nu, mu, chopperDetectorDistance, hdm):
    raise NotImplementedError("Only exists in numba implementation so far.")
