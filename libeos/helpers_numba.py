import numba as nb
import numpy as np

@nb.jit(nb.float64[:](nb.float64[:], nb.float64, nb.float64, nb.float64),
        nopython=True, parallel=True, cache=True)
def merge_frames(tof_e, tofCut, tau, total_offset):
    # fast implementation of the merging function used in file_reader.AmorData.merge_frames
    tof_e_out = np.empty(tof_e.shape, dtype=np.float64)
    dt = (tofCut-tau)
    for ti in nb.prange(tof_e.shape[0]):
        tof_e_out[ti] = ((tof_e[ti]-dt)%tau)+total_offset  # tof shifted to 1 frame
    return tof_e_out

@nb.jit(nb.int64[:](nb.float64[:], nb.uint32[:], nb.int64[:]),
        nopython=True, parallel=True, cache=True)
def extract_walltime(tof_e, dataPacket_p, dataPacketTime_p):
    # assigning every event the wall time of the event packet (absolute time of pulse ?start?)
    totalNumber = np.shape(tof_e)[0]
    wallTime_e = np.empty(totalNumber, dtype=np.int64)
    for i in nb.prange(len(dataPacket_p)-1):
        for j in range(dataPacket_p[i], dataPacket_p[i+1]):
            wallTime_e[j] = dataPacketTime_p[i]
    for j in range(dataPacket_p[-1], totalNumber):
        wallTime_e[j] = dataPacketTime_p[-1]
    return wallTime_e

@nb.jit(nb.types.Tuple((nb.int64[:], nb.float64[:], nb.float64[:], nb.boolean[:]))
                (nb.float64[:, :], nb.uint32[:], nb.int64, nb.int64),
        nopython=True, parallel=True, cache=True)
def filter_project_x(pixelLookUp, pixelID_e, ymin, ymax):
    # project events on z-axis and create filter for events outside of y-range
    events = pixelID_e.shape[0]
    detY_e = np.empty(events, dtype=np.int64)
    detZ_e = np.empty(events, dtype=np.int64)
    detXdist_e = np.empty(events, dtype=np.float64)
    delta_e = np.empty(events, dtype=np.float64)
    mask_e = np.empty(events, dtype=nb.boolean)
    for i in nb.prange(events):
        # resolve pixel ID into y and z indicees, x position and angle
        detY_e[i] = pixelLookUp[pixelID_e[i]-1, 0]
        detZ_e[i] = pixelLookUp[pixelID_e[i]-1, 1]
        detXdist_e[i] = pixelLookUp[pixelID_e[i]-1, 2]
        delta_e[i] = pixelLookUp[pixelID_e[i]-1, 3]
        # define mask and filter y range
        mask_e[i] = (ymin<=detY_e[i]) & (detY_e[i]<=ymax)
    return (detZ_e, detXdist_e, delta_e, mask_e)

@nb.jit(nb.types.Tuple((nb.float64[:], nb.float64[:], nb.boolean[:]))
                (nb.float64[:], nb.float64[:], nb.float64[:], nb.boolean[:],
                 nb.float64, nb.float64, nb.float64, nb.float64, nb.float64, nb.float64),
        nopython=True, parallel=True, cache=True)
def calculate_derived_properties_focussing(tof_e, detXdist_e, delta_e, mask_e,
                                           lmin, lmax, nu, mu, chopperDetectorDistance, hdm):
    events = tof_e.shape[0]
    alphaF_e = np.empty(events, dtype=np.float64)
    lamda_e = np.empty(events, dtype=np.float64)
    qz_e = np.empty(events, dtype=np.float64)
    mask_e_out = np.empty(events, dtype=nb.boolean)
    denom_f1 = 1.e13*hdm
    for i in nb.prange(events):
        lamda_e[i] = denom_f1*tof_e[i]/(chopperDetectorDistance+detXdist_e[i])
        mask_e_out[i] = mask_e[i] & ((lmin<=lamda_e[i]) & (lamda_e[i]<=lmax))
        alphaF_e[i] = nu-mu+delta_e[i]
        qz_e[i] = 4*np.pi*np.sin(np.deg2rad(alphaF_e[i]))/lamda_e[i]
    return lamda_e, qz_e, mask_e
