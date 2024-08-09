"""
Classes describing the AMOR instrument configuration used during reduction.
"""

import logging
import numpy as np

from . import const

class Detector:
    nBlades  = 14  # number of active blades in the detector
    nWires   = 32  # number of wires per blade
    nStripes = 64  # number of stipes per blade
    angle    = np.deg2rad(5.1)  # deg  angle of incidence of the beam on the blades (def: 5.1)
    dZ       = 4.0*np.sin(angle)  # mm  height-distance of neighboring pixels on one blade
    dX       = 4.0*np.cos(angle)  # mm  depth-distance of neighboring pixels on one blace
    bladeZ   = 10.455  # mm  distance between detector blades
    zero     = 0.5*nBlades*bladeZ  # mm  vertical center of the detector
    distance = 4000.  # mm  distance from focal point to leading blade edge

class Grid:

    def __init__(self, qResolution, qzRange):
        self.lamdaCut = const.lamdaCut
        self.dldl = 0.005     # Delta lambda / lambda
        self.qResolution = qResolution
        self.qzRange = qzRange

    def q(self):
        resolutions = [0.005, 0.01, 0.02, 0.025, 0.04, 0.05, 0.1, 1]
        a, b = np.histogram([self.qResolution], bins = resolutions)
        dqdq = np.matmul(b[:-1],a)
        if dqdq != self.qResolution:
            logging.info(f'#   changed resolution to {dqdq}')
        qq = 0.01
        # linear up to qq
        q_grid = np.arange(0, qq, qq*dqdq)
        # exponential from qq on
        q_grid = np.append(q_grid, qq*(1.+dqdq)**np.arange(int(np.log(self.qzRange[1]/qq)/np.log(1+dqdq))))
        q_grid = q_grid[q_grid>=self.qzRange[0]]
        return q_grid

    def lamda(self):
        lamdaMax = 16
        lamdaMin = self.lamdaCut
        lamda_grid = lamdaMin*(1+self.dldl)**np.arange(int(np.log(lamdaMax/lamdaMin)/np.log(1+self.dldl)+1))
        return lamda_grid

    def z(self):
        return np.arange(Detector.nBlades*Detector.nWires+1)

    def lz(self):
        return np.ones(( np.shape(self.lamda()[:-1])[0], np.shape(self.z()[:-1])[0] ))

    def delta(self, detectorDistance):
        # unused for now
        bladeAngle = np.rad2deg( 2. * np.arcsin(0.5*Detector.bladeZ / detectorDistance) )
        blade_grid = np.arctan( np.arange(33) * Detector.dZ / ( detectorDistance + np.arange(33) * Detector.dX) )
        blade_grid = np.rad2deg(blade_grid)
        stepWidth  = blade_grid[1] - blade_grid[0]
        blade_grid = blade_grid - 0.2 * stepWidth

        delta_grid = []
        for b in np.arange(Detector.nBlades-1):
            delta_grid = np.concatenate((delta_grid, blade_grid), axis=None)
            blade_grid = blade_grid + bladeAngle
            delta_grid = delta_grid[delta_grid<blade_grid[0]-0.5*stepWidth]
        delta_grid = np.concatenate((delta_grid, blade_grid), axis=None)

        return -np.flip(delta_grid) + 0.5*Detector.nBlades * bladeAngle
