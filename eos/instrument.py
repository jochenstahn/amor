"""
Classes describing the AMOR instrument configuration used during reduction.
"""

import logging
import numpy as np

from . import const

try:
    from functools import cache
except ImportError:
    # python <3.9
    def cache(func): return func


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

    delta_z: np.ndarray
    pixelLookUp: np.ndarray

    @staticmethod
    def resolve_pixels():
        """
        Determine spatial coordinats and angles from pixel number,
        does only have to be computed once for the detector
        """
        if hasattr(Detector, 'pixelLookUp'):
            return
        nPixel = Detector.nWires * Detector.nStripes * Detector.nBlades
        pixelID = np.arange(nPixel)
        (bladeNr, bPixel) = np.divmod(pixelID, Detector.nWires * Detector.nStripes)
        (bZi, detYi)      = np.divmod(bPixel, Detector.nStripes)                     # z index on blade, y index on detector
        detZi             = bladeNr * Detector.nWires + bZi                          # z index on detector
        detX              = bZi * Detector.dX                                        # x position in detector
        # detZ              = Detector.zero - bladeNr * Detector.bladeZ - bZi * Detector.dZ      # z position on detector
        bladeAngle        = np.rad2deg( 2. * np.arcsin(0.5*Detector.bladeZ / Detector.distance) )
        delta             = (Detector.nBlades/2. - bladeNr) * bladeAngle \
                            - np.rad2deg( np.arctan(bZi*Detector.dZ / ( Detector.distance + bZi * Detector.dX) ) )
        delta_z      = delta[detYi==1]
        pixel_lookup=np.vstack((detYi.T, detZi.T, detX.T, delta.T)).T
        Detector.delta_z = delta_z
        Detector.pixelLookUp = pixel_lookup

# guarantee that pixelLookUp has been computed
Detector.resolve_pixels()

class LZGrid:
    dldl = 0.005  # Delta lambda / lambda

    # as using cahced results, make sure the object is not modified
    @property
    def qResolution(self):
        return self._qResolution
    @property
    def qzRange(self):
        return self._qzRange

    def __init__(self, qResolution, qzRange, lambda_overwrite=None):
        self._qResolution = qResolution
        self._qzRange = qzRange
        if lambda_overwrite is None:
            self.lamdaMax = const.lamdaMax
            self.lamdaCut = const.lamdaCut
        else:
            self.lamdaCut, self.lamdaMax = lambda_overwrite

    @property
    @cache
    def shape(self):
        # gives the shape of the grid, not of the bin-edges
        return (self.lamda().shape[0]-1, self.z().shape[0]-1)

    @cache
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

    @cache
    def lamda(self):
        lamdaMax = self.lamdaMax
        lamdaMin = self.lamdaCut
        lamda_grid = lamdaMin*(1+self.dldl)**np.arange(int(np.log(lamdaMax/lamdaMin)/np.log(1+self.dldl)+1))
        return lamda_grid

    @cache
    def z(self):
        # TODO: shouldn't this be -0.5 to be the edges of each pixel?
        return np.arange(Detector.nBlades*Detector.nWires+1)

    @cache
    def lz(self):
        return np.ones(( self.lamda().shape[0]-1, self.z().shape[0]-1))

    @cache
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
