"""
Classes used to calculate projections/binnings from event data onto given grids.
"""

import logging
from typing import Protocol

import numpy as np
from dataclasses import dataclass

from .event_data_types import EventDatasetProtocol
from .instrument import Detector, LZGrid
from .normalization import LZNormalisation

@dataclass
class ProjectedReflectivity:
    R: np.ndarray
    dR: np.ndarray
    Q: np.ndarray
    dQ: np.ndarray

    @property
    def data(self):
        """
        Return combined data compatible with storing as columns in orso file.
            Q, R, dR, dQ
        """
        return np.array([self.Q, self.R, self.dR, self.dQ]).T

    def data_for_time(self, time):
        tme = np.ones(np.shape(self.Q))*time
        return np.array([self.Q, self.R, self.dR, self.dQ, tme]).T

    def scale(self, factor):
        self.R *= factor
        self.dR *= factor

    def autoscale(self, range):
        filter_q = (range[0]<=self.Q) & (self.Q<=range[1])
        filter_q &= self.dR>0
        if filter_q.sum()>0:
            scale = (self.R[filter_q]/self.dR[filter_q]).sum()/(self.R[filter_q]**2/self.dR[filter_q]).sum()
            self.scale(scale)
            logging.info(f'      scaling factor = {scale}')
            return scale
        else:
            logging.warning('      automatic scaling not possible')
            return 1.0

    def stitch(self, other: 'ProjectedReflectivity'):
        # find scaling factor between two reflectivities at points both are not zero
        filter_q = np.logical_not(np.isnan(other.R*self.R))
        filter_q &= self.R>0
        filter_q &= other.R>0
        R1 = self.R[filter_q]
        dR1 = self.dR[filter_q]
        R2 = other.R[filter_q]
        dR2 = other.dR[filter_q]
        if len(R1)>0:
            scale = (R1**2*R2**2/(dR1**2*dR2**2)).sum() / (R1**3*R2/(dR1**2*dR2**2)).sum()
            self.scale(scale)
            logging.info(f'      scaling factor = {scale}')
            return scale
        else:
            logging.warning('      automatic scaling not possible')
            return 1.0

    def subtract(self, R, dR):
        # subtract another dataset with same q-points
        self.R -= R
        self.dR = np.sqrt(self.dR**2+dR**2)

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        plt.errorbar(self.Q, self.R, xerr=self.dQ, yerr=self.dR, **kwargs)
        plt.yscale('log')
        plt.xlabel('Q / $\\AA^{-1}$')
        plt.ylabel('R')

class ProjectionInterface(Protocol):
    def project(self, dataset: EventDatasetProtocol, monitor: float): ...
    def clear(self): ...
    def plot(self, **kwargs): ...
    def update_plot(self): ...

class LZProjection:
    grid: LZGrid
    lamda: np.ndarray
    alphaF: np.ndarray
    is_normalized: bool

    data: np.recarray

    def __init__(self, tthh: float, grid: LZGrid):
        self.grid = grid
        self.is_normalized = False

        alphaF_z  = tthh + Detector.delta_z
        lamda_l  = self.grid.lamda()
        lamda_c = (lamda_l[:-1]+lamda_l[1:])/2

        lz_shape = self.grid.lz()

        self.lamda  = lz_shape*lamda_c[:, np.newaxis]
        self.alphaF = lz_shape*alphaF_z[np.newaxis, :]
        self.data = np.zeros(self.alphaF.shape, dtype=[
            ('I', np.float64),
            ('mask', bool),
            ('ref', np.float64),
            ('err', np.float64),
            ('res', np.float64),
            ('qz', np.float64),
            ('qx', np.float64),
            ('norm', np.float64),
            ]).view(np.recarray)
        self.data.mask = True
        self.monitor = 0.

    @classmethod
    def from_dataset(cls, dataset: EventDatasetProtocol, grid: LZGrid, has_offspecular=False):
        tthh  = dataset.geometry.nu - dataset.geometry.mu
        output = cls(tthh, grid)
        output.correct_gravity(dataset.geometry.detectorDistance)
        if has_offspecular:
            alphaI_lz = grid.lz()*(dataset.geometry.mu+dataset.geometry.kap+dataset.geometry.kad)
            output.calculate_q(alphaI_lz)
        else:
            output.calculate_q()
        return output

    def correct_gravity(self, detector_distance):
        self.alphaF += np.rad2deg( np.arctan( 3.07e-10 * detector_distance * self.lamda**2 ) )

    def calculate_q(self, alphaI=None):
        if alphaI is None:
            self.data.qz = 4.0*np.pi*np.sin(np.deg2rad(self.alphaF))/self.lamda
            self.data.qx = 0.*self.data.qz
        else:
            self.data.qz = 2.0*np.pi*(np.sin(np.deg2rad(self.alphaF))+np.sin(np.deg2rad(alphaI)))/self.lamda
            self.data.qx = 2.0*np.pi*(np.cos(np.deg2rad(self.alphaF))-np.cos(np.deg2rad(alphaI)))/self.lamda

        if self.data.qz[0,self.data.qz.shape[1]//2]  < 0:
            # assuming a 'measurement from below' when center of detector at negative qz
            self.data.qz *= -1

        self.calculate_q_resolution()

    def calculate_q_resolution(self):
        res_lz    = self.grid.lz() * 0.022**2
        res_lz    = res_lz + (0.008/self.alphaF)**2
        self.data.res    = self.data.qz * np.sqrt(res_lz)

    def apply_theta_filter(self, theta_range):
        # Filters points within theta range
        self.data.mask &= (self.alphaF<theta_range[0])|(self.alphaF>theta_range[1])

    def apply_theta_mask(self, theta_range):
        # Mask points outside theta range
        self.data.mask &= self.alphaF>=theta_range[0]
        self.data.mask &= self.alphaF<=theta_range[1]

    def apply_lamda_mask(self, lamda_range):
        # Mask points outside lambda range
        self.data.mask &= self.lamda>=lamda_range[0]
        self.data.mask &= self.lamda<=lamda_range[1]

    def apply_norm_mask(self, norm: LZNormalisation):
        # Mask points where normliazation is nan
        self.data.mask &= np.logical_not(np.isnan(norm.norm))

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        """
            Project dataset on grid and add to intensity.
            Can be called multiple times to sequentially add events.
        """
        # TODO: maybe move monitor calculation in here instead of reduction?
        e = dataset.data.events
        int_lz, *_  = np.histogram2d(e.lamda, e.detZ, bins = (self.grid.lamda(), self.grid.z()))
        self.data.I += int_lz
        self.monitor += monitor
        # in case the intensity changed one needs to normalize again
        self.is_normalized = False

    def clear(self):
        # empty data
        self.data[:] = 0
        self.data.mask = True

    @property
    def I(self):
        output = self.data.I[:]
        output[np.logical_not(self.data.mask)] = np.nan
        return output / self.monitor

    def calc_error(self):
        # calculate error bars for resulting intensity after normalization
        self.data.err = self.data.ref * np.sqrt( 1./(self.data.I+.1) + 1./self.data.norm )

    def normalize_over_illuminated(self, norm: LZNormalisation):
        """
        Normalize the dataaset and take into account a difference in
        detector angle for measurement and reference.
        """
        norm_lz = norm.norm
        thetaN_z = Detector.delta_z+norm.angle
        thetaN_lz = np.ones_like(norm_lz)*thetaN_z
        thetaN_lz = np.where(np.absolute(thetaN_lz)>5e-3, thetaN_lz, np.nan)
        self.data.mask &=  (np.absolute(thetaN_lz)>5e-3)
        ref_lz = (self.data.I*np.absolute(thetaN_lz))/(norm_lz*np.absolute(self.alphaF))
        ref_lz *= norm.monitor/self.monitor
        ref_lz[np.logical_not(self.data.mask)] = np.nan
        self.data.norm = norm_lz
        self.data.ref = ref_lz
        self.calc_error()
        self.is_normalized = True

    def normalize_no_footprint(self, norm: LZNormalisation):
        norm_lz = norm.norm
        ref_lz = (self.data.I/norm_lz)
        ref_lz *= norm.monitor/self.monitor
        ref_lz[np.logical_not(self.data.mask)] = np.nan
        self.data.norm = norm_lz
        self.data.ref = ref_lz
        self.calc_error()
        self.is_normalized = True

    def scale(self, factor: float):
        if not self.is_normalized:
            raise ValueError("Dataset needs to be normalized, first")
        self.data.ref *= factor
        self.data.err *= factor

    def project_on_qz(self):
        if not self.is_normalized:
            raise ValueError("Dataset needs to be normalized, first")
        q_q       = self.grid.q()
        weights_lzf = self.data.norm[self.data.mask]
        q_lzf = self.data.qz[self.data.mask]
        R_lzf = self.data.ref[self.data.mask]
        dR_lzf = self.data.err[self.data.mask]
        dq_lzf = self.data.res[self.data.mask]

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

        return ProjectedReflectivity(R_q, dR_q, (q_q[1:]+q_q[:-1])/2., dq_q)

    ############## potential speedup not used right now, needs to be tested ####################
    @classmethod
    def histogram2d_lz(cls, lamda_e, detZ_e, bins):
        """
        Perform binning operation equivalent to numpy bin2d for the sepcial case
        of the second dimension using integer positions (pre-defined pixels).
        Based on the devide_bin algorithm below.
        """
        dimension = bins[1].shape[0]-1
        if not (np.array(bins[1])==np.arange(dimension+1)).all():
            raise ValueError("histogram2d_lz requires second bin dimension to be contigous integer range")
        binning = cls.devide_bin(lamda_e, detZ_e.astype(np.int64), bins[0], dimension)
        return np.array(binning), bins[0], bins[1]

    @classmethod
    def devide_bin(cls, lambda_e, position_e, lamda_edges, dimension):
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
        left_list = cls.devide_bin(lambda_e[left_region], position_e[left_region],
                                             lamda_edges[:split_idx+1], dimension)
        right_region = np.logical_not(left_region)
        right_list = cls.devide_bin(lambda_e[right_region], position_e[right_region],
                                              lamda_edges[split_idx:], dimension)
        return left_list+right_list

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        if 'colorbar' in kwargs:
            cmap=True
            del(kwargs['colorbar'])
        else:
            cmap=False

        if self.is_normalized:
            if not 'norm' in kwargs:
                kwargs['norm'] = LogNorm(2e-3, 2.0)
            self._graph = plt.pcolormesh(self.lamda, self.alphaF, self.data.ref, **kwargs)
            if cmap:
                plt.colorbar(label='R')
        else:
            if not 'norm' in kwargs:
                kwargs['norm'] = LogNorm()
            self._graph = plt.pcolormesh(self.lamda, self.alphaF, self.data.I, **kwargs)
            if cmap:
                plt.colorbar(label='I / cpm')
        plt.xlabel('$\\lambda$ / $\\AA$')
        plt.ylabel('$\\Theta$ / °')

    def update_plot(self):
        """
        Inline update of last plot by just updating the data.
        """
        if self.is_normalized:
            self._graph.set_array(self.data.ref)
        else:
            self._graph.set_array(self.data.I)
