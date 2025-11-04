"""
Classes used to calculate projections/binnings from event data onto given grids.
"""

import logging
import warnings
from abc import ABC, abstractmethod
from typing import List, Tuple, Union

import numpy as np
from dataclasses import dataclass

from .event_data_types import EventDatasetProtocol
from .instrument import Detector, LZGrid
from .normalization import LZNormalisation

class ProjectionInterface(ABC):
    @abstractmethod
    def project(self, dataset: EventDatasetProtocol, monitor: float): ...

    @abstractmethod
    def clear(self): ...

    @abstractmethod
    def plot(self, **kwargs): ...

    @abstractmethod
    def update_plot(self): ...

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


class LZProjection(ProjectionInterface):
    grid: LZGrid
    lamda: np.ndarray
    alphaF: np.ndarray
    is_normalized: bool

    data: np.recarray
    _dtype = np.dtype([
            ('I', np.float64),
            ('mask', bool),
            ('ref', np.float64),
            ('err', np.float64),
            ('res', np.float64),
            ('qz', np.float64),
            ('qx', np.float64),
            ('norm', np.float64),
            ])

    def __init__(self, tthh: float, grid: LZGrid):
        self.grid = grid
        self.is_normalized = False

        alphaF_z  = tthh + Detector.delta_z
        lamda_l  = self.grid.lamda()
        lamda_c = (lamda_l[:-1]+lamda_l[1:])/2

        lz_shape = self.grid.lz()

        self.lamda  = lz_shape*lamda_c[:, np.newaxis]
        self.alphaF = lz_shape*alphaF_z[np.newaxis, :]
        self.data = np.zeros(self.alphaF.shape, dtype=self._dtype).view(np.recarray)
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
        self.monitor = 0.

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

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        if 'colorbar' in kwargs:
            cmap=True
            del(kwargs['colorbar'])
        else:
            cmap=False

        if self.is_normalized:
            I = self.data.ref
        else:
            I = self.data.I


        if not 'norm' in kwargs:
            vmin = I[(I>0)].min()
            vmax = np.nanmax(I)
            kwargs['norm'] = LogNorm(vmin, vmax, clip=True)


        # suppress warning for wrongly sorted y-axis pixels (blades overlap)
        with warnings.catch_warnings(action='ignore', category=UserWarning):
            self._graph = plt.pcolormesh(self.lamda, self.alphaF, I, **kwargs)
        if cmap:
            if self.is_normalized:
                plt.colorbar(label='R')
            else:
                plt.colorbar(label='I / cpm')
        plt.xlabel('$\\lambda$ / $\\AA$')
        plt.ylabel('$\\Theta$ / °')
        plt.xlim(self.lamda[0,0], self.lamda[-1,0])
        af = self.alphaF[self.data.mask]
        plt.ylim(af.min(), af.max())
        plt.title('Wavelength vs. Reflection Angle')

        self._graph_axis = plt.gca()
        plt.connect('button_press_event', self.draw_qline)

    def update_plot(self):
        """
        Inline update of previous plot by just updating the data.
        """
        from matplotlib.colors import LogNorm
        if self.is_normalized:
            I = self.data.ref
        else:
            I = self.data.I

        if isinstance(self._graph.norm, LogNorm):
            vmin = I[(I>0)].min()*0.5
        else:
            vmin = 0
        vmax = np.nanmax(I)
        self._graph.set_array(I)
        self._graph.norm.vmin = vmin
        self._graph.norm.vmax = vmax

        if self.is_normalized:
            self._graph.set_array(self.data.ref)
        else:
            self._graph.set_array(self.data.I)

    def draw_qline(self, event):
        if event.inaxes is not self._graph_axis:
            return
        from matplotlib import pyplot as plt
        tbm = self._graph_axis.figure.canvas.manager.toolbar.mode
        if event.button is plt.MouseButton.LEFT and tbm=='':
            slope = event.ydata/event.xdata
            xmax = 12.5
            self._graph_axis.plot([0, xmax], [0, slope*xmax], '-', color='grey')
            self._graph_axis.text(event.xdata, event.ydata, f'q={np.deg2rad(slope)*4.*np.pi:.3f}', backgroundcolor='white')
            plt.draw()
        if event.button is plt.MouseButton.RIGHT and tbm=='':
            for art in list(self._graph_axis.lines)+list(self._graph_axis.texts):
                art.remove()
            plt.draw()

ONLY_MAP = ['colorbar', 'cmap', 'norm']

class ReflectivityProjector(ProjectionInterface):
    lzprojection: LZProjection
    data: ProjectedReflectivity
    # TODO: maybe implement direct 1d projection in here

    def __init__(self, lzprojection, norm):
        self.lzprojection = lzprojection
        self.norm = norm

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        self.lzprojection.project(dataset, monitor)
        self.lzprojection.normalize_over_illuminated(self.norm)
        self.data = self.lzprojection.project_on_qz()

    def clear(self):
        self.lzprojection.clear()

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        for key in ONLY_MAP:
            if key in kwargs: del(kwargs[key])

        self._graph = plt.errorbar(self.data.Q, self.data.R, xerr=self.data.dQ, yerr=self.data.dR, **kwargs)
        self._graph_axis = plt.gca()
        plt.title('Reflectivity (might be improperly normalized)')
        plt.yscale('log')
        plt.xlabel('Q / $\\AA^{-1}$')
        plt.ylabel('R')

    def update_plot(self):
        ln, _, (barsx, barsy) = self._graph

        yerr_top = self.data.R+self.data.dR
        yerr_bot = self.data.R-self.data.dR
        xerr_top = self.data.Q+self.data.dQ
        xerr_bot = self.data.Q-self.data.dQ

        new_segments_x = [np.array([[xt, y], [xb, y]]) for xt, xb, y in zip(xerr_top, xerr_bot, self.data.R)]
        new_segments_y = [np.array([[x, yt], [x, yb]]) for x, yt, yb in zip(self.data.Q, yerr_top, yerr_bot)]
        barsx.set_segments(new_segments_x)
        barsy.set_segments(new_segments_y)

        ln.set_ydata(self.data.R)


class YZProjection(ProjectionInterface):
    y: np.ndarray
    z: np.ndarray

    data: np.recarray
    _dtype = np.dtype([
            ('cts', np.float64),
            ('I', np.float64),
            ('err', np.float64),
            ])

    def __init__(self):
        self.z = np.arange(Detector.nBlades*Detector.nWires+1)-0.5
        self.y = np.arange(Detector.nStripes+1)-0.5
        self.data = np.zeros((self.y.shape[0]-1, self.z.shape[0]-1), dtype=self._dtype).view(np.recarray)
        self.monitor = 0.

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        detYi, detZi, detX, delta = Detector.pixelLookUp[dataset.data.events.pixelID-1].T

        cts , *_ = np.histogram2d(detYi, detZi, bins=(self.y, self.z))
        self.data.cts += cts
        self.monitor += monitor

        self.data.I = self.data.cts / self.monitor
        self.data.err = np.sqrt(self.data.cts) / self.monitor

    def clear(self):
        self.data[:] = 0
        self.monitor = 0.

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        if 'colorbar' in kwargs:
            cmap=True
            del(kwargs['colorbar'])
        else:
            cmap=False

        vmax = self.data.I.max()

        if not 'norm' in kwargs:
            vmin = self.data.I[(self.data.I>0)].min()*0.5
            kwargs['norm'] = LogNorm(vmin, vmax)

        self._graph = plt.pcolormesh(self.y, self.z, self.data.I.T, **kwargs)
        if cmap:
            plt.colorbar(label='I / cpm')

        plt.xlabel('Y')
        plt.ylabel('Z')
        plt.xlim(self.y[0], self.y[-1])
        plt.ylim(self.z[-1], self.z[0])
        plt.title('Horizontal Pixel vs. Vertical Pixel')

        self._graph_axis = plt.gca()
        plt.connect('button_press_event', self.draw_yzcross)

    def update_plot(self):
        """
        Inline update of previous plot by just updating the data.
        """
        from matplotlib.colors import LogNorm
        if isinstance(self._graph.norm, LogNorm):
            vmin = self.data.I[(self.data.I>0)].min()*0.5
        else:
            vmin = 0
        vmax = self.data.I.max()
        self._graph.set_array(self.data.I.T)
        self._graph.norm.vmin = vmin
        self._graph.norm.vmax = vmax

    def draw_yzcross(self, event):
        if event.inaxes is not self._graph_axis:
            return
        from matplotlib import pyplot as plt
        tbm = self._graph_axis.figure.canvas.manager.toolbar.mode
        if event.button is plt.MouseButton.LEFT and tbm=='':
            self._graph_axis.plot([event.xdata, event.xdata], [self.z[0], self.z[-1]], '-', color='grey')
            self._graph_axis.plot([self.y[0], self.y[-1]], [event.ydata, event.ydata], '-', color='grey')
            self._graph_axis.text(event.xdata, event.ydata, f'({event.xdata:.1f}, {event.ydata:.1f})', backgroundcolor='white')
            plt.draw()
        if event.button is plt.MouseButton.RIGHT and tbm=='':
            for art in list(self._graph_axis.lines)+list(self._graph_axis.texts):
                art.remove()
            plt.draw()

class YTProjection(YZProjection):
    theta: np.ndarray

    def __init__(self, tthh: float):
        dd = Detector.delta_z[1]-Detector.delta_z[0]
        delta = np.hstack([Detector.delta_z, Detector.delta_z[-1]+dd])-dd/2.
        self.theta  = tthh + delta
        super().__init__()

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        if 'colorbar' in kwargs:
            cmap=True
            del(kwargs['colorbar'])
        else:
            cmap=False

        if not 'norm' in kwargs:
            kwargs['norm'] = LogNorm()

        self._graph = plt.pcolormesh(self.y, self.theta, self.data.I.T, **kwargs)
        if cmap:
            plt.colorbar(label='I / cpm')

        plt.xlabel('Y')
        plt.ylabel('Theta / °')
        plt.xlim(self.y[0], self.y[-1])
        plt.ylim(self.theta[-1], self.theta[0])
        plt.title('Horizontal Pixel vs. Angle')

        self._graph_axis = plt.gca()
        plt.connect('button_press_event', self.draw_tzcross)

    def draw_tzcross(self, event):
        if event.inaxes is not self._graph_axis:
            return
        from matplotlib import pyplot as plt
        tbm = self._graph_axis.figure.canvas.manager.toolbar.mode
        if event.button is plt.MouseButton.LEFT and tbm=='':
            self._graph_axis.plot([event.xdata, event.xdata], [self.theta[0], self.theta[-1]], '-', color='grey')
            self._graph_axis.plot([self.y[0], self.y[-1]], [event.ydata, event.ydata], '-', color='grey')
            self._graph_axis.text(event.xdata, event.ydata, f'({event.xdata:.1f}, {event.ydata:.1f})', backgroundcolor='white')
            plt.draw()
        if event.button is plt.MouseButton.RIGHT and tbm=='':
            for art in list(self._graph_axis.lines)+list(self._graph_axis.texts):
                art.remove()
            plt.draw()


class TofZProjection(ProjectionInterface):
    tof: np.ndarray
    z: np.ndarray

    data: np.recarray
    _dtype = np.dtype([
            ('cts', np.float64),
            ('I', np.float64),
            ('err', np.float64),
            ])

    def __init__(self, tau, foldback=False, combine=1):
        self.z = np.arange(Detector.nBlades*Detector.nWires+1)-0.5
        if foldback:
            self.tof = np.arange(0, tau, 0.0005*combine)
        else:
            self.tof = np.arange(0, 2*tau, 0.0005*combine)
        self.data = np.zeros((self.tof.shape[0]-1, self.z.shape[0]-1), dtype=self._dtype).view(np.recarray)
        self.monitor = 0.

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        detYi, detZi, detX, delta = Detector.pixelLookUp[dataset.data.events.pixelID-1].T

        cts , *_ = np.histogram2d(dataset.data.events.tof, detZi, bins=(self.tof, self.z))
        self.data.cts += cts
        self.monitor += monitor

        self.data.I = self.data.cts / self.monitor
        self.data.err = np.sqrt(self.data.cts) / self.monitor

    def clear(self):
        self.data[:] = 0
        self.monitor = 0.

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        from matplotlib.colors import LogNorm

        if 'colorbar' in kwargs:
            cmap=True
            del(kwargs['colorbar'])
        else:
            cmap=False

        if not 'norm' in kwargs:
            kwargs['norm'] = LogNorm()

        self._graph = plt.pcolormesh(self.tof*1e3, self.z, self.data.I.T, **kwargs)
        if cmap:
            plt.colorbar(label='I / cpm')

        plt.xlabel('Time of Flight / ms')
        plt.ylabel('Z')
        plt.xlim(self.tof[0]*1e3, self.tof[-1]*1e3)
        plt.ylim(self.z[-1], self.z[0])
        plt.title('Time of Flight vs. Vertical Pixel')

        self._graph_axis = plt.gca()
        plt.connect('button_press_event', self.draw_tzcross)

    def update_plot(self):
        """
        Inline update of previous plot by just updating the data.
        """
        from matplotlib.colors import LogNorm
        if isinstance(self._graph.norm, LogNorm):
            vmin = self.data.I[(self.data.I>0)].min()*0.5
        else:
            vmin = 0
        vmax = self.data.I.max()
        self._graph.set_array(self.data.I.T)
        self._graph.norm.vmin = vmin
        self._graph.norm.vmax = vmax

    def draw_tzcross(self, event):
        if event.inaxes is not self._graph_axis:
            return
        from matplotlib import pyplot as plt
        tbm = self._graph_axis.figure.canvas.manager.toolbar.mode
        if event.button is plt.MouseButton.LEFT and tbm=='':
            self._graph_axis.plot([event.xdata, event.xdata], [self.z[0], self.z[-1]], '-', color='grey')
            self._graph_axis.plot([self.tof[0]*1e3, self.tof[-1]*1e3], [event.ydata, event.ydata], '-', color='grey')
            self._graph_axis.text(event.xdata, event.ydata, f'({event.xdata:.2f}, {event.ydata:.1f})', backgroundcolor='white')
            plt.draw()
        if event.button is plt.MouseButton.RIGHT and tbm=='':
            for art in list(self._graph_axis.lines)+list(self._graph_axis.texts):
                art.remove()
            plt.draw()

class TofProjection(ProjectionInterface):
    tof: np.ndarray

    data: np.recarray
    _dtype = np.dtype([
            ('cts', np.float64),
            ('I', np.float64),
            ('err', np.float64),
            ])

    def __init__(self, tau, foldback=False):
        if foldback:
            self.tof = np.arange(0, tau, 0.0005)
        else:
            self.tof = np.arange(0, 2*tau, 0.0005)
        self.data = np.zeros(self.tof.shape[0]-1, dtype=self._dtype).view(np.recarray)
        self.monitor = 0.

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        cts , *_ = np.histogram(dataset.data.events.tof, bins=self.tof)
        self.data.cts += cts
        self.monitor += monitor

        self.data.I = self.data.cts / self.monitor
        self.data.err = np.sqrt(self.data.cts) / self.monitor

    def clear(self):
        self.data[:] = 0
        self.monitor = 0.

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        for key in ONLY_MAP:
            if key in kwargs: del(kwargs[key])


        self._graph = plt.plot(self.tof[:-1]*1e3, self.data.I, **kwargs)

        plt.xlabel('Time of Flight / ms')
        plt.ylabel('I / cpm')
        plt.xlim(self.tof[0]*1e3, self.tof[-1]*1e3)
        plt.title('Time of Flight')

    def update_plot(self):
        """
        Inline update of previous plot by just updating the data.
        """
        self._graph[0].set_ydata(self.data.I.T)

class LProjection(ProjectionInterface):
    lamda: np.ndarray

    data: np.recarray
    _dtype = np.dtype([
            ('cts', np.float64),
            ('I', np.float64),
            ('err', np.float64),
            ])

    def __init__(self):
        self.lamda = np.linspace(3.0, 12.0, 91)
        self.data = np.zeros(self.lamda.shape[0]-1, dtype=self._dtype).view(np.recarray)
        self.monitor = 0.

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        cts , *_ = np.histogram(dataset.data.events.lamda, bins=self.lamda)
        self.data.cts += cts
        self.monitor += monitor

        self.data.I = self.data.cts / self.monitor
        self.data.err = np.sqrt(self.data.cts) / self.monitor

    def clear(self):
        self.data[:] = 0
        self.monitor = 0.

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        for key in ONLY_MAP:
            if key in kwargs: del(kwargs[key])


        self._graph = plt.plot(self.lamda[:-1], self.data.I, **kwargs)

        plt.xlabel('Wavelength / Angstrom')
        plt.ylabel('I / cpm')
        plt.xlim(self.lamda[0], self.lamda[-1])
        plt.title('Wavelength')

    def update_plot(self):
        """
        Inline update of previous plot by just updating the data.
        """
        self._graph[0].set_ydata(self.data.I.T)

class TProjection(ProjectionInterface):
    theta: np.ndarray
    z: np.ndarray

    data: np.recarray
    _dtype = np.dtype([
            ('cts', np.float64),
            ('I', np.float64),
            ('err', np.float64),
            ])

    def __init__(self, tthh):
        self.z = np.arange(Detector.nBlades*Detector.nWires+1)-0.5
        dd = Detector.delta_z[1]-Detector.delta_z[0]
        delta = np.hstack([Detector.delta_z, Detector.delta_z[-1]+dd])-dd/2.
        self.theta = tthh+delta
        self.data = np.zeros(self.theta.shape[0]-1, dtype=self._dtype).view(np.recarray)
        self.monitor = 0.

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        detYi, detZi, detX, delta = Detector.pixelLookUp[dataset.data.events.pixelID-1].T

        cts , *_ = np.histogram(detZi, bins=self.z)
        self.data.cts += cts
        self.monitor += monitor

        self.data.I = self.data.cts / self.monitor
        self.data.err = np.sqrt(self.data.cts) / self.monitor

    def clear(self):
        self.data[:] = 0
        self.monitor = 0.

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        for key in ONLY_MAP:
            if key in kwargs: del(kwargs[key])


        self._graph = plt.plot(self.theta[:-1], self.data.I, **kwargs)

        plt.xlabel('Reflection Angle / °')
        plt.ylabel('I / cpm')
        plt.xlim(self.theta[-1], self.theta[0])
        plt.title('Theta')

    def update_plot(self):
        """
        Inline update of previous plot by just updating the data.
        """
        self._graph[0].set_ydata(self.data.I.T)


class CombinedProjection(ProjectionInterface):
    """
    Allows to put multiple projections together to conveniently generate combined graphs.
    """
    projections: List[ProjectionInterface]
    projection_placements: List[Union[Tuple[int, int], Tuple[int, int, int, int]]]
    grid_size: Tuple[int, int]


    def __init__(self, grid_rows, grid_cols, projections, projection_placements):
        self.projections = projections
        self.projection_placements = projection_placements
        self.grid_size = grid_rows, grid_cols

    def project(self, dataset: EventDatasetProtocol, monitor: float):
        for pi in self.projections:
            pi.project(dataset, monitor)

    def clear(self):
        for pi in self.projections:
            pi.clear()

    def plot(self, **kwargs):
        from matplotlib import pyplot as plt
        fig = plt.gcf()
        axs = fig.add_gridspec(self.grid_size[0], self.grid_size[1])
        # axs = fig.add_gridspec(self.grid_size[0]+1, self.grid_size[1],
        #                        height_ratios=[1.0 for i in range(self.grid_size[0])]+[0.2])
        self._axes = []
        for pi, placement in zip(self.projections, self.projection_placements):
            if len(placement) == 2:
                ax = fig.add_subplot(axs[placement[0], placement[1]])
            else:
                ax = fig.add_subplot(axs[placement[0]:placement[1], placement[2]:placement[3]])
            pi.plot(**dict(kwargs))
        # Create the RangeSlider
        # from matplotlib.widgets import RangeSlider
        # slider_ax = fig.add_subplot(axs[self.grid_size[0], :])
        # self._slider = RangeSlider(slider_ax, "Plot Range", 0., 1., valinit=(0., 1.))
        # self._slider.on_changed(self.update_range)

    def update_plot(self):
        for pi in self.projections:
            pi.update_plot()

    # def update_range(self, event):
    #     ...