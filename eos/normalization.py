"""
Defines how to normalize a focusing reflectometry dataset by a reference measurement.
"""
import logging
import os
import numpy as np
from typing import List, Optional

from .event_data_types import EventDatasetProtocol
from .header import Header
from .options import NormalisationMethod
from .instrument import Detector, LZGrid


class LZNormalisation:
    file_list = List[str]
    angle: float
    monitor: float
    norm: np.ndarray

    def __init__(self, reference:EventDatasetProtocol, normalisationMethod: NormalisationMethod, grid: LZGrid):
        self.angle = reference.geometry.nu-reference.geometry.mu
        lamda_e = reference.data.events.lamda
        detZ_e = reference.data.events.detZ
        self.monitor = np.sum(reference.data.pulses.monitor)
        norm_lz, _, _ = np.histogram2d(lamda_e, detZ_e, bins=(grid.lamda(), grid.z()))
        norm_lz = np.where(norm_lz>2, norm_lz, np.nan)
        if normalisationMethod==NormalisationMethod.direct_beam:
            self.norm = np.flip(norm_lz, 1)
        else:
            # correct for reference sm reflectivity
            lamda_l = grid.lamda()
            theta_z = self.angle+Detector.delta_z
            lamda_lz = (grid.lz().T*lamda_l[:-1]).T
            theta_lz = grid.lz()*theta_z
            qz_lz = 4.0*np.pi*np.sin(np.deg2rad(theta_lz))/lamda_lz
            # TODO: introduce variable for `m` and propably for the slope
            # Correct reflectivity of m=5 supermirror
            Rsm_lz = np.ones(np.shape(qz_lz))
            Rsm_lz = np.where(qz_lz>0.0217, 1-(qz_lz-0.0217)*(0.0625/0.0217), Rsm_lz)
            Rsm_lz = np.where(qz_lz>0.0217*5, np.nan, Rsm_lz)
            self.norm = norm_lz/Rsm_lz
        self.file_list = [os.path.basename(entry) for entry in reference.file_list]

    @classmethod
    def from_file(cls, filename, check_hash=None) -> Optional['LZNormalisation']:
        self = super().__new__(cls)
        with open(filename, 'rb') as fh:
            hash = str(np.load(fh, allow_pickle=True))
            self.file_list = np.load(fh, allow_pickle=True)
            self.angle = np.load(fh, allow_pickle=True)
            self.norm = np.load(fh, allow_pickle=True)
            self.monitor = np.load(fh, allow_pickle=True)
        if check_hash is not None and hash != check_hash:
            logging.info('    file hash does not match this reduction configuration')
            raise ValueError('file hash does not match this reduction configuration')
        return self

    @classmethod
    def unity(cls, grid:LZGrid) -> 'LZNormalisation':
        logging.warning(f'normalisation is unity')
        self = super().__new__(cls)
        self.norm = grid.lz()
        self.file_list = []
        self.angle = 1.
        self.monitor = 1.
        return self

    @classmethod
    def model(cls, grid:LZGrid) -> 'LZNormalisation':
        # generate a normalization based on angular and wavelength distribution model
        # TODO: add options for sample size for better absolute normalization
        logging.warning(f'normalisation is model')
        self = super().__new__(cls)
        self.angle = 1.0
        self.monitor = 4e6

        lamda_l  = grid.lamda()
        lamda_c = (lamda_l[:-1]+lamda_l[1:])/2

        delta = np.rad2deg(np.arctan2(grid.z(), Detector.distance))/2.0
        delta_c = (delta[:-1]+delta[1:])/2-delta.mean()

        # approximate spectrum by Maxwell-Boltzmann and intensity by linear footprint
        a = 3.8
        Ilambda = np.sqrt(2./np.pi)*lamda_c**2/a**3*np.exp(-lamda_c**2/(2.*a**2))
        Idelta = np.where(abs(delta_c)<0.75, (self.angle-delta_c), np.nan)

        self.norm = 1e6*Ilambda[:, np.newaxis]*Idelta[np.newaxis, :]

        return self

    def safe(self, filename, hash):
        with open(filename, 'wb') as fh:
            np.save(fh, hash, allow_pickle=False)
            np.save(fh, np.array(self.file_list), allow_pickle=False)
            np.save(fh, np.array(self.angle), allow_pickle=False)
            np.save(fh, self.norm, allow_pickle=False)
            np.save(fh, self.monitor, allow_pickle=False)

    def update_header(self, header:Header):
        header.measurement_additional_files = self.file_list
