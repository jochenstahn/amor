"""
Defines how to normalize a focusing reflectometry dataset by a reference measurement.
"""
import logging
import os
import numpy as np
from typing import List


from .event_data_types import EventDatasetProtocol
from .header import Header
from .options import NormalisationMethod
from .instrument import Grid


class Normalisation:
    normFileList = List[str]
    normAngle: float
    normMonitor: float
    norm_lz: np.ndarray

    def __init__(self, reference:EventDatasetProtocol, normalisationMethod: NormalisationMethod, grid: Grid):
        self.normAngle = reference.geometry.nu-reference.geometry.mu
        lamda_e = reference.data.events.lamda
        detZ_e = reference.data.events.detZ
        self.normMonitor = np.sum(reference.data.pulses.monitor)
        norm_lz, _, _ = np.histogram2d(lamda_e, detZ_e, bins=(grid.lamda(), grid.z()))
        norm_lz = np.where(norm_lz>2, norm_lz, np.nan)
        if normalisationMethod==NormalisationMethod.direct_beam:
            self.norm_lz = np.flip(norm_lz, 1)
        else:
            # correct for reference sm reflectivity
            lamda_l = grid.lamda()
            theta_z = self.normAngle+reference.geometry.delta_z
            lamda_lz = (grid.lz().T*lamda_l[:-1]).T
            theta_lz = grid.lz()*theta_z
            qz_lz = 4.0*np.pi*np.sin(np.deg2rad(theta_lz))/lamda_lz
            # TODO: introduce variable for `m` and propably for the slope
            Rsm_lz = np.ones(np.shape(qz_lz))
            Rsm_lz = np.where(qz_lz>0.0217, 1-(qz_lz-0.0217)*(0.0625/0.0217), Rsm_lz)
            Rsm_lz = np.where(qz_lz>0.0217*5, np.nan, Rsm_lz)
            self.norm_lz = norm_lz/Rsm_lz
        self.normFileList = [os.path.basename(entry) for entry in reference.file_list]

    @classmethod
    def from_file(cls, filename) -> 'Normalisation':
        logging.warning(f'normalisation matrix: found and using {filename}')
        self = super().__new__(cls)
        with open(filename, 'rb') as fh:
            self.normFileList = np.load(fh, allow_pickle=True)
            self.normAngle = np.load(fh, allow_pickle=True)
            self.norm_lz = np.load(fh, allow_pickle=True)
            self.normMonitor = np.load(fh, allow_pickle=True)
        return self

    @classmethod
    def unity(cls, grid:Grid) -> 'Normalisation':
        logging.warning(f'normalisation is unity')
        self = super().__new__(cls)
        self.norm_lz = grid.lz()
        self.normFileList = []
        self.normAngle = 1.
        self.normMonitor = 1.
        return self

    def safe(self, filename):
        with open(filename, 'wb') as fh:
            np.save(fh, np.array(self.normFileList), allow_pickle=False)
            np.save(fh, np.array(self.normAngle), allow_pickle=False)
            np.save(fh, self.norm_lz, allow_pickle=False)
            np.save(fh, self.normMonitor, allow_pickle=False)

    def update_header(self, header:Header):
        header.measurement_additional_files = self.normFileList
