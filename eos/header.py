"""
Class to handle Orso header information that changes gradually during the reduction process.
"""

import platform
import sys
from datetime import datetime
from typing import List, Literal

from orsopy import fileio

from . import __version__


class Header:
    """orso compatible output file header content"""
    owner: fileio.Person
    experiment: fileio.Experiment
    sample: fileio.Sample
    measurement_instrument_settings: fileio.InstrumentSettings
    measurement_scheme: Literal["angle- and energy-dispersive", "angle-dispersive", "energy-dispersive"]
    measurement_data_files: List[fileio.File]
    measurement_additional_files: List[fileio.File]


    def __init__(self):
        self.measurement_data_files          = []
        self.measurement_additional_files    = []

        self.reduction = fileio.Reduction(
            software    = fileio.Software('eos', version=__version__),
            call        = 'placeholder',
            computer    = platform.node(),
            timestamp   = datetime.now(),
            creator     = None,
            corrections = ['histogramming in lambda and alpha_f',
                           'gravity'],
            )
    #-------------------------------------------------------------------------------------------------
    def data_source(self):
        return fileio.DataSource(
            self.owner,
            self.experiment,
            self.sample,
            fileio.Measurement(
                instrument_settings = self.measurement_instrument_settings,
                scheme              = self.measurement_scheme,
                data_files          = self.measurement_data_files,
                additional_files    = self.measurement_additional_files,
                ),
        )
    #-------------------------------------------------------------------------------------------------
    def columns(self):
        cols = [
            fileio.Column('Qz', '1/angstrom', 'normal momentum transfer'),
            fileio.Column('R', '', 'specular reflectivity'),
            fileio.ErrorColumn(error_of='R', error_type='uncertainty', distribution='gaussian', value_is='sigma'),
            fileio.ErrorColumn(error_of='Qz', error_type='resolution', distribution='gaussian', value_is='sigma'),
            ]
        return cols
    #-------------------------------------------------------------------------------------------------
    def orso_header(self, columns=None, extra_columns=[]):
        """
        Generate ORSO header from a copy of this class' data.
        """
        ds = fileio.DataSource.from_dict(self.data_source().to_dict())
        red = fileio.Reduction.from_dict(self.reduction.to_dict())
        if columns is None:
            columns = self.columns()
        return fileio.Orso(ds, red, columns+extra_columns)
