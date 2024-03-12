"""
Class to handle Orso header information that changes gradually during the reduction process.
"""

import platform
import sys
from datetime import datetime

from orsopy import fileio

from . import __version__


class Header:
    """orso compatible output file header content"""

    def __init__(self):
        self.owner                           = None
        self.experiment                      = None
        self.sample                          = None
        self.measurement_instrument_settings = None
        self.measurement_scheme              = None
        self.measurement_data_files          = []
        self.measurement_additional_files    = []

        self.reduction = fileio.Reduction(
            software    = fileio.Software('eos', version=__version__),
            call        = self.create_call_string(),
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
    #-------------------------------------------------------------------------------------------------
    def create_call_string(self):
        callString = ' '.join(sys.argv)
        if '-Y' not in callString:
            callString += f' -Y {datetime.now().year}'
        return callString
    #-------------------------------------------------------------------------------------------------
    def call_string():
        base = 'python eos.py '

        inpt = ''
        if clas.year:
            inpt += f' --year {clas.year}'
        else:
            inpt += f' --year {datetime.now().year}'
        if clas.dataPath:
            inpt += f' --dataPath {clas.dataPath}'
        if clas.subtract:
            inpt += f' -subtract {clas.subtract}'
        if clas.normalisationFileIdentifier:
            inpt += f' -r {clas.normalisationFileIdentifier}'
        # get file list somehow
        if ...:
            inpt += f' -n {file_list}'
        else:
            inpt += f' -n {clas.fileIdentifier}'

        otpt = ''
        if outputFormats != 'Rqz.ort':
            otpt =  f" -of  '{outputFormats}'"
        if clas.outputName:
            otpt += f' -o {clas.outputName}'
        else:
            pass
            # default name
            
        mask = ''    
        if clas.yRange:
            mask += f' -y {clas.yRange}'
        if clas.lambdaRange:
            mask += f' -l {clas.lambdaRange}'
        if clas.thetaRange:
            mask += f' -- thetaRange {clas.thetaRange}'
        elif clas.thetaRangeR:
            mask += f' -t {clas.thetaRangeR}'
        if clas.qzRange:
            mask += f' -q {clas.qzRange}'
        if clas.qResolution:
            mask += f' -a {clas.qResolution}'

        para = ''
        if clas.chopperPhase:
            para += f' --chopperPhase {clas.chopperPhase}'
        if clas.chopperPhaseOffset:
            para += f' --chopperPhaseOffset {clas.chopperPhaseOffset}'
        if clas.mu:
            para += f' --mu {clas.mu}'
        elif clas.muOffset:
            para += f' --muOffset {clas.muOffset}'
        if clas.nu:
            para += f' --nu {clas.nu}'

        if clas.sampleModel:
            modl =  f" --sampleModel '{clas.sampleModel}'"

        acts = ''
        if clas.autoscale:
            acts += f' --autoscale {clas.autoscale}'
        if clas.scale:
            acts += f' --scale {clas.scale}'
        if clas.timeSlize:
            acts += f' --timeSlize {clas.timeSlize}'

        mlst = base + '\n' + inpt + '\n' + outp 
        if mask:
            mlst += '\n' + mask
        if para:
            mlst += '\n' + para
        if acts:
            mlst += '\n' + acts
        if mask:
            mlst += '\n' + modl

        return  mlst
