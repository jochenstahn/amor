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

    def __init__(self, config):
        self.owner                           = None
        self.experiment                      = None
        self.sample                          = None
        self.measurement_instrument_settings = None
        self.measurement_scheme              = None
        self.measurement_data_files          = []
        self.measurement_additional_files    = []

        self.reduction = fileio.Reduction(
            software    = fileio.Software('eos', version=__version__),
            call        = self.call_string(config),
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
    def call_string(self, config):
        self.experiment_config = config.experiment
        self.reader_config = config.reader
        self.reduction_config = config.reductoin
        self.output_config = config.output
        
        base = 'python eos.py'

        inpt = ''
        if self.reader_config.year:
            inpt += f' --year {self.reader_config.year}'
        else:
            inpt += f' --year {datetime.now().year}'
        if self.reader_config.dataPath != '.':
            inpt += f' --dataPath {self.reader_config.dataPath}'
        if self.reduction_config.subtract:
            inpt += f' -subtract {self.reduction_config.subtract}'
        if self.reduction_config.normalisationFileIdentifier:
            inpt += f' -r {" ".join(self.reduction_config.normalisationFileIdentifier)}'
        # get file list somehow
        if False:
            pass
            #inpt += f' -n {file_list}'
        else:
            inpt += f' -n {" ".join(self.reduction_config.fileIdentifier)}'

        otpt = ''
        if self.reduction_config.qResolution:
            otpt += f' -q {self.reduction_config.qResolution}'
        if self.output_config.outputFormats != 'Rqz.ort':
            otpt =  f' -of {" ".join(self.output_config.outputFormats)}'
        if self.output_config.outputName:
            otpt += f' -o {self.output_config.outputName}'
        else:
            pass
            # default name
            
        mask = ''    
        if self.experiment_config.yRange != [11, 41]:
            mask += f' -y {" ".join(str(ii) for ii in self.experiment_config.yRange)}'
        if self.experiment_config.lambdaRange!= [2, 15]:
            mask += f' -l {" ".join(str(ff) for ff in self.experiment_config.lambdaRange)}'
        if self.reduction_config.thetaRange != [-12, 12]:
            mask += f' --thetaRange {" ".join(str(ff) for ff in self.reduction_config.thetaRange)}'
        elif self.reduction_config.thetaRangeR != [-12, 12]:
            mask += f' -t {" ".join(str(ff) for ff in self.reduction_config.thetaRangeR)}'
        if self.experiment_config.qzRange!= [0.005, 0.3]:
            mask += f' -q {" ".join(str(ff) for ff in self.experiment_config.qzRange)}'

        para = ''
        if self.experiment_config.chopperPhase != -13.5:
            para += f' --chopperPhase {self.experiment_config.chopperPhase}'
        if self.experiment_config.chopperPhaseOffset != -5:
            para += f' --chopperPhaseOffset {self.experiment_config.chopperPhaseOffset}'
        if self.experiment_config.mu:
            para += f' --mu {self.experiment_config.mu}'
        elif self.experiment_config.muOffset:
            para += f' --muOffset {self.experiment_config.muOffset}'
        if self.experiment_config.nu:
            para += f' --nu {self.experiment_config.nu}'

        modl = ''
        if self.experiment_config.sampleModel:
            modl += f" --sampleModel '{self.experiment_config.sampleModel}'"

        acts = ''
        if self.reduction_config.autoscale:
            acts += f' --autoscale {" ".join(str(ff) for ff in self.reduction_config.autoscale)}'
        if self.reduction_config.scale != [1]:
            acts += f' --scale {self.reduction_config.scale}'
        if self.reduction_config.timeSlize:
            acts += f' --timeSlize {" ".join(str(ff) for ff in self.reduction_config.timeSlize)}'

        #experiment_config = ExperimentConfig(
        #    offSpecular                  = clas.offSpecular,
        #    )

        mlst = base + '  ' + inpt + '  ' + otpt 
        if mask:
            mlst += '  ' + mask
        if para:
            mlst += '  ' + para
        if acts:
            mlst += '  ' + acts
        if modl:
            mlst += '  ' + modl

        return  mlst
