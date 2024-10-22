"""
Classes for stroing various configurations needed for reduction.
"""
from dataclasses import dataclass, field
from typing import Optional, Tuple
from datetime import datetime
from os import path
import numpy as np

import logging

class Defaults:
    # fileIdentifier
    outputPath                  = '.'
    rawPath                     = ['.', path.join('.','raw'), path.join('..','raw'), path.join('..','..','raw')]
    year                        = datetime.now().year
    normalisationFileIdentifier = []
    normalisationMethod         = 'o'
    monitorType                 = 'auto'
    # subtract
    outputName                  = "fromEOS"
    outputFormat                = ['Rqz.ort']
    incidentAngle               = 'alphaF'
    qResolution                 = 0.01
    #timeSlize
    scale                       = [1]
    # autoscale
    lambdaRange                 = [2., 15.]
    thetaRange                  = [-12., 12.]
    thetaRangeR                 = [-0.75, 0.75]
    yRange                      = [11, 41]
    qzRange                     = [0.005, 0.30]
    chopperSpeed                = 500
    chopperPhase                = -13.5
    chopperPhaseOffset          = 7
    muOffset                    = 0
    mu                          = 0
    nu                          = 0
    sampleModel                 = None
    lowCurrentThreshold         = 50
    #
    
    

@dataclass
class ReaderConfig:
    year: int
    rawPath: Tuple[str]
    startTime: Optional[float] = 0

@dataclass
class ExperimentConfig:
    incidentAngle: str 
    chopperPhase: float
    yRange: Tuple[float, float]
    lambdaRange: Tuple[float, float]
    qzRange: Tuple[float, float]
    monitorType: str
    lowCurrentThreshold: float

    sampleModel: Optional[str] = None
    chopperPhaseOffset: float = 0
    mu: Optional[float] = None
    nu: Optional[float] = None
    muOffset: Optional[float] = None

@dataclass
class ReductionConfig:
    normalisationMethod: str
    qResolution: float
    qzRange: Tuple[float, float]
    thetaRange: Tuple[float, float]
    thetaRangeR: Tuple[float, float]

    fileIdentifier: list = field(default_factory=lambda: ["0"])
    scale: list = field(default_factory=lambda: [1]) #per file scaling; if less elements than files use the last one

    autoscale: Optional[Tuple[bool, bool]] = None
    subtract: Optional[str] = None
    normalisationFileIdentifier: Optional[list] = None
    timeSlize: Optional[list] = None

@dataclass
class OutputConfig:
    outputFormats: list
    outputName: str
    outputPath: str

@dataclass
class EOSConfig:
    reader: ReaderConfig
    experiment: ExperimentConfig
    reduction: ReductionConfig
    output: OutputConfig
    
    _call_string_overwrite=None
    
    #@property
    #def call_string(self)->str:
    #    if self._call_string_overwrite:
    #        return self._call_string_overwrite
    #    else:
    #        return self.calculate_call_string()
    
    def call_string(self):
        base = 'python eos.py'
        
        inpt = ''
        if self.reader.year:
            inpt += f' -Y {self.reader.year}'
        else:
            inpt += f' -Y {datetime.now().year}'
        if np.shape(self.reader.rawPath)[0] == 1:
            inpt += f' --rawPath {self.reader.rawPath}'
        if self.reduction.subtract:
            inpt += f' -subtract {self.reduction.subtract}'
        if self.reduction.normalisationFileIdentifier:
            inpt += f' -n {" ".join(self.reduction.normalisationFileIdentifier)}'
        if self.reduction.fileIdentifier:
            inpt += f' -f {" ".join(self.reduction.fileIdentifier)}'

        otpt = ''
        if self.reduction.qResolution:
            otpt += f' -r {self.reduction.qResolution}'
        if self.output.outputPath != '.':
            inpt += f' --outputdPath {self.output.outputPath}'
        if self.output.outputName:
            otpt += f' -o {self.output.outputName}'
        if self.output.outputFormats != ['Rqz.ort']:
            otpt += f' -of {" ".join(self.output.outputFormats)}'
            
        mask = ''    
        if self.experiment.yRange != Defaults.yRange:
            mask += f' -y {" ".join(str(ii) for ii in self.experiment.yRange)}'
        if self.experiment.lambdaRange!= Defaults.lambdaRange:
            mask += f' -l {" ".join(str(ff) for ff in self.experiment.lambdaRange)}'
        if self.reduction.thetaRange != Defaults.thetaRange:
            mask += f' -t {" ".join(str(ff) for ff in self.reduction.thetaRange)}'
        elif self.reduction.thetaRangeR != Defaults.thetaRangeR:
            mask += f' -T {" ".join(str(ff) for ff in self.reduction.thetaRangeR)}'
        if self.experiment.qzRange!= Defaults.qzRange:
            mask += f' -q {" ".join(str(ff) for ff in self.experiment.qzRange)}'

        para = ''
        if self.experiment.chopperPhase != Defaults.chopperPhase:
            para += f' --chopperPhase {self.experiment.chopperPhase}'
        if self.experiment.chopperPhaseOffset != Defaults.chopperPhaseOffset:
            para += f' --chopperPhaseOffset {self.experiment.chopperPhaseOffset}'
        if self.experiment.mu:
            para += f' --mu {self.experiment.mu}'
        elif self.experiment.muOffset:
            para += f' --muOffset {self.experiment.muOffset}'
        if self.experiment.nu:
            para += f' --nu {self.experiment.nu}'

        modl = ''
        if self.experiment.sampleModel:
            modl += f" --sampleModel '{self.experiment.sampleModel}'"

        acts = ''
        if self.reduction.autoscale:
            acts += f' --autoscale {" ".join(str(ff) for ff in self.reduction.autoscale)}'
        if self.reduction.scale != Defaults.scale:
            acts += f' --scale {self.reduction.scale}'
        if self.reduction.timeSlize:
            acts += f' --timeSlize {" ".join(str(ff) for ff in self.reduction.timeSlize)}'

        mlst = base + inpt + otpt 
        if mask:
            mlst += mask
        if para:
            mlst += para
        if acts:
            mlst += acts
        if modl:
            mlst += modl

        if len(mlst) > 70:
            mlst = base + '  ' + inpt + '  ' + otpt 
            if mask:
                mlst += '  ' + mask
            if para:
                mlst += '  ' + para
            if acts:
                mlst += '  ' + acts
            if modl:
                mlst += '  ' + modl

        logging.debug(f'Argument list build in EOSConfig.call_string: {mlst}')
        return  mlst

            
