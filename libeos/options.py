"""
Classes for stroing various configurations needed for reduction.
"""
from dataclasses import dataclass, field
from typing import Optional, Tuple
from datetime import datetime


class Defaults:
    #fileIdentifier
    normalisationFileIdentifier = []
    dataPath                    = '.'
    year                        = datetime.now().year
    #subtract
    outputName                  = "fromEOS"
    outputFormat                = ['Rqz.ort']
    incidentAngle               = 'alphaF'
    qResolution                 = 0.01
    incidentAngle               = 'alphaF'
    #timeSlize
    scale                       = [1]
    #autoscale
    lambdaRange                 = [2., 15.]
    thetaRange                  = [-12., 12.]
    thetaRangeR                 = [-0.7, 0.7]
    yRange                      = [11, 41]
    qzRange                     = [0.005, 0.30]
    chopperSpeed                = 500
    chopperPhase                = -13.5
    chopperPhaseOffset          = -5
    muOffset                    = 0
    mu                          = 0
    nu                          = 0
    sampleModel                 = None
    #
    
    

@dataclass
class ReaderConfig:
    year: int
    dataPath: str
    startTime: Optional[float] = 0

@dataclass
class ExperimentConfig:
    incidentAngle: str 
    chopperPhase: float
    yRange: Tuple[float, float]
    lambdaRange: Tuple[float, float]
    qzRange: Tuple[float, float]

    sampleModel: Optional[str] = None
    chopperPhaseOffset: float = 0
    mu: Optional[float] = None
    nu: Optional[float] = None
    muOffset: Optional[float] = None

@dataclass
class ReductionConfig:
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
        if self.reader_config.year:
            inpt += f' -Y {self.reader_config.year}'
        else:
            inpt += f' -Y {datetime.now().year}'
        if self.reader_config.dataPath != '.':
            inpt += f' --dataPath {self.reader_config.dataPath}'
        if self.reduction_config.subtract:
            inpt += f' -subtract {self.reduction_config.subtract}'
        if self.reduction_config.normalisationFileIdentifier:
            inpt += f' -n {" ".join(self.reduction_config.normalisationFileIdentifier)}'
        if self.reduction_config.fileIdentifier:
            inpt += f' -f {" ".join(self.reduction_config.fileIdentifier)}'

        otpt = ''
        if self.reduction_config.qResolution:
            otpt += f' -r {self.reduction_config.qResolution}'
        if self.output_config.outputName:
            otpt += f' -o {self.output_config.outputName}'
        if self.output_config.outputFormats != ['Rqz.ort']:
            otpt += f' -of {" ".join(self.output_config.outputFormats)}'
            
        mask = ''    
        if self.experiment_config.yRange != Defaults.yRange:
            mask += f' -y {" ".join(str(ii) for ii in self.experiment_config.yRange)}'
        if self.experiment_config.lambdaRange!= Defaults.lambdaRange:
            mask += f' -l {" ".join(str(ff) for ff in self.experiment_config.lambdaRange)}'
        if self.reduction_config.thetaRange != Defaults.thetaRange:
            mask += f' -T {" ".join(str(ff) for ff in self.reduction_config.thetaRange)}'
        elif self.reduction_config.thetaRangeR != Defaults.thetaRangeR:
            mask += f' -t {" ".join(str(ff) for ff in self.reduction_config.thetaRangeR)}'
        if self.experiment_config.qzRange!= Defaults.qzRange:
            mask += f' -q {" ".join(str(ff) for ff in self.experiment_config.qzRange)}'

        para = ''
        if self.experiment_config.chopperPhase != Defaults.chopperPhase:
            para += f' --chopperPhase {self.experiment_config.chopperPhase}'
        if self.experiment_config.chopperPhaseOffset != Defaults.chopperPhaseOffset:
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
        if self.reduction_config.scale != Defaults.scale:
            acts += f' --scale {self.reduction_config.scale}'
        if self.reduction_config.timeSlize:
            acts += f' --timeSlize {" ".join(str(ff) for ff in self.reduction_config.timeSlize)}'

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

        return  mlst

            
