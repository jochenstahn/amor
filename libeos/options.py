"""
Classes for stroing various configurations needed for reduction.
"""
from dataclasses import dataclass, field
from typing import Optional, Tuple

@dataclass
class ReaderConfig:
    year: int
    dataPath: str
    startTime: Optional[float] = 0

@dataclass
class ExperimentConfig:
    chopperPhase: float
    yRange: Tuple[float, float]
    lambdaRange: Tuple[float, float]
    qzRange: Tuple[float, float]

    sampleModel: Optional[str] = None
    chopperPhaseOffset: float = 0.0
    mu: Optional[float] = None
    nu: Optional[float] = None
    muOffset: Optional[float] = None
    offSpecular: bool = False

@dataclass
class ReductionConfig:
    qResolution: float
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
    reductoin: ReductionConfig
    output: OutputConfig