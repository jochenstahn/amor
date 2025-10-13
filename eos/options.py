"""
Classes for stroing various configurations needed for reduction.
"""
import argparse
from dataclasses import dataclass, field, Field, fields, MISSING
from typing import get_args, get_origin, List, Optional, Tuple, Union
from datetime import datetime
from os import path
import numpy as np

import logging


try:
    from enum import StrEnum
except ImportError:
    try:
        # python <3.11 try to use backports
        from backports.strenum import StrEnum
    except ImportError:
        # python <3.10 use Enum instead
        from enum import Enum as StrEnum

@dataclass
class CommandlineParameterConfig:
    argument: str # default parameter for command line resutign ins "--argument"
    add_argument_args: dict # all arguments that will be passed to add_argument method
    short_form: Optional[str] = None
    group: str = 'misc'
    priority: int = 0

    def __gt__(self, other):
        """
        Sort required arguments first, then use priority, then name
        """
        return (not self.add_argument_args.get('required', False), -self.priority, self.argument)>(
            not other.add_argument_args.get('required', False), -other.priority, other.argument)

class ArgParsable:
    def __init_subclass__(cls):
        # create a nice documentation string that takes help into account
        cls.__doc__ = cls.__name__ + " Parameters:\n"
        for key, typ in cls.__annotations__.items():
            if get_origin(typ) is Union and type(None) in get_args(typ):
                optional = True
                typ = get_args(typ)[0]
            else:
                optional = False

            value = getattr(cls, key, None)
            try:
                cls.__doc__ += f"    {key} ({typ.__name__})"
            except AttributeError:
                cls.__doc__ += f"    {key}"
            if isinstance(value, Field):
                if value.default is not MISSING:
                    cls.__doc__ += f" = {value.default}"
                if 'help' in value.metadata:
                    cls.__doc__ += f" - {value.metadata['help']}"
            elif value is not None:
                    cls.__doc__ += f" = {value}"
            if optional:
                cls.__doc__ += " [Optional]"
            cls.__doc__ += "\n"
        return cls

    @classmethod
    def get_commandline_parameters(cls) -> List[CommandlineParameterConfig]:
        """
        Return a list of arguments used in building the command line parameters.

        Union types besides Optional are not supported.
        """
        output = []
        for field in fields(cls):
            args={}
            if field.default is not MISSING:
                args['default'] = field.default
                args['required'] = False
            elif field.default_factory is not MISSING:
                args['default'] = field.default_factory()
                args['required'] = False
            else:
                args['required'] = True
            if get_origin(field.type) is Union and type(None) in get_args(field.type):
                # optional argument
                typ = get_args(field.type)[0]
                del(args['default'])
            else:
                typ = field.type
            if get_origin(typ) is list:
                args['nargs'] = '+'
                typ = get_args(typ)[0]
                if get_origin(typ) is tuple:
                    # tuple of items are put together during evaluation
                    typ = get_args(typ)[0]
            elif get_origin(typ) is tuple:
                args['nargs'] = len(get_args(typ))
                typ = get_args(typ)[0]
            if issubclass(typ, StrEnum):
                args['choices'] = [ci.value for ci in typ]
                if field.default is not MISSING:
                    args['default'] = field.default.value
                typ = str

            if typ is bool:
                args['action'] = 'store_false' if field.default else 'store_true'
            else:
                args['type'] = typ

            if 'help' in field.metadata:
                args['help'] = field.metadata['help']

            output.append(CommandlineParameterConfig(
                    field.name,
                    add_argument_args=args,
                    group=field.metadata.get('group', 'misc'),
                    short_form=field.metadata.get('short', None),
                    priority=field.metadata.get('priority', 0),
                    ))
        return output

    @classmethod
    def get_default(cls, key):
        """
        Return the default argument for an attribute, None if it doesn't exist.
        """
        for field in fields(cls):
            if field.name != key:
                continue
            if field.default is not MISSING:
                return field.default
            elif field.default_factory is not MISSING:
                return field.default_factory()
        return None

    def is_default(self, key):
        value = getattr(self, key)
        return value == self.get_default(key)

    @classmethod
    def from_args(cls, args: argparse.Namespace):
        """
        Create the child class from the command line argument Namespace object.
        All attributes that are not needed for this class are ignored.
        """
        inpargs = {}
        for field in fields(cls):
            value = getattr(args, field.name)
            typ = field.type
            if get_origin(field.type) is Union and type(None) in get_args(field.type):
                # optional argument
                typ = get_args(field.type)[0]
            if get_origin(typ) is list:
                item_typ = get_args(typ)[0]
                if get_origin(item_typ) is tuple:
                    # tuple of items are put together during evaluation
                    tuple_length = len(get_args(item_typ))
                    value = [tuple(value[i*tuple_length+j] for j in range(tuple_length)) for i in range(len(value)//tuple_length)]
            if isinstance(typ, type) and issubclass(typ, StrEnum):
                # convert str to enum
                try:
                    value = typ(value)
                except ValueError:
                    choices = [ci.value for ci in typ]
                    raise ValueError(f"Parameter --{field.name} has to be one of {choices}")

            inpargs[field.name] = value
        return cls(**inpargs)

# definition of command line arguments

@dataclass
class ReaderConfig(ArgParsable):
    year: int = field(
            default=datetime.now().year,
            metadata={
                'short': 'Y', 
                'group': 'input data', 
                'help': 'year the measurement was performed',
                },
            )
    rawPath: List[str] = field(
            default_factory=lambda: ['.', path.join('.','raw'), path.join('..','raw'), path.join('..','..','raw')],
            metadata={
                'short': 'rp',
                'group': 'input data',
                'help': 'search paths for hdf files',
                },
            )
    startTime: Optional[float] = field(
            default = None,
            metadata={         
                'short': 'st',     
                'group': 'data manicure',
                'help': 'set time zero other than the start of the data aquisition',
                },
            )

class IncidentAngle(StrEnum):
    alphaF = 'alphaF'
    mu = 'mu'
    nu = 'nu'

class MonitorType(StrEnum):
    auto = 'a'
    proton_charge = 'p'
    time = 't'
    neutron_monitor = 'n'
    debug = 'x'

MONITOR_UNITS = {
    MonitorType.neutron_monitor: 'cnts',
    MonitorType.proton_charge: 'mC',
    MonitorType.time: 's',
    MonitorType.auto: 'various',
    MonitorType.debug: 'mC',
    }

@dataclass
class ExperimentConfig(ArgParsable):
    chopperPhase: float = field(
            default=0,
            metadata={
                'short': 'cp',
                'group': 'instrument settings',
                'help': 'phase between opening of chopper 1 and closing of chopper 2 window',
                },
            ) 
    chopperPhaseOffset: float = field(
            default=-5,
            metadata={
                'short': 'co',
                'group': 'instrument settings',  
                'help': 'phase between chopper 1 index pulse and closing edge',
                },                               
            )                          
    chopperSpeed: float = field(
            default=500,
            metadata={
                'short': 'cs',
                'group': 'instrument settings',
                'help': 'rotation speed of the chopper disks in rpm',
                },
            )
    yRange: Tuple[int, int] = field(
            default=(18, 48),
            metadata={
                'short': 'y',
                'group': 'region of interest',
                'help': 'horizontal pixel range on the detector to be used',
                },
            )
    lambdaRange: Tuple[float, float] = field(
            default_factory=lambda: [3, 12.5],
            metadata={
                'short': 'l',
                'group': 'region of interest',
                'help': 'wavelength range to be used (in angstrom)',
                },
            )
    lowCurrentThreshold: float = field(
            default=50,
            metadata={
                'short': 'pt',
                'group': 'instrument settings',
                'help': 'proton current below which the events are ignored (per chopper pulse)',
                },
            )

    incidentAngle: IncidentAngle = field(
            default=IncidentAngle.alphaF,
            metadata={
                'short': 'ai',
                'group': 'instrument settings',
                'help': 'calculate alphaI = [alphaF], [mu]+kappa+delta_kappa or ([nu]+kappa+delta_kappa)/2',
                },
            )
    alphaF = 'alphaF'
    sampleModel: Optional[str] = field(
            default=None,
            metadata={
                'short': 'sm',
                'group': 'sample',  
                'help': 'orso type string to describe the sample in one line',
                },                               
            )                          
    mu: Optional[float] = field(
            default=None,
            metadata={
                'short': 'mu',
                'group': 'sample',  
                'help': 'inclination of the sample surface w.r.t. the instrument horizon',
                },                               
            )                          
    nu: Optional[float] = field(
            default=None,
            metadata={
                'short': 'nu',
                'group': 'sample',  
                'help': 'inclination of the detector w.r.t. the instrument horizon',
                },                               
            )                          
    muOffset: Optional[float] = field(
            default=0,
            metadata={
                'short': 'm',
                'group': 'sample',  
                'help': 'correction offset for mu misalignment (mu_real = mu_file + mu_offset)',
                },                               
            )                          
    monitorType: MonitorType = field(
            default=MonitorType.proton_charge,
            metadata={
                'short': 'mt',
                'group': 'instrument settings',
                'help': 'one of [a]uto, [p]rotonCurrent, [t]ime or [n]eutronMonitor',
                },                               
            )                          

class NormalisationMethod(StrEnum):
    direct_beam = 'd'
    over_illuminated = 'o'
    under_illuminated = 'u'

@dataclass
class ReflectivityReductionConfig(ArgParsable):
    fileIdentifier: List[str] = field(
            metadata={
                'short':    'f',
                'priority': 100,
                'group':    'input data',
                'help':     'file number(s) or offset (if < 1)',
                },
            )

    qResolution: float = field(
            default=0.01,
            metadata={
                'short': 'r',
                'group': 'data manicure',
                'help': 'output resolution of q-scale Delta q / q',
                },
            )
    qzRange: Tuple[float, float] = field(
            default_factory=lambda: [0.005, 0.51],
            metadata={
                'short': 'q',
                'group': 'region of interest',
                'help': '?',
                },
            )
    thetaRange: Tuple[float, float] = field(
            default_factory=lambda: [-12., 12.],
            metadata={
                'short': 't',
                'group': 'region of interest',
                'help': 'absolute theta region of interest',
                },
            )
    thetaRangeR: Tuple[float, float] = field(
            default_factory=lambda: [-0.75, 0.75],
            metadata={
                'short': 'T',
                'group': 'region of interest',
                'help': 'theta region of interest w.r.t. beam center',
                },
            )
    thetaFilters: List[Tuple[float, float]] = field(
            default_factory=lambda: [],
            metadata={
                'short': 'TF',
                'group': 'region of interest',
                'help': 'add one or more theta ranges that will be filtered in reduction',
                },
            )
    normalisationMethod: NormalisationMethod = field(
            default=NormalisationMethod.over_illuminated,
            metadata={
                'short': 'nm', 
                'priority': 90,
                'group': 'input data',
                'help': 'normalisation method: [o]verillumination, [u]nderillumination, [d]irect_beam'})
    scale: List[float] = field(
            default_factory=lambda: [1.],
            metadata={
                'short': 's',
                'group': 'data manicure',
                'help': '(list of) scaling factors, if less elements than files use the last one',
                },
            ) 
    autoscale: Tuple[float, float] = field(
           default=None,
           metadata={
               'short': 'S',
               'group': 'data manicure',
               'help': '',
               },
           )
    subtract: Optional[str] = field(
            default=None, 
            metadata={
                'short': 'sub',
                'group': 'input data', 
                'help': 'File with R(q_z) curve to be subtracted (in .Rqz.ort format)'})
    normalisationFileIdentifier: Optional[List[str]] = field(
            default_factory=lambda: [None], 
            metadata={
                'short': 'n', 
                'priority': 90,
                'group': 'input data', 
                'help': 'file number(s) of normalisation measurement'})
    timeSlize: Optional[List[float]] = field(
            default= None,
            metadata={
                'short': 'ts',
                'group': 'region of interest',
                'help': 'time slizing <interval> ,[<start> [,stop]]',
                },
            )


class OutputFomatOption(StrEnum):
    Rqz_ort = "Rqz.ort"
    Rqz_orb = "Rqz.orb"
    Rlt_ort = "Rlt.ort"
    Rlt_orb = "Rlt.orb"
    ort = "ort"
    orb = "orb"
    Rqz = "Rqz"
    Rlt = "Rlt"


class PlotColormaps(StrEnum):
    gist_ncar = "gist_ncar"
    viridis = "viridis"
    inferno = "inferno"
    gist_rainbow = "gist_rainbow"
    nipy_spectral = "nipy_spectral"
    jochen_deluxe = "jochen_deluxe"

@dataclass
class ReflectivityOutputConfig(ArgParsable):
    outputFormats: List[OutputFomatOption] = field(
            default_factory=lambda: ['Rqz.ort'],
            metadata={
                'short': 'of',
                'group': 'output',
                'help': 'one of "Rqz[.ort]", "Rlt[.ort]" or both with "ort"',
                },
            )
    outputName: str = field(
            default='fromEOS',
            metadata={
                'short': 'o',
                'group': 'output',
                'help': '?',
                },
            )
    outputPath: str = field(
            default='.',
            metadata={
                'short': 'op',
                'group': 'output',
                'help': '?',
                },
            )
    plot: bool = field(
            default=False,
            metadata={
                'group': 'output',
                'help': 'show matplotlib graphs of results',
                },
            )

    plot_colormap: PlotColormaps = field(
            default=PlotColormaps.gist_ncar,
            metadata={
                'short': 'pcmap',
                'group': 'output',
                'help': 'matplotlib colormap used in lambda-theta graphs when plotting',
                },
            )

    def _output_format_list(self, outputFormat):
        format_list = []
        if OutputFomatOption.ort in outputFormat\
                or OutputFomatOption.Rqz_ort in outputFormat\
                or OutputFomatOption.Rqz in outputFormat:
            format_list.append(OutputFomatOption.Rqz_ort)
        if OutputFomatOption.ort in outputFormat\
                or OutputFomatOption.Rlt_ort in outputFormat\
                or OutputFomatOption.Rlt in outputFormat:
            format_list.append(OutputFomatOption.Rlt_ort)
        if OutputFomatOption.orb in outputFormat\
                or OutputFomatOption.Rqz_orb in outputFormat\
                or OutputFomatOption.Rqz in outputFormat:
            format_list.append(OutputFomatOption.Rqz_orb)
        if OutputFomatOption.orb in outputFormat\
                or OutputFomatOption.Rlt_orb in outputFormat\
                or OutputFomatOption.Rlt in outputFormat:
            format_list.append(OutputFomatOption.Rlt_orb)
        return sorted(format_list, reverse=True)

    def __post_init__(self):
        self.outputFormats = self._output_format_list(self.outputFormats)


# ===================================

@dataclass
class ReflectivityConfig:
    reader: ReaderConfig
    experiment: ExperimentConfig
    reduction: ReflectivityReductionConfig
    output: ReflectivityOutputConfig
    
    _call_string_overwrite=None
    
    #@property
    #def call_string(self)->str:
    #    if self._call_string_overwrite:
    #        return self._call_string_overwrite
    #    else:
    #        return self.calculate_call_string()
    
    def call_string(self):
        base = 'eos'
        
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

        mask += f' -y {" ".join(str(ii) for ii in self.experiment.yRange)}'
        mask += f' -l {" ".join(str(ff) for ff in self.experiment.lambdaRange)}'
        mask += f' -t {" ".join(str(ff) for ff in self.reduction.thetaRange)}'
        mask += f' -T {" ".join(str(ff) for ff in self.reduction.thetaRangeR)}'
        mask += f' -q {" ".join(str(ff) for ff in self.reduction.qzRange)}'

        para = ''
        # TODO: Check if we want these parameters for defaults
        para += f' --chopperPhase {self.experiment.chopperPhase}'
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
        # TODO: Check if should be shown if not default
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

class E2HPlotSelection(StrEnum):
    All = 'all'
    Raw = 'raw'
    YZ = 'Iyz'
    LT = 'Ilt'
    YT = 'Iyt'
    TZ = 'Itz'
    Q = 'Iq'
    L = 'Il'
    T = 'It'
    ToF = 'tof'


class E2HPlotArguments(StrEnum):
    Default = 'def'
    OutputFile = 'file'
    Logarithmic = 'log'
    Linear = 'lin'

@dataclass
class E2HReductionConfig(ArgParsable):
    fileIdentifier: str = field(
            default='0',
            metadata={
                'short':    'f',
                'priority': 100,
                'group':    'input data',
                'help':     'file number(s) or offset (if < 1), events2histogram only accepts one short code',
                },
            )

    show_plot: bool = field(
            default=False,
            metadata={
                'short': 'sp',
                'group': 'output',
                'help': 'show matplotlib graphs of results',
                },
            )

    plot: E2HPlotSelection = field(
            default=E2HPlotSelection.All,
            metadata={
                'short': 'p',
                'group': 'output',
                'help': 'select what to plot or write',
                },
            )

    kafka: bool = field(
            default=False,
            metadata={
                'group': 'output',
                'help': 'send result to kafka for Nicos',
                },
            )

    plotArgs: E2HPlotArguments = field(
            default=E2HPlotArguments.Default,
            metadata={
                'short': 'pa',
                'group': 'output',
                'help': 'select configuration for plot',
                },
            )

    update: bool = field(
            default=False,
            metadata={
                'short': 'u',
                'group': 'output',
                'help': 'keep running in the background and update plot when file is modified, implies --plot',
                },
            )

    fast: bool = field(
            default=False,
            metadata={
                'group': 'input data',
                'help': 'skip some reduction steps to speed up calculations',
                },
            )

    normalizationModel: bool = field(
            default=False,
            metadata={
                'short': 'nm',
                'group': 'input data',
                'help': 'use model for incoming spectrum and divergence to normalize for reflectivity',
                },
            )

    plot_colormap: PlotColormaps = field(
            default=PlotColormaps.jochen_deluxe,
            metadata={
                'short': 'pcmap',
                'group': 'output',
                'help': 'matplotlib colormap used in lambda-theta graphs when plotting',
                },
            )

    max_events: int = field(
            default = 10_000_000,
            metadata={
                'group': 'input data',
                'help':  'maximum number of events read at once',
                },
            )

    thetaRangeR: Tuple[float, float] = field(
            default_factory=lambda: [-0.75, 0.75],
            metadata={
                'short': 'T',
                'group': 'region of interest',
                'help': 'theta region of interest w.r.t. beam center',
                },
            )

    thetaFilters: List[Tuple[float, float]] = field(
            default_factory=lambda: [],
            metadata={
                'short': 'TF',
                'group': 'region of interest',
                'help': 'add one or more theta ranges that will be filtered in reduction',
                },
            )

    fontsize: float = field(
            default=8.,
            metadata={
                'short': 'pf',
                'group': 'output',
                'help': 'font size for graphs',
                },
            )

@dataclass
class E2HConfig:
    reader: ReaderConfig
    experiment: ExperimentConfig
    reduction: E2HReductionConfig
