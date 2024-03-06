import argparse
from datetime import date

from .logconfig import update_loglevel
from .options import ReaderConfig, EOSConfig, ExperimentConfig, OutputConfig, ReductionConfig


def commandLineArgs():
    """
    Process command line argument.
    The type of the default values is used for conversion and validation.
    """
    msg = "eos reads data from (one or several) raw file(s) of the .hdf format, \
           performs various corrections, conversations and projections and exports\
           the resulting reflectivity in an orso-compatible format."
    clas = argparse.ArgumentParser(description = msg)

    input_data = clas.add_argument_group('input data')
    input_data.add_argument("-n", "--fileIdentifier",
                            default = ['0'],
                            nargs = '+',
                            help = "file number(s) or offset (if negative)")
    input_data.add_argument("-r", "--normalisationFileIdentifier",
                            default = [],
                            nargs = '+',
                            help = "file number(s) of normalisation measurement")
    input_data.add_argument("-d", "--dataPath",
                            type = str,
                            default = '.',
                            help = "relative path to directory with .hdf files")
    input_data.add_argument("-Y", "--year",
                            default = date.today().year,
                            type = int,
                            help = "year the measurement was performed")
    input_data.add_argument("-sub", "--subtract",
                            help = "R(q_z) curve to be subtracted (in .Rqz.ort format)")

    output = clas.add_argument_group('output')
    output.add_argument("-o", "--outputName",
                            default = "fromEOS",
                            help = "output file name (withot suffix)")
    output.add_argument("-of", "--outputFormat",
                            nargs = '+',
                            default = ['Rqz.ort'])
    output.add_argument("--offSpecular",
                            type = bool,
                            default = False,
                            )
    output.add_argument("-a", "--qResolution",
                            default = 0.01,
                            type = float,
                            help = "q_z resolution")
    output.add_argument("-ts", "--timeSlize",
                            nargs = '+',
                            type = float,
                            help = "time slizing <interval> ,[<start> [,stop]]")
    output.add_argument("-s", "--scale",
                            nargs = '+',
                            default = [1],
                            type = float,
                            help = "scaling factor for R(q_z)")
    output.add_argument("-S", "--autoscale",
                            nargs = 2,
                            type = float,
                            help = "scale to 1 in the given q_z range")

    masks = clas.add_argument_group('masks')
    masks.add_argument("-l", "--lambdaRange",
                            default = [2., 15.],
                            nargs = 2,
                            type = float,
                            help = "wavelength range")
    masks.add_argument("-t", "--thetaRange",
                            default = [-12., 12.],
                            nargs = 2,
                            type = float,
                            help = "absolute theta range")
    masks.add_argument("-T", "--thetaRangeR",
                            default = [-12., 12.],
                            nargs = 2,
                            type = float,
                            help = "relative theta range")
    masks.add_argument("-y", "--yRange",
                            default = [11, 41],
                            nargs = 2,
                            type = int,
                            help = "detector y range")
    masks.add_argument("-q", "--qzRange",
                            default = [0.005, 0.30],
                            nargs = 2,
                            type = float,
                            help = "q_z range")

    overwrite = clas.add_argument_group('overwrite')
    overwrite.add_argument("-cs", "--chopperSpeed",
                            type = float,
                            help = "chopper speed in rpm")
    overwrite.add_argument("-cp", "--chopperPhase",
                            default = -13.5,
                            type = float,
                            help = "chopper phase")
    overwrite.add_argument("-co", "--chopperPhaseOffset",
                            default = -5,
                            type = float,
                            help = "phase offset between chopper opening and trigger pulse")
    overwrite.add_argument("-m", "--muOffset",
                            default = 0.,
                            type = float,
                            help = "mu offset")
    overwrite.add_argument("-mu", "--mu",
                            default = 0,
                            type = float,
                            help ="value of mu")
    overwrite.add_argument("-nu", "--nu",
                            default = 0,
                            type = float,
                            help = "value of nu")
    overwrite.add_argument("-sm", "--sampleModel",
                            type = str,
                            help = "1-line orso sample model description")

    misc = clas.add_argument_group('misc')
    misc.add_argument('-v', '--verbose', action='store_true')
    misc.add_argument('-vv', '--debug', action='store_true')

    return clas.parse_args()


def expand_file_list(short_notation):
    """Evaluate string entry for file number lists"""
    #log().debug('Executing get_flist')
    file_list=[]
    for i in short_notation.split(','):
        if '-' in i:
            if ':' in i:
                step = i.split(':', 1)[1]
                file_list += range(int(i.split('-', 1)[0]), int((i.rsplit('-', 1)[1]).split(':', 1)[0])+1, int(step))
            else:
                step = 1
                file_list += range(int(i.split('-', 1)[0]), int(i.split('-', 1)[1])+1, int(step))
        else:
            file_list += [int(i)]

    return sorted(file_list)


def output_format_list(outputFormat):
    format_list = []
    if 'ort' in outputFormat or 'Rqz.ort' in outputFormat or 'Rqz' in outputFormat:
        format_list.append('Rqz.ort')
    if 'ort' in outputFormat or 'Rlt.ort' in outputFormat or 'Rlt' in outputFormat:
        format_list.append('Rlt.ort')
    if 'orb' in outputFormat or 'Rqz.orb' in outputFormat or 'Rqz' in outputFormat:
        format_list.append('Rqz.orb')
    if 'orb' in outputFormat or 'Rlt.orb' in outputFormat or 'Rlt' in outputFormat:
        format_list.append('Rlt.orb')
    return sorted(format_list, reverse=True)

def command_line_options():
    clas   = commandLineArgs()
    update_loglevel(clas.verbose, clas.debug)

    reader_config = ReaderConfig(
        year                         = clas.year,
        dataPath                     = clas.dataPath
        )
    experiment_config = ExperimentConfig(
        sampleModel                  = clas.sampleModel,
        chopperPhase                 = clas.chopperPhase,
        chopperPhaseOffset           = clas.chopperPhaseOffset,
        yRange                       = clas.yRange,
        lambdaRange                  = clas.lambdaRange,
        qzRange                      = clas.qzRange,
        offSpecular                  = clas.offSpecular,
        mu                           = clas.mu,
        nu                           = clas.nu,
        muOffset                     = clas.muOffset
        )
    reduction_config = ReductionConfig(
        qResolution                  = clas.qResolution,
        autoscale                    = clas.autoscale,
        thetaRange                   = clas.thetaRange,
        thetaRangeR                  = clas.thetaRangeR,
        fileIdentifier               = clas.fileIdentifier,
        scale                        = clas.scale,
        subtract                     = clas.subtract,
        normalisationFileIdentifier  = clas.normalisationFileIdentifier,
        timeSlize                    = clas.timeSlize
        )
    output_config = OutputConfig(
        outputFormats                = output_format_list(clas.outputFormat),
        outputName                   = clas.outputName
        )

    return EOSConfig(reader_config, experiment_config, reduction_config, output_config)
