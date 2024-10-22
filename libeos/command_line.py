import argparse

from .logconfig import update_loglevel
from .options import ReaderConfig, EOSConfig, ExperimentConfig, OutputConfig, ReductionConfig, Defaults


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
    input_data.add_argument("-f", "--fileIdentifier",
                            required = True,
                            nargs = '+',
                            help = "file number(s) or offset (if < 1)")
    input_data.add_argument("-n", "--normalisationFileIdentifier",
                            default = Defaults.normalisationFileIdentifier,
                            nargs = '+',
                            help = "file number(s) of normalisation measurement")
    input_data.add_argument("-rp", "--rawPath", 
                            type = str,
                            default = Defaults.rawPath,
                            help = "ath to directory with .hdf files")
    input_data.add_argument("-Y", "--year",
                            default = Defaults.year,
                            type = int,
                            help = "year the measurement was performed")
    input_data.add_argument("-sub", "--subtract",
                            help = "R(q_z) curve to be subtracted (in .Rqz.ort format)")
    input_data.add_argument("-nm", "--normalisationMethod",
                            default = Defaults.normalisationMethod,
                            help = "normalisation method: [o]verillumination, [u]nderillumination, [d]irect_beam")
    input_data.add_argument("-mt", "--monitorType",
                            type = str,
                            default = Defaults.monitorType,
                            help = "one of [p]rotonCurrent, [t]ime or [n]eutronMonitor")

    output = clas.add_argument_group('output')
    output.add_argument("-o", "--outputName",
                            default = Defaults.outputName,
                            help = "output file name (withot suffix)")
    output.add_argument("-op", "--outputPath",
                            type = str,
                            default = Defaults.outputPath,
                            help = "path for output")
    output.add_argument("-of", "--outputFormat",
                            nargs = '+',
                            default = Defaults.outputFormat,
                            help = "one of [Rqz.ort, Rlt.ort]")
    output.add_argument("-ai", "--incidentAngle",
                            type = str,
                            default = Defaults.incidentAngle,
                            help = "calulate alpha_i from [alphaF, mu, nu]",
                            )
    output.add_argument("-r", "--qResolution",
                            default = Defaults.qResolution,
                            type = float,
                            help = "q_z resolution")
    output.add_argument("-ts", "--timeSlize",
                            nargs = '+',
                            type = float,
                            help = "time slizing <interval> ,[<start> [,stop]]")
    output.add_argument("-s", "--scale",
                            nargs = '+',
                            default = Defaults.scale,
                            type = float,
                            help = "scaling factor for R(q_z)")
    output.add_argument("-S", "--autoscale",
                            nargs = 2,
                            type = float,
                            help = "scale to 1 in the given q_z range")

    masks = clas.add_argument_group('masks')
    masks.add_argument("-l", "--lambdaRange",
                            default = Defaults.lambdaRange,
                            nargs = 2,
                            type = float,
                            help = "wavelength range")
    masks.add_argument("-t", "--thetaRange",
                            default = Defaults.thetaRange,
                            nargs = 2,
                            type = float,
                            help = "absolute theta range")
    masks.add_argument("-T", "--thetaRangeR",
                            default = Defaults.thetaRangeR,
                            nargs = 2,
                            type = float,
                            help = "relative theta range")
    masks.add_argument("-y", "--yRange",
                            default = Defaults.yRange,
                            nargs = 2,
                            type = int,
                            help = "detector y range")
    masks.add_argument("-q", "--qzRange",
                            default = Defaults.qzRange,
                            nargs = 2,
                            type = float,
                            help = "q_z range")
    masks.add_argument("-ct", "--lowCurrentThreshold",
                            default = Defaults.lowCurrentThreshold,
                            type = float,
                            help = "proton current threshold for discarding neutron pulses")


    overwrite = clas.add_argument_group('overwrite')
    overwrite.add_argument("-cs", "--chopperSpeed",
                            default = Defaults.chopperSpeed,
                            type = float,
                            help = "chopper speed in rpm")
    overwrite.add_argument("-cp", "--chopperPhase",
                            default = Defaults.chopperPhase,
                            type = float,
                            help = "chopper phase")
    overwrite.add_argument("-co", "--chopperPhaseOffset",
                            default = Defaults.chopperPhaseOffset,
                            type = float,
                            help = "phase offset between chopper opening and trigger pulse")
    overwrite.add_argument("-m", "--muOffset",
                            default = Defaults.muOffset,
                            type = float,
                            help = "mu offset")
    overwrite.add_argument("-mu", "--mu",
                            default = Defaults.mu,
                            type = float,
                            help ="value of mu")
    overwrite.add_argument("-nu", "--nu",
                            default = Defaults.nu,
                            type = float,
                            help = "value of nu")
    overwrite.add_argument("-sm", "--sampleModel",
                            default = Defaults.sampleModel,
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
        rawPath                      = clas.rawPath,
        )
    experiment_config = ExperimentConfig(
        sampleModel                  = clas.sampleModel,
        chopperPhase                 = clas.chopperPhase,
        chopperPhaseOffset           = clas.chopperPhaseOffset,
        yRange                       = clas.yRange,
        lambdaRange                  = clas.lambdaRange,
        qzRange                      = clas.qzRange,
        lowCurrentThreshold          = clas.lowCurrentThreshold,
        incidentAngle                = clas.incidentAngle,
        mu                           = clas.mu,
        nu                           = clas.nu,
        muOffset                     = clas.muOffset,
        monitorType                  = clas.monitorType,
        )
    reduction_config = ReductionConfig(
        qResolution                  = clas.qResolution,
        qzRange                      = clas.qzRange,
        autoscale                    = clas.autoscale,
        thetaRange                   = clas.thetaRange,
        thetaRangeR                  = clas.thetaRangeR,
        fileIdentifier               = clas.fileIdentifier,
        scale                        = clas.scale,
        subtract                     = clas.subtract,
        normalisationFileIdentifier  = clas.normalisationFileIdentifier,
        normalisationMethod          = clas.normalisationMethod,
        timeSlize                    = clas.timeSlize,
        )
    output_config = OutputConfig(
        outputFormats                = output_format_list(clas.outputFormat),
        outputName                   = clas.outputName,
        outputPath                   = clas.outputPath,
        )

    return EOSConfig(reader_config, experiment_config, reduction_config, output_config)
