import argparse

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
    clas = argparse.ArgumentParser(description = msg, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    clas.add_argument('-v', '--verbose', action='count', default=0)

    clas_groups = {}

    all_arguments = []
    for cls in [ReaderConfig, ExperimentConfig, OutputConfig, ReductionConfig]:
        all_arguments += cls.get_commandline_parameters()

    all_arguments.sort() # parameters are sorted alphabetically, unless they have higher priority
    for cpc in all_arguments:
        if not cpc.group in clas_groups:
            clas_groups[cpc.group] = clas.add_argument_group(cpc.group)
        if cpc.short_form:
            clas_groups[cpc.group].add_argument(
                    f'-{cpc.short_form}', f'--{cpc.argument}', **cpc.add_argument_args
                    )
        else:
            clas_groups[cpc.group].add_argument(
                    f'--{cpc.argument}', **cpc.add_argument_args
                    )

    #
    # output = clas.add_argument_group('output')
    # output.add_argument("-o", "--outputName",
    #                         default = Defaults.outputName,
    #                         help = "output file name (withot suffix)")
    # output.add_argument("-op", "--outputPath",
    #                         type = str,
    #                         default = Defaults.outputPath,
    #                         help = "path for output")
    # output.add_argument("-of", "--outputFormat",
    #                         nargs = '+',
    #                         default = Defaults.outputFormat,
    #                         help = "one of [Rqz.ort, Rlt.ort]")
    # output.add_argument("-ai", "--incidentAngle",
    #                         type = str,
    #                         default = Defaults.incidentAngle,
    #                         help = "calulate alpha_i from [alphaF, mu, nu]",
    #                         )
    # output.add_argument("-r", "--qResolution",
    #                         default = Defaults.qResolution,
    #                         type = float,
    #                         help = "q_z resolution")
    # output.add_argument("-ts", "--timeSlize",
    #                         nargs = '+',
    #                         type = float,
    #                         help = "time slizing <interval> ,[<start> [,stop]]")
    # output.add_argument("-s", "--scale",
    #                         nargs = '+',
    #                         default = Defaults.scale,
    #                         type = float,
    #                         help = "scaling factor for R(q_z)")
    # output.add_argument("-S", "--autoscale",
    #                         nargs = 2,
    #                         type = float,
    #                         help = "scale to 1 in the given q_z range")
    #
    # masks = clas.add_argument_group('masks')
    # masks.add_argument("-l", "--lambdaRange",
    #                         default = Defaults.lambdaRange,
    #                         nargs = 2,
    #                         type = float,
    #                         help = "wavelength range")
    # masks.add_argument("-t", "--thetaRange",
    #                         default = Defaults.thetaRange,
    #                         nargs = 2,
    #                         type = float,
    #                         help = "absolute theta range")
    # masks.add_argument("-T", "--thetaRangeR",
    #                         default = Defaults.thetaRangeR,
    #                         nargs = 2,
    #                         type = float,
    #                         help = "relative theta range")
    # masks.add_argument("-y", "--yRange",
    #                         default = Defaults.yRange,
    #                         nargs = 2,
    #                         type = int,
    #                         help = "detector y range")
    # masks.add_argument("-q", "--qzRange",
    #                         default = Defaults.qzRange,
    #                         nargs = 2,
    #                         type = float,
    #                         help = "q_z range")
    # masks.add_argument("-ct", "--lowCurrentThreshold",
    #                         default = Defaults.lowCurrentThreshold,
    #                         type = float,
    #                         help = "proton current threshold for discarding neutron pulses")
    #
    #
    # overwrite = clas.add_argument_group('overwrite')
    # overwrite.add_argument("-cs", "--chopperSpeed",
    #                         default = Defaults.chopperSpeed,
    #                         type = float,
    #                         help = "chopper speed in rpm")
    # overwrite.add_argument("-cp", "--chopperPhase",
    #                         default = Defaults.chopperPhase,
    #                         type = float,
    #                         help = "chopper phase")
    # overwrite.add_argument("-co", "--chopperPhaseOffset",
    #                         default = Defaults.chopperPhaseOffset,
    #                         type = float,
    #                         help = "phase offset between chopper opening and trigger pulse")
    # overwrite.add_argument("-m", "--muOffset",
    #                         default = Defaults.muOffset,
    #                         type = float,
    #                         help = "mu offset")
    # overwrite.add_argument("-mu", "--mu",
    #                         default = Defaults.mu,
    #                         type = float,
    #                         help ="value of mu")
    # overwrite.add_argument("-nu", "--nu",
    #                         default = Defaults.nu,
    #                         type = float,
    #                         help = "value of nu")
    # overwrite.add_argument("-sm", "--sampleModel",
    #                         default = Defaults.sampleModel,
    #                         type = str,
    #                         help = "1-line orso sample model description")

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
    update_loglevel(clas.verbose)

    reader_config = ReaderConfig.from_args(clas)
    experiment_config = ExperimentConfig.from_args(clas)
    reduction_config = ReductionConfig.from_args(clas)
    output_config = OutputConfig.from_args(clas)

    return EOSConfig(reader_config, experiment_config, reduction_config, output_config)
