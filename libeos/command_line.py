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

    return clas.parse_args()



def command_line_options():
    clas   = commandLineArgs()
    update_loglevel(clas.verbose)

    reader_config = ReaderConfig.from_args(clas)
    experiment_config = ExperimentConfig.from_args(clas)
    reduction_config = ReductionConfig.from_args(clas)
    output_config = OutputConfig.from_args(clas)

    return EOSConfig(reader_config, experiment_config, reduction_config, output_config)
