"""
eos reduces measurements performed on Amor@SINQ, PSI

Author: Jochen Stahn (algorithms, python draft),
        Artur Glavic (structuring and optimisation of code)
"""

import logging

# need to do absolute import here as pyinstaller requires it
from eos.options import ReflectivityConfig, ReaderConfig, ExperimentConfig, ReflectivityReductionConfig, ReflectivityOutputConfig
from eos.command_line import commandLineArgs
from eos.logconfig import setup_logging, update_loglevel


def main():
    setup_logging()

    # read command line arguments and generate classes holding configuration parameters
    clas = commandLineArgs([ReaderConfig, ExperimentConfig, ReflectivityReductionConfig, ReflectivityOutputConfig],
                           'eos')
    update_loglevel(clas.verbose)

    reader_config = ReaderConfig.from_args(clas)
    experiment_config = ExperimentConfig.from_args(clas)
    reduction_config = ReflectivityReductionConfig.from_args(clas)
    output_config = ReflectivityOutputConfig.from_args(clas)
    config = ReflectivityConfig(reader_config, experiment_config, reduction_config, output_config)

    logging.warning('######## eos - data reduction for Amor ########')

    # only import heavy module if sufficient command line parameters were provided
    from eos.reduction_reflectivity import ReflectivityReduction
    # Create reducer with these arguments
    reducer = ReflectivityReduction(config)
    # Perform actual reduction
    reducer.reduce()

    logging.info('######## eos - finished ########')

if __name__ == '__main__':
    main()
