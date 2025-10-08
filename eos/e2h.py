"""
events2histogram vizualising data from Amor@SINQ, PSI

Author: Jochen Stahn (algorithms, python draft),
        Artur Glavic (structuring and optimisation of code)
"""
import logging

# need to do absolute import here as pyinstaller requires it
from eos.options import E2HConfig, ReaderConfig, ExperimentConfig, E2HReductionConfig
from eos.command_line import commandLineArgs
from eos.logconfig import setup_logging, update_loglevel


def main():
    setup_logging()
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    # read command line arguments and generate classes holding configuration parameters
    clas = commandLineArgs([ReaderConfig, ExperimentConfig, E2HReductionConfig],
                           'events2histogram')
    update_loglevel(clas.verbose)

    reader_config = ReaderConfig.from_args(clas)
    experiment_config = ExperimentConfig.from_args(clas)
    reduction_config = E2HReductionConfig.from_args(clas)
    config = E2HConfig(reader_config, experiment_config, reduction_config)

    logging.warning('######## events2histogram - data vizualization for Amor ########')
    from eos.reduction_e2h import E2HReduction

    # only import heavy module if sufficient command line parameters were provided
    from eos.reduction_reflectivity import ReflectivityReduction
    # Create reducer with these arguments
    reducer = E2HReduction(config)
    # Perform actual reduction
    reducer.reduce()

    logging.info('######## events2histogram - finished ########')

if __name__ == '__main__':
    main()
