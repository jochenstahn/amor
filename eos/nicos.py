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
                           'amor-nicos')
    update_loglevel(clas.verbose)
    if clas.verbose<2:
        # only log info level in logfile
        logger = logging.getLogger()  # logging.getLogger('quicknxs')
        logger.setLevel(logging.INFO)

    reader_config = ReaderConfig.from_args(clas)
    experiment_config = ExperimentConfig.from_args(clas)
    reduction_config = E2HReductionConfig.from_args(clas)
    config = E2HConfig(reader_config, experiment_config, reduction_config)

    logging.warning('######## amor-nicos - Nicos histogram for Amor ########')
    from eos.reduction_kafka import KafkaReduction

    # only import heavy module if sufficient command line parameters were provided
    from eos.reduction_reflectivity import ReflectivityReduction
    # Create reducer with these arguments
    reducer = KafkaReduction(config)
    # Perform actual reduction
    reducer.reduce()

    logging.info('######## amor-nicos - finished ########')

if __name__ == '__main__':
    main()
