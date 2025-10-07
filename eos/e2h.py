import logging

# need to do absolute import here as pyinstaller requires it
from eos.options import E2HConfig, ReaderConfig, ExperimentConfig
from eos.command_line import commandLineArgs
from eos.logconfig import setup_logging, update_loglevel

def main():
    setup_logging()

    # read command line arguments and generate classes holding configuration parameters
    clas = commandLineArgs([ReaderConfig, ExperimentConfig],
                           'events2histogram')
    update_loglevel(clas.verbose)

    reader_config = ReaderConfig.from_args(clas)
    experiment_config = ExperimentConfig.from_args(clas)
    config = E2HConfig(reader_config, experiment_config, )

    logging.warning('######## events2histogram - data vizualization for Amor ########')

    # only import heavy module if sufficient command line parameters were provided
    from eos.reduction_reflectivity import ReflectivityReduction
    # Create reducer with these arguments
    reducer = ReflectivityReduction(config)
    # Perform actual reduction
    reducer.reduce()

    logging.info('######## events2histogram - finished ########')

if __name__ == '__main__':
    main()
