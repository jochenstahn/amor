import os
import cProfile
from unittest import TestCase
from libeos import options, reduction, logconfig

logconfig.setup_logging()
logconfig.update_loglevel(True, False)

# TODO: add test for new features like proton charge normalization

class FullAmorTest(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pr = cProfile.Profile()

    @classmethod
    def tearDownClass(cls):
        cls.pr.dump_stats("profile_test.prof")

    def setUp(self):
        self.pr.enable()
        self.reader_config = options.ReaderConfig(
                year=2023,
                rawPath=(os.path.join('..', "test_data"),),
                )

    def tearDown(self):
        self.pr.disable()
        for fi in ['test.Rqz.ort', '614.norm']:
            try:
                os.unlink(os.path.join(self.reader_config.rawPath[0], fi))
            except FileNotFoundError:
                pass


    def test_time_slicing(self):
        experiment_config = options.ExperimentConfig(
                chopperPhase=-13.5,
                chopperPhaseOffset=-5,
                monitorType=options.Defaults.monitorType,
                lowCurrentThreshold=options.Defaults.lowCurrentThreshold,
                yRange=(11., 41.),
                lambdaRange=(2., 15.),
                qzRange=(0.005, 0.30),
                incidentAngle=options.Defaults.incidentAngle,
                mu=0,
                nu=0,
                muOffset=0.0,
                sampleModel='air | 10 H2O | D2O'
                )
        reduction_config = options.ReductionConfig(
                normalisationMethod=options.Defaults.normalisationMethod,
                qResolution=0.01,
                qzRange=options.Defaults.qzRange,
                thetaRange=(-12., 12.),
                thetaRangeR=(-12., 12.),
                fileIdentifier=["610"],
                scale=[1],
                normalisationFileIdentifier=[],
                timeSlize=[300.0]
                )
        output_config = options.OutputConfig(
                outputFormats=["Rqz.ort"],
                outputName='test',
                outputPath=os.path.join('..', 'test_results'),
                )
        config=options.EOSConfig(self.reader_config, experiment_config, reduction_config, output_config)
        # run three times to get similar timing to noslicing runs
        reducer = reduction.AmorReduction(config)
        reducer.reduce()
        reducer = reduction.AmorReduction(config)
        reducer.reduce()
        reducer = reduction.AmorReduction(config)
        reducer.reduce()

    def test_noslicing(self):
        experiment_config = options.ExperimentConfig(
                chopperPhase=-13.5,
                chopperPhaseOffset=-5,
                monitorType=options.Defaults.monitorType,
                lowCurrentThreshold=options.Defaults.lowCurrentThreshold,
                yRange=(11., 41.),
                lambdaRange=(2., 15.),
                qzRange=(0.005, 0.30),
                incidentAngle=options.Defaults.incidentAngle,
                mu=0,
                nu=0,
                muOffset=0.0
                )
        reduction_config = options.ReductionConfig(
                qResolution=0.01,
                qzRange=options.Defaults.qzRange,
                normalisationMethod=options.Defaults.normalisationMethod,
                thetaRange=(-12., 12.),
                thetaRangeR=(-12., 12.),
                fileIdentifier=["610", "611", "608,612-613", "609"],
                scale=[1],
                normalisationFileIdentifier=["608"],
                autoscale=(True, True)
                )
        output_config = options.OutputConfig(
                outputFormats=["Rqz.ort"],
                outputName='test',
                outputPath=os.path.join('..', 'test_results'),
                )
        config=options.EOSConfig(self.reader_config, experiment_config, reduction_config, output_config)
        reducer = reduction.AmorReduction(config)
        reducer.reduce()
        # run second time to reuse norm file
        reducer = reduction.AmorReduction(config)
        reducer.reduce()
