import os
import cProfile
from unittest import TestCase
from dataclasses import fields, MISSING
from libeos import options, reduction, logconfig

logconfig.setup_logging()
logconfig.update_loglevel(1)

# TODO: add test for new features like proton charge normalization

class FullAmorTest(TestCase):
    @classmethod
    def setUpClass(cls):
        # generate map for option defaults
        cls._field_defaults = {}
        for opt in [options.ExperimentConfig, options.ReductionConfig, options.OutputConfig]:
            defaults = {}
            for field in fields(opt):
                if field.default not in [None, MISSING]:
                    defaults[field.name] = field.default
                elif field.default_factory not in [None, MISSING]:
                    defaults[field.name] = field.default_factory()
            cls._field_defaults[opt.__name__] = defaults
        cls.pr = cProfile.Profile()

    @classmethod
    def tearDownClass(cls):
        cls.pr.dump_stats("profile_test.prof")

    def setUp(self):
        self.pr.enable()
        self.reader_config = options.ReaderConfig(
                year=2025,
                rawPath=[os.path.join('..', "test_data")],
                )

    def tearDown(self):
        self.pr.disable()
        for fi in ['../test_results/test.Rqz.ort', '../test_results/5952.norm']:
            try:
                os.unlink(os.path.join(self.reader_config.rawPath[0], fi))
            except FileNotFoundError:
                pass

    def test_time_slicing(self):
        experiment_config = options.ExperimentConfig(
                chopperSpeed=self._field_defaults['ExperimentConfig']['chopperSpeed'],
                chopperPhase=-13.5,
                chopperPhaseOffset=-5,
                monitorType=self._field_defaults['ExperimentConfig']['monitorType'],
                lowCurrentThreshold=self._field_defaults['ExperimentConfig']['lowCurrentThreshold'],
                yRange=(11., 41.),
                lambdaRange=(2., 15.),
                incidentAngle=self._field_defaults['ExperimentConfig']['incidentAngle'],
                mu=0,
                nu=0,
                muOffset=0.0,
                sampleModel='air | 10 H2O | D2O'
                )
        reduction_config = options.ReductionConfig(
                normalisationMethod=self._field_defaults['ReductionConfig']['normalisationMethod'],
                qResolution=0.01,
                qzRange=self._field_defaults['ReductionConfig']['qzRange'],
                thetaRange=(-12., 12.),
                thetaRangeR=(-12., 12.),
                fileIdentifier=["5980", "5981", "5983-5985"],
                scale=[1],
                normalisationFileIdentifier=[],
                timeSlize=[300.0]
                )
        output_config = options.OutputConfig(
                outputFormats=[options.OutputFomatOption.Rqz_ort],
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
                chopperSpeed=self._field_defaults['ExperimentConfig']['chopperSpeed'],
                chopperPhase=-13.5,
                chopperPhaseOffset=-5,
                monitorType=self._field_defaults['ExperimentConfig']['monitorType'],
                lowCurrentThreshold=self._field_defaults['ExperimentConfig']['lowCurrentThreshold'],
                yRange=(11., 41.),
                lambdaRange=(2., 15.),
                incidentAngle=self._field_defaults['ExperimentConfig']['incidentAngle'],
                mu=0,
                nu=0,
                muOffset=0.0
                )
        reduction_config = options.ReductionConfig(
                normalisationMethod=self._field_defaults['ReductionConfig']['normalisationMethod'],
                qResolution=0.01,
                qzRange=self._field_defaults['ReductionConfig']['qzRange'],
                thetaRange=(-12., 12.),
                thetaRangeR=(-12., 12.),
                fileIdentifier=["5980"],
                scale=[1],
                normalisationFileIdentifier=["5952"],
                autoscale=(0.0, 0.05),
                )
        output_config = options.OutputConfig(
                outputFormats=[options.OutputFomatOption.Rqz_ort],
                outputName='test',
                outputPath=os.path.join('..', 'test_results'),
                )
        config=options.EOSConfig(self.reader_config, experiment_config, reduction_config, output_config)
        reducer = reduction.AmorReduction(config)
        reducer.reduce()
        # run second time to reuse norm file
        reducer = reduction.AmorReduction(config)
        reducer.reduce()
