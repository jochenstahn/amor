import os
import cProfile
from unittest import TestCase
from dataclasses import fields, MISSING
from eos import options, reduction_reflectivity, logconfig

logconfig.setup_logging()
logconfig.update_loglevel(1)

# TODO: add unit tests for individual parts of reduction

class FullAmorTest(TestCase):
    @classmethod
    def setUpClass(cls):
        # generate map for option defaults
        cls._field_defaults = {}
        for opt in [options.ExperimentConfig, options.ReflectivityReductionConfig, options.ReflectivityOutputConfig]:
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
                rawPath=["test_data"],
                )

    def tearDown(self):
        self.pr.disable()
        for fi in ['test_results/test.Rqz.ort', 'test_results/5952.norm']:
            try:
                os.unlink(fi)
            except FileNotFoundError:
                pass

    def test_time_slicing(self):
        experiment_config = options.ExperimentConfig(
                chopperSpeed=self._field_defaults['ExperimentConfig']['chopperSpeed'],
                chopperPhase=-13.5,
                chopperPhaseOffset=-5,
                monitorType=self._field_defaults['ExperimentConfig']['monitorType'],
                lowCurrentThreshold=self._field_defaults['ExperimentConfig']['lowCurrentThreshold'],
                yRange=(18, 48),
                lambdaRange=(3., 11.5),
                incidentAngle=self._field_defaults['ExperimentConfig']['incidentAngle'],
                mu=0,
                nu=0,
                muOffset=0.0,
                sampleModel='air | 10 H2O | D2O'
                )
        reduction_config = options.ReflectivityReductionConfig(
                normalisationMethod=self._field_defaults['ReflectivityReductionConfig']['normalisationMethod'],
                qResolution=0.01,
                qzRange=self._field_defaults['ReflectivityReductionConfig']['qzRange'],
                thetaRange=(-0.75, 0.75),
                fileIdentifier=["6003-6005"],
                scale=[1],
                normalisationFileIdentifier=[],
                timeSlize=[300.0]
                )
        output_config = options.ReflectivityOutputConfig(
                outputFormats=[options.OutputFomatOption.Rqz_ort],
                outputName='test',
                outputPath='test_results',
                )
        config=options.ReflectivityConfig(self.reader_config, experiment_config, reduction_config, output_config)
        # run three times to get similar timing to noslicing runs
        reducer = reduction_reflectivity.ReflectivityReduction(config)
        reducer.reduce()
        reducer = reduction_reflectivity.ReflectivityReduction(config)
        reducer.reduce()
        reducer = reduction_reflectivity.ReflectivityReduction(config)
        reducer.reduce()

    def test_noslicing(self):
        experiment_config = options.ExperimentConfig(
                chopperSpeed=self._field_defaults['ExperimentConfig']['chopperSpeed'],
                chopperPhase=-13.5,
                chopperPhaseOffset=-5,
                monitorType=self._field_defaults['ExperimentConfig']['monitorType'],
                lowCurrentThreshold=self._field_defaults['ExperimentConfig']['lowCurrentThreshold'],
                yRange=(18, 48),
                lambdaRange=(3., 11.5),
                incidentAngle=self._field_defaults['ExperimentConfig']['incidentAngle'],
                mu=0,
                nu=0,
                muOffset=0.0,
                )
        reduction_config = options.ReflectivityReductionConfig(
                normalisationMethod=self._field_defaults['ReflectivityReductionConfig']['normalisationMethod'],
                qResolution=0.01,
                qzRange=self._field_defaults['ReflectivityReductionConfig']['qzRange'],
                thetaRange=(-0.75, 0.75),
                fileIdentifier=["6003", "6004", "6005"],
                scale=[1],
                normalisationFileIdentifier=["5952"],
                autoscale=(0.0, 0.05),
                )
        output_config = options.ReflectivityOutputConfig(
                outputFormats=[options.OutputFomatOption.Rqz_ort],
                outputName='test',
                outputPath='test_results',
                )
        config=options.ReflectivityConfig(self.reader_config, experiment_config, reduction_config, output_config)
        reducer = reduction_reflectivity.ReflectivityReduction(config)
        reducer.reduce()
        # run second time to reuse norm file
        reducer = reduction_reflectivity.ReflectivityReduction(config)
        reducer.reduce()
