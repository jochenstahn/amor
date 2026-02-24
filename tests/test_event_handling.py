import os
import numpy as np

from unittest import TestCase
from datetime import datetime
from copy import deepcopy

from orsopy.fileio import Person, Experiment, Sample, InstrumentSettings, Value, ValueRange, Polarization

from eos import const
from eos.header import Header
from eos.event_data_types import EVENT_BITMASKS, AmorGeometry, AmorTiming, AmorEventStream, \
    EventDataAction, EventDatasetProtocol, PACKET_TYPE, PC_TYPE, PULSE_TYPE, EVENT_TYPE, append_fields
from eos.event_handling import ApplyPhaseOffset, ApplyParameterOverwrites, CorrectChopperPhase, CorrectSeriesTime, \
    AssociatePulseWithMonitor, FilterMonitorThreshold, FilterStrangeTimes, TofTimeCorrection, ApplyMask
from eos.event_analysis import ExtractWalltime, MergeFrames, AnalyzePixelIDs, CalculateWavelength, CalculateQ, \
    FilterQzRange
from eos.options import MonitorType, IncidentAngle, ExperimentConfig


class MockEventData:
    """
    Simulated dataset to be used with event handling unit tests
    """
    geometry: AmorGeometry
    timing: AmorTiming
    data: AmorEventStream

    def __init__(self):
        self.geometry = AmorGeometry(mu=2.0, nu=1.0, kap=0.5, kad=0.0, div=1.5,
                                     chopperSeparation=1000.0, detectorDistance=4000., chopperDetectorDistance=18842.)
        self.timing = AmorTiming(
                ch1TriggerPhase=-9.1, ch2TriggerPhase=6.75,
                chopperPhase=0.17, chopperSpeed=500., tau=0.06
                )
        self.create_data()

    def create_data(self):
        # list of events, here with random time of fligh and pixel location
        events = np.recarray((10000, ), dtype=EVENT_TYPE)
        events.tof = np.random.uniform(low=0., high=0.12, size=events.shape)
        events.pixelID = np.random.randint(0, 28671, size=events.shape)
        events.mask = 0

        # list of data packates containing previous events
        packets = np.recarray((1000,), dtype=PACKET_TYPE)
        packets.start_index = np.linspace(0, events.shape[0]-1, packets.shape[0], dtype=np.uint32)
        packets.time = np.linspace(1700000000000000000, 1700000000000000000+3_600_000,
                                   packets.shape[0], dtype=np.int64)

        # chopper pulses within the measurement time
        pulses = np.recarray((packets.shape[0],), dtype=PULSE_TYPE)
        pulses.monitor = 1.0
        pulses.time = packets.time

        # proton current information with independent timing
        proton_current = np.recarray((50,), dtype=PC_TYPE)
        proton_current.current = 1500.0
        proton_current[np.random.randint(0, proton_current.shape[0]-1, 10)] = 0. # random time with no current
        proton_current.time = np.linspace(1700000000000000300, 1700000000000000000+3_600_000,
                                   proton_current.shape[0], dtype=np.int64)

        self.data = AmorEventStream(events, packets, pulses, proton_current)
        self.orig_data = deepcopy(self.data)


    def append(self, other):
        raise NotImplementedError("Just for testing, no append")

    def update_header(self, header:Header):
        # update a header with the information read from file
        header.owner = Person(name="test user", affiliation='PSI')
        header.experiment = Experiment(title='test experiment', instrument='amor',
                                       start_date=datetime.now(), probe="neutron")
        header.sample = Sample(name='test sample')
        header.measurement_instrument_settings = InstrumentSettings(incident_angle=Value(1.5, 'deg'),
            wavelength = ValueRange(3.0, 12.5, 'angstrom'),
            polarization = Polarization.unpolarized)

class TestActionClass(TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Create test classes to be used
        """
        class T1(EventDataAction):
            def perform_action(self, event: EventDatasetProtocol):
                event.data.events.mask += 1
        class T2(EventDataAction):
            def perform_action(self, event: EventDatasetProtocol):
                event.data.events.mask += 2
        class T4(EventDataAction):
            def perform_action(self, event: EventDatasetProtocol):
                event.data.events.mask += 4
        cls.T1=T1; cls.T2=T2; cls.T4=T4

        class H1(EventDataAction):
            def perform_action(self, event: EventDatasetProtocol):
                ...
            def update_header(self, header:Header) ->None:
                header.sample.name = 'h1'
        class H2(EventDataAction):
            def perform_action(self, event: EventDatasetProtocol):
                ...
            def update_header(self, header: Header) -> None:
                header.sample.name = 'h2'
        class HN(EventDataAction):
            def __init__(self, n):
                self._n = n
            def perform_action(self, event: EventDatasetProtocol):
                ...
            def update_header(self, header: Header) -> None:
                header.sample.name = self._n
        cls.H1=H1; cls.H2=H2; cls.HN = HN

    def setUp(self):
        self.d = MockEventData()
        self.header = Header()
        self.d.update_header(self.header)

    def test_individual(self):
        t1 = self.T1()
        t2 = self.T2()
        t4 = self.T4()

        np.testing.assert_array_equal(self.d.data.events.mask, 0)
        t1.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, 1)
        t2.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, 3)
        t4.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, 7)

    def test_header(self):
        h1 = self.H1()
        h2 = self.H2()
        h3 = self.HN('h3')
        h4 = self.HN('h4')

        h1.update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h1')
        h2.update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h2')
        h3.update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h3')
        h4.update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h4')

    def test_combination(self):
        t1 = self.T1()
        t2 = self.T2()
        t4 = self.T4()
        t12 = t1 | t2
        t24 = t2 | t4
        t1224 = t1 | t2 | t2 | t4
        t1224b = t12 | t24

        np.testing.assert_array_equal(self.d.data.events.mask, 0)
        t12.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, 3)
        t24.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, 9)

        t1224.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, 18)
        t1224b.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, 27)


    def test_combine_header(self):
        h1 = self.H1()
        h2 = self.H2()
        h3 = self.HN('h3')
        h4 = self.HN('h4')

        (h1|h2).update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h2')
        (h2|h1).update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h1')
        (h3|h4).update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h4')
        (h4|h3).update_header(self.header)
        self.assertEqual(self.header.sample.name, 'h3')

    def test_abstract_misssing(self):
        with self.assertRaises(TypeError):
            class E(EventDataAction):
                ...
            _ = E()

    def test_hash(self):
        """
        Check that hashes of different actions are different but
        instances of same action have same hash
        """
        t1 = self.T1()
        t1b = self.T1()
        t2 = self.T2()
        t4 = self.T4()
        h3 = self.HN('h3')
        h3b = self.HN('h3')
        h4 = self.HN('h4')

        self.assertNotEqual(t1.action_hash(), t2.action_hash())
        self.assertNotEqual(t2.action_hash(), t4.action_hash())
        self.assertNotEqual(t1.action_hash(), t4.action_hash())
        self.assertNotEqual(h3.action_hash(), h4.action_hash())
        self.assertEqual(t1.action_hash(), t1b.action_hash())
        self.assertEqual(h3.action_hash(), h3b.action_hash())


class TestSimpleActions(TestCase):
    def setUp(self):
        self.d = MockEventData()
        self.header = Header()
        self.d.update_header(self.header)

    def test_chopper_phase(self):
        cp = CorrectChopperPhase()
        cp.perform_action(self.d)
        np.testing.assert_array_equal(
                self.d.data.events.tof,
                self.d.orig_data.events.tof+
                self.d.timing.tau*(self.d.timing.ch1TriggerPhase-self.d.timing.chopperPhase/2)/180
                )

    def _extract_walltime(self):
        # Extract wall time for events and orig copy
        wt = ExtractWalltime()
        d = self.d.data
        self.d.data = self.d.orig_data
        wt.perform_action(self.d)
        self.d.data = d
        wt.perform_action(self.d)

    def test_extract_walltime(self):
        self._extract_walltime()
        # wallTime should be always a time present in the packet times
        np.testing.assert_array_equal(np.isin(self.d.data.events.wallTime, self.d.data.packets.time), True)
        # make sure extraction works on both original and copy
        np.testing.assert_array_equal(self.d.data.events.wallTime, self.d.orig_data.events.wallTime)

    def test_series_time(self):
        corr = 100
        ct = CorrectSeriesTime(corr)

        with self.assertRaises(ValueError):
            ct.perform_action(self.d)

        self._extract_walltime()


        ct.perform_action(self.d)
        np.testing.assert_array_equal(
                self.d.data.pulses.time,
                self.d.orig_data.pulses.time-corr
                )
        np.testing.assert_array_equal(
                self.d.data.events.wallTime,
                self.d.orig_data.events.wallTime-corr
                )
        np.testing.assert_array_equal(
                self.d.data.proton_current.time,
                self.d.orig_data.proton_current.time-corr
                )

    def test_associate_monitor(self):
        amPC = AssociatePulseWithMonitor(MonitorType.proton_charge)
        amT = AssociatePulseWithMonitor(MonitorType.time)
        amN = AssociatePulseWithMonitor(MonitorType.neutron_monitor)

        self.d.data.pulses.monitor = 13
        amN.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.pulses.monitor, 1)

        self.d.data.pulses.monitor = 13
        amT.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.pulses.monitor, np.float32(2*self.d.timing.tau))

        self.d.data.pulses.monitor = 13
        amPC.perform_action(self.d)
        pcm = self.d.data.proton_current.current *2*self.d.timing.tau*1e-3
        np.testing.assert_array_equal(np.isin(self.d.data.pulses.monitor, pcm), True)

    def test_filter_monitor_threashold(self):
        amPC = AssociatePulseWithMonitor(MonitorType.proton_charge)
        fmt = amPC | FilterMonitorThreshold(1000.)
        fma = amPC | FilterMonitorThreshold(2000.)
        fm0 = amPC | FilterMonitorThreshold(-1.0)

        with self.assertRaises(ValueError):
            fmt.perform_action(self.d)

        self._extract_walltime()
        fm0.perform_action(self.d)
        self.assertEqual(self.d.data.events.mask.sum(), 0)
        fmt.perform_action(self.d)
        # calculate, which events should have 0 monitor
        zero_times = self.d.data.pulses.time[self.d.data.pulses.monitor==0]
        zero_sum = np.isin(self.d.data.events.wallTime, zero_times).sum()
        self.assertEqual(self.d.data.events.mask.sum(), zero_sum*EVENT_BITMASKS['MonitorThreshold'])
        # filter all events
        self.d.data.events.mask = 0
        fma.perform_action(self.d)
        self.assertEqual(self.d.data.events.mask.sum(), self.d.data.events.shape[0]*EVENT_BITMASKS['MonitorThreshold'])

    def test_filter_strage_times(self):
        st = FilterStrangeTimes()

        st.perform_action(self.d)
        self.assertEqual(self.d.data.events.mask.sum(), 0)

        # half events should be strange times (outside of ToF frame)
        self.d.data.events.tof += self.d.timing.tau
        st.perform_action(self.d)
        self.assertEqual(self.d.data.events.mask.sum(),
                         (self.d.data.events.tof>2*self.d.timing.tau).sum()*EVENT_BITMASKS['StrangeTimes'])

    def test_apply_phase_offset(self):
        action = ApplyPhaseOffset(12.5)
        action.perform_action(self.d)
        self.assertEqual(self.d.timing.ch1TriggerPhase, 12.5)

    def test_apply_parameter_overwrites(self):
        action = ApplyParameterOverwrites(ExperimentConfig(muOffset=0.25, mu=3.5, nu=4.5))
        action.perform_action(self.d)
        self.assertEqual(self.d.geometry.mu, 3.5)
        self.assertEqual(self.d.geometry.nu, 4.5)

        action = ApplyParameterOverwrites(ExperimentConfig(muOffset=0.25))
        action.perform_action(self.d)
        self.assertEqual(self.d.geometry.mu, 3.75)

        action = ApplyParameterOverwrites(ExperimentConfig(sampleModel='air | Si | Fe'))
        action.update_header(self.header)
        self.assertIsNotNone(self.header.sample.model)

    def test_apply_sample_model_file(self):
        if os.path.isfile('test.yaml'):
            os.remove('test.yaml')
        action = ApplyParameterOverwrites(ExperimentConfig(sampleModel='test.yaml'))
        action.update_header(self.header)
        self.assertIsNone(self.header.sample.model)

        with open('test.yaml', 'w') as fh:
            fh.write("""stack: air | Si | Fe""")

        try:
            action = ApplyParameterOverwrites(ExperimentConfig(sampleModel='test.yaml'))
            action.update_header(self.header)
            self.assertEqual(self.header.sample.model.stack, 'air | Si | Fe')
        finally:
            os.remove('test.yaml')

    def test_tof_time_correction(self):
        action = TofTimeCorrection()
        with self.assertRaises(ValueError):
            action.perform_action(self.d)

        new_events = append_fields(self.d.data.events, [('delta', np.float64)])
        new_events.delta = 10.0
        self.d.data.events = new_events
        tof_before = self.d.data.events.tof.copy()
        action.perform_action(self.d)
        np.testing.assert_allclose(
            self.d.data.events.tof,
            tof_before - (10.0 / 180.0) * self.d.timing.tau
        )

        self.d.create_data()
        new_events = append_fields(self.d.data.events, [('delta', np.float64)])
        new_events.delta = 10.0
        self.d.data.events = new_events
        tof_before = self.d.data.events.tof.copy()
        action = TofTimeCorrection(correct_chopper_opening=False)
        action.perform_action(self.d)
        np.testing.assert_allclose(
            self.d.data.events.tof,
            tof_before - (self.d.geometry.kad / 180.0) * self.d.timing.tau
        )

    def test_apply_mask(self):
        self.d.data.events = self.d.data.events[:6].copy()
        self.d.data.events.mask[:] = [0, 1, 2, 3, 4, 5]

        action = ApplyMask()
        action.perform_action(self.d)
        self.assertEqual(self.d.data.events.shape[0], 1)
        self.assertEqual(self.d.data.events.mask[0], 0)

        self.d.create_data()
        self.d.data.events = self.d.data.events[:6].copy()
        self.d.data.events.mask[:] = [0, 1, 2, 3, 4, 5]
        action = ApplyMask(bitmask_filter=EVENT_BITMASKS['MonitorThreshold'])
        action.perform_action(self.d)
        np.testing.assert_array_equal(self.d.data.events.mask, np.array([0, EVENT_BITMASKS['MonitorThreshold']],
                                                                        dtype=np.int32))

    def test_merge_frames(self):
        action = MergeFrames(lamdaCut=0.0)
        action.perform_action(self.d)
        self.assertEqual(self.d.data.events.tof.shape, self.d.orig_data.events.tof.shape)
        np.testing.assert_array_compare(lambda x,y: x<=y, self.d.data.events.tof, self.d.orig_data.events.tof)
        self.assertTrue((-self.d.timing.tau<=self.d.data.events.tof).all())
        np.testing.assert_array_less(self.d.data.events.tof, self.d.timing.tau)

        action = MergeFrames(lamdaCut=2.0)
        self.d.data.events.tof = self.d.orig_data.events.tof[:]
        action.perform_action(self.d)
        tofCut = 2.0*self.d.geometry.chopperDetectorDistance/const.hdm*1e-13
        self.assertTrue((tofCut-self.d.timing.tau<=self.d.data.events.tof).all())
        self.assertTrue((self.d.data.events.tof<=tofCut+self.d.timing.tau).all())

    def test_analyze_pixel_ids(self):
        action = AnalyzePixelIDs((1000, 1001))
        action.perform_action(self.d)
        self.assertIn('detZ', self.d.data.events.dtype.names)
        self.assertIn('detXdist', self.d.data.events.dtype.names)
        self.assertIn('delta', self.d.data.events.dtype.names)
        self.assertEqual(
            np.bitwise_and(self.d.data.events.mask, EVENT_BITMASKS['yRange']).astype(bool).sum(),
            self.d.data.events.shape[0]
        )
        # TODO: maybe add a test actually checking correct detector-id resolution

    def test_calculate_wavelength(self):
        action = CalculateWavelength((3.0, 5.0))
        with self.assertRaises(ValueError):
            action.perform_action(self.d)

        new_events = append_fields(self.d.data.events, [('detXdist', np.float64)])
        new_events.detXdist = 0.0
        self.d.data.events = new_events
        action.perform_action(self.d)
        self.assertIn('lamda', self.d.data.events.dtype.names)
        flt = self.d.data.events.mask!=EVENT_BITMASKS['LamdaRange']
        # check all wavelength in range not filtered
        np.testing.assert_array_less(self.d.data.events.lamda[flt], 5.0)
        np.testing.assert_array_less(3.0, self.d.data.events.lamda[flt])
        # check all wavelength out of range filtered
        flt = self.d.data.events.mask==EVENT_BITMASKS['LamdaRange']
        self.assertTrue(((self.d.data.events.lamda[flt]<3.0)|(self.d.data.events.lamda[flt]>5.0)).all())

    def test_calculate_q(self):
        action = CalculateQ(IncidentAngle.alphaF)
        with self.assertRaises(ValueError):
            action.perform_action(self.d)

        # TODO: add checks for actual resulting values

        new_events = append_fields(self.d.data.events, [('lamda', np.float64), ('delta', np.float64)])
        new_events.lamda = 5.0
        new_events.delta = 0.0
        self.d.data.events = new_events
        action.perform_action(self.d)
        self.assertIn('qz', self.d.data.events.dtype.names)
        self.assertNotIn('qx', self.d.data.events.dtype.names)
        action.update_header(self.header)
        self.assertEqual(self.header.measurement_scheme, 'angle- and energy-dispersive')

        self.d.create_data()
        new_events = append_fields(self.d.data.events, [('lamda', np.float64), ('delta', np.float64)])
        new_events.lamda = 5.0
        new_events.delta = 0.0
        self.d.data.events = new_events
        action = CalculateQ(IncidentAngle.mu)
        action.perform_action(self.d)
        self.assertIn('qz', self.d.data.events.dtype.names)
        self.assertIn('qx', self.d.data.events.dtype.names)
        action.update_header(self.header)
        self.assertEqual(self.header.measurement_scheme, 'energy-dispersive')

        self.d.create_data()
        new_events = append_fields(self.d.data.events, [('lamda', np.float64), ('delta', np.float64)])
        new_events.lamda = 5.0
        new_events.delta = 0.0
        self.d.data.events = new_events
        action = CalculateQ(IncidentAngle.nu)
        action.perform_action(self.d)
        self.assertIn('qz', self.d.data.events.dtype.names)
        self.assertNotIn('qx', self.d.data.events.dtype.names)
        action.update_header(self.header)
        self.assertEqual(self.header.measurement_scheme, 'energy-dispersive')

    def test_filter_qz_range(self):
        action = FilterQzRange((0.1, 0.2))
        with self.assertRaises(ValueError):
            action.perform_action(self.d)

        self.d.data.events = self.d.data.events[:5].copy()
        new_events = append_fields(self.d.data.events, [('qz', np.float64)])
        new_events.qz = np.array([0.05, 0.1, 0.15, 0.2, 0.25])
        self.d.data.events = new_events
        action.perform_action(self.d)
        np.testing.assert_array_equal(
            self.d.data.events.mask,
            np.array([1, 0, 0, 0, 1], dtype=np.int32) * EVENT_BITMASKS['qRange']
        )
