"""
Define an event dataformat that performs reduction actions like wavelength calculation on per-event basis.
With large number of events these actions can be time consuming so they use numba based functions.
"""
import numpy as np
import logging

from typing import Tuple

from . import const
from .event_data_types import EventDataAction, EventDatasetProtocol, append_fields, EVENT_BITMASKS
from .helpers import filter_project_x, merge_frames, extract_walltime
from .instrument import Detector
from .options import IncidentAngle
from .header import Header

class ExtractWalltime(EventDataAction):
    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        wallTime = extract_walltime(dataset.data.events.tof,
                                    dataset.data.packets.start_index,
                                    dataset.data.packets.time)
        logging.debug(f'        expending event stream by wallTime')
        new_events = append_fields(dataset.data.events, [('wallTime', wallTime.dtype)])
        new_events.wallTime = wallTime
        dataset.data.events = new_events

class MergeFrames(EventDataAction):
    def __init__(self, lamdaCut=None):
        self.lamdaCut=lamdaCut

    def perform_action(self, dataset: EventDatasetProtocol)->None:
        if self.lamdaCut is None:
            lamdaCut = const.lamdaCut
        else:
            lamdaCut = self.lamdaCut
        tofCut = lamdaCut*dataset.geometry.chopperDetectorDistance/const.hdm*1e-13
        total_offset = (tofCut +
                        dataset.timing.tau * (dataset.timing.ch1TriggerPhase + dataset.timing.chopperPhase/2)/180)
        dataset.data.events.tof = merge_frames(dataset.data.events.tof, tofCut, dataset.timing.tau, total_offset)


class AnalyzePixelIDs(EventDataAction):
    def __init__(self, yRange: Tuple[int, int]):
        self.yRange = yRange

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        (detZ, detXdist, delta, mask) = filter_project_x(
                Detector.pixelLookUp, d.events.pixelID, self.yRange[0], self.yRange[1]
                )
        ana_events = append_fields(d.events, [
            ('detZ', detZ.dtype), ('detXdist', detXdist.dtype), ('delta', delta.dtype)])
        # add analysis per event
        ana_events.detZ = detZ
        ana_events.detXdist = detXdist
        ana_events.delta = delta
        ana_events.mask += np.logical_not(mask)*EVENT_BITMASKS['yRange']
        d.events = ana_events

class CalculateWavelength(EventDataAction):
    def __init__(self, lambdaRange: Tuple[float, float]):
        self.lambdaRange = lambdaRange

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if not 'detXdist' in dataset.data.events.dtype.names:
            raise ValueError("CalculateWavelength requires dataset with analyzed pixels, perform AnalyzePixelIDs first")

        #lamdaMax = const.lamdaCut+1.e13*dataset.timing.tau*const.hdm/(dataset.geometry.chopperDetectorDistance+124.)

        # lambda
        # TODO: one of the most time consuming actions, could be implemented in numba, instead?
        lamda = (1.e13*const.hdm)*d.events.tof/(dataset.geometry.chopperDetectorDistance+d.events.detXdist)

        final_events = append_fields(d.events, [('lamda', np.float64)])
        # add analysis per event
        final_events.lamda = lamda
        final_events.mask += EVENT_BITMASKS["LamdaRange"]*(
                (self.lambdaRange[0]>lamda) | (lamda>self.lambdaRange[1]))
        d.events = final_events

class CalculateQ(EventDataAction):
    def __init__(self, incidentAngle: IncidentAngle):
        self.incidentAngle = incidentAngle

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if not 'lamda' in dataset.data.events.dtype.names:
            raise ValueError("CalculateQ requires dataset with analyzed wavelength, perform CalculateWavelength first")

        lamda = d.events.lamda

        final_events = append_fields(d.events, [('qz', np.float64)])

        # alpha_f
        # q_z
        if self.incidentAngle == IncidentAngle.alphaF:
            alphaF_e  = dataset.geometry.nu - dataset.geometry.mu + d.events.delta
            final_events.qz = 4*np.pi*(np.sin(np.deg2rad(alphaF_e))/lamda)
        elif self.incidentAngle == IncidentAngle.nu:
            alphaF_e  = (dataset.geometry.nu + d.events.delta + dataset.geometry.kap + dataset.geometry.kad) / 2.
            final_events.qz = 4*np.pi*(np.sin(np.deg2rad(alphaF_e))/lamda)
        else:
            alphaF_e  = dataset.geometry.nu - dataset.geometry.mu + d.events.delta
            alphaI    = dataset.geometry.kap + dataset.geometry.kad + dataset.geometry.mu
            final_events.qz = 2*np.pi * ((np.sin(np.deg2rad(alphaF_e)) + np.sin(np.deg2rad(alphaI)))/lamda)
            final_events = append_fields(final_events, [('qx', np.float64)])
            final_events.qx = 2*np.pi * ((np.cos(np.deg2rad(alphaF_e)) - np.cos(np.deg2rad(alphaI)))/lamda)

        dataset.data.events = final_events

    def update_header(self, header: Header):
        if self.incidentAngle == IncidentAngle.alphaF:
            header.measurement_scheme = 'angle- and energy-dispersive'
        else:
            header.measurement_scheme = 'energy-dispersive'

class FilterQzRange(EventDataAction):
    def __init__(self, qzRange: Tuple[float, float]):
        self.qzRange = qzRange

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if not 'qz' in dataset.data.events.dtype.names:
            raise ValueError("FilterQzRange requires dataset with qz values per events, perform WavelengthAndQ first")

        if self.qzRange[1]<0.5:
            d.events.mask += EVENT_BITMASKS["qRange"]*((self.qzRange[0]>d.events.qz) | (d.events.qz>self.qzRange[1]))
