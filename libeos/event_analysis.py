"""
Define an event dataformat that performs reduction actions like wavelength calculation on per-event basis.
"""
import numpy as np
import logging

from typing import Tuple

from . import const
from .event_data_types import EventDataAction, EventDatasetProtocol, EVENT_TYPE, ANA_EVENT_TYPE, FINAL_EVENT_TYPE
from .helpers import filter_project_x
from .instrument import Detector
from .options import IncidentAngle
from .header import Header


class AnalyzePixelIDs(EventDataAction):
    def __init__(self, yRange: Tuple[int, int]):
        self.yRange = yRange

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if d.events.dtype != EVENT_TYPE:
            raise ValueError("AnalyzeEventData only works on raw AmorEventData, this dataset has already been altered")
        pixelLookUp = self.resolve_pixels()
        # TODO: change numba implementation to use native pixelID type
        (detZ, detXdist, delta, mask) = filter_project_x(
                pixelLookUp, d.events.pixelID.astype(np.int64), self.yRange[0], self.yRange[1]
                )
        ana_events = np.recarray(d.events.shape, dtype=ANA_EVENT_TYPE)
        # copy old data
        for field in d.events.dtype.fields.keys():
            ana_events[field] = d.events[field]
        # add analysis per event
        ana_events.detZ = detZ
        ana_events.detXdist = detXdist
        ana_events.delta = delta
        ana_events.mask = mask
        d.events = ana_events
        dataset.geometry.delta_z = self.delta_z

    def resolve_pixels(self):
        """determine spatial coordinats and angles from pixel number"""
        nPixel = Detector.nWires * Detector.nStripes * Detector.nBlades
        pixelID = np.arange(nPixel)
        (bladeNr, bPixel) = np.divmod(pixelID, Detector.nWires * Detector.nStripes)
        (bZi, detYi)      = np.divmod(bPixel, Detector.nStripes)                     # z index on blade, y index on detector
        detZi             = bladeNr * Detector.nWires + bZi                          # z index on detector
        detX              = bZi * Detector.dX                                        # x position in detector
        # detZ              = Detector.zero - bladeNr * Detector.bladeZ - bZi * Detector.dZ      # z position on detector
        bladeAngle        = np.rad2deg( 2. * np.arcsin(0.5*Detector.bladeZ / Detector.distance) )
        delta             = (Detector.nBlades/2. - bladeNr) * bladeAngle \
                            - np.rad2deg( np.arctan(bZi*Detector.dZ / ( Detector.distance + bZi * Detector.dX) ) )
        self.delta_z      = delta[detYi==1]
        return np.vstack((detYi.T, detZi.T, detX.T, delta.T)).T

class TofTimeCorrection(EventDataAction):
    def __init__(self, correct_chopper_opening: bool = True):
        self.correct_chopper_opening = correct_chopper_opening

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if d.events.dtype != ANA_EVENT_TYPE:
            raise ValueError("TofTimeCorrection requires dataset with analyzed events, perform AnalyzedEventData first")

        if self.correct_chopper_opening:
            d.events.tof -= ( d.events.delta / 180. ) * dataset.timing.tau
        else:
            d.events.tof -= ( dataset.geometry.kad / 180. ) * dataset.timing.tau

class WavelengthAndQ(EventDataAction):
    def __init__(self, lambdaRange: Tuple[float, float], incidentAngle: IncidentAngle):
        self.lambdaRange = lambdaRange
        self.incidentAngle = incidentAngle

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if d.events.dtype != ANA_EVENT_TYPE:
            raise ValueError("WavelengthAndQ requires dataset with analyzed events, perform AnalyzedEventData first")

        self.lamdaMax = const.lamdaCut+1.e13*dataset.timing.tau*const.hdm/(dataset.geometry.chopperDetectorDistance+124.)

        # lambda
        lamda = (1.e13*const.hdm)*d.events.tof/(dataset.geometry.chopperDetectorDistance+d.events.detXdist)

        final_events = np.recarray(d.events.shape, dtype=FINAL_EVENT_TYPE)
        # copy old data
        for field in d.events.dtype.fields.keys():
            final_events[field] = d.events[field]
        # add analysis per event
        final_events.lamda = lamda
        final_events.mask &= (self.lambdaRange[0]<=lamda) & (lamda<=self.lambdaRange[1])

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
        if d.events.dtype != FINAL_EVENT_TYPE:
            raise ValueError("FilterQzRange requires dataset with fully analyzed events, perform WavelengthAndQ first")

        if self.qzRange[1]<0.5:
            d.events.mask &=  (self.qzRange[0]<=d.events.qz) & (d.events.qz<=self.qzRange[1])

class ApplyMask(EventDataAction):
    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if not 'mask' in d.events.dtype.names:
            logging.debug("ApplyMask performed on dataset without mask")
            return

        logging.info(f'      number of events: total = {d.events.shape[0]:7d}, '
                     f'filtered = {np.logical_not(d.events.mask).sum():7d}')
        d.events = d.events[d.events.mask]
