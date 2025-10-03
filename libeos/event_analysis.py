"""
Define an event dataformat that performs reduction actions like wavelength calculation on per-event basis.
"""
import numpy as np
import logging

from typing import Tuple
from numpy.lib.recfunctions import rec_append_fields

from . import const
from .event_data_types import EventDataAction, EventDatasetProtocol, append_fields, EVENT_BITMASKS
from .helpers import filter_project_x
from .instrument import Detector
from .options import IncidentAngle
from .header import Header


class AnalyzePixelIDs(EventDataAction):
    def __init__(self, yRange: Tuple[int, int]):
        self.yRange = yRange

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        delta_z, pixelLookUp = self.resolve_pixels()
        # TODO: change numba implementation to use native pixelID type
        (detZ, detXdist, delta, mask) = filter_project_x(
                pixelLookUp, d.events.pixelID.astype(np.int64), self.yRange[0], self.yRange[1]
                )
        ana_events = append_fields(d.events, [
            ('detZ', detZ.dtype), ('detXdist', detXdist.dtype), ('delta', delta.dtype)])
        # add analysis per event
        ana_events.detZ = detZ
        ana_events.detXdist = detXdist
        ana_events.delta = delta
        ana_events.mask += np.logical_not(mask)*EVENT_BITMASKS['yRange']
        d.events = ana_events
        dataset.geometry.delta_z = delta_z

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
        delta_z      = delta[detYi==1]
        pixel_lookup=np.vstack((detYi.T, detZi.T, detX.T, delta.T)).T
        return delta_z, pixel_lookup

class TofTimeCorrection(EventDataAction):
    def __init__(self, correct_chopper_opening: bool = True):
        self.correct_chopper_opening = correct_chopper_opening

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if self.correct_chopper_opening:
            d.events.tof -= ( d.events.delta / 180. ) * dataset.timing.tau
        else:
            d.events.tof -= ( dataset.geometry.kad / 180. ) * dataset.timing.tau

class CalculateWavelength(EventDataAction):
    def __init__(self, lambdaRange: Tuple[float, float]):
        self.lambdaRange = lambdaRange

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        if not 'detXdist' in dataset.data.events.dtype.names:
            raise ValueError("CalculateWavelength requires dataset with analyzed pixels, perform AnalyzePixelIDs first")

        #lamdaMax = const.lamdaCut+1.e13*dataset.timing.tau*const.hdm/(dataset.geometry.chopperDetectorDistance+124.)

        # lambda
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

class ApplyMask(EventDataAction):
    def __init__(self, bitmask_filter=None):
        self.bitmask_filter = bitmask_filter

    def perform_action(self, dataset: EventDatasetProtocol) ->None:
        d = dataset.data
        logging.info(f'      number of events: total = {d.events.shape[0]:7d}, '
                     f'filtered = {(d.events.mask!=0).sum():7d}')
        if logging.getLogger().level == logging.DEBUG:
            # only run this calculation if debug level is actually active
            filtered_by_mask = {}
            for key, value in EVENT_BITMASKS.items():
                filtered_by_mask[key] = ((d.events.mask & value)!=0).sum()
            logging.debug(f"        Removed by filters: {filtered_by_mask}")
        if self.bitmask_filter is None:
            d.events = d.events[d.events.mask==0]
        else:
            # remove the provided bitmask_filter bits from the events
            # this means that all bits that are set in bitmask_filter will NOT be used to filter events
            fltr = (d.events.mask & (~self.bitmask_filter)) == 0
            d.events = d.events[fltr]
