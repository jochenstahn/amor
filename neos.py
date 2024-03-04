#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
eos reduces measurements performed on Amor@SINQ, PSI

Author: Jochen Stahn

conventions (not strictly followed, yet):
- array names end with the suffix '_x[y]' with the meaning
    _e  = events
    _tof
    _l  = lambda
    _t  = theta
    _z  = detector z
    _lz = (lambda, detector z)
    _q  = q_z
- to come
"""

__version__ = '2.0'
__date__    = '2024-03-01'

import os
import sys
import subprocess
import logging
import argparse
import numpy as np
from datetime import datetime
import h5py
from orsopy import fileio
import platform
#=====================================================================================================
# TODO:
# - calculate resolution using the chopperPhase
# - deal with background correction
# - format of 'call' + add '-Y' if not supplied
#=====================================================================================================
def commandLineArgs():
    '''
    Process command line argument. 
    The type of the default values is used for conversion and validation.
    '''
    msg = "eos reads data from (one or several) raw file(s) of the .hdf format, \
           performs various corrections, conversations and projections and exports\
           the resulting reflectivity in an orso-compatible format."
    clas = argparse.ArgumentParser(description = msg)

    input_data = clas.add_argument_group('input data')
    input_data.add_argument("-n", "--fileIdentifier",       
                            default = ['0'],                               
                            nargs = '+',
                            help = "file number(s) or offset (if negative)")
    input_data.add_argument("-r", "--normalisationFileIdentifier",          
                            default = [],                               
                            nargs = '+',
                            help = "file number(s) of normalisation measurement")
    input_data.add_argument("-d", "--dataPath",  
                            type = str,
                            default = '.',
                            help = "relative path to directory with .hdf files")
    input_data.add_argument("-Y", "--year",
                            default = datetime.year,
                            type = int,
                            help = "year the measurement was performed")
    input_data.add_argument("-sub", "--subtract",
                            help = "R(q_z) curve to be subtracted (in .Rqz.ort format)")

    output = clas.add_argument_group('output')
    output.add_argument("-o", "--outputName",
                            default = "fromEOS",
                            help = "output file name (withot suffix)")
    output.add_argument("-of", "--outputFormat",
                            nargs = '+',
                            default = ['Rqz.ort'])
    output.add_argument("--offSpecular",
                            type = bool,
                            default = False,
                            )
    output.add_argument("-a", "--qResolution",
                            default = 0.01,
                            type = float,
                            help = "q_z resolution")
    output.add_argument("-ts", "--timeSlize",
                            nargs = '+',
                            type = float,
                            help = "time slizing <interval> ,[<start> [,stop]]")
    output.add_argument("-s", "--scale",
                            nargs = '+',
                            default = [1],
                            type = float,
                            help = "scaling factor for R(q_z)")
    output.add_argument("-S", "--autoscale",
                            nargs = 2,
                            type = float,
                            help = "scale to 1 in the given q_z range")

    masks = clas.add_argument_group('masks')
    masks.add_argument("-l", "--lambdaRange",                            
                            default = [2., 15.],   
                            nargs = 2, 
                            type = float,
                            help = "wavelength range")
    masks.add_argument("-t", "--thetaRange",      
                            default = [-12., 12.],   
                            nargs = 2, 
                            type = float,
                            help = "absolute theta range")
    masks.add_argument("-T", "--thetaRangeR",
                            default = [-12., 12.],
                            nargs = 2,
                            type = float,
                            help = "relative theta range")
    masks.add_argument("-y", "--yRange",          
                            default = [11, 41],
                            nargs = 2, 
                            type = int,  
                            help = "detector y range")
    masks.add_argument("-q", "--qzRange",          
                            default = [0.005, 0.30], 
                            nargs = 2, 
                            type = float,
                            help = "q_z range")

    overwrite = clas.add_argument_group('overwrite')
    overwrite.add_argument("-cs", "--chopperSpeed",
                            type = float,
                            help = "chopper speed in rpm")
    overwrite.add_argument("-cp", "--chopperPhase",    
                            default = -13.5,                     
                            type = float,
                            help = "chopper phase")
    overwrite.add_argument("-co", "--chopperPhaseOffset",    
                            default = -5,                     
                            type = float,
                            help = "phase offset between chopper opening and trigger pulse")
    overwrite.add_argument("-m", "--muOffset",        
                            default = 0.,                     
                            type = float,
                            help = "mu offset")
    overwrite.add_argument("-mu", "--mu",                                              
                            default = 0,
                            type = float,
                            help ="value of mu")
    overwrite.add_argument("-nu", "--nu",                                              
                            default = 0,
                            type = float,
                            help = "value of nu")
    overwrite.add_argument("-sm", "--sampleModel",
                            type = str,
                            help = "1-line orso sample model description")
                            
    return clas.parse_args()
#=====================================================================================================
class Defs:
    '''definition of a series of fixed parameters and constants'''
    hdm       = 6.626176e-34/1.674928e-27               # h / m
    lamdaCut  = 2.5                                     # Aa
#=====================================================================================================
class Header:
    '''orso compatible output file header content'''
    
    def __init__(self):
        self.owner                           = None
        self.experiment                      = None
        self.sample                          = None 
        self.measurement_instrument_settings = None
        self.measurement_scheme              = None
        self.measurement_data_files          = []
        self.measurement_additional_files    = []

        self.reduction = fileio.Reduction(
            software    = fileio.Software('eos', version=__version__),
            call        = ' '.join(sys.argv),
            computer    = platform.node(),
            timestamp   = datetime.now(),
            creator     = None, 
            corrections = ['histogramming in lambda and alpha_f', 
                           'gravity'],
            )
    #-------------------------------------------------------------------------------------------------
    def data_source(self):
        return fileio.DataSource(
            self.owner,
            self.experiment,
            self.sample,
            fileio.Measurement(
                instrument_settings = self.measurement_instrument_settings,
                scheme              = self.measurement_scheme,
                data_files          = self.measurement_data_files,
                additional_files    = self.measurement_additional_files,
                ),
        )
    #-------------------------------------------------------------------------------------------------
    def columns(self):
        cols = [
            fileio.Column('Qz', '1/angstrom', 'normal momentum transfer'),
            fileio.Column('R', '', 'specular reflectivity'),
            fileio.ErrorColumn(error_of='R', error_type='uncertainty', distribution='gaussian', value_is='sigma'),
            fileio.ErrorColumn(error_of='Qz', error_type='resolution', distribution='gaussian', value_is='sigma'),
            ]
        return cols

#=====================================================================================================
class AmorData:
    '''read meta-data and event streams from .hdf file(s), apply filters and conversions'''
    #-------------------------------------------------------------------------------------------------
    def __init__(self):
        global startTime
        self.startTime = startTime 
    #-------------------------------------------------------------------------------------------------
    def read_data(self, short_notation, norm=False):
        self.data_file_numbers = self.expand_file_list(short_notation)
        #self.year = clas.year
        self.file_list = []
        for number in self.data_file_numbers:
            self.file_list.append(self.path_generator(number))  
        ## read specific meta data and measurement from first file
        if norm: 
            self.readHeaderInfo = False
        else:
            self.readHeaderInfo = True

        _detZ_e            = []
        _lamda_e           = []
        _wallTime_e        = []
        for file in self.file_list:
            self.read_individual_data(file, norm)
            _detZ_e        = np.append(_detZ_e,        self.detZ_e)
            _lamda_e       = np.append(_lamda_e,       self.lamda_e)
            _wallTime_e    = np.append(_wallTime_e,    self.wallTime_e)
        self.detZ_e        = _detZ_e
        self.lamda_e       = _lamda_e
        self.wallTime_e    = _wallTime_e
             
    #-------------------------------------------------------------------------------------------------
    def path_generator(self, number):
        fileName = f'amor{clas.year}n{number:06d}.hdf'
        if   os.path.exists(f'{clas.dataPath}/{fileName}'):
            path = clas.dataPath
        elif os.path.exists(fileName):
            path = '.'
        elif os.path.exists(f'./raw/{fileName}'):
            path = './raw'
        elif os.path.exists(f'../raw/{fileName}'):
            path = '../raw'
        elif os.path.exists(f'/afs/psi.ch/project/sinqdata/{clas.year}/amor/{int(number/1000)}/{fileName}'):
            path = '/afs/psi.ch/project/sinqdata/{clas.year}/amor/{int(number/1000)}'
        else:
            sys.exit(f'# ERROR: the file {fileName} is nowhere to be found!')
        return f'{path}/{fileName}'
    #-------------------------------------------------------------------------------------------------
    def expand_file_list(self, short_notation):
        '''Evaluate string entry for file number lists'''
        #log().debug('Executing get_flist')
        file_list=[]
        for i in short_notation.split(','):
            if '-' in i:
                if ':' in i:
                    step = i.split(':', 1)[1]
                    file_list += range(int(i.split('-', 1)[0]), int((i.rsplit('-', 1)[1]).split(':', 1)[0])+1, int(step))
                else:
                    step = 1
                    file_list += range(int(i.split('-', 1)[0]), int(i.split('-', 1)[1])+1, int(step))
            else:
                file_list += [int(i)]   
        return sorted(file_list)
    #-------------------------------------------------------------------------------------------------
    def resolve_pixels(self):
        '''determine spatial coordinats and angles from pixel number'''
        det = Detector()
        nPixel = det.nWires * det.nStripes * det.nBlades
        pixelID = np.arange(nPixel)
        (bladeNr, bPixel) = np.divmod(pixelID, det.nWires * det.nStripes)
        (bZi, detYi)      = np.divmod(bPixel, det.nStripes)                     # z index on blade, y index on detector
        detZi             = bladeNr * det.nWires + bZi                          # z index on detector 
        detX              = bZi * det.dX                                        # x position in detector
        # detZ              = det.zero - bladeNr * det.bladeZ - bZi * det.dZ      # z position on detector
        bladeAngle        = np.rad2deg( 2. * np.arcsin(0.5*det.bladeZ / det.distance) )
        delta             = (det.nBlades/2. - bladeNr) * bladeAngle \
                            - np.rad2deg( np.arctan(bZi*det.dZ / ( det.distance + bZi * det.dX) ) )
        self.delta_z      = delta[detYi==1]
        return np.vstack((detYi.T, detZi.T, detX.T, delta.T)).T
        #return matr
    #-------------------------------------------------------------------------------------------------
    def read_individual_data(self, fileName, norm=False):   
        defs = Defs()
        pixelLookUp = self.resolve_pixels()
        
        self.hdf = h5py.File(fileName, 'r', swmr=True)
        
        if self.readHeaderInfo:
            # read general information and first data set
            print(f'#   meta data from: {self.file_list[0]}')
            self.hdf = h5py.File(self.file_list[0], 'r', swmr=True)

            title            = self.hdf['entry1/title'][0].decode('utf-8')
            proposal_id      = self.hdf['entry1/proposal_id'][0].decode('utf-8')
            user_name        = self.hdf['entry1/user/name'][0].decode('utf-8')
            user_affiliation = 'unknown'
            user_email       = self.hdf['entry1/user/email'][0].decode('utf-8')
            user_orcid       = None
            sampleName       = self.hdf['entry1/sample/name'][0].decode('utf-8')
            model            = self.hdf['entry1/sample/model'][0].decode('utf-8')
            instrumentName   = 'Amor'
            source           = self.hdf['entry1/Amor/source/name'][0].decode('utf-8')
            sourceProbe      = 'neutron'
            
            start_time       = self.hdf['entry1/start_time'][0].decode('utf-8')
            start_date       = start_time.split(' ')[0]
               
            if clas.sampleModel:
                model = clas.sampleModel
            else:
                model = None

            # assembling orso header information
            header.owner = fileio.Person(
                name = user_name,
                affiliation = user_affiliation, 
                contact             = user_email,
                )
            if user_orcid:
                header.owner.orcid  = user_orcid
            header.experiment = fileio.Experiment(
                title               = title, 
                instrument          = instrumentName,
                start_date          = start_date, 
                probe               = sourceProbe, 
                facility            = source,
                proposalID          = proposal_id
                )
            header.sample = fileio.Sample(
                name                = sampleName,
                model               = model,
                sample_parameters   = None,
                )
            header.measurement_scheme     = 'angle- and energy-dispersive'

        self.chopperDistance         = float(np.take(self.hdf['entry1/Amor/chopper/pair_separation'], 0))
        self.detectorDistance        = float(np.take(self.hdf['entry1/Amor/detector/transformation/distance'], 0))
        self.chopperDetectorDistance = self.detectorDistance - float(np.take(self.hdf['entry1/Amor/chopper/distance'], 0))
        self.tofCut                  = defs.lamdaCut * self.chopperDetectorDistance / defs.hdm * 1.e-13
    
        print(f'#   data from file: {fileName}')
        try:
            self.mu   = float(np.take(self.hdf['/entry1/Amor/master_parameters/mu/value'], 0))
            self.nu   = float(np.take(self.hdf['/entry1/Amor/master_parameters/nu/value'], 0))
            self.kap  = float(np.take(self.hdf['/entry1/Amor/master_parameters/kap/value'], 0))
            self.kad  = float(np.take(self.hdf['/entry1/Amor/master_parameters/kad/value'], 0))
            self.div  = float(np.take(self.hdf['/entry1/Amor/master_parameters/div/value'], 0))
            self.chopperSpeed = float(np.take(self.hdf['/entry1/Amor/chopper/rotation_speed/value'], 0))
            self.chopperPhase = float(np.take(self.hdf['/entry1/Amor/chopper/phase/value'], 0))
        except(KeyError, IndexError):
            logging.warning("     using parameters from nicos cache")
            year_date = str(self.start_date).replace('-', '/', 1) 
            cachePath = '/home/amor/nicosdata/amor/cache/'
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-mu/{year_date}')).split('\t')[-1]
            self.mu = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-nu/{year_date}')).split('\t')[-1]
            self.nu = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-kap/{year_date}')).split('\t')[-1]
            self.kap = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-kad/{year_date}')).split('\t')[-1]
            self.kad = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-div/{year_date}')).split('\t')[-1]
            self.div = float(value)
            value = str(subprocess.getoutput(f'/usr/bin/grep "value" {cachePath}nicos-ch1_speed/{year_date}')).split('\t')[-1]
            self.chopperSpeed = float(value)
            self.chopperPhase = clas.chopperPhase
        self.tau     = 30. / self.chopperSpeed

        if clas.muOffset:
            self.mu += clas.muOffset
        if clas.mu:
            self.mu = clas.mu
        if clas.nu:
            self.mu = clas.nu

        # TODO:  figure out real stop time....
        self.ctime=(self.hdf['/entry1/Amor/detector/data/event_time_zero'][-1] 
            - self.hdf['/entry1/Amor/detector/data/event_time_zero'][0]) / 1.e9
        fileDate = datetime.fromisoformat( self.hdf['/entry1/start_time'][0].decode('utf-8') )
        
        # add header content
        if self.readHeaderInfo:
            self.readHeaderInfo = False
            header.measurement_instrument_settings = fileio.InstrumentSettings(
                incident_angle = fileio.ValueRange(round(self.mu+self.kap+self.kad-0.5*self.div, 3),
                                                   round(self.mu+self.kap+self.kad+0.5*self.div, 3),
                                                   'deg'), 
                wavelength = fileio.ValueRange(defs.lamdaCut, clas.lambdaRange[1], 'angstrom'),
                polarization = fileio.Polarization.unpolarized,
                )
            header.measurement_instrument_settings.mu = fileio.Value(round(self.mu, 3), 'deg', comment='sample angle to horizon')
            header.measurement_instrument_settings.nu = fileio.Value(round(self.nu, 3), 'deg', comment='detector angle to horizon')
        if norm: 
            header.measurement_additional_files.append(fileio.File(file=fileName.split('/')[-1], timestamp=fileDate))
        else:
            header.measurement_data_files.append(fileio.File(file=fileName.split('/')[-1], timestamp=fileDate))
        print(f'#     mu = {self.mu:6.3f}, nu = {self.nu:6.3f}, kap = {self.kap:6.3f}, kad = {self.kap:6.3f}')
        
        # TODO: should extract monitor from counts or beam current times time 
        self.monitor1 = self.ctime
        self.monitor2=self.monitor1
        
        # read data event streams
        tof_e        = np.array(self.hdf['/entry1/Amor/detector/data/event_time_offset'][:])/1.e9
        pixelID_e    = np.array(self.hdf['/entry1/Amor/detector/data/event_id'][:], dtype=int)
        totalNumber = np.shape(tof_e)[0]

        wallTime_e = np.empty(totalNumber)
        dataPacket_p = np.array(self.hdf['/entry1/Amor/detector/data/event_index'][:], dtype=np.uint64)  
        dataPacketTime_p = np.array(self.hdf['/entry1/Amor/detector/data/event_time_zero'][:], dtype=np.uint64)/1e9
        #for i, index in enumerate(dataPacket_p):
        #    wallTime_e[index:] = dataPacketTime_p[i]         
        for i in range(len(dataPacket_p)-1):
            wallTime_e[dataPacket_p[i]:dataPacket_p[i+1]] = dataPacketTime_p[i]
        wallTime_e[dataPacket_p[-1]:] = dataPacketTime_p[-1]
        if not self.startTime and not norm:
            self.startTime = wallTime_e[0]
        wallTime_e -= self.startTime
        #print(f'wall time from {wallTime_e[0]} to {wallTime_e[-1]}')
 
        # filter 'strange' tof times > 2 tau
        if True:
            filter_e   = (tof_e <= 2*self.tau)
            tof_e      = tof_e[filter_e]
            pixelID_e  = pixelID_e[filter_e]
            wallTime_e = wallTime_e[filter_e]
            if np.shape(filter_e)[0]-np.shape(tof_e)[0] > 0.5 : 
                print(f'#    strange times: {np.shape(filter_e)[0]-np.shape(tof_e)[0]}')
        tof_e        = np.remainder( tof_e - self.tofCut + self.tau, self.tau) + self.tofCut  # tof shifted to 1 frame
        tof_e        = tof_e + self.tau * clas.chopperPhaseOffset / 180. # correction for time offset between chopper pulse and tof zero

        # resolve pixel ID into y and z indicees, x position and angle
        (detY_e, detZ_e, detXdist_e, delta_e) = pixelLookUp[np.int_(pixelID_e)-1,:].T
        
        # define mask and filter y range
        mask_e = (clas.yRange[0] <= detY_e) & (detY_e <= clas.yRange[1])
        
        # correct tof for beam size effect at chopper:  t_cor = (delta / 180 deg) * tau 
        # TODO: check for correctness
        if not clas.offSpecular:
            tof_e    -= ( delta_e / 180. ) * self.tau  
                  
        # lambda
        lamda_e = 1.e13 * tof_e * defs.hdm / (self.chopperDetectorDistance + detXdist_e)
        self.lamdaMax = defs.lamdaCut + 1.e13 * self.tau * defs.hdm / (self.chopperDetectorDistance + 124.)
        mask_e = np.logical_and(mask_e, (clas.lambdaRange[0] <= lamda_e) & (lamda_e <= clas.lambdaRange[1]))
        
        # alpha_f  
        alphaF_e = self.nu - self.mu + delta_e 
        
        # q_z
        if clas.offSpecular:
            alphaI = self.kap + self.kad + self.mu
            qz_e = 2 * np.pi * ( np.sin( np.deg2rad(alphaF_e) ) + np.sin( np.deg2rad( alphaI ) ) ) / lamda_e
            qx_e = 2 * np.pi * ( np.cos( np.deg2rad(alphaF_e) ) - np.cos( np.deg2rad( alphaI ) ) ) / lamda_e 
            header.measurement_scheme     = 'energy-dispersive',
        else:
            qz_e = 4 * np.pi * np.sin( np.deg2rad(alphaF_e) ) / lamda_e
            # qx_e = 0.
            header.measurement_scheme     = 'angle- and energy-dispersive'
        
        # filter q_z range
        if clas.qzRange[1] < 0.3 and not norm:
            mask_e = np.logical_and(mask_e, (clas.qzRange[0] <= qz_e) & (qz_e <= clas.qzRange[1]))
        
        self.detZ_e        = detZ_e[mask_e]
        self.lamda_e       = lamda_e[mask_e]
        self.wallTime_e    = wallTime_e[mask_e]

        print(f'#     number of events: total = {totalNumber:7d}, filtered = {np.shape(self.lamda_e)[0]:7d}')
#=====================================================================================================
class Detector:
    def __init__(self):
        self.nBlades  = 14                                    #     number of active blades in the detector 
        self.nWires   = 32                                    #     number of wires per blade
        self.nStripes = 64                                    #     number of stipes per blade
        self.angle    = np.deg2rad(5.1)                       # deg  angle of incidence of the beam on the blades (def: 5.1)
        self.dZ       = 4.0 * np.sin(self.angle)              # mm  height-distance of neighboring pixels on one blade
        self.dX       = 4.0 * np.cos(self.angle)              # mm  depth-distance of neighboring pixels on one blace
        self.bladeZ   = 10.455                                # mm  distance between detector blades 
        self.zero     = 0.5 * self.nBlades * self.bladeZ      # mm  vertical center of the detector
        self.distance = 4000.                                 # mm  distance from focal point to leading blade edge 
#=====================================================================================================
def expand_file_list(short_notation):
    '''Evaluate string entry for file number lists'''
    #log().debug('Executing get_flist')
    file_list=[]
    for i in short_notation.split(','):
        if '-' in i:
            if ':' in i:
                step = i.split(':', 1)[1]
                file_list += range(int(i.split('-', 1)[0]), int((i.rsplit('-', 1)[1]).split(':', 1)[0])+1, int(step))
            else:
                step = 1
                file_list += range(int(i.split('-', 1)[0]), int(i.split('-', 1)[1])+1, int(step))
        else:
            file_list += [int(i)]
        
    return sorted(file_list)
#=====================================================================================================
def normalisation_map(short_notation):
    fromHDF = AmorData()
    normalisation_list = expand_file_list(short_notation)
    name = str(normalisation_list[0])
    for i in range(1, len(normalisation_list), 1):
        name = f'{name}_{normalisation_list[i]}'
    if os.path.exists(f'{clas.dataPath}/{name}.norm'):
        print(f'# normalisation matrix: found and using {clas.dataPath}/{name}.norm')
        norm_lz = np.loadtxt(f'{clas.dataPath}/{name}.norm')
        fh = open(f'{clas.dataPath}/{name}.norm', 'r')
        fh.readline()
        normFileList = fh.readline().split('[')[1].split(']')[0].replace('\'', '').split(', ')
        normAngle = float(fh.readline().split('= ')[1])
        fh.close()
        for i, entry in enumerate(normFileList):
             normFileList[i] = entry.split('/')[-1]
        header.measurement_additional_files = normFileList
    else:
        print(f'# normalisation matrix: using the files {normalisation_list}')
        fromHDF.read_data(short_notation, norm=True)
        normAngle     = fromHDF.nu - fromHDF.mu 
        lamda_e       = fromHDF.lamda_e
        detZ_e        = fromHDF.detZ_e
        norm_lz, bins_l, bins_z = np.histogram2d(lamda_e, detZ_e, bins = (grid.lamda(), grid.z()))
        norm_lz = np.where(norm_lz>0, norm_lz, np.nan)
        # correct for the SM reflectivity
        lamda_l  = grid.lamda()
        theta_z  = normAngle + fromHDF.delta_z
        lamda_lz = (grid.lz().T*lamda_l[:-1]).T
        theta_lz = grid.lz()*theta_z 
        qz_lz    = 4.0*np.pi * np.sin(np.deg2rad(theta_lz)) / lamda_lz
        Rsm_lz   = np.ones(np.shape(qz_lz))
        Rsm_lz   = np.where(qz_lz>0.0217, 1-(qz_lz-0.0217)*(0.0625/0.0217), Rsm_lz)
        Rsm_lz   = np.where(qz_lz>0.0217*5, np.nan, Rsm_lz)
        norm_lz  = norm_lz / Rsm_lz
        if len(lamda_e) > 1e6:
            head = ('normalisation matrix based on the measurements\n'
                   f'{fromHDF.file_list}\n'
                   f'nu - mu = {normAngle}\n'
                   f'shape= {np.shape(norm_lz)} (lambda, z)\n'
                   f'measured at mu = {fromHDF.mu:6.3f} deg\n'
                   f'N(l_lambda, z) = theta(z) / sum_i=-1..1 I(l_lambda+i, z)')
            head = head.replace('../', '')
            head = head.replace('./', '')
            head = head.replace('raw/', '')
            np.savetxt(f'{clas.dataPath}/{name}.norm', norm_lz, header = head)
        normFileList = fromHDF.file_list
    return norm_lz, normAngle, normFileList
#=====================================================================================================
def output_format_list(outputFormat):
    format_list = []
    if 'ort' in outputFormat or 'Rqz.ort' in outputFormat or 'Rqz' in outputFormat:
        format_list.append('Rqz.ort')
    if 'ort' in outputFormat or 'Rlt.ort' in outputFormat or 'Rlt' in outputFormat:
        format_list.append('Rlt.ort')
    if 'orb' in outputFormat or 'Rqz.orb' in outputFormat or 'Rqz' in outputFormat:
        format_list.append('Rqz.orb')
    if 'orb' in outputFormat or 'Rlt.orb' in outputFormat or 'Rlt' in outputFormat:
        format_list.append('Rlt.orb')

    return sorted(format_list, reverse=True)
#=====================================================================================================
def project_on_lz(fromHDF, norm_lz, normAngle, lamda_e, detZ_e):
    # projection on lambda-z-grid
    lamda_l  = grid.lamda()
    theta_z  = fromHDF.nu - fromHDF.mu + fromHDF.delta_z
    lamda_lz = (grid.lz().T*lamda_l[:-1]).T
    theta_lz = grid.lz()*theta_z 
       
    thetaN_z  = fromHDF.delta_z + normAngle
    thetaN_lz = np.ones(np.shape(norm_lz))*thetaN_z
    thetaN_lz = np.where(np.absolute(thetaN_lz)>5e-3, thetaN_lz, np.nan)

    mask_lz   = np.where(np.isnan(norm_lz), False, True)
    mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(thetaN_lz)>5e-3, True, False))
    mask_lz   = np.logical_and(mask_lz, np.where(np.absolute(theta_lz)>5e-3, True, False))
    if clas.thetaRange[1]<12:
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz >= clas.thetaRange[0], True, False))
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz <= clas.thetaRange[1], True, False))
    if clas.thetaRangeR[1]<12:
      t0 = fromHDF.nu - fromHDF.mu
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz-t0 >= clas.thetaRangeR[0], True, False))
      mask_lz   = np.logical_and(mask_lz, np.where(theta_lz-t0 <= clas.thetaRangeR[1], True, False))
    if clas.lambdaRange[1]<15:
      mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz >= clas.lambdaRange[0], True, False))
      mask_lz   = np.logical_and(mask_lz, np.where(lamda_lz <= clas.lambdaRange[1], True, False))

    #           gravity correction
    #theta_lz += np.rad2deg( np.arctan( 3.07e-10 * (fromHDF.detectorDistance + detXdist_e) * lamda_lz**2 ) )
    theta_lz += np.rad2deg( np.arctan( 3.07e-10 * fromHDF.detectorDistance * lamda_lz**2 ) ) 

    z_z       = enumerate(theta_z)
    qz_lz     = 4.0*np.pi * np.sin(np.deg2rad(theta_lz)) / lamda_lz
    int_lz, bins_l, bins_z  = np.histogram2d(lamda_e, detZ_e, bins = (lamda_l, grid.z()))
    #           cut normalisation sample horizon
    int_lz    = np.where(mask_lz, int_lz, np.nan) 
    thetaF_lz  = np.where(mask_lz, theta_lz, np.nan)

    ref_lz    = (int_lz * np.absolute(thetaN_lz)) / (norm_lz * np.absolute(thetaF_lz))
    err_lz    = ref_lz * np.sqrt( 1/(int_lz+.1) + 1/norm_lz )
        
    res_lz    = np.ones((np.shape(lamda_l[:-1])[0], np.shape(theta_z)[0])) * 0.022**2
    res_lz    = res_lz + (0.008/theta_lz)**2
    res_lz    = qz_lz * np.sqrt(res_lz)

    return qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz
#=====================================================================================================
def project_on_qz(q_lz, R_lz, dR_lz, dq_lz, norm_lz, mask_lz):
    q_q       = grid.q()
    mask_lzf  = mask_lz.flatten()
    q_lzf     = q_lz.flatten()[mask_lzf]
    R_lzf     = R_lz.flatten()[mask_lzf]
    dR_lzf    = dR_lz.flatten()[mask_lzf]
    dq_lzf    = dq_lz.flatten()[mask_lzf]
    norm_lzf  = norm_lz.flatten()[mask_lzf]
        
    N_q       = np.histogram(q_lzf, bins = q_q, weights = norm_lzf )[0]
    N_q       = np.where(N_q > 0, N_q, np.nan)

    R_q       = np.histogram(q_lzf, bins = q_q, weights = norm_lzf * R_lzf )[0]
    R_q       = R_q / N_q
 
    dR_q      = np.histogram(q_lzf, bins = q_q, weights = (norm_lzf * dR_lzf)**2 )[0]
    dR_q      = np.sqrt( dR_q ) / N_q 

    # TODO: different error propagations for dR and dq! 
    N_q       = np.histogram(q_lzf, bins = q_q, weights = norm_lzf**2 )[0]
    N_q       = np.where(N_q > 0, N_q, np.nan)
    dq_q      = np.histogram(q_lzf, bins = q_q, weights = (norm_lzf * dq_lzf)**2 )[0]
    dq_q      = np.sqrt( dq_q / N_q )

    q_q       = 0.5 * (q_q + np.roll(q_q, 1))

    return q_q[1:], R_q, dR_q, dq_q
#=====================================================================================================
def autoscale(q_q, R_q, dR_q, pR_q=[], pdR_q=[]):
    if len(pR_q) == 0:
        filter_q  = np.where((clas.autoscale[0]<=q_q)&(q_q<=clas.autoscale[1]), True, False)
        filter_q  = np.where(dR_q>0, filter_q, False)
        if len(filter_q[filter_q]) > 0:
            scale = np.sum(R_q[filter_q]**2/dR_q[filter_q]) / np.sum(R_q[filter_q]/dR_q[filter_q])
        else:
            print(f'#     automatic scaling not possible')
            scale = 1.
    else:
        filter_q  = np.where(np.isnan(pR_q*R_q), False, True)
        filter_q  = np.where(R_q>0, filter_q, False)
        filter_q  = np.where(pR_q>0, filter_q, False)
        if len(filter_q[filter_q]) > 0:
            scale = np.sum(R_q[filter_q]**3 * pR_q[filter_q] / (dR_q[filter_q]**2 * pdR_q[filter_q]**2)) \
                  / np.sum(R_q[filter_q]**2 * pR_q[filter_q]**2 / (dR_q[filter_q]**2  * pdR_q[filter_q]**2))
        else:
            print(f'#     automatic scaling not possible')
            scale = 1.
    R_q  /= scale
    dR_q /= scale
    #print(f'#     scaling factor = {scale}')
        
    return R_q, dR_q
#=====================================================================================================
def loadRqz(name):

    if os.path.exists(f'{clas.dataPath}/{name}'):
        fileName = f'{clas.dataPath}/{name}'
    elif os.path.exists(f'{clas.dataPath}/{name}.Rqz.ort'):
        fileName = f'{clas.dataPath}/{name}.Rqz.ort'
    else:
        sys.exit(f'### the background file \'{clas.dataPath}/{name}\' does not exist! => stopping')

    q_q, Sq_q, dS_q = np.loadtxt(fileName, usecols=(0, 1, 2), comments='#', unpack=True)

    return q_q, Sq_q, dS_q, fileName
#=====================================================================================================
class Grid:

    def __init__(self):
        self.det  = Detector()
        self.lamdaCut = Defs.lamdaCut
        self.dldl = 0.005     # Delta lambda / lambda

    def q(self):
        resolutions = [0.005, 0.01, 0.02, 0.025, 0.04, 0.05, 0.1, 1]
        a, b = np.histogram([clas.qResolution], bins = resolutions)
        dqdq = np.matmul(b[:-1],a)
        if dqdq != clas.qResolution:
            print(f'#   changed resolution to {dqdq}')
        qq = 0.01 
        # linear up to qq
        q_grid = np.arange(0, qq, qq*dqdq)
        # exponential from qq on
        q_grid = np.append(q_grid, qq*(1.+dqdq)**np.arange(int(np.log(0.3/qq)/np.log(1+dqdq))))
        return q_grid

    def lamda(self):
        lamdaMax = 16
        lamdaMin = self.lamdaCut
        lamda_grid = lamdaMin*(1+self.dldl)**np.arange(int(np.log(lamdaMax/lamdaMin)/np.log(1+self.dldl)+1))
        return lamda_grid

    def z(self):
        return np.arange(self.det.nBlades*self.det.nWires+1)

    def lz(self):
        return np.ones(( np.shape(self.lamda()[:-1])[0], np.shape(self.z()[:-1])[0] ))

    def delta(self, detectorDistance):
        # unused for now
        bladeAngle = np.rad2deg( 2. * np.arcsin(0.5*self.det.bladeZ / detectorDistance) )
        blade_grid = np.arctan( np.arange(33) * self.det.dZ / ( detectorDistance + np.arange(33) * self.det.dX) )
        blade_grid = np.rad2deg(blade_grid)
        stepWidth  = blade_grid[1] - blade_grid[0]
        blade_grid = blade_grid - 0.2 * stepWidth

        delta_grid = []
        for b in np.arange(self.det.nBlades-1):
            delta_grid = np.concatenate((delta_grid, blade_grid), axis=None)
            blade_grid = blade_grid + bladeAngle
            delta_grid = delta_grid[delta_grid<blade_grid[0]-0.5*stepWidth]
        delta_grid = np.concatenate((delta_grid, blade_grid), axis=None)

        return -np.flip(delta_grid) + 0.5*self.det.nBlades * bladeAngle

#=====================================================================================================
def main():
    global startTime, grid, clas, header
    clas   = commandLineArgs()
    grid   = Grid()
    header = Header()
    startTime = 0
    if not os.path.exists(f'{clas.dataPath}'):
        os.system(f'mkdir {clas.dataPath}')
    fromHDF = AmorData()
    
    print('\n######## eos - data reduction for Amor ########')
    
    # load or create normalisation matrix
    if clas.normalisationFileIdentifier:
        normalise = True
        norm_lz, normAngle, normFileList = normalisation_map(clas.normalisationFileIdentifier[0])
        header.reduction.corrections.append('normalisation with \'additional files\'')
    else:
        normalise = False
        norm_lz   = grid.lz()
        normAngle = 1.

        print('# normalisation matrix: none requested')

    # load R(q_z) curve to be subtracted:
    if clas.subtract:
        sq_q, sR_q, sdR_q, sFileName = loadRqz(clas.subtract)
        subtract = True
        print(f'# loaded background file: {sFileName}')
        header.reduction.corrections.append(f'background from \'{sFileName}\' subtracted')
    else:
        subtract = False
        
    # load measurement data and do the reduction
    datasetsRqz = []
    datasetsRlt = []
    for i, short_notation in enumerate(clas.fileIdentifier):
        print('# reading input:')
        header.measurement_data_files = []
        fromHDF.read_data(short_notation)

        if clas.timeSlize:
            wallTime_e = fromHDF.wallTime_e
            columns = header.columns() + [fileio.Column('time', 's', 'time relative to start of measurement series')]
            headerRqz = fileio.Orso(header.data_source(), header.reduction, columns)

            interval = clas.timeSlize[0]
            try:
                start = clas.timeSlize[1]
            except: 
                start = 0 
            try:
                stop  = clas.timeSlize[2]
            except:
                stop  = wallTime_e[-1] 
            for i, time in enumerate(np.arange(start, stop, interval)):
                print(f'#  time slize {i:4d}', end='\r')

                filter_e = np.where((time < wallTime_e) & (wallTime_e < time+interval), True, False)
                lamda_e  = fromHDF.lamda_e[filter_e] 
                detZ_e   = fromHDF.detZ_e[filter_e]
            
                qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz = project_on_lz(fromHDF, norm_lz, normAngle, lamda_e, detZ_e)
                q_q, R_q, dR_q, dq_q = project_on_qz(qz_lz, ref_lz, err_lz, res_lz, norm_lz, mask_lz)

                filter_q = np.where((clas.qzRange[0] < q_q) & (q_q < clas.qzRange[1]), True, False)
                q_q = q_q[filter_q]
                R_q = R_q[filter_q]
                dR_q = dR_q[filter_q]
                dq_q = dq_q[filter_q]

                if clas.autoscale:
                    R_q, dR_q = autoscale(q_q, R_q, dR_q)

                if subtract:
                    if len(q_q) == len(sq_q):
                        R_q  -= sR_q
                        dR_q = np.sqrt( dR_q**2 + sdR_q**2 )
                    else:
                        subtract = False
                        print(f'# background file {sFileName} not compatible with q_z scale ({len(sq_q)} vs. {len(q_q)})')

                tme_q              = np.ones(np.shape(q_q))*time
                data               = np.array([q_q, R_q, dR_q, dq_q, tme_q]).T
                headerRqz.data_set = f'{i}: time = {time:8.1f} s  to {time+interval:8.1f} s'
                orso_data          = fileio.OrsoDataset(headerRqz, data)
                # make a copy of the header for the next iteration
                headerRqz          = fileio.Orso.from_dict(headerRqz.to_dict())
                datasetsRqz.append(orso_data)
            print('')

        else:
            lamda_e  = fromHDF.lamda_e
            detZ_e   = fromHDF.detZ_e

            qz_lz, ref_lz, err_lz, res_lz, lamda_lz, theta_lz, int_lz, mask_lz = project_on_lz(fromHDF, norm_lz, normAngle, lamda_e, detZ_e)
            try:
                ref_lz *= clas.scale[i]
                err_lz *= clas.scale[i]
            except:
                ref_lz *= clas.scale[-1]
                err_lz *= clas.scale[-1]

            if 'Rqz.ort' in output_format_list(clas.outputFormat):
                headerRqz = fileio.Orso(header.data_source(), header.reduction, header.columns())
                headerRqz.data_set = f'Nr {i} : mu = {fromHDF.mu:6.3f} deg'

                # projection on q-grid 
                q_q, R_q, dR_q, dq_q = project_on_qz(qz_lz, ref_lz, err_lz, res_lz, norm_lz, mask_lz)

                filter_q = np.where((clas.qzRange[0] < q_q) & (q_q < clas.qzRange[1]), True, False)
                q_q = q_q[filter_q]
                R_q = R_q[filter_q]
                dR_q = dR_q[filter_q]
                dq_q = dq_q[filter_q]

                if clas.autoscale:
                    if i == 0:
                        R_q, dR_q = autoscale(q_q, R_q, dR_q)
                    else: 
                        pRq_z     = datasetsRqz[i-1].data[:,1]
                        pdRq_z    = datasetsRqz[i-1].data[:,2]
                        R_q, dR_q = autoscale(q_q, R_q, dR_q, pRq_z, pdRq_z) 

                if subtract:
                    if len(q_q) == len(sq_q):
                        R_q  -= sR_q
                        dR_q = np.sqrt( dR_q**2 + sdR_q**2 )
                    else:
                        print(f'# backgroung file {sFileName} not compatible with q_z scale ({len(sq_q)} vs. {len(q_q)})')

                data = np.array([q_q, R_q, dR_q, dq_q]).T
                orso_data       = fileio.OrsoDataset(headerRqz, data)
                headerRqz       = fileio.Orso(**headerRqz.to_dict())
                datasetsRqz.append(orso_data)
   
            if 'Rlt.ort' in output_format_list(clas.outputFormat):
                columns = [
                    fileio.Column('Qz', '1/angstrom', 'normal momentum transfer'),
                    fileio.Column('R', '', 'specular reflectivity'),
                    fileio.ErrorColumn(error_of='R', error_type='uncertainty', value_is='sigma'),
                    fileio.ErrorColumn(error_of='Qz', error_type='resolution', value_is='sigma'),
                    fileio.Column('lambda', 'angstrom', 'wavelength'),
                    fileio.Column('alpha_f', 'deg', 'final angle'),
                    fileio.Column('l', '', 'index of lambda-bin'),
                    fileio.Column('t', '', 'index of theta bin'),
                    fileio.Column('intensity', '', 'filtered neutron events per pixel'),
                    fileio.Column('norm', '', 'normalisation matrix'),
                    fileio.Column('mask', '', 'pixels used for calculating R(q_z)'),
                    ]
                #data_source = fromHDF.data_source
                headerRlt = fileio.Orso(header.data_source, header.reduction, columns)
        
                ts, zs=ref_lz.shape
                lindex_lz=np.tile(np.arange(1, ts+1), (zs, 1)).T
                tindex_lz=np.tile(np.arange(1, zs+1), (ts, 1))
  
                j = 0
                for item in zip(
                        qz_lz.T, 
                        ref_lz.T, 
                        err_lz.T, 
                        res_lz.T,
                        lamda_lz.T, 
                        theta_lz.T, 
                        lindex_lz.T, 
                        tindex_lz.T,
                        int_lz.T,
                        norm_lz.T,
                        np.where(mask_lz, 1, 0).T
                        ):
                    data               = np.array(list(item)).T
                    headerRlt          = fileio.Orso(**headerRlt.to_dict())
                    headerRlt.data_set = f'dataset_{i}_{j+1} : alpha_f = {theta_lz[0,j]:6.3f} deg'
                    orso_data          = fileio.OrsoDataset(headerRlt, data)
                    datasetsRlt.append(orso_data)
                    j += 1

    # output
    print('# output:')

    if 'Rqz.ort' in output_format_list(clas.outputFormat):
        print(f'#   {clas.dataPath}/{clas.outputName}.Rqz.ort')
        theSecondLine = f' {header.experiment.title} | {header.experiment.start_date} | sample {header.sample.name} | R(q_z)'
        fileio.save_orso(datasetsRqz, f'{clas.dataPath}/{clas.outputName}.Rqz.ort', data_separator='\n', comment=theSecondLine)

    if 'Rlt.ort' in output_format_list(clas.outputFormat):
        print(f'#   {clas.dataPath}/{clas.outputName}.Rlt.ort')
        theSecondLine = f' {header.experiment.title} | {header.experiment.start_date} | sample {header.sample.name} | R(lambda, theta)'
        fileio.save_orso(datasetsRlt, f'{clas.dataPath}/{clas.outputName}.Rlt.ort', data_separator='\n', comment=theSecondLine)

    print('')
#=====================================================================================================
if __name__ == '__main__':
    main()
