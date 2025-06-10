__version__ = '2025-06-07'

# essential changes with regard to 2024 version 
# - accepts new hdf structure
# TODO:
# - data path for output
# - solve confusion between negative file number and ranges (check for '-')

import os
import sys
import subprocess
import h5py
import glob
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import signal
import logging
from datetime import datetime
from orsopy import fileio

#==============================================================================
#==============================================================================
class Detector:
    def __init__(self):
        self.nBlades = 14                                    # number of active blades in the detector 
        angle        = np.deg2rad( 5.1 )                     # deg  angle of incidence of the beam on the blades (def: 5.1)
        self.dZ      =  4.0 * np.sin(angle)                  # mm  height-distance of neighboring pixels on one blade
        self.dX      =  4.0 * np.cos(angle)                  # mm  depth-distance of neighboring pixels on one blace
        self.bladeZ  = 10.7                                  # mm  distance between detector blades (consistent with nu!)
        self.zero    =  0.5 * self.nBlades * self.bladeZ     # mm  vertical center of the detector 

#==============================================================================
def pixel2quantity():
    det = Detector()
    nPixel = 64 * 32 * det.nBlades
    pixelID = np.arange(nPixel)
    (bladeNr, bPixel) = np.divmod(pixelID, 64*32)
    (bZ, bY)          = np.divmod(bPixel, 64)
    z                 = det.zero - bladeNr * det.bladeZ - bZ * det.dZ
    x                 = (31 - bZ) * det.dX 
    bladeAngle        = np.rad2deg( 2. * np.arcsin(0.5*det.bladeZ / detectorDistance) )
    delta             = (det.nBlades/2. - bladeNr) * bladeAngle - np.rad2deg( np.arctan(bZ*det.dZ / ( detectorDistance + bZ * det.dX) ) )
    quantity = np.vstack((bladeNr.T, bZ.T, bY.T, delta.T, x.T)).T
    
    return quantity

#==============================================================================
def analyse_ev(event_e, tof_e, yMin, yMax, thetaMin, thetaMax):
    
    data_e = np.zeros((len(event_e), 10), dtype=float)
    
    # data_e column description:
    #    0: wall time / s
    #    1: pixelID
    #    2: blade number
    #    3: z on blade
    #    4: y on blade
    #    5: delta / deg  = angle on detector
    #    6: path within detector / mm
    #    7: lambda / angstrom
    #    8: theta / deg
    #    9: q_z / angstrom^-1
    
    data_e[:,0]   = tof_e[:]
    data_e[:,1]   = event_e[:]

    # filter 'strange' tof times > 2 tau
    if True:
        filter_e = (data_e[:,0] <= 2*tau)
        #print(event_e[~filter_e])
        #print(data_e[~filter_e,0])
        data_e = data_e[filter_e,:]
        if np.shape(filter_e)[0]-np.shape(data_e)[0] > 0.5 and verbous: 
            logging.warning(f'##   strange times: {np.shape(filter_e)[0]-np.shape(data_e)[0]}')

    pixelLookUp = pixel2quantity()
    data_e[:,2:7] = pixelLookUp[np.int_(data_e[:,1])-1,:]

    #================================

    # filter y range
    filter_e = (yMin <= data_e[:,4]) & (data_e[:,4] <= yMax) 
    data_e = data_e[filter_e,:]
    
    # correct tof for beam size effect at chopper 
    # t_cor = (delta / 180 deg) * tau 
    data_e[:,0] -= ( data_e[:,5] / 180. ) * tau
      
    # effective flight path length
    #data_e[:,6] = chopperDetectorDistance + data_e[:,6] 
    
    # calculate lambda
    hdm       = 6.626176e-34/1.674928e-27               # h / m
    data_e[:,7] = 1.e13 * data_e[:,0] * hdm / ( chopperDetectorDistance + data_e[:,6] )
    
    # theta 
    # data_e[:,8] = nu - mu + np.rad2deg( np.arctan2(data_e[:,5], detectorDistance) ) 
    data_e[:,8] = nu - mu + data_e[:,5] 

    # gravity compensation
    data_e[:,8] += np.rad2deg( np.arctan( 3.07e-10 * ( detectorDistance +  data_e[:,6]) * data_e[:,7] * data_e[:,7] ) )

    # filter theta range
    filter_l = (thetaMin <= data_e[:,8]) & (data_e[:,8] <= thetaMax) 
    data_e = data_e[filter_l,:]
    
    # q_z
    if mu > -0.25:
        data_e[:,9] = 4*np.pi * np.sin( np.deg2rad( data_e[:,8] ) ) / data_e[:,7]
    else:
        data_e[:,9] = 4*np.pi * np.sin( np.deg2rad( -data_e[:,8] ) ) / data_e[:,7]
    
    # filter q_z range
    #filter_e = (qMin < data_e[:,I7]) & (data_e[:,7] < qMax)
    #data_e = data_e[filter_e,:]

    return data_e

#==============================================================================
class Meta:
  # AMOR hdf dataset with associated properties from metadata
    def __init__(self, filename):
        self.filename = filename

        fh = h5py.File(filename, 'r', swmr=True)

        # for processing

        self.chopperDistance            = float(np.take(fh['/entry1/Amor/chopper/pair_separation'], 0))  # mm
        # the following is the distance from the sample to the detector entry window, not to the center of rotation
        self.detectorDistance           = float(np.take(fh['/entry1/Amor/detector/transformation/distance'], 0)) # mm
        self.chopperDetectorDistance    = self.detectorDistance - float(np.take(fh['entry1/Amor/chopper/distance'], 0))  # mm

        self.lamdaCut                   = 2.5 # Aa

        startDate                       = str(fh['/entry1/start_time'][0].decode('utf-8'))
        self.startDate                  = datetime.strptime(startDate, '%Y-%m-%d %H:%M:%S')
        startDate                       = datetime.timestamp(self.startDate)
        self.countingTime               = float(np.take(fh['/entry1/Amor/detector/data/event_time_zero'], -1))/1e9 - startDate
                                          # not exact for low rates

        ka0 = 0.245 # given inclination of the beam after the Selene guide

        year_date = str(datetime.today()).split(' ')[0].replace("-", "/", 1)

        # deside from where to take the control paralemters
        try:
            self.mu   = float(np.take(fh['/entry1/Amor/instrument_control_parameters/mu'], 0))
            self.nu   = float(np.take(fh['/entry1/Amor/instrument_control_parameters/nu'], 0))
            self.kap  = float(np.take(fh['/entry1/Amor/instrument_control_parameters/kap'], 0))
            self.kad  = float(np.take(fh['/entry1/Amor/instrument_control_parameters/kad'], 0))
            self.div  = float(np.take(fh['/entry1/Amor/instrument_control_parameters/div'], 0))
            chopperSpeed      = float(np.take(fh['/entry1/Amor/chopper/rotation_speed'], 0))
            chopperPhase      = float(np.take(fh['/entry1/Amor/chopper/phase'], 0))
            ch1TriggerPhase   = float(np.take(fh['/entry1/Amor/chopper/ch1_trigger_phase'], 0))
            polarizationConfigLabel = float(np.take(fh['/entry1/Amor/polarization/configuration/value'], 0))
        except (KeyError, IndexError):
            logging.warning(f"     using parameters from nicos cache")
            cachePath = '/home/amor/cache/'
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-mu/'+year_date)).split('\t')[-1]
            self.mu = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-nu/'+year_date)).split('\t')[-1]
            self.nu = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-kap/'+year_date)).split('\t')[-1]
            self.kap = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-kad/'+year_date)).split('\t')[-1]
            self.kad = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-div/'+year_date)).split('\t')[-1]
            self.div = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-ch1_speed/'+year_date)).split('\t')[-1]
            chopperSpeed = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-chopper_speed'+year_date)).split('\t')[-1]
            chopperPhase = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-chopper_phase'+year_date)).split('\t')[-1]
            ch1TriggerPhase = float(value)
            value = str(subprocess.getoutput('/usr/bin/grep "value" '+cachePath+'nicos-ch1_trigger_phase'+year_date)).split('\t')[-1]

        self.tau = 30. / chopperSpeed
 
        self.chopperTriggerPhase = ch1TriggerPhase + chopperPhase/2

        polarizationConfigs = ['undefined', 'oo', 'po', 'mo', 'op', 'pp', 'mp', 'om', 'pm', 'mm']
        self.polarizationConfig = polarizationConfigs[int(polarizationConfigLabel)]

        # for .ort header

        self.title                      = str(fh['entry1/title'][0].decode('utf-8'))

        self.proposal_id                = str(fh['entry1/proposal_id'][0].decode('utf-8'))

        self.userName                   = str(fh['entry1/user/name'][0].decode('utf-8'))
        try:
          self.userEmail                = str(fh['entry1/user/email'][0].decode('utf-8'))
        except:
          self.userEmail                = None

        self.sampleName                 = str(fh['entry1/sample/name'][0].decode('utf-8'))
        self.sampleModel                = str(fh['entry1/sample/model'][0].decode('utf-8'))

        #!!! self.instrumentName             = str(fh['entry1/Amor/name'][0])
        self.instrumentName             = 'Amor'
        #!!! self.instrumentType             = str(fh['entry1/Amor/type'][0])
        self.source                     = str(fh['entry1/Amor/source/name'][0])
        #!!! self.sourceType                 = str(fh['entry1/Amor/source/type'][0])
        #!!! self.sourceProbe                = str(fh['entry1/Amor/source/probe'][0])
        self.sourceProbe                = 'neutron'

        fh.close()

#==============================================================================
def resolveNumbers(dataPath, ident):
    if ident == '0' or '-' in ident[0]:
        try: 
            nnr = int(ident)
        except:
            logging.error("ERROR: '{}' is no valid file identifier!".format(ident))
        fileNames = glob.glob(dataPath+'/*.hdf')
        fileNames.sort()
        fileName  = fileNames[nnr-1] 
        fileName = fileName.split('/')[-1]
        fileNumber = fileName.split('n')[1].split('.')[0].lstrip('0')
        numberString = fileNumber
        numberList   = [int(fileNumber)]
    elif 'a' in ident[0]:
        fileNumber = fileName.split('n')[1].split('.')[0].lstrip('0')
        numberString = fileNumber
        numberList   = [int(fileNumber)]
    else:
        if 'r' in ident: # substitute 'r' (recent) for the actual number
            fileName = glob.glob(dataPath+'/*.hdf')[-1]
            fileName = fileName.split('/')[-1]
            fileNumber = fileName.split('n')[1].split('.')[0].lstrip('0')
            ident = ident.replace('r', fileNumber, 1)
        numberString = ident
        numberList = get_flist(ident) 

    return numberString, numberList 

#==============================================================================
def fileNameCreator(dataPath, ident):
    clas = commandLineArgs()
    # create path/filename
    ident=str(ident)
    if 'a' in ident:
        fileName = ident
        
    else:
        try: 
            nnr = int(ident)
        except:
            logging.error("ERROR: '{}' is no valid file identifier!".format(ident))

        if nnr <= 0 :
            fileName = glob.glob(dataPath+'/*.hdf')[nnr-1]
            fileName = fileName.split('/')[-1]
        else:
            fileName = f'amor{clas.year}n{ident:>06s}'
            #fileName = 'amor2021n{:>06s}'.format(ident)
        
    fileName = fileName.split('.')[0]
    fileName = fileName+'.hdf'
    fileName = dataPath+fileName
    #filename = '/home/software/kafka-to-nexus/'+filename
    #print(fileName)
    
    fileNumber = fileName.split('n')[1].split('.')[0].lstrip('0')
  
    return fileName, fileNumber

#==============================================================================
def selectTime(timeMin, timeMax, dataPacketTime_p, dataPacket_p, detPixelID_e, tof_e):
    dataPacketTime_p = np.array(dataPacketTime_p)/1e9
    dataPacket_p     = np.array(dataPacket_p)
    detPixelID_e     = np.array(detPixelID_e)
    tof_e            = np.array(tof_e)
    startTime = dataPacketTime_p[0]
    stopTime  = dataPacketTime_p[-1]
    if timeMin > (stopTime-startTime): 
        logging.error('ERROR: time interval [{} : {}] s outside measurement time range [0 : {}] s'.format(timeMin, timeMax, stopTime-startTime))
        sys.exit()
    dataPacket_p     = dataPacket_p[(dataPacketTime_p-startTime)<timeMax]
    dataPacketTime_p = dataPacketTime_p[(dataPacketTime_p-startTime)<timeMax]
    dataPacket_p     = dataPacket_p[(dataPacketTime_p-startTime)>=timeMin]
    dataPacketTime_p = dataPacketTime_p[(dataPacketTime_p-startTime)>=timeMin]
    detPixelID_e     = detPixelID_e[dataPacket_p[0]:dataPacket_p[-1]]
    tof_e            = tof_e[dataPacket_p[0]:dataPacket_p[-1]]
    return dataPacket_p, dataPacketTime_p, detPixelID_e, tof_e, 
    
#==============================================================================
def ort_header(fileName):

    #print(fileName)
    meta = Meta(fileName)
    '''
    Build information object for ORSO file headers.
    '''

    #for fidx in self.get_flist(self.Files):
    #  fdata = data_vault['tmp%s'%fidx]
    #  fdate = datetime.strptime(fdata['timestamp'], '%Y-%m-%d %H:%M:%S')
    #  datafiles.append(fileio.File(file = fdata['srcname'], timestamp = fdate))

    inst = fileio.InstrumentSettings(
        incident_angle = fileio.ValueRange(
            mu + meta.kap + meta.kad - meta.div/2,
            mu + meta.kap + meta.kad + meta.div/2,
            'deg'
        ),
        wavelength = fileio.ValueRange(
            lamdaMin,
            lamdaMax, 
            'angstrom'
        ),
    )
    inst.mu = fileio.Value(mu, 'deg', comment = 'sample angle to horizon')
    inst.nu = fileio.Value(nu, 'deg', comment = 'detector angle to horizon')
    inst.kap = fileio.Value(meta.kap, 'deg', comment = 'nominal beam inclination')
    inst.kad = fileio.Value(meta.kad, 'deg', comment = 'offset of beam inclination')
    inst.div = fileio.Value(meta.div, 'deg', comment = 'incoming beam divergence')
    mess = fileio.Measurement(
        instrument_settings = inst,
        data_files = fileName.split('/')[-1],
        counting_time = fileio.Value(meta.countingTime, 's'),
        scheme = 'angle- and energy-dispersive')
    se_dict = None
    smpl = fileio.Sample(
        name  = meta.sampleName,
        model = meta.sampleModel,
        sample_parameters = se_dict
    )
    experiment = fileio.Experiment(
        title = meta.title,
        instrument = meta.instrumentName,
        #instrument_type = meta.instrumentType,
        start_date = meta.startDate,
        probe = meta.sourceProbe,
        #facility = meta.sourceName,
        #source_type = meta.sourceType,
        proposalID = meta.proposal_id,
    )
    ds = fileio.DataSource(
        fileio.Person(
            meta.userName, 
            None, 
            contact = meta.userEmail
        ),
        experiment,
        smpl,
        mess
    )
    red = fileio.Reduction(
        software = fileio.Software(
             'events2histogram', 
             version = __version__
        ),
        timestamp = datetime.now(),
        creator = None, 
        corrections = None,
        computer = None,
        call = 'python '+' '.join(sys.argv)
    )
    cols = [
        fileio.Column('Qz', '1/angstrom'),
        fileio.Column('R', comment = 'uncorrected intensity'),
        fileio.Column('sR', comment = 'sigma of gaussian probability function'),
        fileio.Column('sQz', '1/angstrom', comment = 'resolution based only on sigma_lambda!')
    ]

    header = fileio.Orso(ds, red, cols)

    return header

#==============================================================================
class PlotSelection:

    # header / meta data
    
    def header(self, filename, mu, nu, totalCounts, countingTime, polarizationConfig):
        number = filename.split('n')[1].split('.')[0].lstrip('0')
        header = "#{}   \u03bc={:>1.2f} \u03bd={:>1.2f} {} {:>10} cts   {:>8.1f} s".format(number, mu+5e-3, nu+5e-3, polarizationConfig, totalCounts, countingTime)
        return header

    def headline(self, numberString, totalCounts, polarizationConfig):
        headLine = "#{}   \u03bc={:>1.2f} \u03bd={:>1.2f}  p={}  {:>12,} cts   {:>8.1f} s".format(numberString, mu+5e-3, nu+5e-3, polarizationConfig, totalCounts, countingTime)
        return headLine

    # grids

    def lamda_grid(self):
        dldl      = 0.005 # Delta lambda / lambda
        if foldback:
            lamda_grid = lamdaMin*(1+dldl)**np.arange(int(np.log(lamdaMax/lamdaMin)/np.log(1+dldl)+1))
        else:
            lamda_grid = np.arange(0.01, 2.*lamdaMax-2.*lamdaMin, 0.1)
        return lamda_grid
    
    def theta_grid(self):
        det = Detector()

        bladeAngle = np.rad2deg( 2. * np.arcsin(0.5*det.bladeZ / detectorDistance) )
        blade_grid = np.arctan( np.arange(33) * det.dZ / ( detectorDistance + np.arange(33) * det.dX) )
        blade_grid = np.rad2deg(blade_grid)
        stepWidth = blade_grid[1] - blade_grid[0]
        blade_grid = blade_grid - 0.2 * stepWidth

        delta_grid = []
        for b in np.arange(det.nBlades-1):
            delta_grid = np.concatenate((delta_grid, blade_grid), axis=None)
            blade_grid = blade_grid + bladeAngle
            delta_grid = delta_grid[delta_grid<blade_grid[0]-0.5*stepWidth]
        delta_grid = np.concatenate((delta_grid, blade_grid), axis=None)

        theta_grid = nu - mu - np.flip(delta_grid) + 0.5*det.nBlades * bladeAngle    
  
        theta_grid = theta_grid[theta_grid>=thetaMin]
        theta_grid = theta_grid[theta_grid<=thetaMax]

        return theta_grid

    def q_grid(self):
        dqdq      = 0.010 # Delta q_z / q_z
        q_grid = qMin*(1.+dqdq)**np.arange(int(np.log(qMax/qMin)/np.log(1+dqdq)))
        return q_grid

    # create PNG with several plots

    def all(self, numberString, arg, data_e, polarizationConfig):
        #cmap='gist_earth'
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        y_grid = np.arange(64)
        I_yt, bins_y, bins_t = np.histogram2d(data_e[:,4], data_e[:,8], bins = (y_grid, self.theta_grid()))
        I_lt, bins_l, bins_t = np.histogram2d(data_e[:,7], data_e[:,8], bins = (self.lamda_grid(), self.theta_grid()))
        I_q, bins_q = np.histogram(data_e[:,9], bins = self.q_grid())
        #q_lim = 4*np.pi*np.array([ max( np.sin(self.theta_grid()[0]*np.pi/180.)/self.lamda_grid()[-1] , 1e-4 ),
        #                           min( np.sin(self.theta_grid()[-1]*np.pi/180.)/self.lamda_grid()[0] , 0.03 )])
        q_lim = np.array([qMin, qMax])
        if arg == 'lin':
            #vmin = min(np.min(I_lt), np.min(I_yt))
            vmin = 0
            vmax = max(5, np.max(I_lt), np.max(I_yt))
        else:
            vmin = 0
            vmax = max(1, np.log(np.max(I_lt)+.1)/np.log(10)*1.05, np.log(np.max(I_yt)+.1)/np.log(10)*1.05)
        # I(y, theta)   
        fig = plt.figure()
        axs = fig.add_gridspec(2,3)
        myt = fig.add_subplot(axs[0,0])
        myt.set_title('detector area')
        myt.set_xlabel('$y ~/~ \\mathrm{bins}$')
        myt.set_ylabel('$\\theta ~/~ \\mathrm{deg}$')
        if arg == 'lin':
            myt.pcolormesh(bins_y, bins_t, I_yt.T, cmap=cmap, vmin=vmin, vmax=vmax)
        else:
            myt.pcolormesh(bins_y, bins_t, (np.log(I_yt + 5.e-1) / np.log(10.)).T, cmap=cmap, vmin=vmin, vmax=vmax)
        # I(lambda, theta)   
        mlt = fig.add_subplot(axs[0,1:])
        mlt.set_title('angle- and energy disperse')
        mlt.set_xlabel('$\\lambda ~/~ \\mathrm{\\AA}$')
        mlt.axes.get_yaxis().set_visible(False)
        if arg == 'lin':
            cb = mlt.pcolormesh(bins_l, bins_t, I_lt.T, cmap=cmap, vmin=vmin, vmax=vmax)
        else:
            cb = mlt.pcolormesh(bins_l, bins_t, (np.log(I_lt + 5.e-1) / np.log(10.)).T, cmap=cmap, vmin=vmin, vmax=vmax)
        # I(q_z)
        lqz = fig.add_subplot(axs[1,:])
        lqz.set_title('$I(q_z)$')
        lqz.set_ylabel('$\\log_{10}(\\mathrm{cnts})$')
        lqz.set_xlabel('$q_z~/~\\mathrm{\\AA}^{-1}$')
        lqz.set_xlim(q_lim)
        if arg == 'lin':
            plt.plot(bins_q[:-1], I_q, color='blue', linewidth=0.5)
        else:
            err_q = np.sqrt(I_q+1)
            low_q = np.where(I_q-err_q>0, I_q-err_q, 0.1)
            plt.fill_between(bins_q[:-1], np.log(low_q)/np.log(10), np.log(I_q+err_q/2)/np.log(10), color='lightgrey')
            plt.plot(bins_q[:-1], np.log(I_q+5e-1)/np.log(10), color='blue', linewidth=0.5)
            lw = I_q[ ((q_lim[0] < bins_q[:-1]) & (bins_q[:-1] < q_lim[1])) ].min() 
            plt.ylim(max(-0.1, np.log(lw+.1)/np.log(10)-0.1), )
        #
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=2.8, c='r')
        fig.colorbar(cb, ax=mlt)
        plt.subplots_adjust(hspace=0.6, wspace=0.1)
        plt.savefig(output, format='png', dpi=150)
        #plt.close()

    # create PNG with one plot

    def Iyz(self, numberString, arg, data_e, polarizationConfig):
        det = Detector()
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        y_grid = np.arange(64)
        z_grid = np.arange(det.nBlades*32)
        I_yz, bins_y, bins_z = np.histogram2d(data_e[:,4], (det.nBlades-data_e[:,2])*32-data_e[:,3], bins = (y_grid, z_grid))
        if arg == 'file':
            print('# y z conts')
            for y in range(len(bins_y)-1):
                for z in range(len(bins_z)-1):
                    print(" %6.3f %6.4f %10.3e" %(bins_y[y], bins_z[z], I_yz[y,z]))
                print("")
            return
        elif arg == 'log':
            vmin = 0
            vmax = max(1, np.log(np.max(I_yz)+.1)/np.log(10)*1.05)
            plt.pcolormesh(bins_y[:],bins_z[:],(np.log(I_yz+6e-1)/np.log(10)).T, cmap=cmap, vmin=vmin, vmax=vmax)
        else:
            plt.pcolormesh(bins_y[:],bins_z[:],I_yz.T, cmap=cmap)
        plt.xlabel('$y ~/~ \\mathrm{bins}$')
        plt.ylabel('$z ~/~ \\mathrm{bins}$')
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.colorbar()
        plt.savefig(output, format='png', dpi=150)
        #plt.close()

    def Ilt(self, numberString, arg, data_e, polarizationConfig) :
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        I_lt, bins_l, bins_t = np.histogram2d(data_e[:,7], data_e[:,8], bins = (self.lamda_grid(), self.theta_grid()))
        if arg == 'file':
            print('# lambda theta conts')
            for l in range(len(bins_l)-1):
                for t in range(len(bins_t)-1):
                    print(" %6.3f %6.4f %10.3e" %(bins_l[l], bins_t[t], I_lt[l,t]/ bins_t[t]))
                print("") 
            return    
        elif arg == 'log':
            vmax = max(1, np.log(np.max(I_lt)+.1)/np.log(10)*1.05 )
            plt.pcolormesh(bins_l, bins_t, (np.log(I_lt+I_lt[I_lt>0].min()/2)/np.log(10.)).T, cmap=cmap, vmin=0, vmax=vmax)
        else :
            vmax = max(np.max(I_lt), 5)
            plt.pcolormesh(bins_l, bins_t, I_lt.T, cmap=cmap, vmin=0, vmax=vmax)
        plt.xlim(0,)
        if np.min(bins_t) > 0.01 :
            plt.ylim(bottom=0)
        else:
            plt.ylim(bottom=np.min(bins_t))
        if np.max(bins_t) < -0.01:
            plt.ylim(top=0)
        else:
            plt.ylim(top=np.max(bins_t))
        plt.xlabel('$\\lambda ~/~ \\mathrm{\\AA}$')
        plt.ylabel('$\\theta ~/~ \\mathrm{deg}$')
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.colorbar()
        plt.savefig(output, format='png', dpi=150)
        #plt.close()
    
    def Itz(self, numberString, arg, data_e, polarizationConfig):
        det = Detector()
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        if foldback:
            time_grid = np.arange(0, tau, 0.0005)
        else:
            time_grid = np.arange(0, 2.*tau, 0.0005)
        z_grid = np.arange(det.nBlades*32+1)
     
        I_tz, bins_t, bins_z = np.histogram2d(data_e[:,0], 32*det.nBlades-data_e[:,2]*32-data_e[:,3], bins = (time_grid, z_grid))
        if arg == 'file':
            print('# time z conts')
            for t in range(len(bins_t)-1):
                for z in range(len(bins_z)-1):
                    print(" %6.3f %6.4f %10.3e" %(bins_t[t], bins_z[z], I_tz[t,z]))
                print("")
            return    
        elif arg == 'log':
            vmax = max(2., np.log(np.max(I_tz)+.1)/np.log(10)*1.05 )
            plt.pcolormesh(bins_t, bins_z, (np.log(I_tz+5.e-1)/np.log(10.)).T, cmap=cmap, vmin=0, vmax=vmax)
        else :
            vmax = max(np.max(I_tz), 5)
            plt.pcolormesh(bins_t, bins_z, I_tz.T, cmap=cmap, vmin=0, vmax=vmax)
        if True:
            plt.xlim(0,)
            plt.ylim(0,)
        plt.xlabel('$t ~/~ \\mathrm{s}$')
        plt.ylabel('$z$ pixel row')
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.colorbar()
        plt.savefig(output, format='png', dpi=150)
        #plt.close()
    
    def Iq(self, numberString, arg, data_e, polarizationConfig):
        I_q, bins_q = np.histogram(data_e[:,9], bins = self.q_grid())
        err_q = np.sqrt(I_q+1)
        #q_lim = 4*np.pi*np.array([ max( np.sin(self.theta_grid()[0]*np.pi/180.)/self.lamda_grid()[-1] , 1e-4 ),
        #                           min( np.sin(self.theta_grid()[-1]*np.pi/180.)/self.lamda_grid()[0] , 0.03 )])
        q_lim = [qMin, qMax]
        print(q_lim)
        if arg == 'file':
            header = '# q counts'
            I_q =  np.vstack((bins_q[:-1], I_q, err_q))
            np.savetxt(sys.stdout, I_q.T, header=header)
            logging.info('     use `-p ort` instead of `-p Iq`!')
            return
        elif arg == 'log':
            low_q = np.where(I_q-err_q>0, I_q-err_q, 0.1)
            plt.fill_between(bins_q[:-1], np.log(low_q)/np.log(10), np.log(I_q+err_q/2)/np.log(10), color='lightgrey')
            plt.plot(bins_q[:-1], np.log(I_q+5e-1)/np.log(10), color='blue', linewidth=0.5)
            lw = I_q[ ((q_lim[0] < bins_q[:-1]) & (bins_q[:-1] < q_lim[1])) ].min() 
            plt.ylim(max(-0.1, np.log(lw+.1)/np.log(10)-0.1), )
        else:
            plt.plot(bins_q[:-1], I_q, color='blue', linewidth=0.5)
        plt.ylabel('$\\log_{10}(\\mathrm{cnts})$')
        plt.xlabel('$q_z ~/~ \\mathrm{\\AA}^{-1}$')
        plt.xlim(q_lim)
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png', dpi=150)
        #plt.close()

    def Il(self, numberString, arg, data_e, polarizationConfig):
        I_l, bins_l = np.histogram(data_e[:,7], bins = self.lamda_grid())
        if arg == 'file':
            header = '# lambda counts'
            I_l = np.vstack((bins_l[:-1], I_l)) 
            np.savetxt(sys.stdout, I_l.T, header=header)
            return
        elif arg == 'lin':
            plt.plot(bins_l[:-1], I_l)
            plt.ylabel('$I ~/~ \\mathrm{cnts}$')
        else:
            plt.plot(bins_l[:-1], np.log(I_l+5.e-1)/np.log(10.))
            plt.ylabel('$\\log_{10} I ~/~ \\mathrm{cnts}$')
        plt.xlabel('$\\lambda ~/~ \\mathrm{\\AA}$')
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png', dpi=150)
        #plt.close()
    
    def It(self, numberString, arg, data_e, polarizationConfig):
        I_t, bins_t = np.histogram(data_e[:,8], bins = self.theta_grid())
        if arg == 'file':
            header = '# 2theta counts'
            I_t = np.vstack((bins_t[:-1], I_t))
            np.savetxt(sys.stdout, I_t.T, header=header)
            return
        else:
            plt.plot( I_t, bins_t[:-1])
        plt.xlabel('$\\mathrm{cnts}$')
        plt.ylabel('$\\theta ~/~ \\mathrm{deg}$')
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png', dpi=150)
        #plt.close()

    def tof(self, numberString, arg, data_e, polarizationConfig):
        if foldback:
            time_grid = np.arange(0, 1.3*tau, 0.0005)
        else:
            time_grid = np.arange(0, 2.*tau, 0.0005)
        I_t, bins_t = np.histogram(data_e[:,0], bins = time_grid)
        if arg == 'file':
            header = '# time counts'
            I_t = np.vstack((bins_t[:-1], I_t))
            np.savetxt(sys.stdout, I_t.T, header=header)
            return
        elif arg == 'lin':
            plt.plot(bins_t[:-1]+tau, I_t)
            plt.plot(bins_t[:-1], I_t)
            plt.plot(bins_t[:-1]+2*tau, I_t)
        else:  
            lI_t = np.log(I_t+5.e-1)/np.log(10.)
            plt.plot(bins_t[:-1]+tau, lI_t)
            plt.plot(bins_t[:-1], lI_t)
            plt.plot(bins_t[:-1]+2*tau, lI_t)
        plt.ylabel('log(counts)')
        plt.xlabel('time / s')
        headline = self.headline(numberString, np.shape(data_e)[0], polarizationConfig)
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png')

    def ort(self, numberString, arg, data_e):
        I_q, bins_q = np.histogram(data_e[:,9], bins = self.q_grid())
        sI_q = np.sqrt(I_q)
        sq_q = bins_q[:-1]*0.022
        I_q =  np.vstack((bins_q[:-1], I_q, sI_q, sq_q))

        datasets = []
        fileNumber = get_flist(numberString)[-1]
        fileName = fileNameCreator(dataPath, fileNumber)[0]
        header = ort_header(fileName)
        orso_data=fileio.OrsoDataset(header, I_q.T)
        datasets.append(orso_data)
        fileio.save_orso(datasets, 'e2h.ort', data_separator='\n')

#==============================================================================
#==============================================================================
def endIt(signal, frame):
    print('\n# e2h life mode stopped\n')
    sys.exit(0)
    
#==============================================================================
def get_flist(flist):
    # resolve short notation of filenumbers into a list of filenumbers
    #     e.g. '3,4-9:2,12-14' -> '3, 4, 6, 8, 12, 13, 14'
    out_list=np.array([], dtype=int)
    if ',' in flist:
         fsublists = flist.split(',')
    else: 
        fsublists = [flist]
    for fsublist in fsublists:
        if '-' in fsublist:
            if ':' in fsublist:
                limits, step = fsublist.split(':', 1)
            else:
                step = 1
                limits = fsublist
            for number in range(int(limits.split('-', 1)[0]), 
                                int(limits.split('-', 1)[1])+1, int(step)):
                out_list = np.append(out_list, number)
        else:
            out_list = np.append(out_list, int(fsublist))

    return out_list

#==============================================================================
def process(dataPath, ident, clas):
    #================================
    # constants
    hdm       = 6.626176e-34/1.674928e-27               # h / m
    #================================
    # instrument specific parameters
    #================================
    global lamdaMin, lamdaMax, qMin, qMax, thetaMin, thetaMax
    # defaults   
    lamdaCut  = 2.5     # Aa  used to reshuffle tof
    # data filtering and folding 
    
    #================================
    if clas.lambdaRange:
        lamdaMin  = clas.lambdaRange[0]
        lamdaMax  = clas.lambdaRange[1]
    else:
        lamdaMin  = lamdaCut 

    if clas.timeIntervalAbs:
        timeMin   = clas.timeIntervalAbs[0]
        timeMax   = clas.timeIntervalAbs[1]
    elif clas.timeIntervalInc:
        timeMin   = clas.timeIntervalInc[0] * clas.timeIntervalInc[1]
        timeMax   = clas.timeIntervalInc[0] * (clas.timeIntervalInc[1] + 1.)
    else:
        timeMin = 0
        timeMax = 0

    chopperPhase = clas.chopperPhase
    tofOffset    = clas.TOFOffset
    thetaMin     = clas.thetaRange[0]
    thetaMax     = clas.thetaRange[1]
    yMin         = clas.yRange[0]
    yMax         = clas.yRange[1]
    qMin         = clas.qRange[0]
    qMax         = clas.qRange[1]
    global foldback
    foldback     = not clas.noTOFCorrection

    #================================
    # find and open input file 
    global ev

    data_eSum = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    sumTime = 0

    numberString, numberList = resolveNumbers(dataPath, ident)

    for number in numberList:
        number= str(number)
    
        filename, fileNumber = fileNameCreator(dataPath, number)

        if verbous:
            logging.info('events2histogram processing file ->\033[1m {} \033[0m<-'.format(fileNumber))
    
        for i in range(6):
            ev = h5py.File(filename, 'r', swmr=True)
            try:
                ev['/entry1/Amor/detector/data/event_time_zero'][-1]
                break 
            except (KeyError, IndexError):
                ev.close()
                if i < 5:
                    if verbous:
                        print("no data yet, retrying ({})    ".format(10-2*i), end='\r')
                    time.sleep(2)
                    continue
                else:
                    if verbous:
                        print("# time-out: no longer waiting for data!\a")
                    return
        
        # get and process data
        meta = Meta(filename)

        global mu, nu, tau

        if clas.mu < 98.: 
            mu = clas.mu 
        else:
            mu = meta.mu + clas.muOffset

        if clas.nu < 98.:
            nu = clas.nu
        else:
            nu = meta.nu

        if clas.chopperSpeed:
            tau = 30./ clas.chopperSpeed
        else:
            tau = meta.tau

        try:
            chopperTriggerPhase
        except NameError:
            chopperTriggerPhase = meta.chopperTriggerPhase


        global countingTime, detectorDistance, chopperDetectorDistance, polarizationConfig
        polarizationConfig = meta.polarizationConfig
        detectorDistance = meta.detectorDistance
        chopperDetectorDistance = meta.chopperDetectorDistance
        countingTime = meta.countingTime

        if verbous:
            logging.info("   mu = {:>4.2f} deg,  nu = {:>4.2f} deg".format(mu, nu))
            logging.info(f'   polarization config: {polarizationConfig}')

        try: lamdaMax
        except NameError: lamdaMax = lamdaMin + tau * hdm/chopperDetectorDistance * 1e13
    
        tofOffset = -tau * chopperPhase / 180.                       # mismatch of chopper pulse and time-zero
        tofCut = lamdaCut * chopperDetectorDistance / hdm * 1.e-13  # tof of frame start

        tof_e  = np.array(ev['/entry1/Amor/detector/data/event_time_offset'][:], dtype=np.uint64)/1.e9 + tofOffset # tof 
    
        detPixelID_e = np.array(ev['/entry1/Amor/detector/data/event_id'][:], dtype=np.uint64)  # pixel index
        #totalCounts = len(detPixelID_e)
    
        dataPacket_p = np.array(ev['/entry1/Amor/detector/data/event_index'][:], dtype=np.uint64)  # data packet index

        if timeMax>0: # pick only the time interval defined by `-i` or `-I`
            dataPacketTime_p = np.array(ev['/entry1/Amor/detector/data/event_time_zero'][:], dtype=np.uint64)  # data packet index
            dataPacket_p, dataPacketTime_p, detPixelID_e, tof_e = selectTime(timeMin, timeMax, dataPacketTime_p, dataPacket_p, detPixelID_e, tof_e)
            countingTime = dataPacketTime_p[-1]-dataPacketTime_p[0]+1

        # command line argument not yet defined
        #filterThreshold = 0
        #if filterThreshold:
        #    detPixelID_e, tof_e = filterTof(detPixelID_e, tof_e, dataPacket_p, filterThreshold)
 
        if foldback:
            tof_e  = np.where(tof_e<tofCut, tof_e+2.*tau, tof_e)
            tof_e  = np.where(tof_e>tau+tofCut, tof_e-tau, tof_e)

        data_e = analyse_ev(detPixelID_e, tof_e, yMin, yMax, thetaMin, thetaMax)

        ev.close()

        data_eSum = np.append(data_eSum, data_e, axis=0)
        sumTime += countingTime

    if verbous:
        logging.info("   total counts = {} in {:6.1f} s".format(np.shape(data_e)[0], sumTime))

    if clas.spy:
        number = filename.split('n')[1].split('.')[0].lstrip('0')
        logging.info('chopper speed={:>4.0f} rpm, phase={:>5.3f} deg, tau={} s'.format(30/tau, chopperTriggerPhase, tau))
        logging.info('nr={}, spn={}, cnt={}, tme={}'.format(number, polarizationConfig, np.shape(data_eSum)[0], sumTime))
        logging.info('mu={:>1.2f}, nu={:>1.2f}, kap={:>1.2f}, kad={:>1.2f}, div={:>1.2f}'.format(mu, nu, kap, kad, div))
        sys.exit()

    #================================
    # plotting data
    plotType = clas.plot[0]
    try:
        arg = clas.plot[1]
    except IndexError:
        arg = 'def'
    plott = PlotSelection()
    #string = f"plott.{plotType} (numberString, '{arg}', data_e)"
    try:
        plotFunction = getattr(plott, plotType)
        plotFunction(numberString, arg, data_e, polarizationConfig)
        plt.close()
    except Exception as e:
        logging.error(f"ERROR: '{plotType}' is no known output format!")
        logging.error(f"       original error: {e}")

#==============================================================================   
def commandLineArgs():    
    msg = "events2histogram reads the eventstream from an hdf raw file and \
           creates various histogrammed outputs or plots."
    clas = argparse.ArgumentParser(description = msg)

    clas.add_argument("-b", "--noTOFCorrection", 
                            action='store_true',                       
                            help ="do not correct tof of seond neutron pulse")
    clas.add_argument("-c", "--chopperSpeed",
                            type=float,
                            help ="chopper speed in rpm")
    clas.add_argument("-d", "--dataPath",  
                            help ="relative path to directory with .hdf files")
    clas.add_argument("-D", "--absDataPath",  
                            help ="absolute path to directory with .hdf files")
    clas.add_argument("-f", "--fileIdent",       
                            default='0',                               
                            help ="file number or offset (if negative)")
    clas.add_argument("-I", "--timeIntervalInc", 
                            nargs=2, 
                            type=float,
                            help ="time interval length and increment")
    clas.add_argument("-i", "--timeIntervalAbs",
                            nargs=2, 
                            type=float,
                            help ="absolute time interval to be processed")
    clas.add_argument("-l", "--lambdaRange",                            
                            nargs=2, 
                            type=float,
                            help ="wavelength range to be used")
    clas.add_argument("-M", "--muOffset",        
                            default=0.,                     
                            type=float,
                            help ="mu offset")
    clas.add_argument("-m", "--mu",                                              
                            default=99.,
                            type=float,
                            help ="value of mu")
    clas.add_argument("-n", "--nu",                                              
                            default=99.,
                            type=float,
                            help ="value of nu")
    clas.add_argument("-P", "--chopperPhase",    
                            default=-7.5,                     
                            type=float,
                            help ="chopper phase offset")
    clas.add_argument("-p", "--plot",            
                            default=['all', 'def'], 
                            nargs='+',         
                            help ="select what to plot or write")
    clas.add_argument("-q", "--qRange",          
                            default=[0.005, 0.30], 
                            nargs=2, 
                            type=float,
                            help ="q_z range")
    clas.add_argument("-r", "--iDonno",          
                            action='store_true',                       
                            help ="no idea")
    clas.add_argument("-s", "--spy",             
                            action='store_true',                       
                            help ="report a few key values, no plotting or writing")
    clas.add_argument("-T", "--TOFOffset",       
                            default=0.0,                    
                            type=float,
                            help ="TOF zero offset")
    clas.add_argument("-t", "--thetaRange",      
                            default=[-12., 12.],   
                            nargs=2, 
                            type=float,
                            help ="theta range to be used")
    clas.add_argument("-u", "--update",          
                            action='store_true', 
                            help ="update output every 5 seconds")
    clas.add_argument("-Y", "--year",
                            default = str(datetime.today()).split('-')[0],
                            help = "year, the measurement was performed")
    clas.add_argument("-y", "--yRange",          
                            default=[0, 63],
                            nargs=2, 
                            type=int,  
                            help ="detector y range to be used")

    return clas.parse_args()

#==============================================================================   
def get_dataPath(clas):
    if clas.dataPath:
        dataPath = clas.dataPath + '/'
        if not os.path.exists(dataPath):
             sys.exit('# *** the directory "'+dataPath+'" does not exist ***')
    if clas.absDataPath:
        dataPath = clas.absDataPath + '/'
    elif os.path.exists('./raw'):
        dataPath  = "./raw/"
    elif os.path.exists('../raw'):
        dataPath  = "../raw/"
    else:
        sys.exit('# *** please provide the path to the .hdf data files (-d <rel path> | -D <abs path>, default is "./raw")')

    return dataPath

#==============================================================================   
def get_directDataPath(clas):
        #dataPath = clas.dataPath + '/'
        year = str(datetime.today()).split('-')[0]
        year_date = str(datetime.today()).split(' ')[0].replace("-", "/", 1)
        pNr  = str(subprocess.getoutput('/usr/bin/grep "proposal\t" /home/amor/nicosdata/amor/cache/nicos-exp/'+year_date)[-1]).split('\'')[1]
        dataPath = '/home/amor/nicosdata/amor/data/'+year+'/'+pNr+'/'
        if not os.path.exists(dataPath):
             sys.exit('# *** the directory "'+dataPath+'" does not exist ***')

        return dataPath

#==============================================================================
def main():
    global verbous, output, dataPath
 
    clas = commandLineArgs()

    if clas.update:
        verbous = False
        logging.basicConfig(level=logging.ERROR, format='# %(message)s')
        logging.error('e2h: recreating "e2h_life.png" every 5 s')
        delay = 5
        output = 'e2h_life.png'
        oldtime = 0
        while True:
            dataPath = get_directDataPath(clas)
            fileName = fileNameCreator(dataPath, 0)[0]
            newtime = os.path.getmtime(fileName)
            if newtime-oldtime:
                print('\r# processing   (press ^C to stop)', end='', flush=True)
                process(dataPath, '0', clas)
                oldtime = newtime 
            for i in range(5):
                print('\r# waiting'+'.'*i+' '*(5-i)+' (press ^C to stop)', end='', flush=True)
                time.sleep(delay/5)
            signal.signal(signal.SIGINT, endIt)
    else:
        dataPath = get_dataPath(clas)
        logging.basicConfig(level=logging.INFO, format='# %(message)s')
        verbous = True
        output = 'e2h.png'
        process(dataPath, clas.fileIdent, clas)

#==============================================================================
#============================================================================== 
if __name__ == "__main__":
    main()     


