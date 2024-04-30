__version__ = '2024-03-30'

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
import logging
from datetime import datetime

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
    dZ                = bladeNr * 32 + bZ
    quantity = np.vstack((dZ.T, bY.T, delta.T, x.T)).T
    
    return quantity

#==============================================================================
def analyse_ev(event_e, tof_e, yMin, yMax, thetaMin, thetaMax):
    
    data_e = np.zeros((len(event_e), 9), dtype=float)
    
    # data_e column description:
    #    0: wall time / s
    #    1: pixelID
    #    2: z on detector
    #    3: y on detector
    #    4: delta / deg  = angle on detector
    #    5: path within detector / mm
    #    6: lambda / angstrom
    #    7: theta / deg
    #    8: q_z / angstrom^-1
    
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
    data_e[:,2:6] = pixelLookUp[np.int_(data_e[:,1])-1,:]

    #================================

    # filter y range
    filter_e = (yMin <= data_e[:,3]) & (data_e[:,3] <= yMax) 
    data_e = data_e[filter_e,:]
    
    # correct tof for beam size effect at chopper 
    data_e[:,0] -= ( data_e[:,4] / 180. ) * tau
      
    # effective flight path length
    #data_e[:,5] = chopperDetectorDistance + data_e[:,5] 
    
    # calculate lambda
    hdm       = 6.626176e-34/1.674928e-27               # h / m
    data_e[:,6] = 1.e13 * data_e[:,0] * hdm / ( chopperDetectorDistance + data_e[:,5] )
    
    # theta 
    data_e[:,7] = nu - mu + data_e[:,4] 

    # gravity compensation
    data_e[:,7] += np.rad2deg( np.arctan( 3.07e-10 * ( detectorDistance +  data_e[:,5]) * data_e[:,6] * data_e[:,6] ) )

    # filter theta range
    filter_l = (thetaMin <= data_e[:,7]) & (data_e[:,7] <= thetaMax) 
    data_e = data_e[filter_l,:]
    
    # q_z
    data_e[:,8] = 4*np.pi * np.sin( np.deg2rad( data_e[:,7] ) ) / data_e[:,6]
    
    # filter q_z range
    #filter_e = (qMin < data_e[:,6]) & (data_e[:,6] < qMax)
    #data_e = data_e[filter_e,:]

    return data_e

#==============================================================================
class Meta:
  # AMOR hdf dataset with associated properties from metadata
    def __init__(self, fileName):
        self.fileName = fileName

        fh = h5py.File(fileName, 'r', swmr=True)

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
            self.mu   = float(np.take(fh['/entry1/Amor/master_parameters/mu/value'], 0))
            self.nu   = float(np.take(fh['/entry1/Amor/master_parameters/nu/value'], 0))
            self.kap  = float(np.take(fh['/entry1/Amor/master_parameters/kap/value'], 0))
            self.kad  = float(np.take(fh['/entry1/Amor/master_parameters/kad/value'], 0))
            self.div  = float(np.take(fh['/entry1/Amor/master_parameters/div/value'], 0))
            chSp      = float(np.take(fh['/entry1/Amor/chopper/rotation_speed/value'], 0))
            self.chPh = float(np.take(fh['/entry1/Amor/chopper/phase/value'], 0))
        except (KeyError, IndexError):
            logging.warning(f"     using parameters from nicos cache")
            #cachePath = '/home/amor/nicosdata/amor/cache/'
            cachePath = '/home/nicos/amorcache/'
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
            chSp = float(value)
            self.chPh = np.nan

        if chSp:
            self.tau = 30. / chSp
        else:
            self.tau = 0
 
        try: # not yet correctly implemented in nexus template
            spin = str(fh['/entry1/polarizer/spin_flipper/spin'][0].decode('utf-8'))
            if spin == "b'p'":
                self.spin = 'p'
            elif spin == "b'm'":
                self.spin = 'm'
            elif spin == "b'up'":
                self.spin = 'p'
            elif spin == "b'down'":
                self.spin = 'm'
            elif spin == '?':
                self.spin = '?'
            else:
                self.spin = 'n'  
        except (KeyError, IndexError):
            self.spin = '?'

        fh.close()

#==============================================================================
def resolveNumber(dataPath, ident):
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
    else:
        fileNumber = ident

    return fileNumber 

#==============================================================================
def fileNameCreator(dataPath, ident):
    clas = commandLineArgs()
    ident=str(ident)
    try: 
        nnr = int(ident)
    except:
        logging.error("ERROR: '{}' is no valid file identifier!".format(ident))

    if nnr <= 0 :
        fileName = glob.glob(dataPath+'/*.hdf')[nnr-1]
        fileName = fileName.split('/')[-1]
    else:
        fileName = f'amor{clas.year}n{ident:>06s}'
        
    fileName = fileName.split('.')[0]
    fileName = fileName+'.hdf'
    fileName = dataPath+fileName
    
    fileNumber = fileName.split('n')[-1].split('.')[0].lstrip('0')
  
    return fileName, fileNumber

#==============================================================================
class PlotSelection:

    def headline(self, fileNumber, totalCounts):
        headLine = "#{}   \u03bc={:>1.2f} \u03bd={:>1.2f} {:>12,} cts   {:>8.1f} s".format(fileNumber, mu+5e-3, nu+5e-3, totalCounts, countingTime)
        return headLine

    # grids

    def y_grid(self):
        y_grid = np.arange(yMin, yMax+1, 1)
        return y_grid 

    def lamda_grid(self):
        dldl      = 0.005 # Delta lambda / lambda
        lMin = max(2, lamdaMin)
        lamda_grid = lMin*(1+dldl)**np.arange(int(np.log(lamdaMax/lMin)/np.log(1+dldl)+1))
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

    def all(self, fileNumber, arg, data_e):
        #cmap='gist_earth'
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        I_yt, bins_y, bins_t = np.histogram2d(data_e[:,3], data_e[:,7], bins = (self.y_grid(), self.theta_grid()))
        I_lt, bins_l, bins_t = np.histogram2d(data_e[:,6], data_e[:,7], bins = (self.lamda_grid(), self.theta_grid()))
        I_q, bins_q = np.histogram(data_e[:,8], bins = self.q_grid())
        q_lim = 4*np.pi*np.array([ max( np.sin(self.theta_grid()[0]*np.pi/180.)/self.lamda_grid()[-1] , 1e-4 ),
                                   min( np.sin(self.theta_grid()[-1]*np.pi/180.)/self.lamda_grid()[0] , 0.03 )])
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
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=2.8, c='r')
        fig.colorbar(cb, ax=mlt)
        plt.subplots_adjust(hspace=0.6, wspace=0.1)
        plt.savefig(output, format='png', dpi=150)

    # create PNG with one plot

    def Iyz(self, fileNumber, arg, data_e):
        det = Detector()
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        z_grid = np.arange(det.nBlades*32)
        I_yz, bins_y, bins_z = np.histogram2d(data_e[:,3], data_e[:,2], bins = (self.y_grid(), z_grid))
        if arg == 'log':
            vmin = 0
            vmax = max(1, np.log(np.max(I_yt)+.1)/np.log(10)*1.05)
            plt.pcolormesh(bins_y[:],bins_z[:],(np.log(I_yz+6e-1)/np.log(10)).T, cmap=cmap, vmin=vmin, vmax=vmax)
        else:
            plt.pcolormesh(bins_y[:],bins_z[:],I_yz.T, cmap=cmap)
        plt.xlabel('$y ~/~ \\mathrm{bins}$')
        plt.ylabel('$z ~/~ \\mathrm{bins}$')
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.colorbar()
        plt.savefig(output, format='png', dpi=150)

    def Ilt(self, fileNumber, arg, data_e) :
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        I_lt, bins_l, bins_t = np.histogram2d(data_e[:,6], data_e[:,7], bins = (self.lamda_grid(), self.theta_grid()))
        if arg == 'log':
            vmax = max(1, np.log(np.max(I_lt)+.1)/np.log(10)*1.05 )
            plt.pcolormesh(bins_l, bins_t, (np.log(I_lt+I_lt[I_lt>0].min()/2)/np.log(10.)).T, cmap=cmap, vmin=0, vmax=vmax)
        else :
            vmax = max(np.max(I_lt), 5)
            plt.pcolormesh(bins_l, bins_t, I_lt.T, cmap=cmap, vmin=0, vmax=vmax)
        plt.xlim(0,)
        #if np.min(bins_t) > 0.01 :
        #    plt.ylim(bottom=0)
        #else:
        #    plt.ylim(bottom=np.min(bins_t))
        #if np.max(bins_t) < -0.01:
        #    plt.ylim(top=0)
        #else:
        #    plt.ylim(top=np.max(bins_t))
        plt.xlim(lamdaMin, lamdaMax)
        plt.ylim(thetaMin, thetaMax)
        plt.xlabel('$\\lambda ~/~ \\mathrm{\\AA}$')
        plt.ylabel('$\\theta ~/~ \\mathrm{deg}$')
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.colorbar()
        plt.savefig(output, format='png', dpi=150)
    
    def Itz(self, fileNumber, arg, data_e):
        det = Detector()
        cmap = mpl.cm.gnuplot(np.arange(256))
        cmap[:1, :] = np.array([256/256, 255/256, 236/256, 1])
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        time_grid = np.arange(0, tau, 0.0005)
        z_grid = np.arange(det.nBlades*32+1)
     
        I_tz, bins_t, bins_z = np.histogram2d(data_e[:,0], data_e[:,2], bins = (time_grid, z_grid))
        if arg == 'log':
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
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.colorbar()
        plt.savefig(output, format='png', dpi=150)
    
    def Iq(self, fileNumber, arg, data_e):
        I_q, bins_q = np.histogram(data_e[:,8], bins = self.q_grid())
        err_q = np.sqrt(I_q+1)
        q_lim = 4*np.pi*np.array([ max( np.sin(self.theta_grid()[0]*np.pi/180.)/self.lamda_grid()[-1] , 1e-4 ),
                                   min( np.sin(self.theta_grid()[-1]*np.pi/180.)/self.lamda_grid()[0] , 0.03 )])
        if arg == 'log':
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
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png', dpi=150)

    def Il(self, fileNumber, arg, data_e):
        I_l, bins_l = np.histogram(data_e[:,6], bins = self.lamda_grid())
        if arg == 'lin':
            plt.plot(bins_l[:-1], I_l)
            plt.ylabel('$I ~/~ \\mathrm{cnts}$')
        else:
            plt.plot(bins_l[:-1], np.log(I_l+5.e-1)/np.log(10.))
            plt.ylabel('$\\log_{10} I ~/~ \\mathrm{cnts}$')
        plt.xlabel('$\\lambda ~/~ \\mathrm{\\AA}$')
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png', dpi=150)
    
    def It(self, fileNumber, arg, data_e):
        I_t, bins_t = np.histogram(data_e[:,7], bins = self.theta_grid())
        plt.plot( I_t, bins_t[:-1])
        plt.xlabel('$\\mathrm{cnts}$')
        plt.ylabel('$\\theta ~/~ \\mathrm{deg}$')
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png', dpi=150)

    def tof(self, fileNumber, arg, data_e):
        time_grid = np.arange(0, 1.3*tau, 0.0005)
        I_t, bins_t = np.histogram(data_e[:,0], bins = time_grid)
        if arg == 'lin':
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
        headline = self.headline(fileNumber, np.shape(data_e)[0])
        plt.title(headline, loc='left', y=1.0, c='r')
        plt.savefig(output, format='png')

#==============================================================================
def process(dataPath, ident, clas):
    #================================
    # constants
    hdm       = 6.626176e-34/1.674928e-27               # h / m
    #================================
    # instrument specific parameters
    #================================
    global lamdaMin, lamdaMax, qMin, qMax, thetaMin, thetaMax, yMin, yMax
    # defaults   
    lamdaCut  = 2.5     # Aa  used to reshuffle tof
    # data filtering and folding 
    
    #================================
    if clas.lambdaRange:
        lamdaMin  = clas.lambdaRange[0]
        lamdaMax  = clas.lambdaRange[1]
    else:
        lamdaMin  = lamdaCut 

    chopperPhase = clas.chopperPhase
    tofOffset    = clas.TOFOffset
    thetaMin     = clas.thetaRange[0]
    thetaMax     = clas.thetaRange[1]
    yMin         = clas.yRange[0]
    yMax         = clas.yRange[1]
    qMin         = clas.qRange[0]
    qMax         = clas.qRange[1]

    #================================
    # find and open input file 
    global ev

    data_eSum = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0]])
    sumTime = 0

    number = resolveNumber(dataPath, ident)
    fileName, fileNumber = fileNameCreator(dataPath, str(number))

    if verbous:
        logging.info('life_histogrammer processing file ->\033[1m {} \033[0m<-'.format(fileNumber))

    for i in range(6):
        ev = h5py.File(fileName, 'r', swmr=True)
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
    meta = Meta(fileName)

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
        chPh
    except NameError:
        chPh = meta.chPh
    spin = meta.spin

    global countingTime, detectorDistance, chopperDetectorDistance
    detectorDistance = meta.detectorDistance
    chopperDetectorDistance = meta.chopperDetectorDistance
    countingTime = meta.countingTime

    if verbous:
        logging.info("   mu = {:>4.2f} deg,  nu = {:>4.2f} deg".format(mu, nu))
        if spin == 'u':
            logging.info('   spin <+|')
        elif spin == 'd':
            logging.info('   spin <-|')

    try: lamdaMax
    except NameError: lamdaMax = lamdaMin + tau * hdm/chopperDetectorDistance * 1e13
    
    tofOffset = tau * chopperPhase / 180.                       # mismatch of chopper pulse and time-zero
    tofCut = lamdaCut * chopperDetectorDistance / hdm * 1.e-13  # tof of frame start

    tof_e  = np.array(ev['/entry1/Amor/detector/data/event_time_offset'][:], dtype=np.uint64)/1.e9 + tofOffset # tof 
    
    detPixelID_e = np.array(ev['/entry1/Amor/detector/data/event_id'][:], dtype=np.uint64)  # pixel index
    
    dataPacket_p = np.array(ev['/entry1/Amor/detector/data/event_index'][:], dtype=np.uint64)  # data packet index

    tof_e  = np.where(tof_e<tofCut, tof_e+2.*tau, tof_e)
    tof_e  = np.where(tof_e>tau+tofCut, tof_e-tau, tof_e)

    data_e = analyse_ev(detPixelID_e, tof_e, yMin, yMax, thetaMin, thetaMax)

    ev.close()

    data_eSum = np.append(data_eSum, data_e, axis=0)
    sumTime += countingTime

    if verbous:
        logging.info("   total counts = {} in {:6.1f} s".format(np.shape(data_e)[0], sumTime))

    #================================
    # plotting data
    plotType = clas.plot[0]
    try:
        arg = clas.plot[1]
    except IndexError:
        arg = 'def'
    plott = PlotSelection()
    try:
        plotFunction = getattr(plott, plotType)
        plotFunction(fileNumber, arg, data_e)
        plt.close()
    except Exception as e:
        logging.error(f"ERROR: '{plotType}' is no known output format!")
        logging.error(f"       original error: {e}")

#==============================================================================   
def commandLineArgs():    
    msg = "events2histogram reads the eventstream from an hdf raw file and \
           creates various histogrammed outputs or plots."
    clas = argparse.ArgumentParser(description = msg)

    clas.add_argument("-c", "--chopperSpeed",
                            type=float,
                            help ="chopper speed in rpm")
    clas.add_argument("-d", "--dataPath",  
                            help ="relative path to directory with .hdf files")
    clas.add_argument("-f", "--fileIdent",       
                            default='0',                               
                            help ="file number or offset (if negative)")
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
                            default=-5.,                     
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
    clas.add_argument("-T", "--TOFOffset",       
                            default=0.0,                    
                            type=float,
                            help ="TOF zero offset")
    clas.add_argument("-t", "--thetaRange",      
                            default=[-12., 12.],   
                            nargs=2, 
                            type=float,
                            help ="theta range to be used")
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
    elif os.path.exists('./raw'):
        dataPath  = "./raw/"
    elif os.path.exists('../raw'):
        dataPath  = "../raw/"
    else:
        sys.exit('# *** please provide the path to the .hdf data files (-d <path>, default is "./raw")')

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

    dataPath = get_dataPath(clas)
    logging.basicConfig(level=logging.INFO, format='# %(message)s')
    verbous = True
    output = 'life_plot.png'
    process(dataPath, clas.fileIdent, clas)

#==============================================================================
#==============================================================================
if __name__ == "__main__":
    main()     


