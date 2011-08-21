#!/usr/bin/python
"""Herschel functions"""

from scipy import constants
import astrocon as ac
import math
import pyfits
import os
from hsso import gildas
from datetime import datetime, timedelta
import numpy as np

freq = {'H2O': [556.9360020]}
beameff = [.75]

def fwhm(freq=freq['H2O'][0]*1e9, diam=3.5, unit='arcsec'):
    """Calculate FWHM
    
    freq: rest frequency (Hz)
    diam: diameter (m)"""
    fwhm = 1.2*constants.c/freq/diam
    if unit=='arcsec': fwhm *= 180/math.pi*3600 # convert rad to arcsec
    return fwhm

def size(arcsec, delta=1): # km
    """Calculate projected beam size
    
    arcsec: angular size
    delta: distance to Earth (AU)"""
    return delta*ac.AU*math.sin(arcsec/3600.*math.pi/180)/1e5

class HifiMap(object):
    """Calculate line intensity map and coordinates"""

    def __init__(self, filename, freq0, sideband='USB', subband=1, correct=True):
        self.filename = filename
        # read HSA FITS file
        hdulist = pyfits.open(filename)
        self.ntables = len(hdulist)
        self.npoints = hdulist[1].data.field('flux_1').shape[0]
        self.longitudes = np.zeros((self.ntables-1, self.npoints))
        self.latitudes = np.zeros((self.ntables-1, self.npoints))
        self.fvals = np.zeros((self.ntables-1, self.npoints))
        for j in range(1, self.ntables):
            for k in range(self.npoints):
                # mid-date observing time in UT
                integration_time = hdulist[j].data.field('integration time')[k]
                # if the integration time is not a scalar
                if not isinstance(integration_time, float):
                    integration_time = integration_time[subband-1]
                timestr = datetime(year=1958,month=1,day=1,hour=0,minute=0,second=0) + \
                        timedelta(microseconds=hdulist[j].data.field('obs time')[k] + \
                        integration_time/2.)
                # interpolate and subtract ra and dec of the comet
                self.longitudes[j-1,k] = hdulist[j].data.field('longitude')[k] 
                self.latitudes[j-1,k] = hdulist[j].data.field('latitude')[k] 
                if correct:
                    self.longitudes[j-1,k] -=  gildas.deltadot(timestr,
                        filename="/home/miguel/HssO/Hartley2/python/horizons.txt", column=2)
                    self.latitudes[j-1,k] -= gildas.deltadot(timestr,
                        filename="/home/miguel/HssO/Hartley2/python/horizons.txt", column=3)
                # read frequency and flux
                freq = hdulist[j].data.field('{0}frequency_{1}'.format(sideband.lower(), subband))[k]
                flux = hdulist[j].data.field('flux_{0}'.format(subband))[k]
                vel = gildas.vel(freq, freq0)
                # subtract baseline
                basep = gildas.basepoly(vel, flux, [1.4, 5], deg=1)
                flux -= basep(vel)
                # define line intensity map
                self.fvals[j-1, k] = gildas.intens(flux, vel, [-.5, .5])[0]

def hifimap(filename, freq0, sideband='USB', subband=1, obsid=1, correct=True):
    """Calculate intensity and coordinates"""
    # read HSA FITS file
    hdulist = pyfits.open(filename)
    ntables = len(hdulist)
    npoints = hdulist[1].data.field('flux_1').shape[0]
    longitudes = np.zeros((ntables-1, npoints))
    latitudes = np.zeros((ntables-1, npoints))
    fvals = np.zeros((ntables-1, npoints))
#     f = open('{0}.txt'.format(obsid), 'w')
    for j in range(1, ntables):
        for k in range(npoints):
            # mid-date observing time in UT
#             if subband > 0:
#                 integration_time = hdulist[j].data.field('integration time')[k][subband-1]
#             else:
            integration_time = hdulist[j].data.field('integration time')[k]
            timestr = datetime(year=1958,month=1,day=1,hour=0,minute=0,second=0) + \
                    timedelta(microseconds=hdulist[j].data.field('obs time')[k] + \
                    integration_time/2.)
            # interpolate and subtract ra and dec of the comet
            longitudes[j-1,k] = hdulist[j].data.field('longitude')[k] 
            latitudes[j-1,k] = hdulist[j].data.field('latitude')[k] 
            if correct:
                longitudes[j-1,k] -=  gildas.deltadot(timestr,
                    filename="/home/miguel/HssO/Hartley2/python/horizons.txt", column=2)
                latitudes[j-1,k] -= gildas.deltadot(timestr,
                    filename="/home/miguel/HssO/Hartley2/python/horizons.txt", column=3)
            # read frequency and flux
            freq = hdulist[j].data.field('{0}frequency_{1}'.format(sideband.lower(), subband))[k]
            flux = hdulist[j].data.field('flux_{0}'.format(subband))[k]
            vel = gildas.vel(freq, freq0)
            # subtract baseline
            basep = gildas.basepoly(vel, flux, [1.4, 5], deg=0)
            flux -= basep(vel)
            # define line intensity map
            fvals[j-1, k] = gildas.intens(flux, vel, [-1.4, 1.4])[0]
#             f.write('{0} {1} {2} {3}\n'.format(timestr.isoformat(), longitudes[j-1,k],
#                 latitudes[j-1,k], fvals[j-1,k]))
#     f.close()
    return fvals, longitudes, latitudes

def hsafits(datadir, obsid, reciever):
    """find level 2 FITS file"""
    dirname = '{0}/{1}/level2/{2}/'.format(datadir, obsid, reciever)
    fitsfile = os.listdir(dirname)[0]
    return dirname+fitsfile
