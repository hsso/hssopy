#!/usr/bin/python
"""Herschel related functions"""

from scipy import constants
import math
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits
import os
from os.path import expanduser
from hsso import gildas
from datetime import datetime, timedelta
import numpy as np

freq = {'H2O': [556.9360020]}
beameff = [.75]

def radec(header):
    """Calculate angular coordinates from FITS keywords"""
    # ra and dec FITS keywords
    naxis1 = header['naxis1'] # number of columns
    naxis2 = header['naxis2'] # number of rows
    crval1 = header['crval1']
    cdelt1 = header['cdelt1']
    crpix1 = header['crpix1']
    crval2 = header['crval2']
    cdelt2 = header['cdelt2']
    crpix2 = header['crpix2']
    # define ra and dec grids
    lon = crval1 + cdelt1*(np.arange(1, naxis1+1) - crpix1)
    lat = crval2 + cdelt2*(np.arange(1, naxis2+1) - crpix2)
    return lon, lat
    
def fwhm(freq=freq['H2O'][0]*1e9, diam=3.5, unit='arcsec'):
    """Calculate FWHM
    
    freq -- rest frequency (Hz)
    diam -- diameter (m)"""
    fwhm = 1.2*constants.c/freq/diam
    if unit=='arcsec': fwhm *= 180/math.pi*3600 # convert rad to arcsec
    return fwhm

def sigma(fwhm):
    return fwhm/(2.*np.sqrt(2*np.log(2)))

def size(arcsec, delta=1):
    """Calculate projected beam size at the comet in km
    
    arcsec -- angular size (")
    delta -- distance to Earth (AU)"""
    return delta*constants.au*math.sin(arcsec/3600.*math.pi/180)/1e3

def angsize(au, delta=1):
    """Calculate angular size
    
    au -- size (AU)
    delta -- distance to Earth (AU)"""
    return au/delta*180/math.pi*3600

def finetime(microseconds):
    timestr = datetime(year=1958,month=1,day=1,hour=0,minute=0,second=0) + \
        timedelta(microseconds=microseconds)
    return timestr

class HifiMap(object):
    """Calculate line intensity map and coordinates"""

    def __init__(self):
        """Initial empty arrays"""
        self.longitudes = np.array([])
        self.latitudes = np.array([])
        self.midtime = np.array([])
        self.start = np.array([])
        self.fvals = np.array([])

    def add(self, filename, freq0, sideband='USB', subband=1):
        """Add data from map"""
        self.filename = filename
        # read HSA FITS file
        hdulist = pyfits.open(filename)
        self.ntables = len(hdulist)
        self.npoints = hdulist[1].data.field('flux_1').shape[0]
        self.sigma = sigma(fwhm(freq=freq0*1e9))
        for hdu in hdulist[1:]:
            # mid-date observing time in UT
            self.midtime = np.append(self.midtime,
                [datetime(year=1958, month=1, day=1) +
                timedelta(microseconds=int(ot)) for ot in
                hdu.data.field('obs time')])
            # interpolate and subtract ra and dec of the comet
            self.longitudes = np.append(self.longitudes,
                    hdu.data.field('longitude'))
            self.latitudes = np.append(self.latitudes,
                    hdu.data.field('latitude'))
            # read frequency and flux
            freq = hdu.data.field('{0}frequency_{1}'.format(
                sideband.lower(), subband))
            flux = hdu.data.field('flux_{0}'.format(subband))
            for f, fl, in zip(freq, flux):
                vel = gildas.vel(f, freq0)
                mask = [np.abs(vel) < 5]
                if not hasattr(self, "spec"):
                    self.spec = fl[mask]
                    self.vel = vel[mask]
                else:
                    self.spec = np.vstack((self.spec, fl[mask]))
                    self.vel = np.vstack((self.vel, vel[mask]))
            # define line intensity map
#             self.fvals = np.append(self.fvals,
#                             gildas.intens(flux, vel, [-.5, 1.])[0])

    def plot(self):
        import matplotlib.pyplot as plt
        for i,j in zip(self.vel, self.spec):
            plt.plot(i, j, drawstyle="steps-mid")
        plt.show()

    def correct(self,
            horizons=expanduser("~/HssO/Hartley2/python/horizons.txt")):
        self.longitudes -= np.vectorize(gildas.deltadot)(self.midtime,
                filename=horizons, column=2)
        self.longitudes *= np.cos(self.latitudes * math.pi / 180.)
        self.latitudes -= np.vectorize(gildas.deltadot)(self.midtime,
                filename=horizons, column=3)
        # convert degrees to arcsec
        self.longitudes *= 3600
        self.latitudes *= 3600

    def grid(self, ncell, beameff=.74):
        """grid the data to a uniform grid using Gaussian kernel weights"""
        from scipy.stats import norm
        import matplotlib.pyplot as plt
        xi = np.linspace(self.longitudes.min(), self.longitudes.max(), ncell)
        yi = np.linspace(self.latitudes.min(), self.latitudes.max(), ncell)
        zi = np.zeros((ncell, ncell))
        vel = np.average(self.vel, axis=0)
        for i in range(ncell):
            for j in range(ncell):
                dist = np.sqrt((self.longitudes - xi[i])**2 +
                        (self.latitudes - yi[j])**2)
                sortidx = np.argsort(dist)
                diff = np.append(dist[sortidx[0]], dist[sortidx[1:]] -
                            dist[sortidx[:-1]])
#                 flux = np.sum(self.spec[sortidx,:]*diff[:,np.newaxis]*
#                         norm.pdf(dist[sortidx], 0, self.sigma)[:,np.newaxis],
#                         axis=0)
                flux = np.average(self.spec,
                        weights=norm.pdf(dist, 0, self.sigma), axis=0)
                basep = gildas.basepoly(vel, flux, [1.4, 5], deg=1)
                flux -= basep(vel)
                zi[i,j] = gildas.intens(flux, vel, [-.1, 2.])[0]
        plt.show()
        zi *= .96/beameff
        return xi, yi, zi

    def subtract(self):
        # subtract baseline
        for vel, flux in zip(self.vel, self.spec):
            basep = gildas.basepoly(vel, flux, [1.4, 5], deg=1)
            flux -= basep(vel)

    def griddata(self, ncell):
        """grid the data to a uniform grid using interpolate.griddata"""
        from scipy import interpolate
        xi = np.linspace(self.longitudes.min(), self.longitudes.max(), ncell)
        yi = np.linspace(self.latitudes.min(), self.latitudes.max(), ncell)
        zi = interpolate.griddata((self.longitudes, self.latitudes),
            self.fvals/.75, (xi[None,:], yi[:,None]), method='cubic')
        return xi, yi, zi

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
                    filename=expanduser("~/HssO/Hartley2/python/horizons.txt"),
                    column=2)
                latitudes[j-1,k] -= gildas.deltadot(timestr,
                    filename=expanduser("~/HssO/Hartley2/python/horizons.txt"),
                    column=3)
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
    """find HSA FITS file"""
#     dirname = '{0}/{1}/level2/{2}/'.format(datadir, obsid, reciever)
#     dirname = os.path.join(datadir, str(obsid), 'level2', reciever)
    dirname = os.path.join(datadir, str(obsid), reciever)
    fitsfile = os.listdir(dirname)[0]
    return os.path.join(dirname, fitsfile)
