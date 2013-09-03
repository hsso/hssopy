#!/usr/bin/python
"""Misc functions"""

import physcon as pc
import astrocon as ac
import math
import pyfits
from os.path import join
from hsso import gildas

freq = {'H2O': [556.9360020]}
beameff = [.75]

def fwhm(freq=freq['H2O'], diam=3.5, unit='arcsec'):
    """calculate FWHM"""
    fwhm = 1.2*pc.c/freq/diam
    if unit=='arcsec': fwhm *= 180/math.pi*3600 # convert rad to arcsec
    return fwhm

def size(arcsec, delta=1): # km
    """calculate projected beam size"""
    return delta*ac.AU*math.sin(arcsec/3600.*math.pi/180)/1e5

def fft(hdulist, sideband, subband):
    "Return frequency, flux and frequency throw"
    for i in hdulist[1].header.ascardlist().keys():
        if hdulist[1].header[i] == 'loThrow':
            throw = hdulist[1].header[i[4:]]
            break
    freq = hdulist[1].data.field('{0}frequency_{1}'.format(sideband.lower(), subband))[0]
    flux = hdulist[1].data.field('flux_{0}'.format(subband))[0]
    return freq, flux, throw

class HIFISpectrum(object):

    def __init__(self, fitsfile, sideband="USB", subband=1, byteswap=True, freq0=556.9359877):
        hdus = pyfits.open(fitsfile)
        self.obsid = hdus[0].header['OBS_ID']
        self.backend = hdus[0].header['META_0']
        self.freq = hdus[1].data.field('{0}frequency_{1}'.format(sideband.lower(),
                                    subband))[0]
        self.flux = hdus[1].data.field('flux_{0}'.format(subband))[0]
        self.vel = gildas.vel(self.freq, freq0)
        for i in hdus[1].header.ascardlist().keys():
            if hdus[1].header[i] == 'loThrow':
                self.throw = hdus[1].header[i[4:]]
                break
        if byteswap:
            self.flux = self.flux.byteswap().newbyteorder('L')/.75
            self.freq = self.freq.byteswap().newbyteorder('L')

    def save(self, datadir='/tmp'):
        """Save spectrum to ASCII file"""
        import numpy as np
        np.savetxt(join(datadir, "{}_{}.dat".format(
                    self.obsid, self.backend)),
                    np.transpose((self.freq, self.vel, self.flux)))
