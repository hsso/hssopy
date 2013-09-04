#!/usr/bin/python
"""Misc functions"""

import physcon as pc
import astrocon as ac
import numpy as np
import pyfits
from os.path import join
from hsso import gildas

freq = {'H2O': [556.9360020]}
beameff = [.75]

def fwhm(freq=freq['H2O'], diam=3.5, unit='arcsec'):
    """calculate FWHM"""
    fwhm = 1.2*pc.c/freq/diam
    if unit=='arcsec': fwhm *= 180/np.pi*3600 # convert rad to arcsec
    return fwhm

def size(arcsec, delta=1): # km
    """calculate projected beam size"""
    return delta*ac.AU*np.sin(arcsec/3600.*np.pi/180)/1e5

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

    def __init__(self, fitsfile, subband=1, byteswap=True, freq0=556.9359877):
        hdus = pyfits.open(fitsfile)
        self.freq0 = freq0
        self.obsid = hdus[0].header['OBS_ID']
        self.backend = hdus[0].header['META_0']
        for i in hdus[1].header.ascardlist().keys():
            if hdus[1].header[i] == 'loThrow':
                self.throw = hdus[1].header[i[4:]]
            elif hdus[1].header[i] == 'sideband':
                self.sideband = hdus[1].header[i[4:]]
        self.freq = hdus[1].data.field('{0}frequency_{1}'.format(self.sideband.lower(),
                                    subband))[0]
        self.flux = hdus[1].data.field('flux_{0}'.format(subband))[0]
        self.vel = gildas.vel(self.freq, freq0)
        if byteswap:
            self.flux = self.flux.byteswap().newbyteorder('L')
            self.freq = self.freq.byteswap().newbyteorder('L')

    def save(self, filename, flux="flux"):
        """Save spectrum to ASCII file"""
        np.savetxt(filename, np.transpose((self.freq, self.vel,
                    self.__getattribute__(flux))))

    def add(self, spectrum):
        if np.all(self.freq == spectrum.freq):
            self.flux += spectrum.flux
            self.flux /= 2
        else:
            freq_list = [self.freq, spectrum.freq]
            flux_list = [self.flux, spectrum.flux]
            self.freq, self.flux = gildas.averagen(freq_list, flux_list, goodval=True)
            self.vel = gildas.vel(self.freq, self.freq0)

    def __add__(self, spectrum):
        return self.add(spectrum)

    def fold(self):
        freq_list = [self.freq, self.freq + self.throw]
        flux_list = [self.flux, -self.flux]
        self.freq, self.flux = gildas.averagen(freq_list, flux_list, goodval=True)
        self.vel = gildas.vel(self.freq, self.freq0)

    def scale(self, vel_lim=None):
        if vel_lim:
            maskvel = np.where((self.vel < vel_lim[1]) & (self.vel > vel_lim[0]))
            self.flux -= np.mean(self.flux[maskvel])
        else:
            self.flux -= np.mean(self.flux)

    def fftbase(self, fftlim, shift=0, linelim=1, baselim=3, plot=False):
        from scipy import fftpack
        self.baseflux = self.flux.copy()
        maskline = np.where(np.abs(self.vel - shift) < linelim)
        maskvel = np.where((np.abs(self.vel - shift) < baselim) &
                            (np.abs(self.vel - shift) > linelim))
        func = np.poly1d(np.polyfit(self.freq[maskvel], self.baseflux[maskvel], 3))
        self.baseflux[maskline] = func(self.freq[maskline])

        # FFT
        sample_freq = fftpack.fftfreq(self.flux.size, d=np.abs(self.freq[0]-self.freq[1]))
        sig_fft = fftpack.fft(self.baseflux)
        sig_fft[np.abs(sample_freq) > fftlim] = 0
        if args.debug:
            pidxs = np.where(sample_freq > 0)
            f = sample_freq[pidxs]
            pgram = np.abs(sig_fft)[pidxs]
            plt.loglog(f, pgram)
            plt.axvline(x=fftlim, linestyle='--')
            plt.show()
        self.baseline = np.real(fftpack.ifft(sig_fft))
        # calibrated flux
        self.fluxcal = self.flux - self.baseline
        self.fluxcal *= 0.96/.75
