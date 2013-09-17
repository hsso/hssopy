#!/usr/bin/python
"""Misc functions"""

from scipy import constants
import astrocon as ac
import numpy as np
import pyfits
from os.path import join
from hsso import gildas

freq = {'H2O': 556.9359877}
beameff = [.75]

def fwhm(freq=freq['H2O'], diam=3.5, unit='arcsec'):
    """calculate FWHM"""
    fwhm = 1.2*constants.c/freq/1e9/diam
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

def hififits(datadir, obsid, backend, pol, sideband):
    import glob
    return glob.glob(
        join(datadir, str(obsid), 'level2',
        '{0}-{1}-{2}'.format(backend, pol, sideband),
        'box_001', '*.fits*'))[0]

class HIFISpectrum(object):

    def __init__(self, fitsfile, subband=1, byteswap=True, freq0=556.9359877):
        from datetime import datetime
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
        self.throwvel = gildas.vel(self.freq0+self.throw, freq0)
        self.ra = hdus[1].data.field('longitude')[0]
        self.dec = hdus[1].data.field('latitude')[0]
        self.integration = hdus[1].data.field('integration time')[0]
        date_obs = hdus[0].header['DATE-OBS']
        date_end = hdus[0].header['DATE-END']
        self.start = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f")
        self.end = datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S.%f")
        exp = self.end - self.start
        self.mid_time = self.start + exp/2
        if byteswap:
            self.flux = self.flux.byteswap().newbyteorder('L')
            self.freq = self.freq.byteswap().newbyteorder('L')

    def add(self, spectrum):
        if np.all(self.freq == spectrum.freq):
            self.flux += spectrum.flux
            self.flux /= 2
        else:
            freq_list = [self.freq, spectrum.freq]
            flux_list = [self.flux, spectrum.flux]
            self.freq, self.flux = gildas.averagen(freq_list, flux_list, goodval=True)
            self.vel = gildas.vel(self.freq, self.freq0)

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

    def fftbase(self, fftlim, line = (0,), shift=0, linelim=1, baselim=3, plot=False):
        from scipy import fftpack
        self.baseflux = self.flux.copy()
        for l in line:
            self.maskline = np.where(np.abs(self.vel - l - shift) < linelim)
            self.maskvel = np.where((np.abs(self.vel - l - shift) < baselim) &
                                (np.abs(self.vel - l- shift) > linelim))
            self.func = np.poly1d(np.polyfit(self.freq[self.maskvel],
                                        self.baseflux[self.maskvel], 3))
            self.baseflux[self.maskline] = self.func(self.freq[self.maskline])

        # FFT
        sample_freq = fftpack.fftfreq(self.flux.size, d=np.abs(self.freq[0]-self.freq[1]))
        sig_fft = fftpack.fft(self.baseflux)
        if plot:
            import matplotlib.pyplot as plt
            pidxs = np.where(sample_freq > 0)
            f = sample_freq[pidxs]
            pgram = np.abs(sig_fft)[pidxs]
            plt.loglog(f, pgram)
            plt.axvline(x=fftlim, linestyle='--')
            plt.show()
        sig_fft[np.abs(sample_freq) > fftlim] = 0
        self.baseline = np.real(fftpack.ifft(sig_fft))
        # calibrated flux
        self.fluxcal = self.flux - self.baseline
        self.fluxcal *= 0.96/.75

    def plot(self, flux="flux", twiny=True, filename=None, lim=None):
        import matplotlib.pyplot as plt
        if lim:
            sl = slice(lim, -lim)
        else:
            sl = slice(0, -1)
        plt.plot(self.freq[sl], self.__getattribute__(flux)[sl],
                drawstyle='steps-mid')
        if flux=="flux": plt.plot(self.freq[sl], self.baseline[sl])
        try:
            plt.plot(self.freq[self.maskvel], self.func(self.freq[self.maskvel]))
        except AttributeError:
            pass
        plt.axvline(x=self.freq0, linestyle='--')
        if hasattr(self, 'throw'): plt.axvline(x=self.freq0-self.throw, linestyle='dotted')
        plt.ylabel('$T_{\mathrm{mB}}$ [K]')
        plt.xlabel(r'$\nu$ [GHz]')
        plt.grid(axis='both')
        plt.autoscale(axis='x', tight=True)
        if twiny:
            ax1 = plt.gca()
            # update xlim of ax2
            ax2 = ax1.twiny()
            x1, x2 = ax1.get_xlim()
            ax2.set_xlim(gildas.vel(x1, self.freq0),
                         gildas.vel(x2, self.freq0))
            plt.xlabel('$v$ [km s$^{-1}$]')
        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def save(self, filename, flux="flux"):
        """Save spectrum to ASCII file"""
        np.savetxt(filename, np.transpose((self.freq, self.vel,
                    self.__getattribute__(flux))))

    def tofits(self, filename, columns=("freq", "fluxcal")):
        """Save spectrum to FITS file"""
        cols = []
        for i in columns:
            cols.append(pyfits.Column(name=i, format='E',
                        array=self.__getattribute__(i)))
        tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
        tbhdu.writeto(filename, clobber=True)

def writeto_fits(filename, columns):
    cols = []
    for i, j in columns.iteritems():
        cols.append(pyfits.Column(name=i, format='E', array=j))
    tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
    tbhdu.writeto(filename)
