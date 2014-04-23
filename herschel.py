#!/usr/bin/python
"""Misc functions"""

from scipy import constants
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
    return delta*constants.au*np.sin(arcsec/3600.*np.pi/180)*1e-3

def fft(hdulist, sideband, subband):
    "Return frequency, flux and frequency throw"
    for i in hdulist[1].header.keys():
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

    def __init__(self, fitsfile, subband=1, byteswap=True, freq0=556.9359877,
            beameff=.75):
        from datetime import datetime
        hdus = pyfits.open(fitsfile)
        self.header = hdus[0].header
        self.freq0 = freq0
        self.obsid = hdus[0].header['OBS_ID']
        self.backend = hdus[0].header['META_0']
        for i in hdus[1].header.keys():
            if hdus[1].header[i] == 'loThrow':
                self.throw = hdus[1].header[i[4:]]
            elif hdus[1].header[i] == 'sideband':
                self.sideband = hdus[1].header[i[4:]]
        self.freq = hdus[1].data.field('{0}frequency_{1}'.format(
                    self.sideband.lower(), subband))[0]
        self.flux = hdus[1].data.field('flux_{0}'.format(subband))[0]
        self.flux /= beameff
        self.throwvel = gildas.vel(self.freq0-self.throw, freq0)
        self.ra = hdus[1].data.field('longitude')[0]
        self.dec = hdus[1].data.field('latitude')[0]
        self.integration = hdus[1].data.field('integration time')[0]
        date_obs = hdus[0].header['DATE-OBS']
        date_end = hdus[0].header['DATE-END']
        self.start = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f")
        self.end = datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S.%f")
        self.exp = self.end - self.start
        self.mid_time = self.start + self.exp/2
        if byteswap:
            self.flux = self.flux.byteswap().newbyteorder('L')
            self.freq = self.freq.byteswap().newbyteorder('L')

    @property
    def vel(self):
        return gildas.vel(self.freq, self.freq0)

    def add(self, spectrum):
        if np.all(self.freq == spectrum.freq):
            self.flux += spectrum.flux
            self.flux /= 2
        else:
            freq_list = [self.freq, spectrum.freq]
            flux_list = [self.flux, spectrum.flux]
            self.freq, self.flux = gildas.averagen(freq_list, flux_list,
                    goodval=True)

    def fold(self):
        freq_list = [self.freq, self.freq + self.throw]
        flux_list = [self.flux, -self.flux]
        self.freq, self.flux = gildas.averagen(freq_list, flux_list,
                goodval=True)

    def resample(self, times=2):
        from scipy.signal import resample
        self.flux, self.freq = resample(self.flux, int(len(self.flux)/times),
                                        t=self.freq)
        if hasattr(self, "fluxcal"):
            self.fluxcal = resample(self.fluxcal,
                                        int(len(self.fluxcal)/times))
        if hasattr(self, "baseline"):
            self.baseline = resample(self.baseline,
                                        int(len(self.baseline)/times))

    def scale(self, vel_lim=None):
        """Scale flux by mean value within vel_lim"""
        if vel_lim:
            maskvel = np.where((self.vel < vel_lim[1]) &
                                (self.vel > vel_lim[0]))
            self.flux -= np.mean(self.flux[maskvel])
        else:
            self.flux -= np.mean(self.flux)

    def mask(self, line, shift, linelim, baselim):
        maskline = np.where(np.abs(self.vel - line - shift) < linelim)
        maskvel = np.where((np.abs(self.vel - line - shift) < baselim) &
                                (np.abs(self.vel - line - shift) > linelim))
        func = np.poly1d(np.polyfit(self.freq[maskvel],
                                self.baseflux[maskvel], 3))
        return maskline, maskvel, func

    def fftbase(self, fftlim, line=(0,), shift=0, linelim=1, baselim=3,
                plot=False, throw=False):
        """Fit baseline using FFT

        Parameters
        ----------
        throw : boolean
            True for unfolded data
        """
        from scipy import fftpack
        self.baseflux = self.flux.copy()
        if line:
            # mask emission line
            self.maskline, self.maskvel, self.func = self.mask(line[0], shift,
                    linelim, baselim)
            self.baseflux[self.maskline] = self.func(self.freq[self.maskline])
            if throw:
                maskline, self.maskvelthrow, self.functh = self.mask(
                        self.throwvel, shift, linelim, baselim)
                self.baseflux[maskline] = self.functh(self.freq[maskline])

        # FFT
        sample_freq = fftpack.fftfreq(self.flux.size,
                    d=np.abs(self.freq[0]-self.freq[1]))
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
        self.intens, self.error = gildas.intens(self.fluxcal, self.vel,
                                                (-linelim, linelim))
        self.vshift, self.vshift_e = gildas.vshift(self.fluxcal, self.vel,
                                                (-linelim, linelim))
        self.snr = self.intens/self.error

    def plot(self, x="freq", y="flux", twiny=False, filename=None, lim=None):
        """Plot spectra

        Parameters
        ----------
        lim: int
            slice range that will be used for plotting
        """
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
        label = { "freq": r'$\nu$ [GHz]', "vel":'$v$ [km s$^{-1}$]' }
        if lim:
            mask = np.where(np.abs(self.__getattribute__(x)) < lim)
        else:
            mask = slice(0, -1)
        plt.plot(self.__getattribute__(x)[mask], self.__getattribute__(y)[mask],
                drawstyle='steps-mid')
        if y=="flux" and hasattr(self, "baseline"):
            plt.plot(self.__getattribute__(x)[mask], self.baseline[mask])
            try:
                plt.plot(self.__getattribute__(x)[self.maskvel],
                    self.func(self.freq[self.maskvel]), 'red')
                plt.plot(self.__getattribute__(x)[self.maskline],
                    self.func(self.freq[self.maskline]), 'yellow')
                plt.plot(self.__getattribute__(x)[self.maskvelthrow],
                    self.functh(self.freq[self.maskvelthrow]), 'red')
            except (AttributeError, IndexError):
                pass
        if x == 'freq':
            plt.axvline(x=self.freq0, linestyle='--')
            if hasattr(self, 'throw'):
                plt.axvline(x=self.freq0-self.throw, linestyle='dotted')
        if y == 'fluxcal': plt.axhline(y=0, linestyle='--')
        if x == 'vel': plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
        plt.ylabel('$T_{\mathrm{mB}}$ [K]')
        plt.xlabel(label[x])
        plt.grid(axis='both')
        if twiny:
            ax1 = plt.gca()
            # update xlim of ax2
            ax2 = ax1.twiny()
            x1, x2 = ax1.get_xlim()
            ax2.set_xlim(gildas.vel(x1, self.freq0),
                         gildas.vel(x2, self.freq0))
            plt.xlabel(label["vel"])
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
