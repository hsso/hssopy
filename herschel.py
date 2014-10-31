#!/usr/bin/python
"""Misc functions"""

from scipy import constants, ndimage
import numpy as np
from astropy.io import fits as pyfits
from os.path import join
from hsso import gildas
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import astropy.wcs as pywcs

freq = {'H2O': 556.9359877}

# http://herschel.esac.esa.int/Docs/HIFI/html/ch05s05.html#table-efficiencies
# beam efficiency v0.1
# 0.754 for 557 GHz
# 0.742 for 988 GHz
# 0.715 for 1113 GHz
# 0.710 for 1669 GHz
# beameff = (0.76, 0.75, 0.75, 0.74, 0.64, 0.72, 0.69)

# http://herschel.esac.esa.int/Docs/TechnicalNotes/HIFI_Beam_Efficiencies_17Nov2010.pdf
# efficiencies v1
beameff = (0.76, 0.76, 0.76, 0.76, 0.66, 0.76, 0.76)
sigma = 3.8e-6 # m

# http://herschel.esac.esa.int/twiki/pub/Public/HifiCalibrationWeb/HifiBeamReleaseNote_Sep2014.pdf
# efficiencies v2
ruze_scaling = {
'1H': (0.649, 0.007, 2.857, 0.004, 0.624, 0.007, 2.278, 0.004),
'1V': (0.632, 0.007, 2.834, 0.003, 0.618, 0.007, 2.248, 0.004),
'2H': (0.644, 0.009, 2.807, 0.003, 0.643, 0.009, 2.219, 0.003),
'2V': (0.666, 0.009, 2.812, 0.003, 0.662, 0.009, 2.224, 0.002),
'3H': (0.633, 0.008, 2.817, 0.005, 0.626, 0.008, 2.227, 0.004),
'3V': (0.667, 0.008, 2.864, 0.006, 0.639, 0.008, 2.287, 0.010),
'4H': (0.645, 0.008, 2.817, 0.005, 0.639, 0.008, 2.227, 0.004),
'4V': (0.659, 0.008, 2.838, 0.005, 0.643, 0.008, 2.252, 0.006),
'5H': (0.548, 0.006, 2.696, 0.006, 0.595, 0.006, 2.068, 0.005),
'5V': (0.564, 0.006, 2.728, 0.006, 0.597, 0.006, 2.114, 0.005),
'6H': (0.605, 0.007, 2.831, 0.007, 0.593, 0.007, 2.243, 0.012),
'6V': (0.615, 0.007, 2.866, 0.007, 0.587, 0.007, 2.291, 0.013),
'7H': (0.585, 0.007, 2.806, 0.009, 0.584, 0.007, 2.219, 0.010),
'7V': (0.622, 0.007, 2.817, 0.010, 0.615, 0.007, 2.227, 0.007),
}

def fwhm(freq=freq['H2O'], diam=3.5, unit='arcsec'):
    """calculate FWHM in radians"""
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
    return glob.glob(
        join(datadir, str(obsid), 'level2',
        '{0}-{1}-{2}'.format(backend, pol, sideband),
        'box_001', '*.fits*'))[0]

def pacsfits(datadir, obsid, band):
    """
    Return FITS file with PACS image

    datadir: str
        full path to data directory
    obsid: int
        Herschel obsID
    band: str
        'b' (blue band) or 'r' (red band)
    """
    return glob.glob(join(datadir, str(obsid), 'level2',
    'HPPPMAP{}'.format(band.upper()), '*fits.gz'))[0]

def ruze(band, pol, wave, release):
    """
    Ruze scaling function
    """
    if release > 1:
        eta_mb, sigma_mb = ruze_scaling['{0}{1}'.format(band, pol)][4::2]
        sigma_mb *= 1e-6 # m
    else:
        eta_mb = beameff[band-1]
        sigma_mb = sigma
    return eta_mb * np.exp(-(4*np.pi*sigma_mb/wave)**2)

def wave(freq):
    """calculate wavelength"""
    return constants.c/freq

class HIPESpectrum(object):

    def __init__(self, hdus, subband=1, freq0=freq['H2O'],
            j=1, k=0):
        """Read data from FITS file produced by HIPE"""
        if not isinstance(hdus, pyfits.HDUList): hdus = pyfits.open(hdus)
        for i in hdus[j].header.keys():
            if 'key.META_' in i:
                if hdus[j].header[i] == 'sideband':
                    self.sideband = hdus[j].header[i[4:]]
        self.freq = hdus[j].data.field('{0}frequency_{1}'.format(
                    self.sideband.lower(), subband))[k]
        self.flux = hdus[j].data.field('flux_{0}'.format(subband))[k]
        self.freq0 = freq0
        self.ra = hdus[j].data.field('longitude')[k]
        self.dec = hdus[j].data.field('latitude')[k]

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
        self.ra = np.average((self.ra, spectrum.ra))
        self.dec = np.average((self.dec, spectrum.dec))

    def mask(self, line, shift, linelim, baselim):
        maskline = np.where(np.abs(self.vel - line - shift) < linelim)
        maskvel = np.where((np.abs(self.vel - line - shift) < baselim) &
                                (np.abs(self.vel - line - shift) > linelim))
        func = np.poly1d(np.polyfit(self.freq[maskvel],
                                self.baseflux[maskvel], 3))
        return maskline, maskvel, func

    def fftbase(self, fftlim, line=(0,), shift=0, linelim=(1,), baselim=(3,),
                plot=False, throw=False):
        """Fit baseline using FFT

        Parameters
        ----------
        throw : boolean
            True for unfolded data
        """
        from scipy import fftpack
        self.baseflux = self.flux.copy()
        for i in range(len(line)):
            # mask emission line
            self.maskline, self.maskvel, self.func = self.mask(line[i], shift,
                    linelim[i], baselim[i])
            self.baseflux[self.maskline] = self.func(self.freq[self.maskline])
            if throw:
                maskline, self.maskvelthrow, self.functh = self.mask(
                        self.throwvel, shift, linelim[i], baselim[i])
                self.baseflux[maskline] = self.functh(self.freq[maskline])

        # FFT
        sample_freq = fftpack.fftfreq(self.flux.size,
                    d=np.abs(self.freq[0]-self.freq[1]))
        sig_fft = fftpack.fft(self.baseflux)
        if plot:
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
        for i in range(len(line)):
            self.intens, self.error = gildas.intens(self.fluxcal,
                    self.vel - line[i] - shift, (-linelim[i], linelim[i]))
            self.vshift, self.vshift_e = gildas.vshift(self.fluxcal,
                self.vel - line[i] - shift, (-linelim[i], linelim[i]))
            self.snr = self.intens/self.error

    def plot(self, x="freq", y="flux", twiny=False, filename=None, lim=None):
        """Plot spectra

        Parameters
        ----------
        lim: int
            slice range that will be used for plotting
        """
        from matplotlib.ticker import MultipleLocator
        label = { "freq": r'$\nu$ [GHz]', "vel":'$v$ [km s$^{-1}$]' }
        if lim:
            mask = np.where(np.abs(self.__getattribute__(x)) < lim)
        else:
            mask = slice(0, -1)
        plt.plot(self.__getattribute__(x)[mask], self.__getattribute__(y)[mask],
                drawstyle='steps-mid')
        if y=="flux" and hasattr(self, "baseline"):
            mask = slice(30, -20)
            plt.plot(self.__getattribute__(x)[mask], self.baseline[mask])
            try:
                plt.scatter(self.__getattribute__(x)[self.maskvel],
                    self.func(self.freq[self.maskvel]), color='red')
                plt.plot(self.__getattribute__(x)[self.maskline],
                    self.func(self.freq[self.maskline]), 'yellow')
                plt.scatter(self.__getattribute__(x)[self.maskvelthrow],
                    self.functh(self.freq[self.maskvelthrow]), color='red')
            except (AttributeError, IndexError):
                pass
        if x == 'freq':
            plt.axvline(x=self.freq0, linestyle='--')
            if hasattr(self, 'throw'):
                plt.axvline(x=self.freq0-self.throw, linestyle='dotted')
        if 'fluxcal' in y: plt.axhline(y=0, linestyle='--')
        if x == 'vel' and np.abs(self.__getattribute__(x)[-1]) < 20:
            plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
        plt.ylabel('$T_{\mathrm{mB}}$ [K]')
        plt.xlabel(label[x])
        plt.grid(axis='both')
        plt.autoscale(axis='x', tight=True)
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
            plt.close()
        else:
            plt.show()

    def save(self, filename, flux="flux"):
        """Save spectrum to ASCII file"""
        np.savetxt(filename, np.transpose((self.freq, self.vel,
                    self.__getattribute__(flux))))

class Spectrum(HIPESpectrum):

    def __init__(self, freq, flux, freq0=freq['H2O']):
        self.freq = freq
        self.freq0 = freq0
        self.flux = flux

class HIFISpectrum(HIPESpectrum):

    def __init__(self, hdus, subband=1, byteswap=True, freq0=freq['H2O'],
            j=1, k=0, beameff=0, forwardeff=.96, release=2, pol='H'):
        """Read data from HIFI FITS data file"""
        if not isinstance(hdus, pyfits.HDUList): hdus = pyfits.open(hdus)
        self.pol = pol
        self.header = hdus[0].header
        self.freq0 = freq0
        self.obsid = hdus[0].header['OBS_ID']
        self.backend = hdus[0].header['META_0']
        for i in hdus[j].header.keys():
            if 'key.META_' in i:
                if hdus[j].header[i] == 'Band':
                    self.band = hdus[j].header[i[4:]]
                    self.bandi = int(self.band[0])
                elif hdus[j].header[i] == 'loThrow':
                    self.throw = hdus[j].header[i[4:]]
                elif hdus[j].header[i] == 'sideband':
                    self.sideband = hdus[j].header[i[4:]]
        self.freq = hdus[j].data.field('{0}frequency_{1}'.format(
                    self.sideband.lower(), subband))[k]
        self.flux = hdus[j].data.field('flux_{0}'.format(subband))[k]
        if beameff:
            self.beameff = beameff
        else:
            self.beameff = ruze(self.bandi, self.pol, wave(freq0*1e9), release)
        self.flux *= forwardeff/self.beameff
        if abs(self.throw) > 0:
            self.throwvel = gildas.vel(self.freq0-self.throw, freq0)
        self.ra = hdus[j].data.field('longitude')[k]
        self.dec = hdus[j].data.field('latitude')[k]
        try:
            self.integration = hdus[j].data.field('integration time')[k]
        except KeyError:
            self.integration = hdus[j].data.field('integrations')[k]
        try:
            self.obs_time = hdus[j].data.field('obs time')[k]
        except KeyError:
            self.obs_time = hdus[j].data.field('obsTime')[k]
        # observing (mid)-time
        self.dt = datetime(year=1958, month=1, day=1, hour=0, minute=0, second=0) \
                    + timedelta(microseconds=int(self.obs_time))
        date_obs = hdus[0].header['DATE-OBS']
        date_end = hdus[0].header['DATE-END']
        self.start = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f")
        self.end = datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S.%f")
        self.exp = self.end - self.start
        self.mid_time = self.start + self.exp/2
        if byteswap:
            self.flux = self.flux.byteswap().newbyteorder('L')
            self.freq = self.freq.byteswap().newbyteorder('L')

    def fold(self):
        if abs(self.throw) > 0:
            freq_list = [self.freq, self.freq + self.throw]
            flux_list = [self.flux, -self.flux]
            self.freq, self.flux = gildas.averagen(freq_list, flux_list,
                    goodval=True)
        else:
            print("WARNING: throw is not defined")

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

    def base(self, shift=0, lim=(2, 10), deg=1):
        """Fit polynomial baseline
        """
        basep = gildas.basepoly(self.vel-shift, self.flux, lim, deg=deg)
        # calibrated flux
        self.fluxcal = self.flux - basep(self.vel)

    def int(self, lim=(-1, 1), rmslim=[2,10]):
        """Line intensity"""
        return gildas.intens(self.fluxcal, self.vel, lim, rmslim)

    def vsh(self, lim=(-1, 1), rmslim=[2,10]):
        """Velocity shift"""
        return gildas.vshift(self.fluxcal, self.vel, lim, rmslim)

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

class Pacsmap(object):
    """
    Read PACS photometry map

    fn: string
        Ephemeris ASCII file retrieved from the HORIZONS system containing the
        following columns:
        0. date time (0 1)
        1. Astrometric RA & DEC (2 3)
        19. Heliocentric range & range-rate (4 5)
        20. Observer range & range-rate (6 7)
        27. Sun-Target radial & -vel pos. angle (8 9)
    """

    def __init__(self, obsid, size=60, zoom=0, comet=True, debug=False,
            fn="horizons.txt"):
        """
        Return patch centered on the nucleus

        Parameters
        ----------
        size: float
            size of patch is (2*size, 2*size) in arcseconds
        """
        if isinstance(obsid, pyfits.HDUList):
            self.hdus = obsid
        else:
            self.fitsfile = obsid
            self.hdus = pyfits.open(self.fitsfile)
        # pixel size
        self.cdelt2 = self.hdus[1].header['CDELT2']*3600

        pmap = self.hdus[1].data
        if comet:
            # calculate comet position at midtime
            date_obs = self.hdus[0].header['DATE-OBS']
            date_end = self.hdus[0].header['DATE-END']
            self.start = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f")
            self.end = datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S.%f")
            mid_time = self.start + (self.end-self.start)/2
            # interpolate ra and dec of the comet (start time)
            ra = gildas.deltadot(self.start, filename=fn, column=2)
            dec = gildas.deltadot(self.start, filename=fn, column=3)
            try:
                assert(np.abs(ra - self.hdus[0].header['RA_NOM']) < 4e-5)
                assert(np.abs(dec - self.hdus[0].header['DEC_NOM']) < 4e-5)
            except:
                print(np.abs(ra - self.hdus[0].header['RA_NOM'])*180/np.pi*3600)
            # calculate direction toward the Sun
            phase_ang = gildas.deltadot(self.start, filename=fn, column=8)
            self.psang = 3*np.pi/2 - phase_ang*np.pi/180
            psamv = gildas.deltadot(self.start, filename=fn, column=9)
            self.psamv = (270-psamv)*np.pi/180
            cos, sin = np.cos(self.psang), np.sin(self.psang)
            # origin coordinate is 0 (Numpy and C standards)
            # admit recognized informal extensions of WCS standard
            wcs = pywcs.WCS(self.hdus[1].header, relax=True)
            cometpix = wcs.wcs_world2pix([(ra, dec)], 0)[0]
            com = [int(round(i)) for i in cometpix]
            sh  = cometpix-com
            # shift array to center on comet nucleus
            pmap = ndimage.interpolation.shift(pmap, sh)
            self.pix = np.abs(self.cdelt2)
            self.fov = int(round(size/self.pix))
            # patch with 2fovx2fov
            self.patch = pmap[com[1]-self.fov:com[1]+self.fov+1,
                              com[0]-self.fov:com[0]+self.fov+1]
        else:
            self.patch = pmap
        if zoom > 1: self.patch = ndimage.zoom(self.patch, zoom, order=2)
        if debug:
            plt.imshow(pmap, origin="lower")
            if comet: plt.scatter(*cometpix)
            plt.grid()
            plt.show()
            plt.close()
            plt.imshow(self.patch, origin="lower", interpolation=None,
                    extent=(size, -size, -size, size))
            plt.colorbar()
            plt.grid()
            plt.show()
            plt.close()

    def center_of_mass(self, size=30, percentile=.992, debug=False):
        """
        Calculate center of mass of brightness distribution

        select the top 0.99% pixels to calculate the center of mass
        """
        fov = int(round(size/self.pix))
        center = self.patch.shape[0]/2
        zoom = self.patch[center-fov:center+fov+1,
                            center-fov:center+fov+1]
        hist, bins = np.histogram(zoom.ravel(), normed=True, bins=100)
        threshold = bins[np.cumsum(hist) * (bins[1] - bins[0]) > percentile][0]
        mpatch = np.ma.masked_less(zoom, threshold)
        mapcom = ndimage.measurements.center_of_mass(mpatch)
        # center of mass
        self.com = [int(round(i)) for i in mapcom]
        # fraction of pixel to com
        self.sh  = np.array(mapcom) - self.com
        self.com = [center - fov -.5 + i for i in mapcom]
        print(center, fov, mapcom, self.com)
        if debug:
            plt.imshow(mpatch, origin="lower")
            plt.scatter(*mapcom[::-1], color='r')
            plt.show()
            plt.close()
            plt.imshow(self.patch, origin="lower")
            mapmax = ndimage.measurements.maximum_position(self.patch)[::-1]
            plt.scatter(*mapmax)
            # plot center-of-mass
            plt.scatter(*self.com[::-1], color='r')
            plt.show()
            plt.close()

    def shift(self, center, size=30):
        """
        Shift array to be centerd at center
        """
        # comet position
        self.comet = [self.fov-center[0], self.fov-center[1]]
        self.fov = int(round(size/self.pix))
        self.patch = self.patch[center[0]-self.fov:center[0]+self.fov+1,
                                center[1]-self.fov:center[1]+self.fov+1]

    def gauss_fit(self):
        """fit 2D Gaussian"""
        pass

    def add(self, pmap):
        """average orthogonal scans"""
        self.patch = np.average((self.patch, pmap.patch), axis=0)

    def radprof(self, center=None, binsize=0, rmax=0):
        """calculate radial profile"""
        y, x = np.indices(self.patch.shape)
        if center:
            i, j = center
        else:
            j, i = np.unravel_index(np.argmax(self.patch), self.patch.shape)
        r = np.sqrt((x-i)**2 + (y-j)**2)
        ind = np.argsort(r.flat)
        r *= self.cdelt2
        sr = r.flat[ind]
        sim = self.patch.flat[ind]
        sr = sr[sim >0]
        sim = sim[sim >0]
        # normalize to Jy arcsec-2
        sim /= self.cdelt2**2
        if binsize:
            sr /= binsize
            ri = sr.astype(np.int16)
            deltar = ri[1:] - ri[:-1]
            rind = np.where(deltar)[0]
            rind = np.concatenate(([0], rind+1, [len(ri)]))
            n = rind[1:] - rind[:-1]
            self.rprof = np.array([np.mean(sim[lo:hi]) for lo,hi in
                                    zip(rind[:-1], rind[1:])])
            self.rprof_e = np.array([np.std(sim[lo:hi]) for lo,hi in
                                    zip(rind[:-1], rind[1:])])
            # reset error values
            self.rprof_e = np.where(self.rprof > self.rprof_e, self.rprof_e,
                                    0.99*self.rprof)
            self.rprof_e = np.where(self.rprof_e > 0., self.rprof_e,
                                    1e-2*self.rprof)
            self.r = binsize*(np.unique(ri)[:] + 0.5)
        else:
            self.r = sr
            self.rprof = sim
            self.rprof_e = np.zeros(len(sr))
        if rmax:
            mask = [self.r < rmax]
            self.r = self.r[mask]
            self.rprof = self.rprof[mask]
            self.rprof_e = self.rprof_e[mask]

    def plot(self, size, levels, filename=None):
        """
        Plot image
        """
        import matplotlib.cm as cm
        extent=(size, -size, -size, size)
        plt.imshow(self.patch, origin="lower", interpolation=None,
                cmap=cm.gist_heat_r, extent=extent)
        plt.colorbar()
        plt.autoscale(False)
        # zpatch = ndimage.zoom(patch, {'blue': 4, 'red': 8}[args.band])
        plt.contour(self.patch, origin="lower", levels=np.power(10, levels),
                extent=extent)
        plt.scatter(0, 0, marker='+', s=50, linewidths=1, color='white')
        plt.xlabel("RA offset [arcsec]")
        plt.ylabel("Dec offset [arcsec]")
        radius = .7*size # arcsec
        cos, sin = np.cos(self.psang), np.sin(self.psang)
        x1, y1 = radius*cos, radius*sin
        plt.arrow(x1*.4, y1*.4, x1*.6, y1*.6, fc='k', ec='k', head_width=1.5,
                head_length=2)
        cos, sin = np.cos(self.psamv), np.sin(self.psamv)
        x1, y1 = radius*cos, radius*sin
        plt.arrow(x1*.4, y1*.4, x1*.6, y1*.6, fc='k', ec='k', head_width=1.5,
                head_length=2, linestyle='dotted')
        plt.annotate(r'$\odot$', xy=(self.psang, 1.15*radius), color='black',
                xycoords='polar', va='center', ha='center')
        plt.title(self.start.strftime("%d/%m/%Y"))
        ax = plt.gca()
        ax.set_axis_off()
        extent = ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
        if filename:
            plt.savefig(filename, bbox_inches=extent)
        else:
            plt.show()
        plt.close()

    def tofits(self, filename):
        """
        Save image to FITS file
        """
        hdu = pyfits.PrimaryHDU(self.patch)
        hdu.writeto(filename, clobber=True)
