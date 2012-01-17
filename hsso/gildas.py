#/usr/bin/python
"""CLASS tools"""

import numpy as np
from scipy import interpolate
from scipy import constants
from scipy.integrate import simps
from matplotlib.dates import date2num
import matplotlib.pyplot as plt
from datetime import datetime

def frac_day(dt):
    """Calculate fractional day

    Parameters
    ----------
    dt : datetime.datetime object
    """
    return dt.day + (dt-datetime.datetime(dt.year, dt.month, dt.day))/datetime.timedelta(1)

def gauss(sigma, mu=0):
    return lambda x: np.exp(-(x-mu)**2/(2.*sigma**2))/sigma/np.sqrt(2*np.pi)

def stitchn(freq, flux, goodval=False):
    """stitch subbands

    Parameters
    ----------
    freq : list of arrays
        input frequencies
    flux : list of arrays
        input fluxes
    """
    wave = np.concatenate([i for i in freq ])
    _flux = np.concatenate([i for i in flux ])
    # sort frequencies
    sortval = np.argsort(wave)
    stitchedFreq = wave[sortval]
    stitchedFlux = _flux[sortval]
    if goodval: # return finite fluxes
        goodval = np.isfinite(stitchedFlux)
        return stitchedFreq[goodval], stitchedFlux[goodval]
    else:
        return stitchedFreq, stitchedFlux

def stitch(freq1, flux1, freq2, flux2, goodval=False):
    """stitch subbands"""
    # sort frequencies
    sortval1 = np.argsort(freq1)
    sortval2 = np.argsort(freq2)
    if freq1[sortval1][0] < freq2[sortval2][0]:
        idx = np.where(freq2[sortval2] > freq1[sortval1][-1])
        stitchedFreq = np.append(freq1[sortval1], freq2[sortval2][idx])
        stitchedFlux = np.append(flux1[sortval1], flux2[sortval2][idx])
    else:
        raise Exception
    if goodval: # return finite fluxes
        goodval = np.isfinite(stitchedFlux)
        return stitchedFreq[goodval], stitchedFlux[goodval]
    else:
        return stitchedFreq, stitchedFlux

def average(freqh, fluxh, freqv, fluxv, goodval=False):
    # sort frequencies
    sortval = np.argsort(freqv)
    f = interpolate.interp1d(freqv[sortval], fluxv[sortval], bounds_error=False)
    newfluxv = f(freqh)
    fluxav = (fluxh+newfluxv)/2.
    if goodval:
        # return finite fluxes
        goodval = np.isfinite(fluxav)
        return freqh[goodval], fluxav[goodval]
    else:
        return freqh, fluxav

def averagen(wave, flux, goodval=False):
    """Average spectra
    
    Parameters
    ----------
    wave : list of arrays
        input frequencies
    flux : list of arrays 
        input fluxes
    goodval : bool, optional
        return finite values
    """
    fluxav = 0.
    for i in range(len(wave)):
        # sort frequencies
        sortval = np.argsort(wave[i])
        f = interpolate.interp1d(wave[i][sortval], flux[i][sortval], bounds_error=False)
        fluxav += f(wave[0])
    # average total flux
    fluxav /= len(flux)
    if goodval:
        # return finite fluxes
        goodval = np.isfinite(fluxav)
        return wave[0][goodval], fluxav[goodval]
    else:
        return wave[0], fluxav

def basepoly(v, flux, lim, deg=2, debug=False):
    """Calculate baseline"""
    mask = np.where((np.abs(v) < lim[1]) & (np.abs(v) > lim[0]))
    func = np.poly1d(np.polyfit(v[mask], flux[mask], deg))
    if debug:
        plt.plot(v[mask], flux[mask])
        plt.plot(v[mask], func(v[mask]))
        plt.show()
        plt.close()
    return func

def vel(freq, freq0, deltadot=0.):
    """return velocity scale
    
    deltadot: km/s"""
    return constants.c*1e-3 * (freq0-freq)/freq0 - deltadot

def freq(vel, freq0, deltadot=0.):
    """return frequency scale"""
    return freq0*(1. - (vel + deltadot)/constants.c/1e-3)

def intens(flux, vel, lim=[-1.2, 1.2], rmslim=[2,5]):
    """return intensity in K km/s with statistical error

    Parameters
    ----------
    flux : 1-D array
        Array containing flux data
    vel : 1-D array
        Array containing velocities data

    Returns
    -------
    intensity, rms error : float, float
    """
    # sort velocities
    sortval = np.argsort(vel)
    idx = np.where((vel[sortval] >= lim[0]) & (vel[sortval] <= lim[1]))
    delv = np.average(vel[sortval][idx][1:] - vel[sortval][idx][:-1])
    n = len(flux[sortval][idx])
    stderr = np.sqrt(n) * delv * rms(flux, vel, rmslim)
    return simps(flux[sortval][idx], vel[sortval][idx]), stderr

def rms(flux, vel, lim=[2, 6]):
    """return rms between [xi,xo]"""
    mask = np.where((vel > lim[0]) & (vel < lim[1]))
    return np.std(flux[mask])

def vshift(flux, vel, lim=[-1.2, 1.2]):
    """Calculate velocity offset as weighted average

    Returns
    -------
    out: tuple of arrays
        Returns a tuple with velocity offset and error arrays in m/s
    """
    mask = np.where((vel > lim[0]) & (vel < lim[1]))
    n = len(flux[mask])
    stderr = np.sqrt(np.sum(vel[mask]**2)*np.sum(flux[mask])**2 +
            n*np.dot(vel[mask], flux[mask])**2 -
            np.sum(vel[mask])*np.sum(flux[mask])*np.dot(vel[mask], flux[mask]))/ \
            np.sum(flux[mask])**2 * rms(flux, vel, [2, 20])
    return np.average(vel[mask], weights=flux[mask])*1e3, stderr*1e3

def deltadot(middate, filename="/home/miguel/HssO/Wild2/horizons.txt", column=5):
    """calculate a quantity from JPL horizons at mid-point of the exposures

    input:
        middate: date of the observations obtained from the FITS header
        filename: tabular data from horizons
        columns: quantity to interpolate
    """
    # read year, month, day and hour minute from the table
    ymd, hm = np.loadtxt(filename, dtype='S', usecols=(0,1), unpack=True)
    # read quantity from JPL Horizons ephemris file
    vdot = np.loadtxt(filename, usecols=(column,), unpack=True)
    datenum = np.zeros(len(ymd))
    for i in range(len(ymd)):
        datenum[i] = date2num(datetime.strptime(ymd[i]+hm[i], "%Y-%b-%d%H:%M"))
    # 1D interpolation
    return interpolate.interp1d(datenum, vdot)(date2num(middate))

def wcs(cdelt1, cdelt2, crpix1, crpix2, crval1, crval2):
    """return world coordinate system"""
    RA = crval1 + cdelt1*(np.arange(1, 4) - crpix1)
    dec = crval2 + cdelt2*(np.arange(1, 4) - crpix2)
    return RA, dec

def movaver(x, y, window_len=4):
    """Moving average and resampling"""
    if x[1] < x[0]:
        x = x[::-1]
        y = y[::-1]
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    w = np.ones(window_len, 'd')
    s = np.convolve(w/w.sum(), y, mode='same')
    fint = interpolate.interp1d(x, s)
    skip = window_len/2 + 1
    newx = np.linspace(x[skip], x[-skip], x.shape[0]/window_len)
    return newx, fint(newx)

def movav(x, y, window_len=4):
    """Moving average and resampling"""
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if x.size%window_len:
        newx = x[:-(x.size%window_len)].reshape(x.size/window_len, window_len)
        newy = y[:-(x.size%window_len)].reshape(x.size/window_len, window_len)
    else:
        newx = x.reshape(x.size/window_len, window_len)
        newy = y.reshape(x.size/window_len, window_len)
    return np.average(newx, axis=1), np.average(newy, axis=1)

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def RAtodeg(pos):
    return (pos[0]+pos[1]/60.+pos[2]/3600.)*180./12.

def dectodeg(pos):
    return pos[0]+pos[1]/60.+pos[2]/3600.
