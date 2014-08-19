import numpy as np

def weight_time(header):
    """weight Time, for Time*abs(Fres)/(Tsys**2)"""
    return header['EXPOSURE']*np.abs(header['FRES'])/header['TSYS']**2 * 1e6

def frequency(header, imagfreq=False):
    """define frequency array

    Parameters
    ----------
    header: dict
        header output from the read_class function.
    imagfreq : bool
        create an array with the image frequency.

    Returns
    -------
    out : array
        frequency array of evenly spaced values in GHz.
    """
    rest_frequency = header['RESTF']
    nchan = header['NCHAN']
    voff = header['VOFF']
    foff = header['FOFF']
    fres = header['FRES']
    refchan = header['RCHAN']
    imfreq = header['IMAGE']
    doppler = header['DOPPLER']
    if not imagfreq:
        freq =  rest_frequency + (np.arange(1, nchan+1) - refchan) * fres
    else:
        freq =  imfreq - (np.arange(1, nchan+1) - refchan) * fres
    return freq*1e-3

def spectra_mask(indexes, source, line, tel):
    """select spectra indexes matching filters"""
    mask = [index for (index, d) in enumerate(indexes) if
            source  in d['XSOURC'] and
            line in d["XLINE"] and
            'SMT-{0}'.format(tel) in d['XTEL']]
#             tel in d['XTEL']]
#             d['XTEL'] == 'SMT-10M-{0}'.format(tel)]
    return mask

def pgram_peaks(freq, flux, f, num):
    from scipy.signal import lombscargle
    normval = freq.shape[0]
    pgram = lombscargle(freq, flux, f)/normval
    pgram = np.sqrt(4*(pgram/normval))
    # find local maxima not at the edge of the periodogram
    maxmask = np.r_[False, pgram[1:] > pgram[:-1]] &\
                np.r_[pgram[:-1] > pgram[1:], False]
    sortarg = np.argsort(pgram[maxmask])
    peak_freqs = f[maxmask][sortarg[-num:]]
    peak_flux = pgram[maxmask][sortarg[-num:]]
    return pgram, peak_freqs, peak_flux

def linfunc(a, x, peak_freqs):
    """Target function

    Parameters
    ----------
    a: array
        Coefficients of sin and cos functions
    x: array
        Frequency values
    peak_freqs: array
        Peak frequency values
    """
    num = len(peak_freqs)
    sinwave = [a[i]*np.sin(peak_freqs[i]*x) for i in range(len(peak_freqs))]
    coswave = [a[num+i]*np.cos(peak_freqs[i]*x) for i in range(len(peak_freqs))]
    return np.column_stack(sinwave+coswave)

def fitfunc(p, x, peak_freqs):
    """Target function"""
    sinwave = [p[2*i]*np.sin(peak_freqs[i]*x) + p[2*i+1]*np.cos(peak_freqs[i]*x)
                for i in range(len(peak_freqs))]
    return np.sum(sinwave, axis=0)

def fft_peaks(freq, flux, num):
    """Calculate power spectrum using FFT"""
    sample_freq = fftpack.fftfreq(flux.size, d=np.abs(freq[0]-freq[1]))
    sig_fft = fftpack.fft(flux)
    # Because the resulting power is symmetric, only the positive part of the
    # spectrum needs to be used for finding the frequency
    pidxs = np.where(sample_freq > 0)
    f = sample_freq[pidxs]
    pgram = np.abs(sig_fft)[pidxs]
    # find local maxima not at the edge of the periodogram
    maxmask = np.r_[False, pgram[1:] > pgram[:-1]] &\
                np.r_[pgram[:-1] > pgram[1:], False]
    sortarg = np.argsort(pgram[maxmask])
    peak_freqs = f[maxmask][sortarg[-num:]]
    peak_flux = pgram[maxmask][sortarg[-num:]]
    return f, pgram, peak_freqs, peak_flux

def astroML_peaks(freq, flux, f, num):
    from astroML.time_series import lomb_scargle
    dy = np.ones(len(flux))*1e-3
    pgram = lomb_scargle(freq, flux, dy, 2*np.pi*f, generalized=False)
    # find local maxima not at the edge of the periodogram
    maxmask = np.r_[False, pgram[1:] > pgram[:-1]] &\
                np.r_[pgram[:-1] > pgram[1:], False]
    sortarg = np.argsort(pgram[maxmask])
    peak_freqs = f[maxmask][sortarg[-num:]]
    peak_flux = pgram[maxmask][sortarg[-num:]]
    return f, pgram, peak_freqs, peak_flux
