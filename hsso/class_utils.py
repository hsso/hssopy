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
