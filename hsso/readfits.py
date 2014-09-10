#/usr/bin/python
"""Module to read HIPE fits files"""

import pyfits
import numpy as np
import physcon as pc

class HipeFits:
    """read fits output from HiClass"""
    def __init__(self, filename, j=1):
        self.filename = filename
        self.hdulist = pyfits.open(filename)
        self.cols = self.hdulist[j].columns
        self.hdr = self.hdulist[j].header
        if self.hdulist[j].header.has_key('telescop'):
            self.telescope = self.hdulist[j].header['telescop']
        elif self.hdulist[j].columns.names.__contains__('TELESCOP'):
            self.telescope = self.hdulist[j].data.field('telescop')[0]
        self.crpix1 = self.hdulist[j].header['crpix1']
        if self.crpix1 == 0: self.crpix1 = 1133
        self.cdelt1 = self.hdulist[j].header['cdelt1']*1e-9 # step
        if self.hdr.has_key('restfreq'):
            self.restfreq = self.hdr['restfreq']*1e-9
        self.badval = self.hdulist[j].header['blank']
        if self.hdulist[j].header['deltaf1']*1e-9 > 0:
            self.throw = self.hdulist[j].header['deltaf1']*1e-9
        else:
            self.throw = self.hdulist[j].header['deltaf2']*1e-9
#         self.dateobs = self.hdulist[j].header['date-obs']

    def close(self):
        self.hdulist.close()

    def getFlux(self, k=0, j=1):
        """read flux"""
        # remove bad channels
        x = self.hdulist[j].data.field('data')[k]
        indices = np.where(x != self.badval)
        return x[indices]

    def getFreq(self, k=0, j=1):
        """calculate frequency scale"""
        x = self.hdulist[j].data.field('data')[k]
        if self.cols.names.__contains__('RESTFREQ'):
            freq = self.hdulist[j].data.field('restfreq')[k]*1e-9 + \
                self.cdelt1*(np.arange(1, len(x)+1) - self.crpix1)
        else:
            freq = self.restfreq + self.cdelt1*(np.arange(1,
                len(x)+1) - self.crpix1)
        x = self.hdulist[j].data.field('data')[k]
        indices = np.where(x != self.badval)
        return freq[indices]

def telescop(hdulist, j, i):
    """Return telescope name"""
    if 'telescope' in hdulist[j].header:
        return hdulist[j].header['telescop']
    elif hdulist[j].columns.names.__contains__('TELESCOP'):
        return hdulist[j].data.field('telescop')[i]

def read(filename, rec="WBS", sb="USB", pol="H"):
    """Read FITS file"""
    hdulist = pyfits.open(filename)
    for j in range(1, hdulist[0].header['dsets___']+1):
        # find sideband
        if hdulist[j].header['line'].find(sb) > 0:
            # loop over subbands
            for i in range(hdulist[j].data.field('data').shape[0]):
                tel = telescop(hdulist, j, i)
                # find receiver
                if not tel: tel= "-"+rec[0]+pol
                if tel.find("-"+rec[0]+pol) >= 0:
                    # flux for subband
                    x = hdulist[j].data.field('data')[i]
                    # rest frequency in GHz
#                     restfreq = hdulist[j].data.field('restfreq')[i]*1e-9
                    restfreq = hdulist[j].header['restfreq']*1e-9
                    if tel.find("-H") >= 0: # HRS
                        crpix1 = hdulist[j].header['crpix1']
                        cdelt1 = hdulist[j].header['cdelt1']*1e-9 # step
#                         cdelt1 = hdulist[j].data.field('cdelt1')[i]*1e-9
                    elif tel.find("-W") >= 0: # WBS
                        crpix1 = hdulist[j].header['crpix1']
#             crpix1 = hdulist[j].data.field('crpix1')[i]
                        cdelt1 = hdulist[j].header['cdelt1']*1e-9 # step
                    # remove bad channels
                    badval = hdulist[j].header['blank']
                    indices = np.where(x != badval)
                    # calculate frequency scale
                    freq = restfreq + cdelt1*(np.arange(1, len(x)+1) - crpix1)
                    # short good values
                    goodval = np.argsort(freq[indices])
                    throw = hdulist[j].header['deltaf1']*1e-9
                    return freq[goodval], x[goodval], throw
