#!/usr/bin/python
"""
Plot fits CLASS files
"""

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import physcon as pc
import mars
from scipy import interpolate

# rest frequency
freq0 = 1876.2275000
# Mars velocity wrt Herschel
dvel = 15.9234171 # km/s

freqarr = np.arange(638.041, 726.029, 5e-4) # GHz
avflux = np.zeros(len(freqarr))
obsid = 1342194545
hdulist = pyfits.open('fits/'+str(obsid)+'_5fhf_level2.fits')
for j in range(1, hdulist[0].header['dsets___']+1):
    # skip LSB
    if hdulist[j].header['line'].find('LSB') > 0: continue
    # loop over subbands
    for i in range(hdulist[j].data.field('data').shape[0]):
        tel = hdulist[j].data.field('telescop')[i]
        # skip HRS
        if tel.find("-H") > 0: continue
#         if tel.find("-WV") > 0: continue
        # flux for subband
        x = hdulist[j].data.field('data')[i]
        # rest frequency in GHz
        restfreq = hdulist[j].data.field('restfreq')[i]*1e-9
        if tel.find("-H") > 0: # HRS
            crpix1 = hdulist[j].header['crpix1']
            cdelt1 = hdulist[j].data.field('cdelt1')[i]*1e-9
        elif tel.find("-W") > 0: # WBS
            crpix1 = hdulist[j].header['crpix1']
#             crpix1 = hdulist[j].data.field('crpix1')[i]
            cdelt1 = hdulist[j].header['cdelt1']*1e-9 # step
        # remove bad channels
        badval = hdulist[j].header['blank']
        indices = np.where(x != badval)
        # calculate frequency scale
        freq = restfreq + cdelt1*(np.arange(1, len(x)+1) - crpix1)
        print(j, tel, restfreq, cdelt1, crpix1, freq[0], freq[-1])
        # interpolate to freqarr
        goodval = np.argsort(freq[indices])
        f = interpolate.interp1d(freq[goodval], x[goodval])
        # find freqarr range
        ind = np.where((freqarr > freq[goodval][0]) & (freqarr < freq[goodval][-1]))
        newflux = f(freqarr[ind])
        # zero values of avflux
        tmp = avflux[ind]
        zeroind = np.where(tmp == 0)
        tmp[zeroind] = newflux[zeroind]
        # average subbands
        avind = np.where(tmp != 0)
        tmp[avind] = (tmp[avind] + newflux[avind])/2.
        avflux[ind] = tmp
        # switch to velocity scale
        vel = pc.c*1e-3 * (freq-freq0)/freq0 + dvel + mars.obsids[obsid][0]
#         plt.plot(freq[indices], x[indices], drawstyle='steps-mid')
#         plt.plot(freqarr[ind], newflux, drawstyle='steps-mid')
#         plt.plot(vel[indices], x[indices], drawstyle='steps-mid')
# plt.xlim([-10, 10])
# LSB frequency
freqarr = np.linspace(626.052, 714.0485, len(freqarr)) # GHz
freqarr *= (1. + (dvel + mars.obsids[obsid][0])/pc.c*1e3)
plt.plot(freqarr, avflux, drawstyle='steps-mid')
plt.ylim([9, 16])
plt.ylabel("T$_{mB}$ [K]")
plt.xlabel("freq [GHz]")
# plt.show()
plt.savefig("figures/%s-WBS-LSB.png" % str(obsid))
plt.close()
