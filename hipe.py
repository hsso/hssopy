#!/usr/bin/python

import pyfits
import argparse
import glob
import numpy as np
from os.path import join
import matplotlib.pyplot as plt
import jpl

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--backend', default=('WBS',), nargs='+',
                    choices=('HRS', 'WBS'))
parser.add_argument('--pol', default=('H',), nargs='+', choices=('H', 'V'))
parser.add_argument('--subband', default=range(1,5), type=int, nargs='+',
                    choices=range(1, 5))
parser.add_argument('--sideband', default='USB', choices=('USB', 'LSB'))
parser.add_argument('--datadir', default='./')
parser.add_argument('-o', '--obsid', default="")
parser.add_argument('-f', '--freq', default=None, type=float)
args = parser.parse_args()

if __name__ == "__main__":
    for be in args.backend:
        for p in args.pol:
            hdulist = pyfits.open( glob.glob(
                       join(args.datadir, args.obsid, 'level2',
                       '{0}-{1}-{2}'.format(be, p, args.sideband),
                       'box_001', '*.fits*'))[0])
            for i in hdulist[1].header.ascardlist().keys():
                if i.find('META_') > 0 and hdulist[1].header[i] == 'loThrow':
                    throw = hdulist[1].header[i[4:]]
            print throw
            for subband in args.subband:
                if 'flux_{0}'.format(subband) in hdulist[1].data.names:
                    freq = hdulist[1].data.field('{0}frequency_{1}'.format(
                            args.sideband.lower(), subband))[0]
                    flux = hdulist[1].data.field('flux_{0}'.format(subband))[0]
                    plt.plot(freq, flux, drawstyle='steps-mid',
                             label='{0}-{1} subband {2}'.format(be, p, subband))
    if args.freq:
        plt.axvline(x=args.freq, linestyle='--')
        plt.axvline(x=args.freq-throw, linestyle=':')
    # Read catalog data file
    x0, x1 = plt.gca().get_xlim()
    images = jpl.JPLmol(32003).trans
    img = images[(x0*1e3 < images['FREQ']) & (images['FREQ'] < x1*1e3) &
                (images['LGINT'] > -3.4)
                ]
    for j in img:
        print j['FREQ']*1e-3, j['LGINT']
        plt.axvline(x=j['FREQ']*1e-3, linestyle='--')
    plt.legend()
    plt.show()
