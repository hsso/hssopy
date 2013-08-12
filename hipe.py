#!/usr/bin/python

import pyfits
import argparse
import glob
import numpy as np
from os.path import join
import matplotlib.pyplot as plt

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--backend', default='WBS')
parser.add_argument('--pol', default='H')
parser.add_argument('--subband', default=range(1,5), type=int, nargs='+')
parser.add_argument('--sideband', default='USB')
parser.add_argument('--datadir', default='./')
parser.add_argument('-o', '--obsid', default="75")
parser.add_argument('-f', '--freq', default=None, type=float)
args = parser.parse_args()

hdulist = pyfits.open( glob.glob(
        join(args.datadir, args.obsid, 'level2',
        '{0}-{1}-{2}'.format(args.backend, args.pol, args.sideband.upper()),
        'box_001', '*.fits*'))[0])
for subband in args.subband:
    if 'flux_{0}'.format(subband) in hdulist[1].data.names:
        freq = hdulist[1].data.field('{0}frequency_{1}'.format(args.sideband.lower(), subband))[0]
        flux = hdulist[1].data.field('flux_{0}'.format(subband))[0]
        plt.plot(freq, flux, label='subband {0}'.format(subband), drawstyle='steps-mid')
if args.freq: plt.axvline(x=args.freq, linestyle='--')
plt.legend()
plt.show()
