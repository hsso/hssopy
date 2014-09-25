#!/usr/bin/python
"""
Quick view of PACS data
"""

from herschel import pacsfits
from astropy.io import fits as pyfits
import argparse
import aplpy

parser = argparse.ArgumentParser()
parser.add_argument('obsid', default="")
parser.add_argument('-b', '--band', default="blue", choices=("blue", "red"),
                help="PACS band")
parser.add_argument('--datadir', default='./')
args = parser.parse_args()

fitsfile = pacsfits(args.datadir, args.obsid, args.band[0])
hdus = pyfits.open(fitsfile)
gc = aplpy.FITSFigure(fitsfile)
gc.show_colorscale()
ra = hdus[0].header['RA_NOM']
dec = hdus[0].header['DEC_NOM']
gc.show_markers(ra, dec)
# Center at the comet position
gc.recenter(ra, dec, width=0.02,height=0.02)
gc.save('myfirstplot.png')
