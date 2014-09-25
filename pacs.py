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
gc.show_markers(hdus[0].header['RA_NOM'], hdus[0].header['DEC_NOM'])
gc.save('myfirstplot.png')
