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
ra = hdus[0].header['RA_NOM']
dec = hdus[0].header['DEC_NOM']
# Center at the comet position
# gc.recenter(ra, dec, width=0.02,height=0.02)
gc.show_markers(ra, dec, marker="+", color="white")
gc.show_colorscale(pmin=0.6)
# gc.show_colorscale(stretch="log", vmid=-0.0001)
# gc.show_colorscale(stretch="arcsinh")
gc.show_contour(colors="white")
gc.save('{}_{}.png'.format(args.obsid, args.band))
