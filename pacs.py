#!/usr/bin/python
"""
Quick view of PACS data
"""

from herschel import pacsfits
import argparse
import aplpy

parser = argparse.ArgumentParser()
parser.add_argument('obsid', default="")
parser.add_argument('-b', '--band', default="blue", choices=("blue", "red"),
                help="PACS band")
parser.add_argument('--datadir', default='./')
args = parser.parse_args()

fitsfile = pacsfits(args.datadir, obsid, args.band[0])
gc = aplpy.FITSFigure(fitsfile)
gc.show_colorscale()
