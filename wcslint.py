#!/usr/bin/python

import sys
from astropy.io import fits as pyfits

hdus = pyfits.open(sys.argv[1], mode='update')

for i in hdus:
#     if 'CUNIT1' in i.header and 'Degrees' in i.header['CUNIT1']:
    i.header['CUNIT1'] = 'deg'
    i.header['CUNIT2'] = 'deg'
hdus.flush()
hdus.close()
