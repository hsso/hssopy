#!/usr/bin/python

import pyfits, sys

hdus = pyfits.open(sys.argv[1], mode='update')

for i in hdus:
#     if 'CUNIT1' in i.header and 'Degrees' in i.header['CUNIT1']:
    i.header.update('CUNIT1', 'deg')
    i.header.update('CUNIT2', 'deg')
hdus.flush()
