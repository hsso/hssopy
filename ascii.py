#!/usr/bin/python

import os
import pyfits
import numpy as np

hsadir = '/home/miguel/herschel/Jupiter'
asciidir = '/home/miguel/herschel/ascii'

# for obsid in [1342212195, 1342212184, 1342212192]:
for obsid in [1342212200]:
    for rec in ['HRS']:
        for pol in ['H', 'V']:
            for sb in ['USB', 'LSB']:
                obsdir = os.path.join(hsadir, str(obsid), 'level2',
                    '{0}-{1}-{2}'.format(rec, pol, sb), 'box_001')
                fitsfile = os.listdir(obsdir)[0]
                hdulist = pyfits.open(os.path.join(obsdir, fitsfile))
                print obsdir, fitsfile
                for i in range(1,5):
                    if 'flux_{0}'.format(i) in hdulist[1].columns.names:
                        freq = hdulist[1].data.field('{0}frequency_{1}'.format(sb.lower(), i))[0]
                        flux = hdulist[1].data.field('flux_{0}'.format(i))[0]
                        np.savetxt(os.path.join(asciidir,'{0}_{1}-{2}-{3}_{4}.txt'.format(obsid,
                            rec, pol, sb, i)), np.transpose((freq, flux)))
