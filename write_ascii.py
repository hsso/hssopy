#!/usr/bin/python

import os
import pyfits
import numpy as np

# hsadir = '/mnt/herschel/data/scratch/deval/Saturn/Saturn-1'
# asciidir = '/mnt/herschel/data/scratch/deval/Saturn/ascii/'

def convert_ascii(obsid, hsadir, asciidir):
    for rec in ['HRS', 'WBS']:
        for pol in ['H', 'V']:
            for sb in ['USB', 'LSB']:
                obsdir = os.path.join(hsadir, str(obsid), 'level2',
                    '{0}-{1}-{2}'.format(rec, pol, sb), 'box_001')
                fitsfile = os.listdir(obsdir)[0]
                hdulist = pyfits.open(os.path.join(obsdir, fitsfile))
                print obsdir, fitsfile
                convert_one(hdulist, obsid, rec, pol, sb, asciidir)

def convert_one(hdulist, obsid, rec, pol, sb, asciidir, suffix=''):
    for j in range(1, len(hdulist)):
        for i in range(1,5):
            if 'flux_{0}'.format(i) in hdulist[j].columns.names:
                print '{0}frequency_{1}'.format(sb.lower(), i)
                freq = hdulist[j].data.field('{0}frequency_{1}'.format(sb.lower(), i))[0]
                flux = hdulist[j].data.field('flux_{0}'.format(i))[0]
                outfile = os.path.join(asciidir,
                    '{0}_{1}-{2}-{3}_{4:02}_{5}{6}.txt'.format(obsid, rec, pol,
                        sb, j, i, suffix))
                np.savetxt(outfile, np.transpose((freq, flux)))
                # Set proper permissions
                os.chmod(outfile, 0644)
