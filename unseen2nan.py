import numpy as np
import healpy as hp
import sys

inmap=sys.argv[1]
inmask=sys.argv[2]

map=hp.read_map(inmap,nest=True)
mask=hp.read_map(inmask,nest=True)



mask[np.where(map==0.0)]=0.0

outmap='gal_kq85.fits'

hp.write_map(outmap,mask, nest=True)
