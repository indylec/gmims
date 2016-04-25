#Code to apply a given healpix mask to an unmasked healpix map; both mask and map must have the same resolution. It is assumed that masked pixels have value 0 in the mask file, these are set to the hp.UNSEEN value in the map we want to mask. Also assumed that both maps are in the NEST ordering scheme

import numpy as np
import healpy as hp
import sys

inmap=sys.argv[1]
binno=sys.argv[2]
inmask=sys.argv[3]
#mask_name=sys.argv[3]

map=hp.read_map(inmap, nest=True)
mask=hp.read_map(inmask, nest=True)

map[np.where(mask==0.0)]=hp.UNSEEN
map[np.where(map==0.0)]=hp.UNSEEN

#outfile=inmap.rsplit('.',1)[0]+'_'+mask_name+'.fits'

outfile='g'+binno+'g_512.fits'

hp.write_map(outfile, map,nest=True )
