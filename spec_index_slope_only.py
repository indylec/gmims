import numpy as np
import healpy as hp
from astropy.io import ascii
import matplotlib.pyplot as plt
import sys

coeffs_in=sys.argv[1]
low_name=sys.argv[2]
hi_name=sys.argv[3]
low_freq=float(sys.argv[4])
hi_freq=float(sys.argv[5])
nside=int(sys.argv[6])

data=ascii.read(coeffs_in, guess=False, delimiter=' ')


npix=hp.nside2npix(nside)

pixels=data['pixelno']
#print pixels
slope=data['slope']
#b=data['offset']

good_regions=np.shape(data['pixelno'])
goodno=int(good_regions[0])

slope_betamap=np.ones(npix)*hp.UNSEEN

for iii in range(goodno):
    slope_betamap[[pixels[iii]]]=np.log10(slope[iii])/np.log10(hi_freq/low_freq)
    print slope[iii],np.log10(slope[iii])/np.log10(hi_freq/low_freq)

slope_betamap_out=low_name+'-'+hi_name+'_spec_index_coeff.fits'

hp.write_map(slope_betamap_out,slope_betamap,nest=True)

