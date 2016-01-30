import numpy as np
import healpy as hp
from astropy.io import ascii
import sys
import matplotlib.pyplot as plt



#Read in the slopes from the text file, convert to spectral index and turn into a healpix map.

coeffs_in=sys.argv[1]
hi_in=sys.argv[2]
lo_in=sys.argv[3]

data=ascii.read(coeffs_in, guess=False, delimiter=' ')

g0_freq=1300855232.
g8_freq=1689233472.
haslam_freq=408000000.

npix=hp.nside2npix(4)

pixels=data['pixelno']
slope=data['slope']

si_map=np.zeros(npix)

if hi_in=='g0':
    nu2=g0_freq
elif hi_in=='g8':
    nu2=g8_freq

if lo_in=='g0':
    nu1=g0_freq
elif lo_in=='haslam':
    nu1=haslam_freq
    

spec_index=np.log10(slope)/np.log10(nu2/nu1)
print spec_index.shape

print pixels.shape

for i in range(pixels.shape[0]):
    
    print pixels[i],slope[i],spec_index[i]
    si_map[pixels[i]]=spec_index[i]

#hp.mollview(si_map,nest=True)
#plt.show()
