import numpy as np
import healpy as hp
from astropy.io import ascii
import sys
import matplotlib.pyplot as plt



#Compute classic spectral index, not from T-T plots


low_in=sys.argv[1]
hi_in=sys.argv[2]
low_name=sys.argv[3]
hi_name=sys.argv[4]
lowfreq=float(sys.argv[5])
hifreq=float(sys.argv[6])




#g0_freq=1300855232.
#g8_freq=1689233472.
#haslam_freq=408000000.
#stockert_freq=1420000000.

hi=hp.read_map(hi_in,nest=True)
low=hp.read_map(low_in,nest=True)

hi[np.where(hi==hp.UNSEEN)]=np.nan
low[np.where(low==hp.UNSEEN)]=np.nan

#if hi_name=='g0g':
#    hifreq=g0_freq
#elif hi_name=='g8g':
#    hifreq=g8_freq
#elif hi_name=='stcorr':
#    hifreq=stockert_freq

#if low_name=='g0g':
#    lowfreq=g0_freq
#elif low_name=='g8g':
#    lowfreq=g8_freq
#elif low_name=='hcorr':
#    lowfreq=haslam_freq
#elif low_name=='stcorr':
#    lowfreq=stockert_freq

s=hp.nside2npix(512)

betamap=np.empty(s)

for i in range(s):
    betamap[i]=np.log10(hi[i]/low[i])/np.log10(hifreq/lowfreq)

hp.mollview(betamap,nest=True, min=-6., max= 1., title='Spectral index, {0}-{1}'.format(low_name,hi_name))

plt.show()

outname=low_name+'-'+hi_name+'_spec_index_pix.fits'
hp.write_map(outname, betamap,nest=True)
