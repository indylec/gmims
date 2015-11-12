import numpy as np
import healpy as hp
import theil_sen as ts
import matplotlib.pyplot as plt
from astropy.io import ascii

gmims0_freq=1310025600.0
haslam_freq=408000000.0

denominator=np.log(gmims0_freq/haslam_freq)

gmims0=hp.read_map('TP_000-360_0_healpix_res9_masked.fits',nest=True)
gmims_neg=np.where(gmims0<0)
gmims0[gmims_neg]=hp.UNSEEN
gmims_unseen=np.where(gmims0==hp.UNSEEN)
gmims0[gmims_unseen]=float('NaN')
gmims0=gmims0+0.3935betamapbetam

haslam=hp.read_map('final_haslam_masked.fits',nest=True)
haslam_unseen=np.where(haslam==hp.UNSEEN)
haslam[haslam_unseen]=float('NaN')

betamap=np.zeros(np.size(gmims0))

print "Calculating spectral index for each pixel" 

for i in range (0,np.size(gmims0)):
    betamap[i]=(np.log(gmims0[i]/haslam[i]))/denominator

hp.mollview(betamap,nest=True)

plt.show()





