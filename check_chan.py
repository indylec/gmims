#Plots mean and median of each channel

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys


map_in=sys.argv[1]
stokes_in=sys.argv[2]

map=fits.getdata(map_in)
header=fits.getheader(map_in)

nchan=map.shape[0]
ny=map.shape[1]
nx=map.shape[2]

npix=nx*ny

f0=header['CRVAL3']
df=header['CDELT3']
refpix=header['CRPIX3']

means=np.empty(nchan)
medians=np.empty(nchan)

chans=np.arange(nchan)

for i in range(nchan):
    means[i]=np.mean(map[i,:,:])
    medians[i]=np.median(map[i,:,:])

fig=plt.figure(figsize=(10,6))
ax=fig.add_subplot(111)

ax.plot(chans,means,'k-', label="Mean value")
ax.plot(chans,medians, 'k--',label="Median value")

#add bin limits here, hard-coded

bins_list=np.asarray([205,410,615,820,1025,1230,1435,1640,1845,2047])

for i in range(bins_list.shape[0]):
    ax.axvline(x=bins_list[i], color='k',linestyle='dotted',alpha=0.8)
    ax.text(bins_list[i]-105,0.02,str(i))
    

ax.set_xlim([np.amin(chans),np.amax(chans)])

ax.set_xlabel("Channel number")
ax.set_ylabel("Brightness temp. (K)")

#if stokes_in=='TP':
#   texty=0.1
#elif stokes_in=='Q':
#    texty=-0.02
#elif stokes_in=='U':
#    texty=-0.02
#else:
#    texty=0.05




ax2=ax.twiny()

freqs=((np.arange(nchan)-refpix)*df+f0)/1E9

ax2.plot(freqs,means,color="blue",marker="None", ls="None")
ax2.plot(freqs,medians,color="blue",marker="None", ls="None")
ax2.set_xlabel("Frequency (GHz)")
ax2.set_xlim([np.amin(freqs),np.amax(freqs)])

ax.legend(loc=0)

plt.show()

plt.savefig("plots/means_medians_"+stokes_in+"_bins_corrected.pdf", dpi=200, bbox_inches="tight")








