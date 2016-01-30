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

ax.plot(chans,means,'-k', label="Mean value")
ax.plot(chans,medians, '-b',label="Median value")

ax.set_xlabel("Channel number")
ax.set_ylabel("Brightness temp. (K)")

if stokes_in=='TP':
   texty=0.1
elif stokes_in=='Q':
    texty=-0.02
elif stokes_in=='U':
    texty=-0.02
else:
    texty=0.05

ax.text(2300,texty,stokes_in)


ax2=ax.twiny()

freqs=((np.arange(nchan)-refpix)*df+f0)/1E9

ax2.plot(freqs,means,color="black",marker="None", ls="None")
ax2.plot(freqs,medians,color="black",marker="None", ls="None")
ax2.set_xlabel("Frequency (GHz)")

ax.legend(loc=0)

#plt.show()

plt.savefig("plots/means_medians_"+stokes_in+".pdf", dpi=200, bbox_inches="tight")








