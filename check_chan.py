#Python version of check_chan.pro - looks at proportion of good pixels in each channel of the map

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys


map_in=sys.argv[1]

map=fits.getdata(map_in)
header=fits.getheader(map_in)

nchan=map.shape[0]
ny=map.shape[1]
nx=map.shape[2]

npix=nx*ny

f0=header['CRVAL3']
df=header['CDELT3']
refpix=header['CRPIX3']

ngood=np.empty(nchan)

for i in range nchan:
    ngood[i]=np.count_nonzero(~np.isnan(map[i,:,:]))

fgood=ngood.astpye('f')/npix

fg95=fgood[np.where(fgood>0.95)]

fig=plt.figure(6,4)
ax=fig.subplot(111)

ax.plot(nchan,fgood,'-k')
ax.plot(nchan[np.where(fgood>0.95], fg95, '-r')

plt.show()

#ax2=ax.twiny()







