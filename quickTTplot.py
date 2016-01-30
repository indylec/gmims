import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import ascii

pixelno=int(sys.argv[1])
coeffs_in=sys.argv[2]
hi_in=sys.argv[3]
lo_in=sys.argv[4]

coeffs=ascii.read(coeffs_in, guess=False, delimiter=' ')

a=coeffs['slope']
b=coeffs['offset']

print a[pixelno],b[pixelno]
#a_inv=coeffs['a_inv']
#b_inv=coeffs['b_inv']

if np.isnan(a[pixelno] or b[pixelno] or a_inv[pixelno] or b_inv[pixelno]):
    sys.exit("The chosen pixel is blank")
    


gmims0=hp.read_map(hi_in,nest=True)
gmims_unseen=np.where(gmims0==hp.UNSEEN)
gmims0[gmims_unseen]=float('NaN')


haslam=hp.read_map(lo_in,nest=True)
haslam_unseen=np.where(haslam==hp.UNSEEN)
haslam[haslam_unseen]=float('NaN')


#Order the map arrays into rows each containing the number of pixels in one nside=4 pixel
gmims0_div=np.reshape(gmims0,(192,16384))
#gmims8_div=np.reshape(gmims8,(192,4096))
haslam_div=np.reshape(haslam,(192,16384))


axislimits1=[np.nanmin(haslam_div[pixelno,:])-0.1,np.nanmax(haslam_div[pixelno,:])+0.1,np.nanmin(gmims0_div[pixelno,:])-0.1,np.nanmax(gmims0_div[pixelno,:])+0.1]

axislimits2=[np.nanmin(gmims0_div[pixelno,:])-0.1,np.nanmax(gmims0_div[pixelno,:])+0.1,np.nanmin(haslam_div[pixelno,:])-0.1,np.nanmax(haslam_div[pixelno,:])+0.1]
#print axislimits


x1=np.linspace(10,60,220)
x2=np.arange(-2,4,0.1)
pix=np.arange(0,3145728,1.0)
pix_div=np.reshape(pix,(192,16384))

line=x1*a[pixelno]+b[pixelno]

#print x1[100],line[100]


plt.title('TT plot for pixel number '+str(int(pixelno)))

plt.axis(axislimits1)
plt.scatter(haslam_div[pixelno,:],gmims0_div[pixelno,:],c=pix_div[pixelno,:],marker='+',s=200, cmap='hsv',linewidths=4)
plt.plot(x1,a[pixelno]*x1+b[pixelno],'k-')

cbar=plt.colorbar()
cbar.set_label('HEALPix nested pixel number')
plt.xlabel('Haslam amplitude (K)')
plt.ylabel('GMIMS amplitude (K)')

#plt.subplot(122)
#plt.axis(axislimits2)
#plt.scatter(gmims0_div[pixelno,:],haslam_div[pixelno,:],c=pix_div[pixelno,:],marker='+',s=200,cmap='hsv',linewidths=4)
#plt.plot(x2,a_inv[pixelno]*x2+b_inv[pixelno],'k-',x2,x2/a[pixelno]-b[pixelno]/a[pixelno],'k--')
#plt.colorbar()
#plt.xlabel('GMIMS amplitude (K)')
#plt.ylabel('Haslam amplitude (K)')

#plt.legend(('pixels','linear regression fit','inverted linear regression fit'),loc=0)


plt.show()
