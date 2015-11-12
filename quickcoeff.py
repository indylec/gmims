import sys
import numpy as np
import healpy as hp
import theil_sen_uncertain as ts
import matplotlib.pyplot as plt
from astropy.io import ascii

pixno=float(sys.argv[1])

gmims0=hp.read_map('TP_000-360_0_healpix_res9_masked.fits',nest=True)
gmims_unseen=np.where(gmims0==hp.UNSEEN)
gmims0[gmims_unseen]=float('NaN')


haslam=hp.read_map('final_haslam_masked.fits',nest=True)
haslam_unseen=np.where(haslam==hp.UNSEEN)
haslam[haslam_unseen]=float('NaN')

#Order the map arrays into rows each containing the number of pixels in one nside=4 pixel
gmims0_div=np.reshape(gmims0,(192,16384))
#gmims8_div=np.reshape(gmims8,(192,4096))
haslam_div=np.reshape(haslam,(192,16384))

#find pixels which are good for both maps
gmims_good=np.where(~np.isnan(gmims0_div[pixno,:]))
haslam_good=np.where(~np.isnan(haslam_div[pixno,:]))

good_intersect=np.intersect1d(np.asarray(gmims_good),np.asarray(haslam_good))


print 'Estimating regression coefficients...'
 #evaluate a and b for each pair of haslam, gmims0 regions and each pair of gmims0,gmims8 regions using the theil-sen estimator
a_b=ts.theil_sen(haslam_div[pixno,good_intersect],gmims0_div[pixno,good_intersect])
#a_b_invert=ts.theil_sen(gmims0_div[pixno,good_intersect],haslam_div[pixno,good_intersect])

a=a_b[0]
b=a_b[1]
intercept_1=-b/a
#a_inv=a_b_invert[0]
#b_inv=a_b_invert[1]
#intercept_2=-b_inv/a_inv
    #a_b_2[i,:]=ts.theil_sen(gmims0[i,:],gmims8[i,:])
print '...done!'


print 'slope = ', a
print '95 range for slope = ',a_b[2],' to ',a_b[3]
print 'y-intercept = ', b
print '95 range for intercept = ',a_b[4],' to ',a_b[5]

f=open('offset_and_slopes.txt', 'a')

written=(str(pixno)+' '+str(a)+' '+str(b)+' '+str(a_b[2])+' '+str(a_b[3])+' '+str(a_b[4])+' '+str(a_b[5])+'\n')

f.write(written)

f.close()

#print 'linear coeff = ', a
#print 'y-intercept = ', b
#print 'x-intercept = ', intercept_1

#print 'inverse linear coeff = ', a_inv
#print 'inverse y-intercept = ', b_inv
#print 'inverse x-intercept = ', intercept_2

axislimits1=[np.nanmin(haslam_div[pixno,:])-0.1,np.nanmax(haslam_div[pixno,:])+0.1,np.nanmin(gmims0_div[pixno,:])-0.1,np.nanmax(gmims0_div[pixno,:])+0.1]

axislimits2=[np.nanmin(gmims0_div[pixno,:])-0.1,np.nanmax(gmims0_div[pixno,:])+0.1,np.nanmin(haslam_div[pixno,:])-0.1,np.nanmax(haslam_div[pixno,:])+0.1]
#print axislimits


x1=np.arange(0,90,0.5)
x2=np.arange(-2,4,0.1)
pix=np.arange(0,3145728,1.0)
pix_div=np.reshape(pix,(192,16384))

plt.figure(1)
plt.title('TT plot for pixel number '+str(int(pixno)))
plt.suptitle('TT plot for pixel number '+str(int(pixno)))


plt.axis(axislimits1)
plt.scatter(haslam_div[pixno,:],gmims0_div[pixno,:],c=pix_div[pixno,:],marker='+',s=200, cmap='hsv',linewidths=4)
plt.plot(x1,a*x1+b,'k-')
plt.colorbar()
plt.xlabel('Haslam amplitude (K)')
plt.ylabel('GMIMS amplitude (K)')

#plt.subplot(122)
#plt.axis(axislimits2)
#plt.scatter(gmims0_div[pixno,:],haslam_div[pixno,:],c=pix_div[pixno,:],marker='+',s=200,cmap='hsv',linewidths=4)
#plt.plot(x2,a_inv*x2+b_inv,'k--',x2,x2/a-b/a,'k-')
#plt.colorbar()
#plt.xlabel('GMIMS amplitude (K)')
#plt.ylabel('Haslam amplitude (K)')

#plt.legend(('pixels','linear regression fit','inverted linear regression fit'),loc=0)


plt.show()

