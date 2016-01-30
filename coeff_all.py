import sys
import numpy as np
import healpy as hp
#import theilsen_bootstrap as ts
import theil_sen_uncertain as ts
import matplotlib.pyplot as plt
from astropy.io import ascii

#pixno=float(sys.argv[1])
lo_in=sys.argv[1]
hi_in=sys.argv[2]
lo_name=sys.argv[3]
hi_name=sys.argv[4]
#mask_name=sys.argv[5]

gmims0=hp.read_map(hi_in,nest=True)
gmims_unseen=np.where(gmims0==hp.UNSEEN)
gmims0[gmims_unseen]=float('NaN')
#if hi_name=='stockert':
 #   gmims0-=2.8

haslam=hp.read_map(lo_in,nest=True)
haslam_unseen=np.where(haslam==hp.UNSEEN)
haslam[haslam_unseen]=float('NaN')
#if lo_name=='stockert':
 #   haslam -= 2.8
#elif lo_name=='haslam':
 #   haslam -= 5.8

#Order the map arrays into rows each containing the number of pixels in one nside=4 pixel
gmims0_div=np.reshape(gmims0,(192,16384))
#gmims8_div=np.reshape(gmims8,(192,4096))
haslam_div=np.reshape(haslam,(192,16384))

#find pixels which are good for both maps

text_name=lo_name+'-'+hi_name+'_coeffs.txt'

print 'Estimating regression coefficients...'

f=open(text_name, 'a')
f.write('pixelno slope offset goodpix \n')
f.close()

for i in range (192):
 #evaluate a and b for each pair of haslam, gmims0 regions and each pair of gmims0,gmims8 regions using the theil-sen estimator
    #print 'Current pixel: ',i
    gmims_good=np.where(~np.isnan(gmims0_div[i,:]))
    haslam_good=np.where(~np.isnan(haslam_div[i,:]))

    good_intersect=np.intersect1d(np.asarray(gmims_good),np.asarray(haslam_good))
    #print 'Number of good pixels: ',good_intersect.shape[0]

    if good_intersect.shape[0]>0:
        a_b=ts.theil_sen(haslam_div[i,good_intersect],gmims0_div[i,good_intersect])
        a=a_b[0]
        b=a_b[1]
        #intercept_1=-b/a

        #print 'slope = ', a
        #print '95 range for slope = ',a_b[2],' to ',a_b[3]
        #print 'y-intercept = ', b
        #print '95 range for intercept = ',a_b[4],' to ',a_b[5]

        f=open(text_name, 'a')
        
        written=(str(i)+' '+str(a)+' '+str(b)+' '+str(good_intersect.shape[0])+'\n')
        f.write(written)
        f.close()
        #print 'Values written to file.'

    else:
        #print 'There are no good pixels in this subpixel!'
        f=open(text_name, 'a')
        written_blank=('# No good values for this subpixel.'+'\n')
        f.write(written_blank)
        f.close()
        
        
#print 'linear coeff = ', a
#print 'y-intercept = ', b
#print 'x-intercept = ', intercept_1

#print 'inverse linear coeff = ', a_inv
#print 'inverse y-intercept = ', b_inv
#print 'inverse x-intercept = ', intercept_2

#axislimits1=[np.nanmin(haslam_div[pixno,:])-0.1,np.nanmax(haslam_div[pixno,:])+0.1,np.nanmin(gmims0_div[pixno,:])-0.1,np.nanmax(gmims0_div[pixno,:])+0.1]

#axislimits2=[np.nanmin(gmims0_div[pixno,:])-0.1,np.nanmax(gmims0_div[pixno,:])+0.1,np.nanmin(haslam_div[pixno,:])-0.1,np.nanmax(haslam_div[pixno,:])+0.1]
#print axislimits


#x1=np.arange(0,90,0.5)
#x2=np.arange(-2,4,0.1)
#pix=np.arange(0,3145728,1.0)
#pix_div=np.reshape(pix,(192,16384))

#plt.figure(1)
#plt.title('TT plot for pixel number '+str(int(pixno)))
#plt.suptitle('TT plot for pixel number '+str(int(pixno)))


#plt.axis(axislimits1)
#plt.scatter(haslam_div[pixno,:],gmims0_div[pixno,:],c=pix_div[pixno,:],marker='+',s=200, cmap='hsv',linewidths=4)
#plt.plot(x1,a*x1+b,'k-')
#plt.colorbar()
#plt.xlabel('Haslam amplitude (K)')
#plt.ylabel('GMIMS amplitude (K)')

#plt.subplot(122)
#plt.axis(axislimits2)
#plt.scatter(gmims0_div[pixno,:],haslam_div[pixno,:],c=pix_div[pixno,:],marker='+',s=200,cmap='hsv',linewidths=4)
#plt.plot(x2,a_inv*x2+b_inv,'k--',x2,x2/a-b/a,'k-')
#plt.colorbar()
#plt.xlabel('GMIMS amplitude (K)')
#plt.ylabel('Haslam amplitude (K)')

#plt.legend(('pixels','linear regression fit','inverted linear regression fit'),loc=0)


#plt.show()

