import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import ascii
import theil_sen_uncertain as ts
from scipy import stats

pixelno=int(sys.argv[3])
#coeffs_in=sys.argv[2]
hi_in=sys.argv[1]
lo_in=sys.argv[2]

nside_big=int(sys.argv[4])
nside_small=int(sys.argv[5])

correct=int(sys.argv[6])

subpixel_size=int(hp.nside2npix(nside_big)/hp.nside2npix(nside_small))
print subpixel_size
subpixels=int(hp.nside2npix(nside_small))
print subpixels

#coeffs=ascii.read(coeffs_in, guess=False, delimiter=' ')

#a=coeffs['slope']
#b=coeffs['offset']

#print a[pixelno],b[pixelno]
#a_inv=coeffs['a_inv']
#b_inv=coeffs['b_inv']

#if np.isnan(a[pixelno] or b[pixelno] or a_inv[pixelno] or b_inv[pixelno]):
#    sys.exit("The chosen pixel is blank")
    

gmims0=hp.read_map(hi_in,nest=True)
if correct:
    gmims0-=2.8
gmims_unseen=np.where(gmims0==hp.UNSEEN)
gmims0[gmims_unseen]=float('NaN')



haslam=hp.read_map(lo_in,nest=True)
if correct:
    haslam-=2.7
haslam_unseen=np.where(haslam==hp.UNSEEN)
haslam[haslam_unseen]=float('NaN')

gmimsview=np.ones(hp.nside2npix(nside_big))*hp.UNSEEN
gmimsview[subpixel_size*pixelno:subpixel_size*(pixelno+1)]=gmims0[subpixel_size*pixelno:subpixel_size*(pixelno+1)]
#hp.mollview(gmims0,nest=True,min=0,max=6)
#plt.show()
#hp.mollview(gmimsview,nest=True)
#plt.show()

haslamview=np.ones(hp.nside2npix(nside_big))*hp.UNSEEN
haslamview[subpixel_size*pixelno:subpixel_size*(pixelno+1)]=haslam[subpixel_size*pixelno:subpixel_size*(pixelno+1)]
hp.mollview(haslam,nest=True,min=10.,max=100.)
#plt.show()
hp.mollview(haslamview,nest=True)
hp.mollview(gmimsview,nest=True)
plt.show()

#Order the map arrays into rows each containing the number of pixels in one nside=4 pixel
gmims0_div=np.reshape(gmims0,(subpixels,subpixel_size))
#gmims8_div=np.reshape(gmims8,(192,4096))
haslam_div=np.reshape(haslam,(subpixels,subpixel_size))

gmims_good=np.where(~np.isnan(gmims0_div[pixelno,:]))
haslam_good=np.where(~np.isnan(haslam_div[pixelno,:]))

good_intersect=np.intersect1d(np.asarray(gmims_good),np.asarray(haslam_good))

print 'good pixels in this subpixel: {0}'.format(good_intersect.shape[0])

if good_intersect.shape[0]>60:
        a_b=ts.theil_sen(haslam_div[pixelno,good_intersect],gmims0_div[pixelno,good_intersect],sample=False)
        a=a_b[0]
        b=a_b[1]
        slope_up=a_b[2]
        slope_low=a_b[3]

else:
    sys.exit('Not enough good pixels in the selected region!')

a_ols,b_ols,r,p,std = stats.linregress(haslam_div[pixelno,good_intersect],gmims0_div[pixelno,good_intersect])

print 'a = {0}, b = {1}, slope_low = {2}, slope_up = {3}'.format(a,b,slope_low,slope_up)

print 'Normal linear regression: a = {0}, b = {1}'.format(a_ols,b_ols)

print -b/a

slope_si=np.log10(a)/np.log10(1420000000./408000000.)

pixwise_si=np.log10(gmims0_div[pixelno,:]/haslam_div[pixelno,:])/np.log10(1420000000./408000000.)

pixwise_mean=np.nanmean(pixwise_si)

a_from_pixwise=(1420000000./408000000.)**pixwise_mean

beta_up=(1./np.log(10))*(slope_up-a)/a

beta_low=(1./np.log(10))*(a-slope_low)/a

print "Slope spectral index: {0} + {1} - {2}".format(slope_si,beta_up,beta_low)
print "Average spectral index in subpixel: {0}".format(pixwise_mean)


axislimits1=[np.nanmin(haslam_div[pixelno,:])-0.1,np.nanmax(haslam_div[pixelno,:])+0.1,np.nanmin(gmims0_div[pixelno,:])-0.1,np.nanmax(gmims0_div[pixelno,:])+0.1]

axislimits2=[np.nanmin(gmims0_div[pixelno,:])-0.1,np.nanmax(gmims0_div[pixelno,:])+0.1,np.nanmin(haslam_div[pixelno,:])-0.1,np.nanmax(haslam_div[pixelno,:])+0.1]
#print axislimits

axislimits3=[0,100,0,6]

x1=np.linspace(0,60,240)
x2=np.arange(2,6,0.01)
pix=np.arange(0,hp.nside2npix(nside_big),1.0)
pix_div=np.reshape(pix,(subpixels,subpixel_size))

#line=x1*a+b

#print x1[100],line[100]


plt.title('TT plot for pixel number '+str(int(pixelno)))

plt.axis(axislimits3)
plt.scatter(haslam_div[pixelno,:],gmims0_div[pixelno,:],c=pixwise_si,marker='.',s=50, cmap='hsv',linewidths=0)#,c=pix_div[pixelno,:])
plt.plot(x1,a*x1+b,'k-')
plt.plot(x1,a_ols*x1+b_ols,'k--')
#plt.plot(x1,0.33*x1+2.7,'b-')
plt.plot(x1,a_from_pixwise*x1,'k-.')

cbar=plt.colorbar()
cbar.set_label('Pixelwise spectral index')
plt.xlabel('Haslam amplitude (K)')
plt.ylabel('Stockert amplitude (K)')

#plt.subplot(122)
#plt.axis(axislimits2)
#plt.scatter(gmims0_div[pixelno,:],haslam_div[pixelno,:],c=pix_div[pixelno,:],marker='+',s=200,cmap='hsv',linewidths=4)
#plt.plot(x2,a_inv[pixelno]*x2+b_inv[pixelno],'k-',x2,x2/a[pixelno]-b[pixelno]/a[pixelno],'k--')
#plt.colorbar()
#plt.xlabel('GMIMS amplitude (K)')
#plt.ylabel('Haslam amplitude (K)')

#plt.legend(('pixels','linear regression fit','inverted linear regression fit'),loc=0)


plt.show()
