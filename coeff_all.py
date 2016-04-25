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

nside_big=int(sys.argv[5])
nside_small=int(sys.argv[6])

npix_big=hp.nside2npix(nside_big)

subpixel_size=int(hp.nside2npix(nside_big)/hp.nside2npix(nside_small))
print subpixel_size
subpixels=int(hp.nside2npix(nside_small))
print subpixels
#mask_name=sys.argv[5]

def a2t(freq):
    h=6.626068E-34
    k_B=1.3806503E-23
    T_0=2.726

    x=h*freq/(k_B*T_0)
    
    return (np.exp(x)-1.)**2/(x**2*np.exp(x))

gmims0=hp.read_map(hi_in,nest=True)
gmims0=hp.ud_grade(gmims0,nside_out=256,pess=True,order_in='NESTED', order_out='NESTED')
gmims_unseen=np.where(gmims0==hp.UNSEEN)
gmims0[gmims_unseen]=float('NaN')
if hi_name=='stockert':
    gmims0*= 1.55
    gmims0-=2.8
#gmims0/=a2t(1420000000.)

haslam=hp.read_map(lo_in,nest=True)
haslam=hp.ud_grade(haslam,nside_out=256,pess=True, order_in='NESTED',order_out='NESTED')
haslam_unseen=np.where(haslam==hp.UNSEEN)
haslam[haslam_unseen]=float('NaN')
#haslam/=a2t(408000000.)
if lo_name=='stockert':
    haslam *= 1.55
    haslam -= 2.8
elif lo_name=='haslam':
    haslam -= 2.7
elif lo_name=='haslam-wdipole':
    
    t_matrix=np.empty((npix_big,4))
    t_matrix[:,0]=1.0

#set up the template matrix, Nx4 matrix

    for i in range (npix_big):
        t_matrix[i,1:]=hp.pix2vec(nside_big,i,nest=True)

    w_coeffs=np.asarray([8.9,3.2,0.7,-0.8])

    wdipole=np.matmul(t_matrix,w_coeffs)

    hp.mollview(wdipole, nest=True)

    plt.show()

    #exit()
    haslam -= wdipole
    

#Order the map arrays into rows each containing the number of pixels in one nside=16 pixel
gmims0_div=np.reshape(gmims0,(subpixels,subpixel_size))
#gmims8_div=np.reshape(gmims8,(192,4096))
haslam_div=np.reshape(haslam,(subpixels,subpixel_size))

#find pixels which are good for both maps

text_name=lo_name+'-'+hi_name+'_coeffs.txt'

print 'Estimating regression coefficients...'

f=open(text_name, 'a')
f.write('pixelno slope offset slope_up slope_low goodpix \n')
f.close()

for i in range (hp.nside2npix(nside_small)):
 #evaluate a and b for each pair of haslam, gmims0 regions and each pair of gmims0,gmims8 regions using the theil-sen estimator
    #print 'Current pixel: ',i
    gmims_good=np.where(~np.isnan(gmims0_div[i,:]))
    haslam_good=np.where(~np.isnan(haslam_div[i,:]))

    good_intersect=np.intersect1d(np.asarray(gmims_good),np.asarray(haslam_good))
    #print 'Number of good pixels: ',good_intersect.shape[0]

    if good_intersect.shape[0]>64:
        a_b=ts.theil_sen(haslam_div[i,good_intersect],gmims0_div[i,good_intersect],sample=False)
        a=a_b[0]
        b=a_b[1]
        slope_up=a_b[2]
        slope_low=a_b[3]
        #intercept_1=-b/a

        #print 'slope = ', a
        #print '95 range for slope = ',a_b[2],' to ',a_b[3]
        #print 'y-intercept = ', b
        #print '95 range for intercept = ',a_b[4],' to ',a_b[5]

        f=open(text_name, 'a')
        
        written=(str(i)+' '+str(a)+' '+str(b)+' '+str(slope_up)+' '+str(slope_low)+' '+str(good_intersect.shape[0])+'\n')
        f.write(written)
        f.close()
        #print 'Values written to file.'

    else:
        #print 'There are less than 10 good pixels in this subpixel!'
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

