import numpy as np
import healpy as hp
import theil_sen as ts
import matplotlib.pyplot as plt
from astropy.io import ascii

gmims0_freq=1310025600.0
gmims8_freq=1718926336.0
haslam_freq=408000000.0

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


#initialise the arrays of regression coefficients
a_b=np.zeros((192,2))
a_b_invert=np.zeros((192,2))


print 'Estimating regression coefficients...'
#evaluate a and b for each pair of haslam, gmims0 regions and each pair of gmims0,gmims8 regions using the theil-sen estimator
for i in range(192):
    a_b[i,:]=ts.theil_sen(haslam_div[i,:],gmims0_div[i,:])
    a_b_invert[i,:]=ts.theil_sen(gmims0_div[i,:],haslam_div[i,:])
    #a_b_2[i,:]=ts.theil_sen(gmims0[i,:],gmims8[i,:])
print '...done!'

print a_b
print a_b_invert

#separate the a (slope) and b(offset) values
a=a_b[:,0]
b=a_b[:,1]
intercept_1=-b/a
a_inv=a_b_invert[:,0]
b_inv=a_b_invert[:,1]
intercept_2=-b_inv/a_inv

#ascii.write([a,b,a_inv,b_inv],'coeffs.txt',names=['a','b','a_inv','b_inv'])

#plot the slope and offset for each subpixel
plt.figure(1)
plt.subplot(211)
plt.plot(b,'k-+',lw=3.0)
plt.plot(intercept_2,'r--^')
plt.subplot(212)
plt.plot(b_inv,'k-+',lw=3.0)
plt.plot(intercept_1,'r--^')

plt.show()
#set up the overdetermined system to get out the offsets on each map
#a_matrix=np.matrix(np.vstack((-a,np.ones(15))).T)
#b_matrix=np.transpose(np.matrix(b)) 

#solve the normal equation
#m_hat=np.linalg.inv((a_matrix.T)*a_matrix)*(a_matrix.T)*b_matrix

#calculate the spectral index for each region
#beta=np.log10(a)/np.log10(gmims0_freq/haslam_freq)

#pad the spectral index array with healpix 'missing' values
#missing=np.ones(192-np.shape(beta)[0])*hp.UNSEEN
#beta_map=np.concatenate((beta,missing))

#print 'Writing spectral index map...'
#write the spectral indices to a healpix map of nside=4
#hp.write_map('spec_index_1.fits',beta_map,nest=True,coord='C')
#print '...done!'

#print the offsets:

#print 'm1 = {!s}, m2={!s}'.format(m_hat[0],m_hat[1])



