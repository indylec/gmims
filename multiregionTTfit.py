import numpy as np
import healpy as hp
import theil_sen as ts
import matplotlib.pyplot as plt

gmims0_freq=1310025600.0
gmims8_freq=1718926336.0
haslam_freq=408000000.0

gmims0=hp.read_map('TP_000-360_0_healpix_res8_smoothed.fits',nest=True)
#gmims8=hp.read_map('TP_000-360_8_healpix_res8_smoothed.fits')
haslam=hp.read_map('haslam_nres8.dg.fits',nest=True)

#Order the map arrays into rows each containing the number of pixels in one nside=4 pixel
gmims0_div=np.reshape(gmims0,(192,4096))
#gmims8_div=np.reshape(gmims8,(192,4096))
haslam_div=np.reshape(haslam,(192,4096))

#Trim the blank southern part of the gmims maps, do the same for Haslam
gmims0_div=gmims0_div[0:151,:]
#gmims8_div=gmims8_div[0:151,:]
haslam_div=haslam_div[0:151,:]

#initialise the arrays of regression coefficients
a_b_1=np.zeros((151,2))
a_b_2=np.zeros((151,2))

print 'Estimating regression coefficients...'
#evaluate a and b for each pair of haslam, gmims0 regions and each pair of gmims0,gmims8 regions using the theil-sen estimator
for i in range(151):
    a_b_1[i,:]=ts.theil_sen(haslam_div[i,:],gmims0_div[i,:])
    #a_b_2[i,:]=ts.theil_sen(gmims0[i,:],gmims8[i,:])
print '...done!'

#separate the a (slope) and b(offset) values
a=a_b_1[:,0]
b=a_b_1[:,1]

#set up the overdetermined system to get out the offsets on each map
a_matrix=np.matrix(np.vstack((-a,np.ones(151))).T)
b_matrix=np.transpose(np.matrix(b)) 

#solve the normal equation
m_hat=np.linalg.inv(np.transpose(a_matrix)*a_matrix)*np.transpose(a_matrix)*b_matrix

#calculate the spectral index for each region
beta=np.log10(a)/np.log10(gmims0_freq/haslam_freq)

#pad the spectral index array with healpix 'missing' values
missing=np.ones(192-np.shape(beta)[0])*hp.UNSEEN
beta_map=np.concatenate((beta,missing))

print 'Writing spectral index map...'
#write the spectral indices to a healpix map of nside=4
hp.write_map('haslam.gmims0_beta.fits',beta_map,nest=True,coord='C')
print '...done!'

#print the offsets:

print 'm1 = {!s}, m2={!s}'.format(m_hat[0],m_hat[1])



