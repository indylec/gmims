import numpy as np
import healpy as hp
import theil_sen_uncertain as ts
import matplotlib.pyplot as plt

gmims0=hp.read_map('TP_000-360_0_healpix_res9_masked.fits',nest=True)
gmims_unseen=np.where(gmims0==hp.UNSEEN)
gmims0[gmims_unseen]=float('NaN')


haslam=hp.read_map('final_haslam_masked.fits',nest=True)
haslam_unseen=np.where(haslam==hp.UNSEEN)
haslam[haslam_unseen]=float('NaN')

#Order the map arrays into rows each containing the number of pixels in one nside=4 pixel
gmims0_div=np.reshape(gmims0,(192,16384))
haslam_div=np.reshape(haslam,(192,16384))


#get the 'good' fit pixels
#indices=np.genfromtxt('pixlist2.txt',delimiter='\n')

#regions=np.shape(indices)[0]

#initialise the arrays of regression coefficients
a_b_1=np.zeros((192,2))
a_b_2=np.zeros((192,2))




print 'Estimating regression coefficients...'
#evaluate a and b for each pair of haslam, gmims0 regions and each pair of gmims0,gmims8 regions using the theil-sen estimator
for i in range (0,192):
     a_b_1[i,:]=ts.theil_sen(haslam_div[i,:],gmims0_div[i,:])
     #a_b_2[i,:]=ts.theil_sen(gmims0_div[i,:],haslam_div[i,:])
print '...done!'

#a_b_1=a_b_1[~np.isnan(a_b_1).any(1)]
#a_b_2=a_b_2[~np.isnan(a_b_2).any(1)]


#print 'a_b_1:',a_b_1
#print 'a_b_2:',a_b_2

#separate the a (slope) and b(offset) values
a=a_b_1[:,0]
b=a_b_1[:,1]
#a_inv=a_b_2[:,0]
#b_inv=a_b_2[:,1]

print 'Using ',np.shape(a)[0],' regions'

print 'Calculating overall offsets (monopole) using overdetermined system...'
#set up the overdetermined system to get out the offsets on each map
a_matrix=np.matrix(np.vstack((-a,np.ones(np.shape(a)[0]))).T)
b_matrix=np.transpose(np.matrix(b))

a_inv_matrix=np.matrix(np.vstack((-a_inv,np.ones(np.shape(a_inv)))).T)
b_inv_matrix=np.transpose(np.matrix(b_inv))

#solve the normal equation
m_hat=np.linalg.inv(np.transpose(a_matrix)*a_matrix)*np.transpose(a_matrix)*b_matrix

m_hat_inv=np.linalg.inv(np.transpose(a_inv_matrix)*a_inv_matrix)*np.transpose(a_inv_matrix)*b_inv_matrix

print '...done!'

print 'Using GMIMS vs Haslam:'
print 'm1 = {!s}, m2={!s}'.format(m_hat[0],m_hat[1])

print 'Using Haslam vs GMIMS:'
print 'm1 = {!s}, m2={!s}'.format(m_hat_inv[0],m_hat_inv[1])
