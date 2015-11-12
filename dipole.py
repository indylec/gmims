import numpy as np
import healpy as hp
from astropy.io import *
from math import *

data=ascii.read('offsets_slopes_100_nospur.txt', names=['pixno','slope','intercept','c95sl','c95sh','c95il','c95ih'])


npix=hp.nside2npix(512)

pixels=data['pixno']
slope=data['slope']
b=data['intercept']

good_regions=np.shape(data['pixno'])

ipix=np.arange(npix)
theta,phi=hp.pix2ang(512,ipix,nest=True)

t_matrix=np.zeros((npix,4))

for i in range (npix):

    t_matrix[i,:]=[1,cos(phi[i])*sin(theta[i]),sin(theta[i])*sin(phi[i]),cos(theta[i])]



t0=t_matrix[:,0]
t1=t_matrix[:,1]
t2=t_matrix[:,2]
t3=t_matrix[:,3]

t0=np.reshape(t0,(192,16384))
t1=np.reshape(t1,(192,16384))
t2=np.reshape(t2,(192,16384))
t3=np.reshape(t3,(192,16384))



t0_i=np.zeros(good_regions)
t1_i=np.zeros(good_regions)
t2_i=np.zeros(good_regions)
t3_i=np.zeros(good_regions)

for ii in range(good_regions[0]):

   t0_i[ii]=np.mean(t0[pixels[ii],:])
   t1_i[ii]=np.mean(t1[pixels[ii],:])
   t2_i[ii]=np.mean(t2[pixels[ii],:])
   t3_i[ii]=np.mean(t3[pixels[ii],:])


a_matrix=np.vstack([t0_i,t1_i,t2_i,t3_i]).T


intercept,res,rank,sing=np.linalg.lstsq(a_matrix,b)


print "The x-vector is: ",intercept

gmims,gheader=hp.read_map('gmims_masked_final_spur.fits',nest=True,h=True)

for i in range (npix):
    gmims[i]=gmims[i]-(intercept[0]+cos(phi[i])*sin(theta[i])*intercept[1]+sin(theta[i])*sin(phi[i])*intercept[2]+cos(theta[i])*intercept[3])

hp.mollview(gmims,nest=True)

hp.write_map('gmims_dipole_corrected_spur.fits',gmims,nest=True)
