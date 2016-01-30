import numpy as np
import healpy as hp
from astropy.io import ascii
import matplotlib.pyplot as plt
import sys


coeffs_in=sys.argv[1]
gmims_in=sys.argv[2]


data=ascii.read(coeffs_in, guess=False, delimiter=' ')


npix=hp.nside2npix(512)

pixels=data['pixelno']
slope=data['slope']
b=data['offset']

print b.shape

good_regions=np.shape(data['pixelno'])

ipix=np.arange(npix)
theta,phi=hp.pix2ang(512,ipix,nest=True)

t_matrix=np.zeros((npix,4))

for i in range (npix):

    t_matrix[i,:]=[1,np.cos(phi[i])*np.sin(theta[i]),np.sin(theta[i])*np.sin(phi[i]),np.cos(theta[i])]



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

print a_matrix.shape


x,res,rank,sing=np.linalg.lstsq(a_matrix,b)


print "The x-vector is: "
print x

gmims,gheader=hp.read_map(gmims_in,nest=True,h=True)

fit=np.empty(gmims.shape[0])

fit=(x[0]+np.cos(phi)*np.sin(theta)*x[1]+np.sin(theta)*np.sin(phi)*x[2]+np.cos(theta)*x[3])

gmims_corrected=gmims-fit

corrected_out=coeffs_in.rsplit('.',1)[0]+'_corrected.fits'
corrected_pdf=gmims_in.rsplit('.',1)[0]+'_corrected.pdf'

fit_out=coeffs_in.rsplit('.',1)[0]+'_fit.fits'
fit_pdf=gmims_in.rsplit('.',1)[0]+'_fit.pdf'

fig1=plt.figure(figsize=(8,5),dpi=150)
hp.mollview(fit,nest=True,fig=fig1.number, xsize=8*150,min=-1.,max=1.)
#plt.savefig(fit_pdf,dpi=150, bbox_inches="tight")
#plt.show()

fig2=plt.figure(figsize=(8,5),dpi=150)
hp.mollview(gmims_corrected,nest=True,fig=fig2.number, xsize=8*150,min=-1.,max=40.)
#plt.savefig(corrected_pdf,dpi=150, bbox_inches="tight")
#plt.show()

hp.write_map(corrected_out,gmims_corrected,nest=True)
hp.write_map(fit_out, fit, nest=True)

coeffs_out=coeffs_in.rsplit('.',1)[0]+'_dipole.txt'

f=open(coeffs_out,'a')
f.write(str(x))
f.close()

