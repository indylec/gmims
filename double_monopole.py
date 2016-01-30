import numpy as np
import healpy as hp
from astropy.io import ascii
import matplotlib.pyplot as plt
import sys

coeffs_in=sys.argv[1]
m1_name=sys.argv[2]
m2_name=sys.argv[3]

data=ascii.read(coeffs_in, guess=False, delimiter=' ')

npix=hp.nside2npix(512)

pixels=data['pixelno']
slope=data['slope']
b=data['offset']

slope=slope[np.where(~np.isnan(slope))]
b=b[np.where(~np.isnan(b))]

good_regions=np.shape(b)
goodno=int(good_regions[0])

print goodno

a_matrix_monopole=np.vstack([-slope,np.ones(goodno)]).T

a=a_matrix_monopole
aa_inv=np.linalg.inv(np.matmul(a.T,a))
atb=np.matmul(a.T,b.T)
x_norm=np.matmul(aa_inv,atb)

print x_norm

#x_lin,res,rank,sing=np.linalg.lstsq(a_matrix_monopole,b)

#print x_lin

res_norm=np.matmul(a,x_norm.T)-b

#print res_norm

monopole_out=coeffs_in.rsplit('.',1)[0]+'_double_monopole.txt'

f=open(monopole_out,'a')
f.write('################ MONOPOLES ###################\n')
#f.write('#First row is '+low_name+', second is '+hi_name+'\n')
f.write('#m1 is {0}, m2 is {1} \n'.format(m1_name,m2_name))
f.write('m1 m2\n')
f.write(str(x_norm[0])+' '+str(x_norm[1])+'\n')
#f.write(str(x_lin[0])+' '+str(x_lin[1])+'\n')#f.write(str(m20)+' '+str(m21)+' '+str(m22)+' '+str(m23)+' '+str(m20_err)+' '+str(m21_err)+' '+str(m22_err)+' '+str(m23_err)+'\n')
f.write('##################################################################')
f.close()

#exit()
