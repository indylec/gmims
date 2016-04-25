import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

nside_big=256
npix_big=hp.nside2npix(nside_big)

t_matrix=np.empty((npix_big,4))
t_matrix[:,0]=1.0

#set up the template matrix, Nx4 matrix

for i in range (npix_big):
    t_matrix[i,1:]=hp.pix2vec(nside_big,i,nest=True)

w_coeffs=np.asarray([8.9,3.2,0.7,-0.8])

wdipole=np.matmul(t_matrix,w_coeffs)

hp.mollview(wdipole, nest=True)

plt.show()
