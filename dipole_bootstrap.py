#implements the correct, full version of the dipole fiting in Wehus et al. 2014, outputs two sets of monopole/dipole coefficients, for map m1 and m2 where m2 is at the higher frequecy

import numpy as np
import healpy as hp
from astropy.io import ascii
import matplotlib.pyplot as plt
import sys


coeffs_in=sys.argv[1]
low_in=sys.argv[2]
hi_in=sys.argv[3]
low_name=sys.argv[4]
hi_name=sys.argv[5]




data=ascii.read(coeffs_in, guess=False, delimiter=' ')


npix=hp.nside2npix(512)

pixels=data['pixelno']
#print pixels
slope=data['slope']
b=data['offset']

print b.shape

good_regions=np.shape(data['pixelno'])
goodno=int(good_regions[0])

ipix=np.arange(npix)
theta,phi=hp.pix2ang(512,ipix,nest=True)

t_matrix=np.zeros((npix,4))

#set up the template matrix, Nx8 matrix

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

for ii in range(goodno):
    #print ii, pixels[ii]
    t0_i[ii]=np.mean(t0[int(pixels[ii]),:])
    t1_i[ii]=np.mean(t1[int(pixels[ii]),:])
    t2_i[ii]=np.mean(t2[int(pixels[ii]),:])
    t3_i[ii]=np.mean(t3[int(pixels[ii]),:])

#print t0_i
#print b

a_matrix_monopole=np.vstack([-slope*t0_i,t0_i]).T

a=a_matrix_monopole
aa_inv=np.linalg.inv(np.matmul(a.T,a))
atb=np.matmul(a.T,b.T)
x=np.matmul(aa_inv,atb)

#x,res,rank,sing=np.linalg.lstsq(a_matrix_monopole,b)

print x
#print res

res=np.matmul(a,x.T)-b

print res

exit()



#a_matrix_1=np.vstack([t0_i,t1_i,t2_i,t3_i])
#a_matrix_2=np.vstack([-slope*t0_i,-slope*t1_i,-slope*t2_i,-slope*t3_i])
a_matrix_3=np.vstack([-slope*t0_i,-slope*t1_i,-slope*t2_i,-slope*t3_i,t0_i,t1_i,t2_i,t3_i]).T

#print a_matrix_1.shape
#print a_matrix_2.shape
print a_matrix_3.shape

#x,res,rank,sing=np.linalg.lstsq(a_matrix_1,b)
#x,res,rank,sing=np.linalg.lstsq(a_matrix_2,b)
x,res,rank,sing=np.linalg.lstsq(a_matrix_3,b)


print "The x-vector is: "
print x

#Now start the bootstrap to get the uncertainties, 100 iterations
x_bootstrap_8=np.empty((100,8))
#x_bootstrap_4=np.empty((100,4))

for j in range(100):
    #print 'Loop {0}'.format(j)
    
    rand_indices=np.random.randint(0,100,int(good_regions[0]))

    a_matrix_bootstrap_8=np.empty((goodno,8))
    #a_matrix_bootstrap_4=np.empty((goodno,4))
    b_bootstrap=np.empty(goodno)
    
    for k in range(goodno):
        a_matrix_bootstrap_8[k,:]=a_matrix_3[rand_indices[k],:]
        #a_matrix_bootstrap_4[k,:]=a_matrix_2[rand_indices[k],:]
        #a_matrix_bootstrap_4[k,:]=a_matrix_1[rand_indices[k],:]
        b_bootstrap[k]=b[rand_indices[k]]

    x_bootstrap_8[j,:], rb, rab, sb = np.linalg.lstsq(a_matrix_bootstrap_8,b_bootstrap)
    #x_bootstrap_4[j,:], rb, rab, sb = np.linalg.lstsq(a_matrix_bootstrap_4,b_bootstrap)

m10=x[0]
#m10_err=np.std(x_bootstrap_4[:,0])
m10_err=np.std(x_bootstrap_8[:,0])
m11=x[1]
#m11_err=np.std(x_bootstrap_4[:,1])
m11_err=np.std(x_bootstrap_8[:,1])
m12=x[2]
#m12_err=np.std(x_bootstrap_4[:,2])
m12_err=np.std(x_bootstrap_8[:,2])
m13=x[3]
#m13_err=np.std(x_bootstrap_4[:,3])
m13_err=np.std(x_bootstrap_8[:,3])

m20=x[4]
m20_err=np.std(x_bootstrap_8[:,4])
m21=x[5]
m21_err=np.std(x_bootstrap_8[:,5])
m22=x[6]
m22_err=np.std(x_bootstrap_8[:,6])
m23=x[7]
m23_err=np.std(x_bootstrap_8[:,7])


#generate fits, corrected maps, and output everything

low,lowheader=hp.read_map(low_in,nest=True,h=True)
if low_name=='haslam':
    low -= 5.8
elif low_name=='stockert':
    low -= 2.8

hi,hiheader=hp.read_map(hi_in,nest=True,h=True)
if hi_name=='stockert':
    hi-=2.8

fit_low=(m10+np.cos(phi)*np.sin(theta)*m11+np.sin(theta)*np.sin(phi)*m12+np.cos(theta)*m13)

fit_hi=(m20+np.cos(phi)*np.sin(theta)*m21+np.sin(theta)*np.sin(phi)*m22+np.cos(theta)*m23)

low_corrected=low-fit_low
hi_corrected=hi-fit_hi

low_corrected_out=low_name+'_'+coeffs_in.rsplit('.',1)[0]+'_simple_corrected.fits'
hi_corrected_out=hi_name+'_'+coeffs_in.rsplit('.',1)[0]+'_corrected.fits'
#low_corrected_pdf=gmims_in.rsplit('.',1)[0]+'_corrected.pdf'

low_fit_out=low_name+'_'+coeffs_in.rsplit('.',1)[0]+'simple_fit.fits'
hi_fit_out=hi_name+'_'+coeffs_in.rsplit('.',1)[0]+'_fit.fits'
#fit_pdf=gmims_in.rsplit('.',1)[0]+'_fit.pdf'

#fig1=plt.figure(figsize=(8,5),dpi=150)
#hp.mollview(fit,nest=True,fig=fig1.number, xsize=8*150,min=-1.,max=1.)
#plt.savefig(fit_pdf,dpi=150, bbox_inches="tight")
#plt.show()

#fig2=plt.figure(figsize=(8,5),dpi=150)
#hp.mollview(gmims_corrected,nest=True,fig=fig2.number, xsize=8*150,min=-1.,max=40.)
#plt.savefig(corrected_pdf,dpi=150, bbox_inches="tight")
#plt.show()

hp.write_map(low_corrected_out,low_corrected,nest=True)
hp.write_map(hi_corrected_out,hi_corrected,nest=True)
hp.write_map(low_fit_out, fit_low, nest=True)
hp.write_map(hi_fit_out, fit_hi, nest=True)

coeffs_out=coeffs_in.rsplit('.',1)[0]+'_simple_dipoles.txt'

f=open(coeffs_out,'a')
f.write('################ MONOPOLE + DIPOLE COEFFICIENTS###################\n')
f.write('#First row is '+low_name+', second is '+hi_name+'\n')
#f.write('#Coefficients for {0} \n'.format(low_name))
f.write('x0 x1 x2 x3 x0err x1err x2err x3err \n')
f.write(str(m10)+' '+str(m11)+' '+str(m12)+' '+str(m13)+' '+str(m10_err)+' '+str(m11_err)+' '+str(m12_err)+' '+str(m13_err)+'\n')
f.write(str(m20)+' '+str(m21)+' '+str(m22)+' '+str(m23)+' '+str(m20_err)+' '+str(m21_err)+' '+str(m22_err)+' '+str(m23_err)+'\n')
f.write('##################################################################')
f.close()
