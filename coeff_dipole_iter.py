import numpy as np
import healpy as hp
from scipy.stats import theilslopes
from astropy.io import ascii
import sys

#get input arguments: lo, hi map names and resolutions
lo_in=sys.argv[1]
hi_in=sys.argv[2]
lo_name=sys.argv[3]
hi_name=sys.argv[4]
nside_map=int(sys.argv[5])
nside_subregion=int(sys.argv[6])
lo_freq=float(sys.argv[7])
hi_freq=float(sys.argv[8])

#compute the relevant healpix parameters

npix_map=hp.nside2npix(nside_map)
npix_subregion=hp.nside2npix(nside_subregion)

subpixel_size=int(npix_map/npix_subregion)
subpixels=npix_subregion

#output names
out=lo_name+'_'+hi_name

#set up the template matrix, N_map x 4 matrix

## t_matrix_map=np.empty((npix_map,4))
## t_matrix_map[:,0]=1.0

## for i in range (npix_map):
##         t_matrix_map[i,1:]=hp.pix2vec(nside_map,i,nest=True)

t_matrix_subregion=np.empty((npix_subregion,4))
t_matrix_subregion[:,0]=1.0

for i in range (npix_subregion):
        t_matrix_subregion[i,1:]=hp.pix2vec(nside_subregion,i,nest=True)

t_matrix_map=np.empty((npix_map,4))
t_matrix_map[:,0]=1.0

for i in range (npix_map):
        t_matrix_map[i,1:]=hp.pix2vec(nside_map,i,nest=True)

#hard-code wehus dipole coeffs

w_coeffs=np.asarray([8.9,3.2,0.7,-0.8])

#read in the data

hi_data=hp.read_map(hi_in,nest=True)
hi_data=hp.ud_grade(hi_data,nside_out=256,pess=True,order_in='NESTED', order_out='NESTED')
hi_unseen=np.where(hi_data==hp.UNSEEN)
hi_data[hi_unseen]=float('NaN')
if hi_name=='stockert':
    hi_data*= 1.55
    hi_data-=2.8

lo_data=hp.read_map(lo_in,nest=True)
lo_data=hp.ud_grade(lo_data,nside_out=256,pess=True, order_in='NESTED',order_out='NESTED')
lo_unseen=np.where(lo_data==hp.UNSEEN)
lo_data[lo_unseen]=float('NaN')

if lo_name=='stockert':
    lo_data *= 1.55
    lo_data -= 2.8
elif lo_name=='haslam':
    lo_data -= 2.7

nsim=101   

iter=10


#start the bootstrap loop
md_tot=np.zeros((nsim,4))

for sim in range(nsim):

    print 'sim ', sim

    hi_iter=hi_data
    lo_iter=lo_data

    #Start the iteration loop

    for i in range(iter):

        print 'iter ', i

        #Set up the coefficient list

        temp_coeffs=np.empty((subpixels,4))
        md_temp=np.zeros(4)

        #Order the map arrays into rows each containing the number of pixels in one nside=16 pixel

        hi_div=np.reshape(hi_iter,(subpixels,subpixel_size))
        lo_div=np.reshape(lo_iter,(subpixels,subpixel_size))

        print 'Estimating regression coefficients...'

        for i in range (subpixels):
            hi_good=np.where(~np.isnan(hi_div[i,:]))
            lo_good=np.where(~np.isnan(lo_div[i,:]))
            good_intersect=np.intersect1d(np.asarray(hi_good),np.asarray(lo_good))

            if good_intersect.shape[0]>64:
                init_coeffs[i,:]=theilslopes(hi_div[i,good_intersect],lo_div[i,good_intersect])

            else:
                temp_coeffs[i,:]=np.nan


        #retain the good pixels only (do this for first iteration)

        slope_temp=temp_coeffs[:,0]
        b_temp=temp_coeffs[:,1]

        #spec index prior
        if i == 0:
            prior=np.empty(2)
            prior[0]=-4.0
            prior[1]=-2.0
            prior=(hi_freq/lo_freq)**prior

            goodslopes=np.where(np.logical_or(slope>=prior[0],slope<=prior[1]))

            pixels=np.arange(subpixels)
            pixels=pixels[goodslopes]

            good_regions=int(pixels.shape[0])

            if sim != 0:
                pixels=np.random.randint(0,100,good_regions)

            
        b=b_temp[pixels]
        slope=slope_temp[pixels]

        A=t_matrix_subregion[pixels,:]

        md_temp,resx,rank,sing=np.linalg.lstsq(A,b)

        md_tot[sim,:] += md_temp

        if np.all(md_temp/md_tot[sim,:] < 0.01):
            if sim == 0:
                hi_corrected_out=lo_name+'_'+hi_name+'_corrected_map.fits'
                hp.write_map(hi_corrected_out,hi_iter,nest=True)

                slope_betamap=np.ones(npix_region)*hp.UNSEEN
                slope_betamap[pixels]=np.log10(slope)/np.log10(hi_freq/lo_freq)
                hp.write_map(out+'_slope_betamap.fits',slope_betamap,nest=True)

                betamap=np.empty(npix_map)
                betamap=np.log10(hi_iter/lo_iter)/np.log10(hi_freq/lo_freq)
                hp.write_map(out+'_pix_betamap.fits',betamap,nest=True)
                


                continue

        else:
            hi_iter-=np.matmul(t_matrix_map,md_temp)

        

#Output
    
m10=md_tot[0,0]
m10_err=np.std(x_bootstrap_4[:,0])

m11=x[0,1]
m11_err=np.std(md_tot[:,1])

m12=x[0,2]
m12_err=np.std(md_tot[:,2])

m13=x[0,3]
m13_err=np.std(md_tot[:,3])

coeffs_out=out+'_iter_bootstrap_mdtot.txt'
       
f=open(coeffs_out,'a')
f.write('################ MONOPOLE + DIPOLE COEFFICIENTS###################\n')
#f.write('#First row is '+low_name+', second is '+hi_name+'\n')
f.write('#Coefficients for gmims'# {0} \n'.format(low_name))
f.write('x0 x1 x2 x3 x0err x1err x2err x3err \n')
f.write(str(m10)+' '+str(m11)+' '+str(m12)+' '+str(m13)+' '+str(m10_err)+' '+str(m11_err)+' '+str(m12_err)+' '+str(m13_err)+'\n')
#f.write(str(m20)+' '+str(m21)+' '+str(m22)+' '+str(m23)+' '+str(m20_err)+' '+str(m21_err)+' '+str(m22_err)+' '+str(m23_err)+'\n')
f.write('##################################################################')
f.close()

hi_fit=np.matmul(t_matrix_map,md_tot[0,:])
hp.write_map(out+'_fit.fits', fit_hi, nest=True)

   
    
