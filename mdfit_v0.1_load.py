import sys
import numpy as np
import healpy as hp

nu_lo=408000000.
nu_hi=1420000000.

nside_map=256
npix_map=hp.nside2npix(nside_map)
nside_regions=4
npix=hp.nside2npix(nside_regions)
npix_sub=npix_map/npix

params=np.load("md_params.npy")
A_pos=np.load("md_A_pos.npy")
b_pos=np.load("md_b_pos.npy")
npos=np.load("md_npos.npy")


#make matrix template
T=np.empty((npix,4))
T[:,0]=1.0
for i in range(npix):
    T[i,1:]=hp.pix2vec(nside_regions,i,nest=True)

#find valid pixels

prior=np.empty(2)
prior[0]=-4.0
prior[1]=-2.0
prior=(nu_hi/nu_lo)**prior

print prior


pix=np.arange(npix)
for i in range (npix):
    if params[1,i]<prior[0] or params[1,i] >prior[1]:
        params[1:,i]=hp.UNSEEN
pix=pix[np.where(params[1,:]!=hp.UNSEEN)]
numval=pix.size

#print numval

slope=params[1,:]
slope=slope[pix]
T=T[pix,:]

offset=params[2,:]
b=offset[pix]

A=np.empty((numval,8))

for i in range(4):
    A[:,i]=-slope*T[:,i]
    
for i in range(4):
    A[:,i+4]=T[:,i]

print A.shape

x,resx,rank,sing=np.linalg.lstsq(A,b)

print x

res=np.matmul(A,x.T)-b

print res.size

mu_res=np.mean(res)
print mu_res
sigma_res=np.sqrt(np.var(res))
#print res.size

#4-sigma cut:

for i in range(numval):
    if np.abs(res[i]-mu_res)/sigma_res > 4.:
        params[1:,pix[i]]=hp.UNSEEN
        A[i,:]=0.
        b[i]=0.

np.save("md_A",A)
np.save("md_b",b)
np.save("md_x",x)
np.save("md_pix",pix)

md_tot=np.zeros(8)

def dist(md):

    for i in range (npos):
        if b_pos[i]-np.sum(A_pos[i,:]*(md_tot+md)) < 0.0:
            return 1.E30

    distance=np.sum((np.matmul(A,md)-b)**2)
    return distance

def mcmc_search(md):

    n=int(md.size)
    m=1000000
    mcmc_rms=1.E-3
    dp=np.ones(n)*mcmc_rms

    p0=md
    d0=dist(p0)

    #while d0 == 1.E30:
        #for

    p=p0
    accept=np.zeros(n)
    reject=np.zeros(n)

    for i in range(1,m):
        for j in range(n):
            p[j]=p0[j]+dp[j]*np.random.randn()

            d=dist(p)

            if d<d0:
                p0=p
                d0=d
                accept[j]=accept[j]+1
                
            else:
                reject[j]=reject[j]+1
                

        if np.mod(i,1E4) == 0:
            print "i = ",i
            print "p0 = ", p0
            print "d0 = ",d0
            if all(accept==0):
                if all(dp>0.01*mcmc_rms):
                    dp=0.5*dp
                else:
                    exit()
            accept[:]=0
            reject[:]=0
    md=p0

    return md


md=np.zeros(8)
md=mcmc_search(md)

print "mc_md = ",mc_md

    

