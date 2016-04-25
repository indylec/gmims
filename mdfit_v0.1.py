import sys
import numpy as np
import healpy as hp
from scipy.stats.mstats import theilslopes
from astropy.io import ascii

lo_in=sys.argv[1]
hi_in=sys.argv[2]
#mask_in=sys.argv[3]

nu_lo=408000000.
nu_hi=1420000000.

#text_name=lo_name+'-'+hi_name+'_coeffs.txt'

nside_map=256
npix_map=hp.nside2npix(nside_map)
nside_regions=4
npix=hp.nside2npix(nside_regions)
npix_sub=npix_map/npix

#mask_map=hp.read_map(mask_in,nest=True)

map_lo=hp.read_map(lo_in,nest=True)
map_lo_ma=np.ma.masked_where(map_lo==hp.UNSEEN,map_lo)
#map_lo_ma=np.ma.masked_where(mask_map==1.0,map_lo)

#print npix_map
#print map_lo_ma.compressed().size

map_hi=hp.read_map(hi_in,nest=True)
map_hi_ma=np.ma.masked_where(map_hi==hp.UNSEEN,map_hi)
#map_hi_ma=np.ma.masked_where(mask_map==1.0,map_lhi)

#print npix_map
#print map_hi_ma.compressed().size

print 'Finding coldest pixels...'

npos=0.
b_pos=np.zeros(1000)
A_pos=np.zeros((1000,8))

md=np.asarray([8.9,  3.2,  0.7, -0.8])

map_pos=map_lo
map_pos[map_pos==hp.UNSEEN]=np.nan
while any(~np.isnan(map_pos)):
    myloc=np.nanargmin(map_pos)
    b_pos[npos]=map_pos[myloc]
    A_pos[npos,0]=1.0
    vec=hp.pix2vec(nside_map,myloc,nest=True)
    A_pos[npos,1:4]=vec
    #print npos, b_pos[npos], np.matmul(A_pos[npos,:4],md)
    npos+=1
    pixlist=hp.query_disc(nside_map,vec,10.*np.pi/180.,nest=True)
    map_pos[pixlist]=np.nan

print npos

#exit()

map_pos=map_hi
map_pos[map_pos==hp.UNSEEN]=np.nan
while any(~np.isnan(map_pos)):
    myloc=np.nanargmin(map_pos)
    b_pos[npos]=map_pos[myloc]
    A_pos[npos,4]=1.0
    vec=hp.pix2vec(nside_map,myloc,nest=True)
    A_pos[npos,5:8]=vec
    #print b_pos[npos], A_pos[npos], npos
    npos+=1
    pixlist=hp.query_disc(nside_map,vec,10.*np.pi/180.,nest=True)
    map_pos[pixlist]=np.nan

print npos

#exit()

#get slopes and offsets

print 'Estimating regression coefficients...'

params=np.empty((5,npix))

for i in range(npix):
    #print "Working on pixel "+str(i)
    params[0,i]=i
    x=map_lo_ma[i*npix_sub:(i+1)*npix_sub]
    y=map_hi_ma[i*npix_sub:(i+1)*npix_sub]
    if (x.compressed().size == y.compressed().size) & (x.compressed().size > 0.25*npix_sub):
        params[1:,i]=theilslopes(y,x)
    else:
        params[1:,i]=hp.UNSEEN

np.save("md_params",params)
np.save("md_A_pos",A_pos)
np.save("md_b_pos",b_pos)
np.save("md_npos",npos)

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

print numval

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


#print A.shape

x,resx,rank,sing=np.linalg.lstsq(A,b)

print "x= ",x

res=np.matmul(A,x.T)-b

mu_res=np.mean(res)
sigma_res=np.sqrt(np.var(res))

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

print "Starting mcmc search..."

md_tot=np.zeros(8)

def dist(md):

    for i in range (int(npos)):
        if b_pos[i]-np.sum(A_pos[i,:]*(md_tot+md)) < 0.0:
            return 1.E30

    distance=np.sum((np.matmul(A,md)-b)**2)
    return distance

def mcmc_search(md):

    n=int(md.size)
    m=1000000
    mcmc_rms=1.E-2
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
            print "accept = ",accept
            print "reject = ",reject
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

#exit()
    






