
import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from scipy import ndimage
#from scipy import interpolate as ip
from astropy.io import fits

#take command line input of fits file name
input_file=sys.argv[1]
nside=int(sys.argv[2])
zero_middle=int(sys.argv[3])
scale=int(sys.argv[4])
stokes=sys.argv[5]
bin=int(sys.argv[6])
#gal=int(sys.argv[5])


## def getrdvec(header):
##     naxis1=header['naxis1']
##     naxis2=header['naxis2']
##     crval1=header['crval1']
##     crpix1=header['crpix1']
##     cdelt1=header['cdelt1']
##     crval2=header['crval2']
##     crpix2=header['crpix2']
##     cdelt2=header['cdelt2']
    
##     ra=(arange(naxis1) + 1 - crpix1)*cdelt1 + crval1
##     dec=(arange(naxis2) + 1 - crpix2)*cdelt2 + crval2

##     #Create 2d theta & phi arrays (in colatitude and longitude, radians)

##     theta=np.pi/2.0-(dec*np.pi/180.0)
##     phi=ra*2.*np.pi/360.0+2.*np.pi

##     theta=np.tile(np.reshape(theta,(naxis2,1)),(1,naxis1))
##     phi=np.tile(phi,(naxis2,1))

##     rdvec=hp.pixelfunc.ang2vec(theta,phi)


##     return rdvec

    

    
    
def regrid(nside,inmap):
#def regrid (nside, inmap, channel):
    #regrids the given channel of the GMIMS map into HEALPix format using interpolation

    npix=hp.pixelfunc.nside2npix(nside)

    hdulist=fits.open(inmap)
    header0=hdulist[0].header
    tt=hdulist[0].data

    naxis1=header0['naxis1']
    naxis2=header0['naxis2']
    crval1=header0['crval1']
    crpix1=header0['crpix1']
    cdelt1=header0['cdelt1']
    crval2=header0['crval2']
    crpix2=header0['crpix2']
    cdelt2=header0['cdelt2']
    
    #if header0['NAXIS3']==3:
    #    ttonechan=tt[channel,:,:]
    #    print 'Regridding channel number {!s}.'.format(channel)
    
    #else:
    if len(tt.shape)==3:
        ttonechan=tt[0,:,:]
    else:
        ttonechan=tt
    if scale:
        ttonechan=ttonechan/1000. #stockert only

   

     #Pad the ends for better interpolation with 10 pixels on either side
    npad=10
    ttonechan=np.concatenate((ttonechan[:,-1-npad:-1],ttonechan,ttonechan[:,1:1+npad]),axis=1)
    ttlen=ttonechan.shape[1]
    #print 'Padded array is {!s} pixels wide'.format(ttlen)

     #Define HEALPix directions
    #x,y,z=hp.pixelfunc.pix2vec(nside,np.arange(npix),nest=False)
    #rdvecs=np.column_stack((x,y,z))

    #get HEALPix grid in 2D angles
    theta,phi=hp.pixelfunc.pix2ang(nside,np.arange(npix))
    dec=(np.pi*0.5-theta)*180.0/np.pi
    ra=phi*180.0/np.pi
    if zero_middle:
        ra[np.where(ra>180.)]-=360. #this is used for maps which define their coords as 180 to -180

    #print ra
    #print dec

    #print 'cdelt1:',cdelt1

    #get pixel values to interpolate at
    rapix = (ra - crval1)/cdelt1+crpix1-1+npad
    if np.amin(rapix)<npad:
        rapix[np.where(rapix<npad)] += naxis1 -1
    
    decpix= (dec-crval2)/cdelt2 + crpix2-1

    #print 'first 100 values of rapix:',rapix[0:100]
    #print 'first 100 values of decpix:',decpix[0:100]
    #print 'ttonechan.shape is', ttonechan.shape

    print "Interpolating..."
    #interpolate the tt map onto the healpix grid!
    #coords=np.vstack((rapix,decpix)).T
    #x=np.arange(ttlen)
    #y=np.arange(naxis2)

    #    x2d=np.tile(x,(naxis2,1))
    #    y2d=np.tile(np.reshape(y,(naxis2,1)),(1,ttlen))
    #
    #xycoords=np.vstack((x2d.ravel(),y2d.ravel())).T

    #temptt=ip.griddata(xycoords,ttonechan.ravel(),coords,method='cubic')
   
    temptt=ndimage.map_coordinates(ttonechan,[decpix,rapix],order=3,cval=-1.6375e+30,prefilter=False)

    #f=ip.interp2d(x,y,ttonechan,kind='cubic',copy=True,bounds_error=False,fill_value=-1.6375e+30)
    #temptt=f(decpix,rapix)
    

    print "...done!"

    

    #print 'temptt.shape is', temptt.shape
    #print 'first 500 values of temptt:',temptt[0:500]

    #r=hp.rotator.Rotator(coord=['G','C'])
    #theta_gal,phi_gal=r(theta,phi)
    #galpixnos=hp.pixelfunc.ang2pix(nside,theta_gal,phi_gal)

    return temptt#,galpixnos



map=regrid(nside,input_file) #galpixnos
#map,galpixnos=regrid(nside,input_file, bin)

#if gal:
map_nest=hp.pixelfunc.reorder(map,r2n=True)
   

#else:

    #map_gal=map[galpixnos] 
    #map_gal_nest=hp.pixelfunc.reorder(map_gal,r2n=True)
    

 



#outfile=stokes+'cube_bin_'+str(bin)+'_nest_'+str(nside)+'.fits'
outfile=stokes+'_bin'+str(bin)+'_nest_'+str(nside)+'.fits'
#outfile=input_file.rsplit('.',1)[0]+'_healpix_nest_512.fits'
#if gal:
 #   hp.mollview(map_nest,nest=True,min=3.0, max=4.0)
 #   plt.show()
 #   hp.write_map(outfile,map_nest,nest=True)

#else:
hp.mollview(map_nest,nest=True)#,min=3.0, max=4.0)
plt.show()
hp.write_map(outfile,map_nest,nest=True)

    
print 'Wrote map to ',outfile




    

   

   

    

    

    


    
    

    
    
    
