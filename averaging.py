#Somewhat different to averaging.pro -- this code averages the GMIMS data into the desired number of bins, keeping the number of channels per bin as close to equal as possible or discarding channels if you so choose. Borrows from stack.py programs used with galfacts data


import sys
import numpy as np 
from astropy.io import fits


filein=sys.argv[1]
nobins=int(sys.argv[2])
stokes=sys.argv[3]
discard=int(sys.argv[4])


data=fits.getdata(filein)
header=fits.getheader(filein)

nx=header['NAXIS1']
ny=header['NAXIS2']
nchan=header['NAXIS3']

f0=header['CRVAL3']
df=header['CDELT3']
freq=np.arange(nchan)*df+f0

st_binwidth=nchan/nobins
rmdr=nchan%nobins

#If we discard channels, get rid of half the remainder at the front and half at the back; if not, add 1 channel to each bin until there are none left.

#if discard == True:
#    first=rmdr/2
    
#    for i in range(nobins):
        
#        schan=first+i*binwidth
#        echan=schan+binwidth
#        mapmean=np.mean(data[schan:echan,:,:])
#        freqmean=np.mean(freq[schan:echan])

#        newhead=header.copy()
        
 #       newhead['CRVAL3']=(str(freqmean))
 #       newhead['CRPIX3']=(1.0)
 #       newhead['BUNIT']=('K')
 #       newhead['OBJECT']=('GMIMS {0} bin {1}'.format(stokes,str(i)))
 #       newhead['ORIGINAL']=(filein)
 #       newhead['COMMENT']='{0} channels in original'.format(nchan)
 #       newhead['COMMENT']='Bin starts at channel {0}'.format(schan)
 #       newhead['COMMENT']='Bin ends at channel {0}'.format(echan)

 #       fits.writeto(stokes+'_bin_'+str(i),mapmean,newhead)

if discard==True:
    schan=rmdr/2
    print 'Discard is ON'
    outname=stokes+'_bin_discard'+str(i)+'.fits'
else:
    schan=0
    print 'Discard is OFF'
    
    
for i in range(nobins):
    print 'averaging into bin {0}'.format(i+1)
    if discard==False and i < rmdr:
        binwidth = st_binwidth+1
    else:
        binwidth=st_binwidth

    #schan=first+i*binwidth
    echan=schan+binwidth
    print 'binwidth is {0}'.format(binwidth)
    print 'start channel is {0}, end channel is {1}'.format(schan,echan)
    mapmean=np.mean(data[schan:echan,:,:],axis=0)
    freqmean=np.mean(freq[schan:echan])
    print 'Average frequency is {0} Hz'.format(freqmean)

    newhead=header.copy()
        
    del newhead['CRVAL3']
    del newhead['CRPIX3']
    del newhead['CDELT3']
    del newhead['CTYPE3']
    newhead['BUNIT']=('K')
    newhead['OBJECT']=('GMIMS {0} bin {1}'.format(stokes,str(i)))
    newhead['ORIGINAL']=(filein)
    newhead['COMMENT']='Frequency is {0} Hz'.format(freqmean)
    newhead['COMMENT']='{0} channels in original'.format(nchan)
    newhead['COMMENT']='Bin starts at channel {0}'.format(schan)
    newhead['COMMENT']='Bin ends at channel {0}'.format(echan)

    schan=echan

    if discard==True:
        outname=stokes+'_bin_discard_'+str(i)+'.fits'
    if discard==False:
        outname=stokes+'_bin_no_discard_'+str(i)+'.fits'
    fits.writeto(outname,mapmean,newhead)
    print 'binned map written to '+outname

            
        

        

    
