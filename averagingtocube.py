#Extension of averaging.py, lets you bin gmims cube to a given number of bins but outputs this as a cube rather than individual maps 

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

outcube=np.empty((nobins,ny,nx))
freqmean=np.empty(nobins)
schan_arr=np.empty(nobins)
echan_arr=np.empty(nobins)

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
    outcube[i,:,:]=np.mean(data[schan:echan,:,:],axis=0)
    freqmean[i]=np.mean(freq[schan:echan])
    #print 'Average frequency is {0} Hz'.format(freqmean)
    schan_arr[i]=schan
    echan_arr[i]=echan
    schan=echan

newhead=header.copy()
        
    
newhead['BUNIT']=('K')
newhead['OBJECT']=('GMIMS {0} binned cube'.format(stokes))
newhead['ORIGINAL']=(filein)
#newhead['COMMENT']='Frequency is {0} Hz'.format(freqmean)
newhead['COMMENT']='{0} channels in original'.format(nchan)
#newhead['COMMENT']='Bin starts at channel {0}'.format(schan)
#newhead['COMMENT']='Bin ends at channel {0}'.format(echan)
newhead['COMMENT']='Bin/frequency lookup table in extension'
    
#create primary HDU and HDUlist
binnedhdu=fits.PrimaryHDU(outcube,header=newhead)
binnedhdulist=fits.HDUList([binnedhdu])

bin_numbers=np.arange(nobins)+1

#create binary table columns containing info about each plane of the cube
col1=fits.Column(name='bin',format='d',array=bin_numbers)
col2=fits.Column(name='start_channel',format='d', array=schan_arr)
col3=fits.Column(name='end_channel',format='d', array=echan_arr)
col4=fits.Column(name='mean_frequency',format='E',array=freqmean)

cols=fits.ColDefs([col1,col2,col3,col4])
bintablehdu=fits.BinTableHDU.from_columns(cols)

binnedhdulist.append(bintablehdu)




if discard==True:
    outname=stokes+'_cube_bin_discard.fits'
if discard==False:
    outname=stokes+'_cube_bin_no_discard.fits'
    
binnedhdulist.writeto(outname)
