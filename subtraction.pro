pro subtraction,infile,goodfraction,endchans,x1,x2,y1,y2,outfile,subtract_map

;********************************************************************
; Routine that flags out channels with less than goodfraction of good
; data, then finds the median value within a box delimited by x1,x2,y1,y2
; It then subtracts this value from the whole map. , averages the
; 2000-channel data cube down to around 10 channels of equal size and
; subtracts these averages from the corresponding channel blocks in
; the initial map. 

; The outputs are two fits files: one containing the averaged
; channels, an error map produced using the std.dev of the
; binning (as an image extension) and a binary table extension
; specifying the start and end channels of the chunks going into each
; bin; the other contains the resulting 'subtracted' map which gives
; an idea of the variations in zero-level across the frequency
;



;********************************************************************


;Read in the file and determine basic parameters

map=readfits(infile,header)

mapsize = SIZE(map)

nx = mapsize[1]
ny = mapsize[2]
nchan = mapsize[3]

npix = nx * ny

ngood = LONARR(nchan)

;Calculate the appropriate frequency range from the header values

f0=sxpar(header,'CRVAL3')
df=sxpar(header,'CDELT3')
refpix=sxpar(header, 'CRPIX3')

freq=(findgen(nchan)-refpix)*df+f0

;determine the good-to-total pixel ratio for each channel

FOR ichan = 0, nchan-1 DO BEGIN 
   void = WHERE(finite(map[*,*,ichan]), ngood_count)
   ngood[ichan]=ngood_count
ENDFOR

fgood= FLOAT(ngood)/npix
fgood_frac= WHERE(fgood GE goodfraction,goodno) ;this is the 'good' subset of channels

fg=fgood[fgood_frac]

good_freq=freq[fgood_frac]

good_map=map[*,*,fgood_frac]

good_medians=dblarr(goodno)

;find the pixels which are 'good' across the whole range of frequencies

good_intersect=where(finite(good_map[x1:x2,y1:y2,0]),mc,complement=bad_union)

FOR ii=0,goodno-1 DO BEGIN

   good_indices=where(finite(good_map[x1:x2,y1:y2,ii]),tmc,complement=bad_indices)
   good_intersect=setintersection(good_indices,good_intersect)
   
ENDFOR

;Take the median for each channel

FOR gi=0,goodno-1 DO BEGIN

   good_map_i=good_map[x1:x2,y1:y2,gi]
   good_medians[gi]=median(good_map_i[good_intersect])

ENDFOR

;Subtract the median from the whole map

median_cube=rebin(reform(good_medians,1,1,goodno),nx,ny,goodno)

goodmap=good_map-median_cube

;******************************************************************************

;figure out where the 'gap' is

FOR ii=1,goodno-1 DO BEGIN

   IF fgood_frac[ii]-fgood_frac[ii-1] GE 200 THEN BEGIN 

      startgap=ii-1
      endgap=ii

   ENDIF

ENDFOR

print,FORMAT='("Gap begins at channel ",I0," and ends at channel ",I0)',fgood_frac[startgap],fgood_frac[endgap]

;split the good channels up into two chunks on either side of the gap

good_indices_1=fgood_frac[0:startgap]
good_indices_2=fgood_frac[endgap:goodno-1]

size1=size(good_indices_1)
size2=size(good_indices_2)

print,FORMAT='("First chunk contains ",I0," channels; second chunk contains ",I0," channels")',size1[1],size2[1]

bin_size=goodno/endchans ;this is the nominal bin size that will be adjusted according to the number of channels

print, FORMAT='("Nominal bin size is ",I0," channels per bin, assuming ",I0," good channels and ",I0," bins")',bin_size,goodno,endchans
;************************************************************

;Prepare the binary table output
xtable={bin:0,stchan:0,echan:0,mfr:0L}
xtable_arr=replicate(xtable,20)

;The final subtracted output


sub=fltarr(nx,ny,goodno,/nozero)


;now we deal with the chunks individually 

chunk1_bins=size1[1]/bin_size

r1=size1[1] mod bin_size

binsize1=bin_size

while r1 gt chunk1_bins do begin

   binsize1+=(r1/chunk1_bins)
   r1=r1 mod chunk1_bins

endwhile

print, FORMAT='("There are ",I0," bins in the first chunk. Current bin size is ",I0," ; remainder is ",I0)',chunk1_bins,binsize1,r1

;r1 is now smaller than the number of bins and we will just add one channel to
;the first n=r1 bins while averaging
 

binned1=fltarr(nx,ny,chunk1_bins,/nozero)
mfr1=fltarr(chunk1_bins,/nozero)
errbin1=fltarr(nx,ny,chunk1_bins,/nozero)

k1=0
k2=binsize1

FOR b1=0,chunk1_bins-1 DO BEGIN

   IF r1 gt 0 THEN BEGIN
      binned1[*,*,b1]=mean(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)
      mfr1[b1]= mean(good_freq[k1:k2])
      errbin1[*,*,b1]=stddev(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)

      sub[*,*,k1:k2]=goodmap[*,*,k1:k2]-rebin(binned1[*,*,b1],nx,ny,binsize1+1)
      
      
      r1--

      xtable.bin=b1+1
      xtable.stchan=fgood_frac[k1]
      xtable.echan=fgood_frac[k2-1]
      xtable.mfr=mfr1[b1]

      xtable_arr[b1]=xtable
      

      IF b1 LT chunk1_bins-1 then begin
         k1=k2+1
         k2=k2+1+binsize1
      ENDIF

      

   ENDIF ELSE BEGIN
      
      binned1[*,*,b1]=mean(goodmap[*,*,k1:k2-1],dimension=3,/DOUBLE,/NAN)
      mfr1[b1]= mean(good_freq[k1:k2-1])
      errbin1[*,*,b1]=stddev(goodmap[*,*,k1:k2-1],dimension=3,/DOUBLE,/NAN)

      sub[*,*,k1:k2-1]=goodmap[*,*,k1:k2-1]-rebin(binned1[*,*,b1],nx,ny,binsize1)


      xtable.bin=b1+1
      xtable.stchan=fgood_frac[k1]
      xtable.echan=fgood_frac[k2-1]
      xtable.mfr=mfr1[b1]

      xtable_arr[b1]=xtable
      
      IF b1 LT chunk1_bins-1 THEN BEGIN
         k1=k2
         k2=k2+binsize1
      ENDIF

      


   ENDELSE

ENDFOR

;Check that we have reached the end of chunk 1

print, FORMAT='("k2 is ",I0,"; end of first chunk is ",I0)',k2-1,startgap

;Start on chunk 2

chunk2_bins=size2[1]/bin_size

r2=size2[1] mod bin_size

binsize2=bin_size



while r2 gt chunk2_bins do begin

   binsize2+=(r2/chunk2_bins)
   r2=r2 mod chunk2_bins

endwhile

print, FORMAT='("There are ",I0," bins in the second chunk. Current bin size is ",I0," ; remainder is ",I0)',chunk2_bins,binsize2,r2

;r2 is now smaller than the number of bins and we will just add one channel to
;the first n=r2 bins while averaging
 

binned2=fltarr(nx,ny,chunk2_bins,/nozero)
mfr2=fltarr(chunk2_bins,/nozero)
errbin2=fltarr(nx,ny,chunk2_bins,/nozero)

k1=endgap
k2=endgap+binsize2

FOR b2=0,chunk2_bins-1 DO BEGIN

   IF r2 gt 0 THEN BEGIN
      binned2[*,*,b2]=mean(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)
      mfr2[b2]= mean(good_freq[k1:k2])
      errbin1[*,*,b2]=stddev(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)

      sub[*,*,k1:k2]=goodmap[*,*,k1:k2]-rebin(binned2[*,*,b2],nx,ny,binsize2+1)



      r2--

      xtable.bin=chunk1_bins+b2+1
      xtable.stchan=fgood_frac[k1]
      xtable.echan=fgood_frac[k2]
      xtable.mfr=mfr2[b2]

      xtable_arr[chunk1_bins+b2]=xtable
     
      IF b2 LT chunk2_bins-1 then begin
         k1=k2+1
         k2=k2+1+binsize2
      ENDIF

   ENDIF ELSE BEGIN
      
      binned2[*,*,b2]=mean(goodmap[*,*,k1:k2-1],dimension=3,/DOUBLE,/NAN)

      mfr2[b2]= mean(good_freq[k1:k2-1])
      errbin1[*,*,b2]=stddev(goodmap[*,*,k1:k2-1],dimension=3,/DOUBLE,/NAN)
      
      sub[*,*,k1:k2-1]=goodmap[*,*,k1:k2-1]-rebin(binned2[*,*,b2],nx,ny,binsize2)
      
      
      r2--

      xtable.bin=chunk1_bins+b2+1
      xtable.stchan=fgood_frac[k1]
      xtable.echan=fgood_frac[k2-1]
      xtable.mfr=mfr2[b2]

      xtable_arr[chunk1_bins+b2]=xtable
      

      IF b2 LT chunk2_bins-1 THEN BEGIN
         k1=k2
         k2=k2+binsize2
      ENDIF


   ENDELSE

ENDFOR

print, FORMAT='("k2 is ",I0,"; end of second chunk is ",I0)',k2,goodno-1






;***********************************************************************************************
;***********************************************************************************************
;                                 OUTPUT ROUTINES
;***********************************************************************************************
;***********************************************************************************************




binned_final=[[[binned1]],[[binned2]]]
error_final=[[[errbin1]],[[errbin2]]]
xtable_arr=xtable_arr[0:(chunk1_bins+chunk2_bins-1)]

sub_index=findgen(goodno)


header_sub=header




sxaddpar, header, 'EXTEND', 'T'
sxaddpar, header, 'CRPIX3', '1.0', '  Reference pixel',/PDU
sxaddpar, header, 'CDELT3', '1.0', '  Pixel size in world coordinate units',/PDU
sxaddpar, header, 'CRVAL3','1.0', '  Reference pixel value',/PDU
sxaddpar, header, 'CTYPE3', 'CHANNEL', '  3rd type axis',/PDU
sxaddpar, header, 'ORIGINAL',infile,'This file was reduced from the original file shown here',/PDU

sxaddpar,header_sub,'COMMENT','The channels shown here are channels with more than 95% "good" pixels and can be separated by gaps'

print, 'Writing new extended file to ',outfile

mwrfits,binned_final,outfile,header
mwrfits,error_final,outfile
mwrfits,xtable_arr,outfile

print, 'Writing subtraction map to ',subtract_map
mwrfits, sub,subtract_map,header_sub

fxbhmake,header_sub_ext,1,'SUBEXT'
fxbaddcol,col1,header_sub_ext,sub_index,'New channel number'
fxbaddcol,col2,header_sub_ext,fgood_frac,'Original channel number'
fxbaddcol,col3,header_sub_ext,good_freq,'Frequency (Hz)'

fxbcreate, unit, subtract_map,header_sub_ext

fxbwrite,unit,sub_index,col1,1
fxbwrite,unit,fgood_frac,col2,1
fxbwrite,unit,good_freq,col3,1

fxbfinish,unit

fxbcreate,unit2,outfile,header_sub_ext,extensionno

print,FORMAT='("The details of the original and new channels & frequency were written as extension number",I0,"to",A0)', extensionno,outfile

fxbwrite,unit2,sub_index,col1,1
fxbwrite,unit2,fgood_frac,col2,1
fxbwrite,unit2,good_freq,col3,1

fxbfinish,unit2

print,'Done.'


END

