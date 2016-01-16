pro averaging, infile, goodfraction,endchans, outfile

;********************************************************************
; Routine to average the 2000-channel GMIMS cube to 'endchans'
; channels.
; Outputs fits file with given number of channels, variance map of the
; binning and a binary extension table with mean frequency per binned
; channel and range of channels. Automatically removes channels with less 
; than goodfraction fraction of good pixels (say 95%)
;

map=mrdfits(infile)

mapsize = SIZE(map)

nx = mapsize[1]
ny = mapsize[2]
nchan = mapsize[3]

npix = nx * ny

ngood = LONARR(nchan)

;determine the good-to-total pixel ratio for each channel

FOR ichan = 0, nchan-1 DO BEGIN 
   void = WHERE(finite(map[*,*,ichan]), ngood_count)
   ngood[ichan]=ngood_count
ENDFOR

fgood= FLOAT(ngood)/npix
fgood_frac= WHERE(fgood GE goodfraction,goodno) ;this is the 'good' subset of channels

goodmap=map[*,*,fgood_frac]

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



k1=0
k2=binsize1

FOR b1=0,chunk1_bins-1 DO BEGIN

   IF r1 gt 0 THEN BEGIN
      binned1[*,*,b1]=mean(goodmap[*,*,k1:k2],dimension=3)
      r1--
      IF b1 LT chunk1_bins-1 then begin
         k1=k2+1
         k2=k2+1+binsize1
      ENDIF

   ENDIF ELSE BEGIN
      
      binned1[*,*,b1]=mean(goodmap[*,*,k1:k2-1],dimension=3)
      IF b1 LT chunk1_bins-1 THEN BEGIN
         k1=k2
         k2=k2+binsize1
      ENDIF


   ENDELSE

ENDFOR

;Check that we have reached the end of chunk 1

print, FORMAT='("k2 is ",I0,"; end of first chunk is ",I0)',k2,startgap

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

k1=endgap
k2=endgap+binsize2

FOR b2=0,chunk2_bins-1 DO BEGIN

   IF r2 gt 0 THEN BEGIN
      binned2[*,*,b2]=mean(goodmap[*,*,k1:k2],dimension=3)
      r2--
      IF b2 LT chunk2_bins-1 then begin
         k1=k2+1
         k2=k2+1+binsize2
      ENDIF

   ENDIF ELSE BEGIN
      
      binned2[*,*,b2]=mean(goodmap[*,*,k1:k2-1],dimension=3)
      IF b2 LT chunk2_bins-1 THEN BEGIN
         k1=k2
         k2=k2+binsize2
      ENDIF


   ENDELSE

ENDFOR


print, FORMAT='("k2 is ",I0,"; end of second chunk is ",I0)',k2,goodno-1

binned_final=[[[binned1]],[[binned2]]]

mwrfits,binned_final,outfile

END






      




      
      
