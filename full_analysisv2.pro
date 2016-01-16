
pro full_analysisv2,infile,goodfraction,endchans,x1,x2,y1,y2,averagefile,outplot,rippleplot

;********************************************************************
; 



;********************************************************************


;Read in the file and determine basic parameters

map=readfits(infile,header)

mapsize = SIZE(map)

name=strtrim(strmid(infile,8,7),2)

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
raw_medians=dblarr(goodno)
onesub_medians=raw_medians

;find the pixels which are 'good' across the whole range of
;frequencies for the defined 'empty box'

good_intersect=where(finite(good_map[x1:x2,y1:y2,0]),mc,complement=bad_union)

FOR ii=0,goodno-1 DO BEGIN

   good_indices=where(finite(good_map[x1:x2,y1:y2,ii]),tmc,complement=bad_indices)
   good_intersect=setintersection(good_indices,good_intersect)
   
ENDFOR

;'raw' intersection

raw_good_intersect=where(finite(good_map[*,*,0]))

FOR ii=0,goodno-1 DO BEGIN

   raw_good_indices=where(finite(good_map[*,*,ii]))
   raw_good_intersect=setintersection(raw_good_indices,raw_good_intersect)
   
ENDFOR

;Take the median for each channel

FOR gi=0,goodno-1 DO BEGIN

   raw_map_i=good_map[*,*,gi]
   raw_medians[gi]=median(raw_map_i[raw_good_intersect])
   good_map_i=good_map[x1:x2,y1:y2,gi]
   good_medians[gi]=median(good_map_i[good_intersect])

ENDFOR

;Subtract the median from the whole map ---------

median_cube_1=rebin(reform(good_medians,1,1,goodno),nx,ny,goodno)

goodmap=good_map-median_cube_1

FOR gi=0,goodno-1 DO BEGIN

   onesub_map_i=goodmap[*,*,gi]
   onesub_medians[gi]=median(onesub_map_i[raw_good_intersect])

ENDFOR
   
;********************   BINNING ALGORITHM *************************************
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

;********AVERAGING************************************************************

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
start1=fltarr(chunk1_bins)
end1=fltarr(chunk1_bins)

k1=0
k2=binsize1

FOR b1=0,chunk1_bins-1 DO BEGIN

   IF r1 gt 0 THEN BEGIN
      binned1[*,*,b1]=mean(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)
      mfr1[b1]= mean(good_freq[k1:k2])
      errbin1[*,*,b1]=stddev(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)
      start1[b1]=k1
      end1[b1]=k2
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
      start1[b1]=k1
      end1[b1]=k2
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
start2=fltarr(chunk2_bins)
end2=fltarr(chunk2_bins)

k1=endgap
k2=endgap+binsize2

FOR b2=0,chunk2_bins-1 DO BEGIN

   IF r2 gt 0 THEN BEGIN
      binned2[*,*,b2]=mean(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)
      mfr2[b2]= mean(good_freq[k1:k2])
      errbin1[*,*,b2]=stddev(goodmap[*,*,k1:k2],dimension=3,/DOUBLE,/NAN)
      start2[b2]=k1
      end2[b2]=k2
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
      start2[b2]=k1
      end2[b2]=k2
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

IF k2 NE goodno-1 THEN BEGIN
   k2=goodno-1
   end2[chunk2_bins-1]=goodno-1
ENDIF


binned_final=[[[binned1]],[[binned2]]]
error_final=[[[errbin1]],[[errbin2]]]
start_indices=[start1,start2]
end_indices=[end1,end2]
xtable_arr=xtable_arr[0:(chunk1_bins+chunk2_bins-1)]



;******************************* EXTRACT 'RIPPLE' FROM SUBTRACTED MAP
;************************

sub_size=size(sub)

sub_index=indgen(sub_size[3])

bins=chunk1_bins+chunk2_bins

good_sub_medians=fltarr(sub_size[3])
new_medians=good_sub_medians

good_sub_intersect=where(finite(sub[*,*,0]),mc,complement=bad_union)


FOR ii=0,sub_size[3]-1 DO BEGIN

   good_sub_indices=where(finite(sub[*,*,ii]),stmc,complement=bad_sub_indices)
   good_sub_intersect=setintersection(good_sub_indices,good_sub_intersect)
   
ENDFOR


FOR gi=0,sub_size[3]-1 DO BEGIN

   good_sub_map_i=sub[*,*,gi]
   good_sub_medians[gi]=median(good_sub_map_i[good_sub_intersect])

ENDFOR

;; FOR mi=0,mapsize[3]-1 DO BEGIN

;;  good_means[mi]=mean(map[*,*,mi],/double,/nan)

;; ENDFOR 

new_map=fltarr(sub_size[1],sub_size[2],sub_size[3])
median_cube_2=new_map
new_average=fltarr(sub_size[1],sub_size[2],bins)
new_err=new_average
ripple=fltarr(sub_size[3])



FOR ii=0, bins-1 DO BEGIN

chunksize=(end_indices[ii]-start_indices[ii])+1
chunk=goodmap[*,*,start_indices[ii]:end_indices[ii]]
median_chunk=good_sub_medians[start_indices[ii]:end_indices[ii]]
smoothed_m=ts_smooth(median_chunk,5)
ripple[start_indices[ii]:end_indices[ii]]=smoothed_m
median_cube_chunk=rebin(reform(smoothed_m,1,1,chunksize),sub_size[1],sub_size[2],chunksize)
new_chunk=chunk-median_cube_chunk ;;;;; at this point, also save each 'median_cube_chunk' into a larger cube which has the same dimesions as the full map and then output this for use on other maps
new_map[*,*,start_indices[ii]:end_indices[ii]]=new_chunk
median_cube_2[*,*,start_indices[ii]:end_indices[ii]]=median_cube_chunk

new_average[*,*,ii]=mean(new_chunk,dimension=3,/DOUBLE,/NAN)
new_err[*,*,ii]=stddev(new_chunk,dimension=3,/DOUBLE,/NAN)

ENDFOR


FOR gi=0,sub_size[3]-1 DO BEGIN

   new_map_i=new_map[*,*,gi]
   new_medians[gi]=median(new_map_i[raw_good_intersect])

ENDFOR

;**********************************************************************************************
;                                       PLOTTING BIT
;******************************************************************************************

;!PLOT.MULTI=[0,3,2]

set_plot, 'ps'
device,file=outplot
device,/color

plot,fgood_frac,onesub_medians,yrange=[-0.1,0.4],psym=1,title=name+' Raw(red), 1 sub (black) & 2 subs(green)',xtitle='Frequency channel',ytitle='Median value (K)'
oplot,fgood_frac,raw_medians,psym=1,col=1
oplot, fgood_frac,new_medians,psym=1,col=2

for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].stchan,!y.crange,color=3
for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].echan,!y.crange,color=3
;plots,60,0.08,psym=1,col=1
;xyouts,80,0.08,name,col=1

device,/close

print, 'Median plot written to ',outplot 

device,file=rippleplot
device,/color

plot,fgood_frac,good_sub_medians,psym=1,title='Median brightness of background map',xtitle='frequency channel',ytitle='Median brightness (K)'

device, /close

set_plot,'x'

;set_plot, 'ps'
;device,file=outplot2
;device,/color

;plot,fgood_frac,new_medians,psym=1,title='Ripple vs final median',xtitle='Frequency channel',ytitle=';Median value'
;oplot,fgood_frac,ripple,psym=1,col=1
;oplot, fgood_frac,onesub_medians,psym=1,col=2

;for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].stchan,!y.crange,color=3
;for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].echan,!y.crange,color=3
;plots,60,0.08,psym=1,col=1
;xyouts,80,0.08,name,col=1

;device,/close
;set_plot,'x'

;print, 'Ripple comparison plot written to ',outplot2 

;set_plot, 'ps'
;device,file=outplot3
;device,/color

;plot,fgood_frac,good_medians,psym=1,yrange=[-0.7,0.7],title='Raw medians vs box medians',xtitle='Frequency channel',ytitle='Median value'
;name='Channels with 95% or more good pixels'
;oplot,fgood_frac,raw_medians,psym=1,col=2
;oplot, fgood_frac,onesub_medians,psym=1,col=2

;for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].stchan,!y.crange,color=0
;for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].echan,!y.crange,color=0
;plots,60,0.08,psym=1,col=1
;xyouts,80,0.08,name,col=1

;device,/close
;set_plot,'x'

;print, 'Median plot written to ',outplot 


;***********************************************************************************************
;***********************************************************************************************
;                                 OUTPUT ROUTINES
;***********************************************************************************************
;***********************************************************************************************

sxaddpar, header, 'CRPIX1',0.0,'Reference pixel'
sxaddpar,header,'CRVAL1',0.0,'Reference pixel value'
sxaddpar, header, 'EXTEND', 'T'
sxaddpar, header, 'CRPIX3', '1.0', '  Reference pixel',/PDU
sxaddpar, header, 'CDELT3', '1.0', '  Pixel size in world coordinate units',/PDU
sxaddpar, header, 'CRVAL3','1.0', '  Reference pixel value',/PDU
sxaddpar, header, 'CTYPE3', 'CHANNEL', '  3rd type axis',/PDU
sxaddpar, header, 'ORIGINAL',infile,'This file was reduced from the original file shown here',/PDU

;sxaddpar,header_sub,'COMMENT','The channels shown here are channels with more than 95% "good" pixels and can be separated by gaps'

print, 'Writing new extended file to ',averagefile



mwrfits,new_average,averagefile,header

fxbhmake,header_ext_2,1,'SUBEXT'
fxbaddcol,col1,header_ext_2,sub_index,'New channel number'
fxbaddcol,col2,header_ext_2,fgood_frac,'Original channel number'
fxbaddcol,col3,header_ext_2,good_freq,'Frequency (Hz)'
fxbaddcol,col4,header_ext_2,start_indices,'Channels where a bin begins'
fxbaddcol,col5,header_ext_2,end_indices,'Channels where a bin ends'

fxbcreate, unit, averagefile,header_ext_2

fxbwrite,unit,sub_index,col1,1
fxbwrite,unit,fgood_frac,col2,1
fxbwrite,unit,good_freq,col3,1
fxbwrite,unit,start_indices,col4,1
fxbwrite,unit,end_indices,col5,1


fxbfinish,unit

;mwrfits,new_err,averagefile

mwrfits,median_cube_1,averagefile ;cube used for the 'first subtraction'
mwrfits,median_cube_2,averagefile ;cube used for the 'second subtraction'
mwrfits,xtable_arr,averagefile


print,'Done.'


END
