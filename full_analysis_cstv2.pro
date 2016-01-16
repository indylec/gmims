pro full_analysis_cstv2,infile,reffile,averagefile

;********************************************************************
; 



;********************************************************************

;********************************************************************
;********Read in the file and determine basic parameters*************
;********************************************************************
map=readfits(infile,header)

map=map[0:360,*,*]

median_cube_1=readfits(reffile,EXTEN_NO=2)
median_cube_2=readfits(reffile,EXTEN_NO=3)

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

;Read in the 'good' indices and bin limits from the first binning 

;reftable_1=mrdfits(reffile,2,/use_colnum)
reftable_2=mrdfits(reffile,1)

fgood_frac=reftable_2.(1)
good_freq=reftable_2.(2)
start_indices=reftable_2.(3)
end_indices=reftable_2.(4)

gsize=size(fgood_frac)
goodno=gsize[1]

good_map=map[*,*,fgood_frac]

good_medians=dblarr(goodno)
raw_medians=dblarr(goodno)
onesub_medians=raw_medians



;; ;find the pixels which are 'good' across the whole range of frequencies

;; good_intersect=where(finite(good_map[x1:x2,y1:y2,0]),mc,complement=bad_union)

;; FOR ii=0,goodno-1 DO BEGIN

;;    good_indices=where(finite(good_map[x1:x2,y1:y2,ii]),tmc,complement=bad_indices)
;;    good_intersect=setintersection(good_indices,good_intersect)
   
;; ENDFOR

;; ;'raw' intersection

;; raw_good_intersect=where(finite(good_map[*,*,0]))

;; FOR ii=0,goodno-1 DO BEGIN

;;    raw_good_indices=where(finite(good_map[*,*,ii]))
;;    raw_good_intersect=setintersection(raw_good_indices,raw_good_intersect)
   
;; ENDFOR

;; ;Take the median for each channel

;; FOR gi=0,goodno-1 DO BEGIN

;;    raw_map_i=good_map[*,*,gi]
;;    raw_medians[gi]=median(raw_map_i[raw_good_intersect])   
;;    good_map_i=good_map[x1:x2,y1:y2,gi]
;;    good_medians[gi]=median(good_map_i[good_intersect])

;; ENDFOR

;Subtract the median from the whole map

;median_cube=rebin(reform(good_medians,1,1,goodno),nx,ny,goodno)

goodmap=good_map-median_cube_1

;FOR gi=0,goodno-1 DO BEGIN

 ;  onesub_map_i=goodmap[*,*,gi]
  ; onesub_medians[gi]=median(onesub_map_i[raw_good_intersect])

;ENDFOR

bsize=size(start_indices)

bins=bsize[1]

;print, bins

;binned=fltarr(nx,ny,bins)

;print, size(binned)
;errbin=binned
;mfr=fltarr(bins)

;sub=fltarr(nx,ny,goodno,/nozero)

;xtable={bin:0,stchan:0,echan:0,mfr:0L}
;xtable_arr=replicate(xtable,bins)



;**************FIRST AVERAGING AND SUBTRACTION***************************

;FOR ii=0,bins-1 DO BEGIN

;print, FORMAT='("Averaging bin number ",I0,".")',ii

;binsize=end_indices[ii]-start_indices[ii]

;print, FORMAT='("This bin contains ",I0," channels.")',binsize

;print,"debug check: this is the SIZE arary for the 'binned' array"
;print,size(binned)

;binned[*,*,ii]=mean(goodmap[*,*,start_indices[ii]:end_indices[ii]],dimension=3,/double,/nan)

;mfr[ii]=mean(good_freq[start_indices[ii]:end_indices[ii]])

;errbin[*,*,ii]=stddev(goodmap[*,*,start_indices[ii]:end_indices[ii]],dimension=3,/double,/nan)

;sub[*,*,start_indices[ii]:end_indices[ii]]=goodmap[*,*,start_indices[ii]:end_indices[ii]]-rebin(binned[*,*,ii],nx,ny,binsize+1)

;xtable_arr[ii].bin=[ii+1]
;xtable_arr[ii].stchan=fgood_frac[start_indices[ii]]
;xtable_arr[ii].echan=fgood_frac[end_indices[ii]]
;xtable_arr[ii].mfr=mfr[ii]


;ENDFOR






;******************************* EXTRACT 'RIPPLE' FROM SUBTRACTED MAP
;************************

;; sub_size=size(sub)

;; sub_index=indgen(sub_size[3])

;; good_sub_medians=fltarr(sub_size[3])
;; new_medians=good_sub_medians

;; good_sub_intersect=where(finite(sub[*,*,0]),mc,complement=bad_union)


;; FOR ii=0,sub_size[3]-1 DO BEGIN

;;    good_sub_indices=where(finite(sub[*,*,ii]),stmc,complement=bad_sub_indices)
;;    good_sub_intersect=setintersection(good_sub_indices,good_sub_intersect)
   
;; ENDFOR


;; FOR gi=0,sub_size[3]-1 DO BEGIN

;;    good_sub_map_i=sub[*,*,gi]
;;    good_sub_medians[gi]=median(good_sub_map_i[good_sub_intersect])

;;  ENDFOR

;; FOR mi=0,mapsize[3]-1 DO BEGIN

;;  good_means[mi]=mean(map[*,*,mi],/double,/nan)

;; ENDFOR 

new_map=fltarr(nx,ny,goodno)
sm_map=new_map
new_average=fltarr(nx,ny,bins)
new_err=new_average
ripple=fltarr(goodno)

;; FOR ii=0, bins-1 DO BEGIN

;; chunksize=(end_indices[ii]-start_indices[ii])+1
;; chunk=goodmap[*,*,start_indices[ii]:end_indices[ii]]
;; median_chunk=good_sub_medians[start_indices[ii]:end_indices[ii]]
;; smoothed_m=ts_smooth(median_chunk,5)
;; median_cube_chunk=rebin(reform(smoothed_m,1,1,chunksize),sub_size[1],sub_size[2],chunksize)
;; new_chunk=chunk-median_cube_chunk
;; new_map[*,*,start_indices[ii]:end_indices[ii]]=new_chunk
;; sm_map[*,*,start_indices[ii]:end_indices[ii]]=median_cube_chunk

;; new_average[*,*,ii]=mean(new_chunk,dimension=3,/DOUBLE,/NAN)
;; new_err[*,*,ii]=stddev(new_chunk,dimension=3,/DOUBLE,/NAN)

;; ENDFOR

new_map=goodmap-median_cube_2

FOR ii=0,bins-1 DO BEGIN

new_chunk=new_map[*,*,start_indices[ii]:end_indices[ii]]
new_average[*,*,ii]=mean(new_chunk,dimension=3,/double,/nan)
new_err[*,*,ii]=stddev(new_chunk,dimension=3,/double,/nan)


ENDFOR


;; FOR gi=0,sub_size[3]-1 DO BEGIN

;;    new_map_i=new_map[*,*,gi]
;;    new_medians[gi]=median(new_map_i[raw_good_intersect])

;; ENDFOR

;**********************************************************************************************
;                                       PLOTTING BIT
;******************************************************************************************

;; !PLOT.MULTI=[0,3,2]

;; set_plot, 'ps'
;; device,file=outplot
;; device,/color

;; plot,fgood_frac,onesub_medians,yrange=[0.1,1.0],psym=1,title=name+' Raw(red), 1 sub (black) & 2 subs(green)',xtitle='Frequency channel',ytitle='Median value (K)'
;; oplot,fgood_frac,raw_medians,psym=1,col=1
;; oplot, fgood_frac,new_medians,psym=1,col=2

;; for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].stchan,!y.crange,color=3
;; for pp=0,bins-1 do oplot,fltarr(2)+xtable_arr[pp].echan,!y.crange,color=3

;; device, /close
;; set_plot,'x'

;; print, 'Median plot written to ',outplot


;***********************************************************************************************
;***********************************************************************************************
;                                 OUTPUT ROUTINES
;***********************************************************************************************
;***********************************************************************************************


sxaddpar, header, 'EXTEND', 'T'
sxaddpar, header, 'CRPIX3', '1.0', '  Reference pixel',/PDU
sxaddpar, header, 'CDELT3', '1.0', '  Pixel size in world coordinate units',/PDU
sxaddpar, header, 'CRVAL3','1.0', '  Reference pixel value',/PDU
sxaddpar, header, 'CTYPE3', 'CHANNEL', '  3rd type axis',/PDU
sxaddpar, header, 'ORIGINAL',infile,'This file was reduced from the original file shown here',/PDU

;sxaddpar,header_sub,'COMMENT','The channels shown here are channels with more than 95% "good" pixels and can be separated by gaps'

print, 'Writing new extended file to ',averagefile

mwrfits,new_average,averagefile,header
mwrfits,new_err,averagefile
;mwrfits,xtable_arr,averagefile

;; fxbhmake,header_ext_2,1,'SUBEXT'
;; fxbaddcol,col1,header_ext_2,sub_index,'New channel number'
;; fxbaddcol,col2,header_ext_2,fgood_frac,'Original channel number'
;; fxbaddcol,col3,header_ext_2,good_freq,'Frequency (Hz)'
;; fxbaddcol,col4,header_ext_2,start_indices,'Channels where a bin begins'
;; fxbaddcol,col5,header_ext_2,end_indices,'Channels where a bin ends'

;; fxbcreate, unit, averagefile,header_ext_2

;; fxbwrite,unit,sub_index,col1,1
;; fxbwrite,unit,fgood_frac,col2,1
;; fxbwrite,unit,good_freq,col3,1
;; fxbwrite,unit,start_indices,col4,1
;; fxbwrite,unit,end_indices,col5,1


;; fxbfinish,unit

print,'Done.'


END
