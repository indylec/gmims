pro sub2, rawfile,averagefile,x1,x2,y1,y2,finalout
  
;**************************************************************************
;description goes here
;
;***************************************************************************

; Data input routines: map, binary extensions

map=mrdfits(rawfile,0,header)


avgtable=mrdfits(averagefile,2,/use_colnum) ; contains bin number, start channel, end channel and mean frequency
indextable=mrdfits(averagefile,3)  ; contains indices of the good, subtracted map, corresponding indices of the raw map, and the  frequencie of the good channels

avgtablesize=size(avgtable)

bins=avgtablesize[1]

fgood_frac=indextable.(1)
good_freq=indextable.(2)

startchans=fltarr(bins)
endchans=fltarr(bins)
start_indices=fltarr(bins)
end_indices=fltarr(bins)

nostart=n_elements(startchans)
noend=n_elements(endchans)

for ii=0,bins-1 do begin
startchans[ii]=avgtable[ii].c2
endchans[ii]=avgtable[ii].c3
start_indices[ii]=WHERE(fgood_frac EQ avgtable[ii].c2)
end_indices[ii]=WHERE(fgood_frac EQ avgtable[ii].c3)

endfor

good_map=map[*,*,fgood_frac]

mapsize=size(good_map)

good_medians=dblarr(mapsize[3])

new_map=fltarr(mapsize[1],mapsize[2],mapsize[3])

sm_map=new_map

good_means=good_medians

new_medians=good_medians


good_intersect=where(finite(map[*,*,0]),mc,complement=bad_union)


FOR ii=0,mapsize[3]-1 DO BEGIN

   good_indices=where(finite(map[*,*,ii]),tmc,complement=bad_indices)
   good_intersect=setintersection(good_indices,good_intersect)
   
ENDFOR


FOR gi=0,mapsize[3]-1 DO BEGIN

   good_map_i=map[x1:x2,y1:y2,gi]
   good_medians[gi]=median(good_map_i[good_intersect])

ENDFOR


median_cube=rebin(reform(good_medians,1,1,goodno),nx,ny,goodno)

goodmap=good_map-median_cube

;; FOR mi=0,mapsize[3]-1 DO BEGIN

;;  good_means[mi]=mean(map[*,*,mi],/double,/nan)

;; ENDFOR 

;;;;;;;;;; MUST REWRITE THE NEXT BIT TO JUST SUBTRACT MOVING AVERAGE
;;;;;;;;;; FROM THE 'GOODMAP' ACCORDING TO CHANNEL

;;;;;;;;;; FINAL BIT IS TO AVERAGE ALONG THE SAME CHUNKS AS BEFORE


;;;;;;;;;; WHEN THIS IS A MODULE OF THE CODE WANT TO OUTPUT COMPARISON
;;;;;;;;;; BETWEEN 1ST AVERAGING AND 2ND AVERAGING

FOR ii=0, bins-1 DO BEGIN

chunksize=(end_indices[ii]-start_indices[ii])+1
chunk=goodmap[*,*,start_indices[ii]:end_indices[ii]]
median_chunk=good_medians[start_indices[ii]:end_indices[ii]]
smoothed_m=ts_smooth(median_chunk,5)
3d_median=rebin(reform(smoothed_m,1,1,chunksize),mapsize[1],mapsize[2],chunksize)
new_chunk=chunk-3d_median
new_map[*,*,start_indices[ii]:end_indices[ii]]=new_chunk
sm_map[*,*,start_indices[ii]:end_indices[ii]]=3d_median

ENDFOR


FOR gi=0,mapsize[3]-1 DO BEGIN

   new_map_i=new_map[*,*,gi]
   new_medians[gi]=median(new_map_i[good_intersect])

ENDFOR



;!PLOT.MULTI=[0,3,2]

set_plot, 'ps'
device,file=outfile
device,/color

plot,fgood_frac,good_medians,psym=1,xtitle='Frequency channel',ytitle='Median value'
name='Channels with 95% or more good pixels'
oplot,fgood_frac,new_medians,psym=1,col=1

for pp=0,nostart-1 do oplot,fltarr(2)+startchans[pp],!y.crange,color=2
for pp=0,noend-1 do oplot,fltarr(2)+endchans[pp],!y.crange,color=3
;plots,60,0.08,psym=1,col=1
;xyouts,80,0.08,name,col=1

device,/close
set_plot,'x'


sxdelpar,header,'ORIGINAL'
sxdelpar,header,'COMMENT'
sxaddpar,header,'COMMENT','Each channel is comprised of identical pixels which each have the value of the moving average of the median over 5 channels'

mwrfits,sm_map,smoothout,header



END
