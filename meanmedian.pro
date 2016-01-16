pro meanmedian, infile,averagefile,outfile


map=mrdfits(infile,0,header)



avgtable=mrdfits(averagefile,2,/use_colnum)
subtable=mrdfits(infile,1)

tablesize=size(avgtable)

bins=tablesize[1]

fgood_frac=subtable.(1)
good_freq=subtable.(2)

startchans=fltarr(bins)
endchans=fltarr(bins)

nostart=n_elements(startchans)
noend=n_elements(endchans)

for ii=0,bins-1 do begin
startchans[ii]=avgtable[ii].c2
endchans[ii]=avgtable[ii].c3

endfor

mapsize=size(map)

good_medians=dblarr(mapsize[3])

good_means=good_medians


good_intersect=where(finite(map[*,*,0]),mc,complement=bad_union)


FOR ii=0,mapsize[3]-1 DO BEGIN

   good_indices=where(finite(map[*,*,ii]),tmc,complement=bad_indices)
   good_intersect=setintersection(good_indices,good_intersect)
   
ENDFOR


FOR gi=0,mapsize[3]-1 DO BEGIN

   good_map_i=map[*,*,gi]
   good_medians[gi]=median(good_map_i[good_intersect])

ENDFOR

FOR mi=0,mapsize[3]-1 DO BEGIN

 good_means[mi]=mean(map[*,*,mi],/double,/nan)

ENDFOR



;!PLOT.MULTI=[0,3,2]

set_plot, 'ps'
device,file=outfile
device,/color

plot,fgood_frac,good_medians,psym=1,xtitle='Frequency channel',ytitle='Median value'
name='Channels with 95% or more good pixels'
;oplot,fgood_frac,good_means,psym=1,col=1

for pp=0,nostart-1 do oplot,fltarr(2)+startchans[pp],!y.crange,color=2
for pp=0,noend-1 do oplot,fltarr(2)+endchans[pp],!y.crange,color=3
;plots,60,0.08,psym=1,col=1
;xyouts,80,0.08,name,col=1

device,/close
set_plot,'x'



END
