pro medianplot, infile, x1,x2,y1,y2,outfile


map=mrdfits(infile,0,header)

map_box=map[x1:x2,y1:y2,*]

mapsize=size(map)


nx = mapsize[1]
ny = mapsize[2]
nchans=mapsize[3]

npix = nx * ny

ngood = LONARR(nchans)

FOR ichan = 0, nchans-1 DO BEGIN 
   void = WHERE(finite(map[*,*,ichan]), ngood_count)
   ngood[ichan]=ngood_count
ENDFOR

fgood = FLOAT(ngood)/ npix

print, fgood[500:510]

fgood95=where(fgood GE 0.95,good95)
fg95=fgood[fgood95]

good_map=map_box[*,*,fgood95]

print, good95

good_medians=dblarr(good95)

good_intersect=where(finite(good_map[*,*,0]),mc,complement=bad_union)


FOR ii=0,good95-1 DO BEGIN

   good_indices=where(finite(good_map[*,*,ii]),tmc,complement=bad_indices)
   good_intersect=setintersection(good_indices,good_intersect)
   
ENDFOR


FOR gi=0,good95-1 DO BEGIN

   good_map_i=good_map[*,*,gi]
   good_medians[gi]=median(good_map_i[good_intersect])

ENDFOR




;!PLOT.MULTI=[0,3,2]

set_plot, 'ps'
device,file=outfile
device,/color

plot,fgood95,good_medians,psym=1,title=strmid(infile,5,10)+' corners at '+strtrim(string(x1),2)+','+strtrim(string(y1),2)+' and '+strtrim(string(x2),2)+','+strtrim(string(y2),2),xtitle='Frequency channel',ytitle='Median value'
name='Channels with 95% or more good pixels'
;oplot,fgood95,good_medians,psym=1,col=1
;plots,60,0.08,psym=1,col=1
;xyouts,80,0.08,name,col=1

device,/close
set_plot,'x'



END
