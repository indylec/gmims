pro mediancalc, in,igood,ngood,xstart,xend,ystart,yend,good_medians




map_box=in[xstart:xend,ystart:yend,*]



good_map=map_box[*,*,igood]

good_medians=dblarr(ngood)

good_intersect=where(finite(good_map[*,*,0]),mc,complement=bad_union)


FOR ii=0,ngood-1 DO BEGIN

   good_indices=where(finite(good_map[*,*,ii]),tmc,complement=bad_indices)
   good_intersect=setintersection(good_indices,good_intersect)
   
ENDFOR


FOR gi=0,ngood-1 DO BEGIN

   good_map_i=good_map[*,*,gi]
   good_medians[gi]=median(good_map_i[good_intersect])

ENDFOR

;out=good_medians

END


pro median_multi,infile,outfile,x1,x2,x3,x4,x5,x6,y1,y2,y3,y4

map=mrdfits(infile,0,header)

mapsize=size(map)


nx = mapsize[1]
ny = mapsize[2]
nchans=mapsize[3]

npix = nx * ny

num_good = LONARR(nchans)

FOR ichan = 0, nchans-1 DO BEGIN 
   void = WHERE(finite(map[*,*,ichan]), ngood_count)
   num_good[ichan]=ngood_count
ENDFOR

fgood = FLOAT(num_good)/ npix

fgood95=where(fgood GE 0.95,good95)
fg95=fgood[fgood95]

plotarray=dblarr(good95,6)

mediancalc,map,fgood95,good95,x1,x2,y3,y4,plotarray[*,0]

mediancalc,map,fgood95,good95,x3,x4,y3,y4,plotarray[*,1]

mediancalc,map,fgood95,good95,x5,x6,y3,y4,plotarray[*,2]

mediancalc,map,fgood95,good95,x1,x2,y1,y2,plotarray[*,3]

mediancalc,map,fgood95,good95,x3,x4,y1,y2,plotarray[*,4]

mediancalc,map,fgood95,good95,x5,x6,y1,y2,plotarray[*,5]





set_plot, 'ps'
device,file=outfile
device,/color
!P.MULTI=[0,3,2]

plot,fgood95,plotarray[*,0],psym=1,title=strmid(infile,5,10)+' corners at '+strtrim(string(x1),2)+','+strtrim(string(y3),2)+' and '+strtrim(string(x2),2)+','+strtrim(string(y4),2),xtitle='Frequency channel',ytitle='Median value',yrange=[-0.4,0.1]
;name='Channels with 95% or more good pixels'
;oplot,fgood95,good_medians,psym=1,col=1
;plots,60,0.08,psym=1,col=1
;xyouts,80,0.08,name,col=1

plot,fgood95,plotarray[*,1],psym=1,title=strmid(infile,5,10)+' corners at '+strtrim(string(x3),2)+','+strtrim(string(y3),2)+' and '+strtrim(string(x4),2)+','+strtrim(string(y4),2),xtitle='Frequency channel',ytitle='Median value',yrange=[-0.4,0.1]

plot,fgood95,plotarray[*,2],psym=1,title=strmid(infile,5,10)+' corners at '+strtrim(string(x5),2)+','+strtrim(string(y3),2)+' and '+strtrim(string(x6),2)+','+strtrim(string(y4),2),xtitle='Frequency channel',ytitle='Median value',yrange=[-0.4,0.1]

plot,fgood95,plotarray[*,3],psym=1,title=strmid(infile,5,10)+' corners at '+strtrim(string(x1),2)+','+strtrim(string(y1),2)+' and '+strtrim(string(x2),2)+','+strtrim(string(y2),2),xtitle='Frequency channel',ytitle='Median value',yrange=[-0.4,0.1]

plot,fgood95,plotarray[*,4],psym=1,title=strmid(infile,5,10)+' corners at '+strtrim(string(x3),2)+','+strtrim(string(y1),2)+' and '+strtrim(string(x4),2)+','+strtrim(string(y2),2),xtitle='Frequency channel',ytitle='Median value',yrange=[-0.4,0.1]

plot,fgood95,plotarray[*,5],psym=1,title=strmid(infile,5,10)+' corners at '+strtrim(string(x5),2)+','+strtrim(string(y1),2)+' and '+strtrim(string(x6),2)+','+strtrim(string(y2),2),xtitle='Frequency channel',ytitle='Median value',yrange=[-0.4,0.1]

device,/close
set_plot,'x'



END
