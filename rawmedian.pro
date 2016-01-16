
pro rawmedian,rawfile,binfile,outplot


map=mrdfits(rawfile,0,header)

f0=sxpar(header,'CRVAL3')
df=sxpar(header,'CDELT3')
refpix=sxpar(header, 'CRPIX3')

map_info=mrdfits(binfile,1)
xtable_arr=mrdfits(binfile,4)

startchans=map_info.(3)
endchans=map_info.(4)

bins=size(startchans)

fgood_frac=map_info.(1)

good_map=map[*,*,fgood_frac]

goodno=size(fgood_frac,/n_elements)
print,goodno
medians=dblarr(goodno)


raw_good_intersect=where(finite(good_map[*,*,0]))

FOR ii=0,goodno-1 DO BEGIN

   raw_good_indices=where(finite(good_map[*,*,ii]))
   raw_good_intersect=setintersection(raw_good_indices,raw_good_intersect)
   
ENDFOR

FOR gi=0,goodno-1 DO BEGIN

   raw_map_i=good_map[*,*,gi]
   medians[gi]=median(raw_map_i[raw_good_intersect])

ENDFOR


set_plot,'ps'
device,file=outplot
device,/color

plot, fgood_frac,medians,psym=1,subtitle='Median values per channel for TP_000-072',xtitle='frequency channel',ytitle='Median value (K)',xstyle=8,ymargin=[6,4]
axis,xaxis=1,xrange=(!X.CRANGE-refpix)*df+f0,xstyle=1,xtitle= 'Frequency (Hz)',xcharsize=0.9
for pp=0,bins[1]-1 do oplot,fltarr(2)+xtable_arr[pp].stchan,!y.crange,color=1,linestyle=2
for pp=0,bins[1]-1 do oplot,fltarr(2)+xtable_arr[pp].echan,!y.crange,color=1,linestyle=2

device,/close

set_plot,'x'


end
