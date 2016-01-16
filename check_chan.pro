PRO check_chan, file

map=readfits(file,header)

mapsize = SIZE(map)

nx = mapsize[1]
ny = mapsize[2]
nchan = mapsize[3]

npix = nx * ny

ngood = LONARR(nchan)

f0=sxpar(header,'CRVAL3')
df=sxpar(header,'CDELT3')
refpix=sxpar(header, 'CRPIX3')


FOR ichan = 0, nchan-1 DO BEGIN 
   void = WHERE(finite(map[*,*,ichan]), ngood_count)
   ngood[ichan]=ngood_count
ENDFOR


fgood = FLOAT(ngood)/ npix
fgood80=where(fgood GE 0.8, good80)
fgood90=where(fgood GE 0.9, good90)
fgood95=where(fgood GE 0.95,good95)
fgood98=where(fgood GE 0.98,good98)

fg95=fgood[fgood95]


print,FORMAT='( "80% of good points: ",I0,"/",I0," channels")',good80,nchan
print,FORMAT='( "90% of good points: ",I0,"/",I0," channels")',good90,nchan
print,FORMAT='( "95% of good points: ",I0,"/",I0," channels")',good95,nchan
print,FORMAT='( "98% of good points: ",I0,"/",I0," channels")',good98,nchan

filename='good_channels_plot_'+strmid(file,5,10)+'_2.ps'
set_plot, 'ps'
device,file=filename
device,/color
Plot, fgood, PSYM = 1,subtitle=strmid(file,5,10),xtitle='Frequency channel',ytitle='Fraction of good pixels',xstyle=8
name='Channels with 95% or more good pixels'
axis,xaxis=1,xrange=(!X.CRANGE-refpix)*df+f0,xstyle=1,xtitle= 'Frequency (Hz)',xcharsize=0.9

oplot,fgood95,fg95,psym=1,col=1
plots,60,0.405,psym=1,col=1
xyouts,80,0.4,name,col=1

device, /close
set_plot,'x'

END

