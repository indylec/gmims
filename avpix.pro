pro avpix

read_fits_map,'gmims_dipole_corrected_spur.fits',map,ordering='nested'
badval=where(map eq !healpix.bad_value)
map[badval]=!values.f_nan

map_div=reform(map,16384,192)


average=fltarr(192)

for i=0,191 do begin

   average[i]=mean(map_div[*,i],/nan)

endfor

mollview,average,/asinh,/nest,png='coarse_corrected_map.png'

end
