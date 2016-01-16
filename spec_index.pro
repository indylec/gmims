pro spec_index

gmims_freq=1310025600.0
haslam_freq=408000000.0


denominator=alog(gmims_freq/haslam_freq)

read_fits_map,'gmims_dipole_corrected.fits',gmims,ordering=nested

s=size(gmims)

gmims_neg=where(gmims LT 0)
gmims[gmims_neg]=!healpix.bad_value
gmims_bad=where(gmims eq !healpix.bad_value)
gmims[gmims_bad]=!values.f_nan
;gmims+=0.289447

read_fits_map,'./haslam/final_haslam_masked.fits',haslam,ordering=nested
haslam_bad=where(haslam eq !healpix.bad_value)
haslam[haslam_bad]=!values.f_nan

betamap=make_array(s[1])

for i=0,s[1]-1 do begin

betamap[i]=(alog(gmims[i]/haslam[i]))/denominator

endfor

mollview,betamap,/nest,/asinh,max=-2,png='spec_index_map_dipolcorr.png'

end
