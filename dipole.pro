pro dipole,a_b_in

data=fltarr(7,192)
openr,lun,a_b_in,/get_lun
readf,lun,data
close,/all

good_pixels=where(finite(data[1,*]),good_count)
data=data[*,good_pixels]

npix=nside2npix(512)

pixels=data[0,*]
slope=data[1,*]
b=data[2,*]
sloperr=data[4,*]-data[3,*]
offseterr=data[6,*]-data[5,*]
zeroerr=where(offseterr ne 0,count)
print,count
offseterr=offseterr[zeroerr]
b=b[zeroerr]
b_red=b/offseterr
pixels=pixels[zeroerr]

;help,offseterr
;print,offseterr
;print,b[0]

ipix=lindgen(npix)
pix2ang_nest,512,ipix,theta,phi

t_matrix=fltarr(npix,4)

t0 = replicate(1,npix)
t1 = cos(phi)*sin(theta)
t2 = sin(phi)*sin(theta)
t3 = cos(theta)
;for i=0,npix-1 do begin
;
;t_matrix[*,i]=[1,cos(phi[i])*sin(theta[i]),sin(theta[i])*sin(phi[i]),cos(theta[i])]
;
;endfor

;print,size(t_matrix)

;t0=t_matrix[0,*]
;t1=t_matrix[1,*]
;t2=t_matrix[2,*]
;t3=t_matrix[3,*]

t0=reform(t0,16384,192)
t1=reform(t1,16384,192)
t2=reform(t2,16384,192)
t3=reform(t3,16384,192)

t0_i=fltarr(count)
t1_i=fltarr(count)
t2_i=fltarr(count)
t3_i=fltarr(count)

for ii=0,count-1 do begin

   t0_i[ii]=mean(t0[*,pixels[ii]])/offseterr[ii]
   t1_i[ii]=mean(t1[*,pixels[ii]])/offseterr[ii]
   t2_i[ii]=mean(t2[*,pixels[ii]])/offseterr[ii]
   t3_i[ii]=mean(t3[*,pixels[ii]])/offseterr[ii]
endfor

t0_i=reform(t0_i,1,count)
t1_i=reform(t1_i,1,count)
t2_i=reform(t2_i,1,count)
t3_i=reform(t3_i,1,count)

help,t0_i
help,t1_i
help,t2_i
help,t3_i

;slope_col=rebin(slope,4,155)
a_matrix=[t0_i,t1_i,t2_i,t3_i]
b_red=reform(b_red,1,count)
;a_matrix_1=-slope_col*a_matrix
a_matrix_t=transpose(a_matrix)
help,a_matrix
help,a_matrix_t
help,b_red
cov=la_invert(a_matrix_t##a_matrix)
help,cov
err=diag_matrix(cov)
help,err
err=sqrt(err)

;print,a_matrix
;print,size(a_matrix_1)
;a_matrix=reform(a_matrix,4,192)
;a=[a_matrix_1,a_matrix]

;print,size(a)

x=la_least_squares(a_matrix,b_red,residual=resid, status=status,method=2)


print,"The x-vector is: ",x

print,"The uncertainties are:",err

chisq=total((b_red-(x[0]*t0_i+x[1]*t1_i+x[2]*t2_i+x[3]*t3_i))^2)

print,"residual= ",resid

print,"Chi squared = ",chisq

print,"DoF = ",count-4

print,"reduced chi-squared= ",chisq/(count-4)

print,b[0]

read_fits_map,'../fits_files/healpix/gmims_final_masked_nospur.fits',gmims,hdr,exthdr
fit =  x[0]*t0 + x[1]*t1 + x[2]*t2 + x[3]*t3
fit = reform(fit,npix)
mollview, fit, /nest, min=-1.5, max=1.5,units='K',title='Dipole fit';,png='masked_fit_v2.png'
offsets = fltarr(192)
offsets[*] = !values.f_nan
offsets[pixels] = b
mollview, offsets, /nest, min=-1.5, max=1.5,units='K',title='Offset map';,png='masked_offset_v2.png'
;for i=0,npix-1 do begin
;
;   gmims[i]=gmims[i]-(x[0]+cos(phi[i])*sin(theta[i])*x[1]+sin(theta[i])*sin(phi[i])*x[2]+cos(theta[i])*x[3])
;
;endfor

mollview,gmims-fit,/nest,min = -1, max = 2,units='K',title='Corrected GMIMS map';,png='masked_corrected_v2.png'

;write_fits_map,'gmims_masked_corrected_v2.fits',gmims-fit,exthdr,ordering=nested

end
