function slope,x_1,x_2,y_1,y_2

result = (y_2-y_1)/(x_2-x_1)
return,result

end



function theilsen,x,y

sx=n_elements(x)
sy=n_elements(y)

if sx ne sy then begin
print,'x and y need to have the same length'
   return, -1
endif

n=float(sx)

print,n,' good pixels.'

;print, 'Beginning bootstrapping...'

;slope_bootstrap=make_array(niter,/nozero)
;intercept_bootstrap=make_array(niter,/nozero)
combno=double(n*(n-1)*0.5)
print,combno
slopes=make_array(combno,/nozero)
i1=fix(n*randomu(seed,combno))
i2=fix(n*randomu(seed,combno))
nn=where(abs(x[i2]-x[i1]) ne 0,count_zero)
i1=i1[nn]
i2=i2[nn]
slopes = slope(x[i1],x[i2],y[i1],y[i2])
sorted_slopes=slopes[sort(slopes)]
   
median_slope=median(sorted_slopes)

intercepts=make_array(n,/nozero)
upper=intercepts
lower=intercepts

   for c=0,n-1 do begin
      intercepts[c]=y[c]-median_slope*x[c]
      intercept=median(intercepts)
   endfor

z=double(1.96*sqrt(n*(n-1)*(2*n+5)/18.))
;print,z
ru=long64(round((combno+z)/2+1.))
rl=long64(round((combno-z)/2.))
;print,rl,ru
slope_l=sorted_slopes[rl]
slope_u=sorted_slopes[ru]

for c=0,n-1 do begin
      upper[c]=y[c]-slope_u*x[c]
      off_up=median(upper)
   endfor

for c=0,n-1 do begin
      lower[c]=y[c]-slope_l*x[c]
      off_low=median(lower)
   endfor


  ts=[median_slope,intercept,slope_l,slope_u,off_low,off_up]

return,ts

end



pro offsetsv3, gmims_in,haslam_in

print,'Reading maps and preparing data...'
read_fits_map,gmims_in,gmims,ordering='nested'
read_fits_map,haslam_in,haslam,ordering='nested'

haslam=haslam[*,0]
;haslam*=0.001
haslam-=2.725

gmims_div=reform(gmims,16384,192)
haslam_div=reform(haslam,16384,192)

print,'...done!'

for i=0,192-1 do begin

   print,'Working on pixel',i
   
   gmims_badval=where(gmims_div[*,i] eq !healpix.bad_value)
   gmims_div[gmims_badval,i]=!values.f_nan

   haslam_badval=where(haslam_div[*,i] eq !healpix.bad_value)
   haslam_div[haslam_badval,i]=!values.f_nan

   gmims_good=where(finite(gmims_div[*,i]),count1)
   print,count1,' good gmims pixels'
   haslam_good=where(finite(haslam_div[*,i]),count2)
   print,count2,' good haslam pixels'
   
   if (count1 gt 0) and (count2 gt 0) then begin

      good_intersect=setintersection(gmims_good,haslam_good)

      print,n_elements(good_intersect),' good pixels in both maps'
      
      print, 'Pixel is good - beginning TS regression'
      a_b=theilsen(haslam_div[good_intersect,i],gmims_div[good_intersect,i])

      openw,lun,'offsets_slopes_masked_no_boot_v3.txt',/append,/get_lun
      printf,lun,[i,a_b]
      close,lun
      free_lun,lun
      
   endif else begin

      print,'No good data in this subpixel'
      openw,lun,'offsets_slopes_masked_no_boot_v3.txt',/append,/get_lun
      empty=make_array(6,value=!values.f_nan)
      printf,lun,i,empty
      close,lun
      free_lun,lun
  endelse

endfor

end
