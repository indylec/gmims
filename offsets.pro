function slope,x_1,x_2,y_1,y_2

result = (1-2*(x_1 gt x_2))*((y_2-y_1)/abs(x_2-x_1))
return,result

end



function theilsen,x,y,nsamples,niter

sx=n_elements(x)
sy=n_elements(y)

if sx ne sy then begin
   print,'x and y need to have the same length'
   return, -1
endif

n=sx

print,n,' good pixels.'

print, 'Beginning bootstrapping...'

slope_bootstrap=make_array(niter,/nozero)
intercept_bootstrap=make_array(niter,/nozero)

for ii=0,niter-1 do begin

   print, 'bootstrap loop', ii

   if n LT 100 then begin
      ix=sort(x)
      combno=(n*(n-1)*0.5)
      slopes=make_array(combno,/nozero)
      pair=combigen(n,2)

      for ind=0,combno-1 do begin
         i=ix[pair[ind,0]]
         j=ix[pair[ind,1]]
         slopes[ind] = slope(x[i],x[j],y[i],y[j])
         endfor
      
   endif else begin
      i1=fix(n*randomu(seed,nsamples))
      i2=fix(n*randomu(seed,nsamples))
      nn=where(abs(x[i2]-x[i1]) ne 0,count_zero)
      ;print,count_zero
      i1=i1[nn]
      i2=i2[nn]
      slopes = slope(x[i1],x[i2],y[i1],y[i2])
   endelse
   loop_slope=median(slopes)

   intercepts=make_array(n,/nozero)

   for c=0,n-1 do begin
      intercepts[c]=y[c]-loop_slope*x[c]
      intercept=median(intercepts)
   endfor

   slope_bootstrap[ii]=loop_slope
   intercept_bootstrap[ii]=intercept

endfor

print,'...bootstrapping done for this pixel!'
  

final_slope=median(slope_bootstrap)
c95=prank(slope_bootstrap,[5,95])
final_intercept=median(intercept_bootstrap)
c95i=prank(intercept_bootstrap,[5,95])



ts=[final_slope,final_intercept,c95[0],c95[1],c95i[0],c95i[1]]

return,ts

end



pro offsets, gmims_in,haslam_in

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
      a_b=theilsen(haslam_div[good_intersect,i],gmims_div[good_intersect,i],1e7,250)

      openw,lun,'offsets_slopes_masked_250_v2.txt',/append,/get_lun
      printf,lun,[i,a_b]
      close,lun
      free_lun,lun
      
   endif else begin

      print,'No good data in this subpixel'
      openw,lun,'offsets_slopes_masked_250_v2.txt',/append,/get_lun
      empty=make_array(6,value=!values.f_nan)
      printf,lun,i,empty
      close,lun
      free_lun,lun
      endelse

endfor

end

