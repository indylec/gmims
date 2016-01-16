PRO TTPLOT,file1,file2,outfile,subpixel

read_fits_map,file1,map1
read_fits_map,file2,map2

map1_div=reform(map1,192,16384)
map2_div=reform(map2,192,16384)

subpix1=reform(map1_div[subpix,*],16384)
subpix2=reform(map2_div[subpix,*],16384)

coeff1=robust_line(subpix1,subpix2,/bisect)
coeff2=robust_line(subpix2,subpix1,/bisect)

set_plot,'ps'
device,file=outfile
device,/color


plot,subpix1,subpix2,psym=1,xrange=[0,100000],xtitle='haslam (K)',ytitle='gmims (K)'



END
