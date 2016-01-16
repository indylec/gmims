pro offsetread,pixel,slope,offset

pixel=fltarr(192)
slope=pixel
offset=pixel

openr,lun,'offsets_slopes.txt',/get_lun

for i=0,191 do begin

   readf,lun,a,b,c
   pixel[i]=a
   slope[i]=b
   offset[i]=c
   readf,lun

endfor

close,lun


end
