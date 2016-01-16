pro display_offsets,a_b_in

data=fltarr(7,149)
openr,lun,a_b_in,/get_lun
readf,lun,data
close,/all

pixels=data[0,*]
b=data[2,*]

array=fltarr(192)

for i=0,148 do begin

array[pixels[i]]=b[i]

endfor

array[where(array eq 0)]=!healpix.bad_value

mollview,array,/nest,min=-1,max=2


end



