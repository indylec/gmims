pro monopole,a_b_in

data=fltarr(7,192)
openr,lun,a_b_in,/get_lun
readf,lun,data
close,/all

npix=nside2npix(512)

a=data[1,*]
b=data[2,*]
ones=fltarr(192)+1.0
ones=rebin(ones,1,192)

a=[-a,ones]

print,a
print,b

x=la_least_squares(a,b)


print,"The x-vector is: ",x


end
