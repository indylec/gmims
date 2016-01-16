PRO tableread,in,out


fxbopen,unit,in,2

fxbreadm,unit,[1,2,3,4],bin,startchan,endchan,mfr

s=size(bin)

openw,10,out

printf,10,'Bin & ',' Start Channel &',' End Channel &',' Mean frequency \\'

for ii=0,8 do begin
 printf,10,bin[ii],' & ',startchan[ii],' & ',endchan[ii],' & ',mfr[ii],'\\'

endfor
close,10

END
