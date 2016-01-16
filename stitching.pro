pro stitching, infile,outfile


OPENR, lun, infile, /GET_LUN
; Read one line at a time, saving the result into array
files = ''
line = ''
WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  files = [files, line] & $
ENDWHILE
; Close the file and free the file unit
FREE_LUN, lun

;print, files

filesize=size(files)
nofiles=filesize[1]

files=files[1:nofiles-1]

print,FORMAT='("Reading in ",A0)',files[0]
map0=readfits(files[0],header)
print,FORMAT='("Reading in ",A0)',files[1]
map1=readfits(files[1])
print,FORMAT='("Reading in ",A0)',files[2]
map2=readfits(files[2])
print,FORMAT='("Reading in ",A0)',files[3]
map3=readfits(files[3])
print,FORMAT='("Reading in ",A0)',files[4]
map4=readfits(files[4])

map1=map1[0:359,*,*]
map2=map2[0:359,*,*]
map3=map3[0:359,*,*]
map4=map4[1:359,*,*]

final_map=[map4,map3,map2,map1,map0]

SXADDPAR,header,'CRVAL1',0.0000,'  Reference pixel value'
SXADDPAR,header,'CRPIX1',0.00,'  Reference pixel'
SXDELPAR,header,'ORIGINAL'
SXADDPAR,header,'COMMENT','This map has been stitched together from  5 reduced blocks'

writefits,outfile,final_map,header



END
