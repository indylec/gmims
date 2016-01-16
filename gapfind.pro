FUNCTION GAPFIND, map, good_ratio

mapsize = SIZE(map)

nx = mapsize[1]
ny = mapsize[2]
nchan = mapsize[3]

npix = nx * ny

ngood = LONARR(nchan)

FOR ichan = 0, nchan-1 DO BEGIN 
   void = WHERE(finite(map[*,*,ichan]), ngood_count)
   ngood[ichan]=ngood_count
ENDFOR

fgood= FLOAT(ngood)/npix
fgood_frac= WHERE(fgood GE good_ratio,goodno)


FOR ii=1,goodno-1 DO BEGIN

   IF fgood_frac[ii]-fgood_frac[ii-1] GE 200 THEN BEGIN 

      startgap=ii-1
      endgap=ii

   ENDIF

ENDFOR

good_chunk_1=fgood_frac[0:startgap]
good_chunk_2=fgood_frac[endgap:goodno-1]


   
 RETURN,
END
