FUNCTION NSIDE2RES, nside
ns = nside
nres = 0
WHILE ns GT 1 DO BEGIN
    IF ns LT 1 THEN BEGIN
        PRINT, 'NSIDE2RES: Illegal NSIDE:',nside
        nres = -1
        BREAK
    ENDIF
    ns = ns / 2
    nres = nres + 1
ENDWHILE
RETURN, nres
END

PRO getcoord, crval1,crpix1,cdelt1,crval2,crpix2,cdelt2,head
;
; Gets coordinate info from old-fashioned fits map
;
rdcard, head, 'CRVAL1  ', crval1
rdcard, head, 'CRPIX1  ', crpix1
rdcard, head, 'CDELT1  ', cdelt1
rdcard, head, 'CRVAL2  ', crval2
rdcard, head, 'CRPIX2  ', crpix2
rdcard, head, 'CDELT2  ', cdelt2

END

PRO rdcard, head, name, var
; Reads data from a fits header card
ip = WHERE(STRCMP(head, name, 8))
ncard = N_ELEMENTS(ip)

CASE ncard OF
    1: READS, STRMID(head[ip],11), var
    0: PRINT, 'RDCARD: '+name+' card not found' 
    ELSE: BEGIN
        PRINT, 'RDCARD: multiple ',+name+' cards found --- Using first'
        READS, STRMID(head[ip[0]],11), var
    END
ENDCASE
END


PRO regrid_GMIMS, nside, channel
;
; Regrids the DRAO polarization survey into HEALPix format
;

prog = 'regrid_GMIMS'
version='0.1'

npix = NSIDE2NPIX(nside)
; Read in raw data 2-D arrays

dir = '/scratch/cosmo4_1/jpl/lowfreq/'
READ_FITS_MAP,'TP_000-360_cstv2.fits',tt,headt

tt=tt[*,*,channel]

; Pad at ends for better interpolation
npad = 10
s = SIZE(tt)
last = s[1] - 1
tt = [tt[last-npad:last-1,*],tt,tt[1:npad,*]]
ttlen = s[1]

; Define HEALPix pixel directions
PIX2VEC_NEST, nside, lindgen(npix), rdvecs

VEC2ANG, rdvecs, dec, ra, /astro

; Interpolate values at HEALPix pixels from the 2-D grids:

getcoord, crval1,crpix1,cdelt1,crval2,crpix2,cdelt2,headt
rapix =  (ra  - crval1)/cdelt1 + crpix1 - 1 + npad
IF (MIN(rapix) LT npad) THEN rapix[WHERE(rapix LT npad)] += ttlen - 1
decpix = (dec - crval2)/cdelt2 + crpix2 - 1

temptt = INTERPOLATE(tt,rapix,decpix,cubic=-0.5,missing=!values.F_NAN)
tt = 0
; Adjust zero-level:

;temptt = temptt - 3000.



; Write out

header = [ $
"TTYPE1  = 'TEMPERATURE '       /label for field 1                 ", $
"TFORM1  = 'E       '           /data format of field: 4-byte REAL ", $
"TUNIT1  = 'K      '           /physical unit of field 1          ", $
"HISTORY  "+prog+": Version "+version, $
"HISTORY  "+prog+": Map is channel"+string(channel)+"of binned GMIMS map" $
]

outfil = 'TP_000-360_'+strtrim(string(channel),1)+'_healpix_res'+STRING(nside2res(nside),FORMAT='(I1)')+'.fits'
WRITE_FITS_MAP, outfil,temptt,header, Coordsys='C',Ordering='nested'

QUIT:
END
