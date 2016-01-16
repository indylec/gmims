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

FUNCTION getrdvec, headt, stt 
;
; Creates a set of 3D position vectors corresponding to the pixels in
; a 2-D sky map with FITS header headt, described by SIZE array stt.
crval1 = 0d0
crval2 = 0d0
getcoord, crval1,crpix1,cdelt1,crval2,crpix2,cdelt2,headt
ra  = (LINDGEN(stt[1]) + 1 - crpix1)*cdelt1 + crval1
dec = (LINDGEN(stt[2]) + 1 - crpix2)*cdelt2 + crval2

; Create 2d theta and phi arrays:
ra  = REBIN(ra,stt[1],stt[2])
dec = REFORM(dec,1,stt[2])
dec = REBIN(dec,stt[1],stt[2])

ANG2VEC, dec, ra, rdvec, /ASTRO

RETURN, rdvec
END

PRO map2hp, nside, lbvec, tt, ttv0, htt, tvar
;
; Bins a map array into HEALPix. Inputs:
;
; nside : size parameter of output HEALpix map
; lbvec : array 3-D position vectors for each map pixel, 
; tt    : input 2D map
; ttv0  : variance per input pixel.
; htt   : output Healpix map
; tvar  : output variance map

; Find corresponding pixels:
VEC2PIX_NEST, nside, lbvec, rawpix

good = WHERE(FINITE(tt))
goodpix = rawpix[good]
goodtt = tt[good]
; Find original pixels that correspond to each Healpix pixel:
count = HISTOGRAM(goodpix, OMIN=om, REVERSE_INDICES=ri)

; Reject cases with < 2 good original pixels)
mask = WHERE(count gt 1) 
npix = NSIDE2NPIX(nside)
pix = LINDGEN(npix) + om
pix = pix[mask]
PRINT, 'Min, Max pixel worth binning:', MINMAX(pix)
PRINT, 'Number of pixels worth binning:', N_ELEMENTS(pix)

count = count[mask]
PRINT, 'Maximum number of input pixels per output pixel:', MAX(count)
PRINT, 'Mean number', MEAN(count)

i1 = ri[mask]
i2 = ri[mask+1] - 1
htt = REPLICATE(0./0.,npix)  ; Default is NaN = blank
tvar = fltarr(npix)
FOR j=0L,N_ELEMENTS(count)-1 DO BEGIN
    htt[pix[j]]  = MEAN(goodtt[ri[i1[j]:i2[j]]])
    tvar[pix[j]] = count[j]
ENDFOR

END

PRO regrd_drao, nside
;
; Regrids the DRAO polarization survey into HEALPix format
;

prog = 'regrd_drao'
version = '0.5'
d2r = !dpi/180d0
ttv0  = 20.^2  ; Variance in Stockert data (mK^2)
polv0 = 12.^2  ; Variance in DRAO data (mK^2)

; Euler matrix for Galactic to B1950 coords (from NOD2 header):
matrix =     [[-0.0669887394152D0, 0.492728466075D0, -0.867600811151D0], $
              [-0.872755765852D0, -0.45034695802D0,  -0.188374601723D0], $
              [-0.483538914632D0,  0.744584633283D0,  0.460199784784D0]]  

npix = NSIDE2NPIX(nside)
; Read in raw data 2-D arrays

dir = '/scratch/cosmo4_1/jpl/lowfreq/'
READ_FITS_MAP, dir+'Stockert_bin.fits',tt,headt

; Set blanks to NAN
bad = WHERE(tt EQ -32000)
tt[bad] = 0./0.

; Pad at ends for better interpolation
npad = 10
s = SIZE(tt)
last = s[1] - 1
tt = [tt[last-npad:last-1,*],tt,tt[1:npad,*]]
ttlen = s[1]

READ_FITS_MAP, dir+'26m_survey_Q.fits',q,headq
READ_FITS_MAP, dir+'26m_survey_U.fits',u,headu  ; NB IAU sign convention!

s = SIZE(q)
rdvecs = getrdvec(headq,s)

last = s[1] - 1
qqlen = s[1]
height = s[2]

; Replace north pole strip by its mean:
q[*,s[2]-1] = mean(q[*,s[2]-1])
u[*,s[2]-1] = mean(u[*,s[2]-1])

qu = [[[q]], [[u]]]
q = 0
u = 0

; Interpolate into Healpix pixels:

; tidy up: match at RA = 0 = 360
seam = 0.5*(qu[0,*,*] + qu[last,*,*])
qu = [qu[last-npad:last-1,*,*],seam,qu[1:last-1,*,*],seam,qu[1:npad,*,*]]

; Define HEALPix pixel directions
PIX2VEC_NEST, nside, lindgen(npix), lbvecs

; Stockert Survey is in B1950 coords

rdvecs = ROTATE_COORD(lbvecs, Euler = matrix)
VEC2ANG, rdvecs, dec, ra, /astro

; Interpolate values at HEALPix pixels from the 2-D grids:

getcoord, crval1,crpix1,cdelt1,crval2,crpix2,cdelt2,headt
rapix =  (ra  - crval1)/cdelt1 + crpix1 - 1 + npad
IF (MIN(rapix) LT npad) THEN rapix[WHERE(rapix LT npad)] += ttlen - 1
decpix = (dec - crval2)/cdelt2 + crpix2 - 1

temptt = INTERPOLATE(tt,rapix,decpix,cubic=-0.5,missing=!values.F_NAN)
tt = 0
; Adjust zero-level:

temptt = temptt - 3000.

; DRAO survey is in J2000 coords:

rdvecs = ROTATE_COORD(lbvecs, in='G', out='C')
VEC2ANG, rdvecs, dec, ra, /astro
getcoord, crval1,crpix1,cdelt1,crval2,crpix2,cdelt2,headq
rapix =  (ra  - crval1)/cdelt1 + crpix1 - 1 + npad
IF (MIN(rapix) LT npad) THEN rapix[WHERE(rapix LT npad)] += qulen - 1
decpix = (dec - crval2)/cdelt2 + crpix2 - 1

qu_heal = FLTARR(npix,2)

; The following assignments use a FORTRAN-like quirk of IDL, namely
;                  x[0,etc] = array 
; does the same as x[*,etc] = array
; but much faster.
;
qu_heal[0,0] =   INTERPOLATE(qu[*,*,0],rapix,decpix,cubic=-0.5,missing=0./0.)
qu_heal[0,1] = - INTERPOLATE(qu[*,*,1],rapix,decpix,cubic=-0.5,missing=0./0.)
; NB Switch to COSMO sign convention!
qu = 0

void = ROTATE_COORD(rdvecs,in='C',out='G',stokes=qu_heal)
void = 0

qu_heal[0,1] = - qu_heal[*,1]  ; Restore IAU sign convention

liqu = FLOAT( [[temptt], [qu_heal]] )
qu_heal = 0
temptt = 0

; Write out

header = [ $
"TTYPE1  = 'TEMPERATURE '       /label for field 1                 ", $
"TFORM1  = 'E       '           /data format of field: 4-byte REAL ", $
"TUNIT1  = 'mK      '           /physical unit of field 1          ", $
"TTYPE2  = 'Q_POLARIZATION '    /label for field 2                 ", $
"TFORM2  = 'E       '           /data format of field: 4-byte REAL ", $
"TUNIT2  = 'mK      '            /physical unit of field 2         ", $
"TTYPE3  = 'U_POLARIZATION '    /                                  ", $
"TFORM3  = 'E       '           /                                  ", $
"TUNIT3  = 'mK      '           /                                  ", $
"POLCCONV= 'IAU     '           / Stokes (Q,U) signs follow IAU convention", $
"HISTORY  "+prog+": Version "+version, $
"HISTORY  "+prog+": T Map is Stockert 1420 MHz Survey", $
"HISTORY  "+prog+": Q,U maps are Penticton 1400 MHz Survey", $
"HISTORY  "+prog+": Baseline of 3000 mK subtracted from original T map" $
]

outfil = 'L_iqumap_r'+STRING(nside2res(nside),FORMAT='(I1)')+'.fits'
WRITE_TQU, outfil,liqu, Coordsys='G',/NESTED,Xhdr=header

QUIT:
END
