;+
; NAME:
;   mirim_xytov2v3
;
; PURPOSE:
;   Convert MIRI imager SCA x,y pixel locations to JWST v2,v3 locations
;
; CALLING SEQUENCE:
;   mirim_xytov2v3,x,y,v2,v3,filter,[refdir=,xan=xan,yan=yan]
;
; INPUTS:
;   x       - SCA pixel x coordinate (? indexed)
;   y       - SCA pixel y coordinate (? index)
;   filter - Filter name (e.g., 'F770W')
;
; OPTIONAL INPUTS:
;   refdir - Root directory for distortion files
;
; OUTPUT:
;   v2     - JWST v2 coordinate in arcsec
;   v3     - JWST v3 coordinate in arcsec
;
; OPTIONAL OUTPUT:
;   xan    - JWST xan coordinate in arcsec
;   yan    - JWST yan coordinate in arcsec
;
; COMMENTS:
;   Works with CDP6 delivery files.  Inverse function is mirim_v2v3toxy
;   Based on the IDL example code provided by Alistair Glasse in
;   CDP-6.
;
;   Note that both input and output can be vectors of numbers.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   08-Sep-2016  Written by David Law (dlaw@stsci.edu)
;   17-Oct-2016  Input/output v2/v3 in arcsec (D. Law)
;-
;------------------------------------------------------------------------------


pro mirim_xytov2v3,xpixel,ypixel,v2,v3,filter,refdir=refdir,xan=xan,yan=yan

if (~keyword_set(refdir)) then $
  refdir=concat_dir(ml_getenv('JWSTTOOLS_DIR'),'mirim/distfiles/cdp6/')

reffile='MIRI_FM_MIRIMAGE_DISTORTION_06.03.00.fits'
reffile=concat_dir(refdir,reffile)
; Read global header
hdr=headfits(reffile)

; AI-matrix
AI=mrdfits(reffile,'AI matrix')
; BI-matrix
BI=mrdfits(reffile,'BI matrix')
; TI-matrix
TI=mrdfits(reffile,'TI matrix')
; MI-matrix
MI=mrdfits(reffile,'MI matrix')
; Boresight offsets
boresight=mrdfits(reffile,'boresight offsets')

; How many input points?
npoints=n_elements(xpixel)
if (n_elements(ypixel) ne npoints) then begin
  splog,'Input xpixel,ypixel array size mismatch!'
  return
endif
UnitVec=fltarr(npoints)+1.0

; What is the boresight index in the table?
indx=where(boresight.filter eq filter)
if (indx eq -1) then begin
  splog,'Bad boresight filter!'
  return
endif

; Transform to SCA pixel position without boresight offsets
SCA = fltarr(npoints,3)
SCA[*,0]=xpixel-boresight[indx].COL_OFFSET
SCA[*,1]=ypixel-boresight[indx].ROW_OFFSET
SCA[*,2]=UnitVec

; Transform to Detector Focal Plane coordinates
DFP = MI ## SCA
DFP_X = DFP[*,0]
DFP_Y = DFP[*,1]

; Transform to Entract Focal Plane coordinates
Xout = [[UnitVec],[DFP_X],[DFP_X^2],[DFP_X^3],[DFP_X^4]]
Yout = [[UnitVec],[DFP_Y],[DFP_Y^2],[DFP_Y^3],[DFP_Y^4]]
Xin=fltarr(npoints)
Yin=fltarr(npoints)
for i=0,npoints-1 do begin
  Xin[i]=transpose(Yout[i,*]) ## AI ## Xout[i,*]
  Yin[i]=transpose(Yout[i,*]) ## BI ## Xout[i,*]
endfor
EFP=[[Xin],[Yin],[UnitVec]]

; Transform to XAN, YAN
JWST_XYAN = TI ## EFP

; Transform to V2,V3
JWST=fltarr(npoints,3)
JWST[*,0] = JWST_XYAN[*,0]
JWST[*,1] = -JWST_XYAN[*,1]-7.8
JWST[*,2] = JWST_XYAN[*,2]

; Output vectors in units of arcsec
v2=JWST[*,0]*60.
v3=JWST[*,1]*60.
xan=JWST_XYAN[*,0]*60.
yan=JWST_XYAN[*,1]*60.

return
end

