;+
; NAME:
;   mirim_v2v3toxy
;
; PURPOSE:
;   Convert JWST v2,v3 locations to MIRI Imager SCA x,y pixel locations
;
; CALLING SEQUENCE:
;   mirim_v2v3toxy,v2,v3,x,y,filter,[refdir=,/xan]
;
; INPUTS:
;   v2     - JWST v2 coordinate in arcminutes
;   v3     - JWST v3 coordinate in arcminutes
;   filter - Filter name (e.g., 'F770W')
;
; OPTIONAL INPUTS:
;   refdir - Root directory for distortion files
;   /xan   - Specifies that input coordinates are actually xan,yan frame
;
; OUTPUT:
;   x       - SCA pixel x coordinate (? indexed)
;   y       - SCA pixel y coordinate (? index)
;
; OPTIONAL OUTPUT:
;
; COMMENTS:
;   Works with CDP6 delivery files.  Inverse function is mirim_xytov2v3.
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
;-
;------------------------------------------------------------------------------

pro mirim_v2v3toxy,v2,v3,xpixel,ypixel,filter,refdir=refdir,xan=xan

if (~keyword_set(refdir)) then $
  refdir=concat_dir(ml_getenv('JWSTTOOLS_DIR'),'mirim/distfiles/cdp6/')

reffile='MIRI_FM_MIRIMAGE_DISTORTION_06.03.00.fits'
reffile=concat_dir(refdir,reffile)
; Read global header
hdr=headfits(reffile)

; A-matrix
A=mrdfits(reffile,'A matrix')
; B-matrix
B=mrdfits(reffile,'B matrix')
; T-matrix
T=mrdfits(reffile,'T matrix')
; M-matrix
M=mrdfits(reffile,'M matrix')
; Boresight offsets
boresight=mrdfits(reffile,'boresight offsets')

; How many input points?
npoints=n_elements(v2)
if (n_elements(v3) ne npoints) then begin
  splog,'Input v2,v3 array size mismatch!'
  return
endif
UnitVec=fltarr(npoints)+1.0

; Set up JWST V2,V3 arrays
JWST=fltarr(npoints,3)
; If inputs were really in XAN, YAN:
if (keyword_set(xan)) then begin
  JWST[*,0]=v2
  JWST[*,1]=-v3-7.8
  JWST[*,2]=UnitVec
; Otherwise inputs really were v2,v3
endif else begin
  JWST[*,0]=v2
  JWST[*,1]=v3
  JWST[*,2]=UnitVec
endelse

; Convert to XAN, YAN
JWST_XYAN=fltarr(npoints,3)
JWST_XYAN[*,0] = JWST[*,0]
JWST_XYAN[*,1] = -JWST[*,1]-7.8
JWST_XYAN[*,2] = JWST[*,2]

; Compute MIRI Entrance Focal Plane coordinates
EFP=T ## JWST_XYAN
; Components
EFP_X=EFP[*,0]
EFP_Y=EFP[*,1]

; Transform from EFP to MIRI Detector Focal Plane
Xin = [[UnitVec], [EFP_X], [EFP_X^2], [EFP_X^3], [EFP_X^4]]
Yin = [[UnitVec], [EFP_Y], [EFP_Y^2], [EFP_Y^3], [EFP_Y^4]]
Xout=fltarr(npoints)
Yout=fltarr(npoints)
for i=0,npoints-1 do begin
  Xout[i]=transpose(Yin[i,*]) ## A ## Xin[i,*]
  Yout[i]=transpose(Yin[i,*]) ## B ## Xin[i,*]
endfor
DFP=[[Xout],[Yout],[UnitVec]]

; Transform to SCA pixel position
SCA = M ## DFP

; Add boresight offset
; What is the boresight index in the table?
indx=where(boresight.filter eq filter)
if (indx eq -1) then begin
  splog,'Bad boresight filter!'
  return
endif
SCA[*,0] += boresight[indx].COL_OFFSET
SCA[*,1] += boresight[indx].ROW_OFFSET

; Output vectors
xpixel=SCA[*,0]
ypixel=SCA[*,1]

return
end
