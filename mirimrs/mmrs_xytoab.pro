;+
; NAME:
;   mmrs_xytoab
;
; PURPOSE:
;   Convert MRS detector coordinates to MRS local alpha,beta coordinates
;
; CALLING SEQUENCE:
;   mmrs_xytoab,x,y,a,b,channel,[slicenum=,slicename=,refdir=]
;
; INPUTS:
;   x      - X coordinate in 0-indexed pixels
;   y      - Y coordinate in 0-indexed pixels
;   channel - channel name (e.g, '1A')
;
; OPTIONAL INPUTS:
;   refdir - Root directory for distortion files
;
; OUTPUT:
;   a       - Alpha coordinate in arcsec
;   b       - Beta coordinate in arcsec
;
; OPTIONAL OUTPUT:
;   slicenum - Slice number (e.g., 11)
;   slicename - Slice name (e.g., 211A for ch2, slice 11, sub-band A)
;
; COMMENTS:
;   Works with CDP4 delivery files.  Inverse function is mmrs_abtoxy.pro
;   Not all input x,y can actually map to alpha,beta because some pixels
;   fall between slices.  alpha,beta for these are set to -999.
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
;   30-July-2015  Written by David Law (dlaw@stsci.edu)
;-
;------------------------------------------------------------------------------


pro mmrs_xytoab,x,y,a,b,channel,slicenum=slicenum,slicename=slicename,refdir=refdir

if (~keyword_set(refdir)) then $
  refdir='/Users/dlaw/jwst/mirimrs/distfiles/cdp4/'

; Strip input channel into components, e.g.
; if channel='1A' then
; ch=1 and sband='A'
ch=fix(strmid(channel,0,1))
sband=strmid(channel,1,1)

; Ensure we're not using integer inputs
xdbl=double(x)
ydbl=double(y)

; Determine input reference FITS file
case channel of
  '1A': reffile='MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_04.02.00.fits'
  '1B': reffile='MIRI_FM_MIRIFUSHORT_12MEDIUM_DISTORTION_04.02.00.fits'
  '1C': reffile='MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_04.02.00.fits'
  '2A': reffile='MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_04.02.00.fits'
  '2B': reffile='MIRI_FM_MIRIFUSHORT_12MEDIUM_DISTORTION_04.02.00.fits'
  '2C': reffile='MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_04.02.00.fits'
  '3A': reffile='MIRI_FM_MIRIFULONG_34SHORT_DISTORTION_04.02.00.fits'
  '3B': reffile='MIRI_FM_MIRIFULONG_34MEDIUM_DISTORTION_04.02.00.fits'
  '3C': reffile='MIRI_FM_MIRIFULONG_34LONG_DISTORTION_04.02.00.fits'
  '4A': reffile='MIRI_FM_MIRIFULONG_34SHORT_DISTORTION_04.02.00.fits'
  '4B': reffile='MIRI_FM_MIRIFULONG_34MEDIUM_DISTORTION_04.02.00.fits'
  '4C': reffile='MIRI_FM_MIRIFULONG_34LONG_DISTORTION_04.02.00.fits'
  else: begin
    print,'Invalid band'
    return
    end
endcase
reffile=concat_dir(refdir,reffile)

; Read global header
hdr=headfits(reffile)
; Get beta zeropoint and spacing from header
beta0=fxpar(hdr,'B_ZERO'+strcompress(string(ch),/remove_all))
dbeta=fxpar(hdr,'B_DEL'+strcompress(string(ch),/remove_all))

d2c_alpha=mrdfits(reffile,'Alpha-CH'+strcompress(string(ch),/remove_all))
d2c_slice=mrdfits(reffile,'Slice-Number')

; Define slice for these pixels
slicenum=fix(d2c_slice[x,y])-ch*100
slicename=strcompress(string(slicenum+ch*100)+sband,/remove_all)

; Define an index0 where slice number is physical
; (i.e., not between slices)
index0=where((slicenum gt 0)and(slicenum le 99),nindex0)

; Initialize a,b to -999.
; (they will only be set to something else if the pixel
; lands on a valid slice)
a=replicate(-999.,n_elements(slicenum))
b=replicate(-999.,n_elements(slicenum))

; Define beta for these pixels
if (nindex0 gt 0) then $
  b[index0]=beta0+(slicenum[index0]-1.)*dbeta

; Define alpha for these pixels
for k=0,nindex0-1 do begin
  thiscoeff=d2c_alpha[slicenum[index0[k]]-1]; Alpha coefficients for this slice
  thisalpha=0.
  for i=0,4 do begin
    for j=0,4 do begin
      coind=1+(i*5)+j
      thisalpha += thiscoeff.(coind)*(((xdbl[index0[k]]-thiscoeff.(0))^j)   * (ydbl[index0[k]]^i))
      endfor
    endfor
    a[index0[k]]=thisalpha
endfor

return
end
