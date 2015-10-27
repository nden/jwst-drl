;+
; NAME:
;   mmrs_v2v3toab
;
; PURPOSE:
;   Convert JWST v2,v3 coordinates to MRS local alpha,beta coordinates
;
; CALLING SEQUENCE:
;   mmrs_v2v3toab,v2,v3,a,b,channel,[refdir=]
;
; INPUTS:
;   v2      - V2 coordinate in arcmin
;   v3      - V3 coordinate in arcmin
;   channel - channel name (e.g, '1A')
;
; OPTIONAL INPUTS:
;   refdir - Root directory for distortion files
;
; OUTPUT:
;   a       - Alpha coordinate in arcsec
;   b       - Beta coordinate in arcsec
;
; COMMENTS:
;   Works with CDP4 delivery files.  Inverse function is mmrs_abtov2v3.pro
;
; EXAMPLES:
;
; BUGS:
;   Current the beta values returned by this function are WRONG.
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   30-July-2015  Written by David Law (dlaw@stsci.edu)
;-
;------------------------------------------------------------------------------

pro mmrs_v2v3toab,v2,v3,a,b,channel,refdir=refdir

if (~keyword_set(refdir)) then $
  refdir=concat_dir(ml_getenv('JWSTTOOLS_DIR'),'mirimrs/distfiles/cdp5b/')

; Strip input channel into components, e.g.
; if channel='1A' then
; ch=1 and sband='A'
ch=fix(strmid(channel,0,1))
sband=strmid(channel,1,1)

; Ensure we're not using integer inputs
v2dbl=double(v2)
v3dbl=double(v3)

; Determine input reference FITS file
case channel of
  '1A': reffile='MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_5B.02.00.fits'
  '1B': reffile='MIRI_FM_MIRIFUSHORT_12MEDIUM_DISTORTION_5B.02.00.fits'
  '1C': reffile='MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_5B.02.00.fits'
  '2A': reffile='MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_5B.02.00.fits'
  '2B': reffile='MIRI_FM_MIRIFUSHORT_12MEDIUM_DISTORTION_5B.02.00.fits'
  '2C': reffile='MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_5B.02.00.fits'
  '3A': reffile='MIRI_FM_MIRIFULONG_34SHORT_DISTORTION_5B.02.00.fits'
  '3B': reffile='MIRI_FM_MIRIFULONG_34MEDIUM_DISTORTION_5B.02.00.fits'
  '3C': reffile='MIRI_FM_MIRIFULONG_34LONG_DISTORTION_5B.02.00.fits'
  '4A': reffile='MIRI_FM_MIRIFULONG_34SHORT_DISTORTION_5B.02.00.fits'
  '4B': reffile='MIRI_FM_MIRIFULONG_34MEDIUM_DISTORTION_5B.02.00.fits'
  '4C': reffile='MIRI_FM_MIRIFULONG_34LONG_DISTORTION_5B.02.00.fits'
  else: begin
    print,'Invalid band'
    return
    end
endcase
reffile=concat_dir(refdir,reffile)

; Read v2,v3 -> alpha,beta table
convtable=mrdfits(reffile,'V2/V3->al,be')
; Determine which rows we need
alindex=where(strcompress(convtable.(0),/remove_all) eq 'U_CH'+channel+',al')
beindex=where(strcompress(convtable.(0),/remove_all) eq 'U_CH'+channel+',be')
if ((alindex lt 0)or(beindex lt 0)) then exit
; Trim to relevant al, be rows for this channel
conv_al=convtable[alindex]
conv_be=convtable[beindex]

a=conv_al.(1)+conv_al.(2)*v2dbl + $
      conv_al.(3)*v3dbl + conv_al.(4)*v2dbl*v3dbl
b=conv_be.(1)+conv_be.(2)*v2dbl + $
      conv_be.(3)*v3dbl + conv_be.(4)*v2dbl*v3dbl
stop
return
end
