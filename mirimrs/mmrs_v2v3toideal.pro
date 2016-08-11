;+
; NAME:
;   mmrs_v2v3toideal
;
; PURPOSE:
;   Convert V2, V3 coordinates to MRS Ideal coordinates
;
; CALLING SEQUENCE:
;   mmrs_
;
; INPUTS:
;   v2      - V2 coordinate in arcmin
;   v3      - V3 coordinate in arcmin
;
; OUTPUT:
;   XIdl    - XIdeal coordinate in arcsec
;   YIdl    - YIdeal coordinate in arcsec
;
; COMMENTS:
;   Does not do tangent-plane correction for Ideal frame; this is
;   minute effect
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
;   09-Aug-2016  Written by David Law (dlaw@stsci.edu)
;-
;------------------------------------------------------------------------------

pro mmrs_v2v3toideal,v2,v3,xidl,yidl

; Define V2REF, V3REF at middle of 1A field (alpha=beta=0.0)
V2REF = -8.3942412d ; In arcmin
V3REF = -5.3123744d ; In arcmin

xidl=-(v2-V2REF)*60.
yidl=(v3-V3REF)*60.

return
end
