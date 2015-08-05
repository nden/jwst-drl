;+
; NAME:
;   mmrs_siaf
;
; PURPOSE:
;   Figure out the siaf locations for MRS slices based on reference files
;
; CALLING SEQUENCE:
;   mrs_siaf,channel
;
; INPUTS:
;   channel - channel name (e.g, '1A')
;
; OPTIONAL INPUTS:
;   rootdir - Root directory for distortion files
;
; OUTPUT:
;   siaf_[channel].txt  Slicer corner coordinates
;   siaf_[channel]ab.ps  Plot of alpha,beta corner coordinates
;   siaf_[channel]v2v3.ps  Plot of v2,v3 corner coordinates
;
; COMMENTS:
;   Works with CDP4 delivery files.  Re
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

pro mmrs_siaf,channel,rootdir=rootdir

if (~keyword_set(rootdir)) then $
  rootdir='~/jwst/trunk/mirimrs/distfiles/cdp4/'

; Strip input channel into components, e.g.
; if channel='1A' then
; ch=1 and sband='A'
channel=strupcase(channel)
ch=fix(strmid(channel,0,1))
sband=strmid(channel,1,1)

; Determine output files to put the results in
outfile='siaf_'+channel+'.txt'
openw,lun,outfile,/get_lun,width=250

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
reffile=concat_dir(rootdir,reffile)

; Read global header
hdr=headfits(reffile)
; Get beta zeropoint and spacing from header
beta0=fxpar(hdr,'B_ZERO'+strcompress(string(ch),/remove_all))
dbeta=fxpar(hdr,'B_DEL'+strcompress(string(ch),/remove_all))

; Read FoV alpha boundaries
extname='FoV-CH'+strcompress(string(ch),/remove_all)
alphalimits=mrdfits(reffile,extname)

; Read alpha,beta -> v2,v3 table
convtable=mrdfits(reffile,'al,be->V2/V3')
; Determine which rows we need
v2index=where(strcompress(convtable.(0),/remove_all) eq 'T_CH'+channel+',V2')
v3index=where(strcompress(convtable.(0),/remove_all) eq 'T_CH'+channel+',V3')
if ((v2index lt 0)or(v3index lt 0)) then exit
; Trim to relevant v2, v3 rows for this channel
conv_v2=convtable[v2index]
conv_v3=convtable[v3index]

; Determine number of slices
nslices=n_elements(alphalimits)
; Create a 1-indexed vector of slice numbers and slice names
; (the names will be of the form 112A for ch 1, slice 12,
; band A)
slicenum=indgen(nslices)+1
slicename=string(ch*100+slicenum)+sband

; Figure out beta boundaries of each slice
beta1=beta0+(slicenum-0.5)*dbeta; Lower bound
beta2=beta1+dbeta; Upper bound

; Convert from our list of maximum and minimum alpha,beta
; to actual corner coordinates for each slice
alpha_corners=fltarr(4,nslices)
beta_corners=fltarr(4,nslices)
; Order is lower-left, upper-left, upper-right, lower-right
for i=0,nslices-1 do begin
  alpha_corners[0,i]=alphalimits[i].(0)
  alpha_corners[1,i]=alphalimits[i].(0)
  alpha_corners[2,i]=alphalimits[i].(1)
  alpha_corners[3,i]=alphalimits[i].(1)
  beta_corners[0,i]=beta1[i]
  beta_corners[1,i]=beta2[i]
  beta_corners[2,i]=beta2[i]
  beta_corners[3,i]=beta1[i]
endfor

; Compute corner coordinates for an inscribed footprint
inscr_alpha=fltarr(4)
inscr_beta=fltarr(4)
inscr_alpha[0]=max(alpha_corners[0,*])
inscr_alpha[1]=max(alpha_corners[0,*])
inscr_alpha[2]=min(alpha_corners[2,*])
inscr_alpha[3]=min(alpha_corners[3,*])
inscr_beta[0]=min(beta_corners)
inscr_beta[1]=max(beta_corners)
inscr_beta[2]=max(beta_corners)
inscr_beta[3]=min(beta_corners)

; Convert to v2,v3 corner coordinates
v2_corners=fltarr(4,nslices)
v3_corners=fltarr(4,nslices)
for i=0,nslices-1 do begin
  for j=0,3 do begin
    v2_corners[j,i]=conv_v2.(1)+conv_v2.(2)*alpha_corners[j,i] + $
      conv_v2.(3)*beta_corners[j,i] + conv_v2.(4)*alpha_corners[j,i]*beta_corners[j,i]
   v3_corners[j,i]=conv_v3.(1)+conv_v3.(2)*alpha_corners[j,i] + $
      conv_v3.(3)*beta_corners[j,i] + conv_v3.(4)*alpha_corners[j,i]*beta_corners[j,i]
  endfor
endfor

; Convert corners of inscribed box
inscr_v2=conv_v2.(1)+conv_v2.(2)*inscr_alpha + $
      conv_v2.(3)*inscr_beta + conv_v2.(4)*inscr_alpha*inscr_beta
inscr_v3=conv_v3.(1)+conv_v3.(2)*inscr_alpha + $
      conv_v3.(3)*inscr_beta + conv_v3.(4)*inscr_alpha*inscr_beta

; Print all of the corner coordinates to a file
printf,lun,'# SliceName SliceNum a_ll b_ll v2_ll v3_ll a_ul b_ul v2_ul v3_ul a_ur b_ur v2_ur v3_ur a_lr b_lr v2_lr v3_lr'
printf,lun,channel,'   -1',inscr_alpha[0],inscr_beta[0],inscr_v2[0],inscr_v3[0],$
  inscr_alpha[1],inscr_beta[1],inscr_v2[1],inscr_v3[1],$
  inscr_alpha[2],inscr_beta[2],inscr_v2[2],inscr_v3[2],$
  inscr_alpha[3],inscr_beta[3],inscr_v2[3],inscr_v3[3]
for i=0,nslices-1 do begin
  printf,lun,slicename[i],slicenum[i],$
    alpha_corners[0,i],beta_corners[0,i],v2_corners[0,i],v3_corners[0,i],$
    alpha_corners[1,i],beta_corners[1,i],v2_corners[1,i],v3_corners[1,i],$
    alpha_corners[2,i],beta_corners[2,i],v2_corners[2,i],v3_corners[2,i],$
    alpha_corners[3,i],beta_corners[3,i],v2_corners[3,i],v3_corners[3,i]
endfor

; Plot the corners in alpha,beta for this subband
plotname='siaf_'+channel+'ab.ps'
set_plot,'ps'
device,filename=plotname,/color
loadct,39
plot,alpha_corners[*,0],beta_corners[*,0],/nodata,xrange=[min(alpha_corners),max(alpha_corners)],yrange=[min(beta_corners),max(beta_corners)],xstyle=1,ystyle=1,xtitle='Alpha',ytitle='Beta',xcharsize=1.3,ycharsize=1.3,title=channel
seed=56
colors=randomu(seed,nslices)*250
for i=0,nslices-1 do begin
  oplot,[alpha_corners[*,i],alpha_corners[0,i]],[beta_corners[*,i],beta_corners[0,i]],color=colors[i]
endfor
oplot,inscr_alpha,inscr_beta
device,/close

; Plot the corners in v2,v3
plotname='siaf_'+channel+'v2v3.ps'
device,filename=plotname,/color
loadct,39
plot,v2_corners[*,0],v3_corners[*,0],/nodata,xrange=[min(v2_corners),max(v2_corners)],yrange=[min(v3_corners),max(v3_corners)],xstyle=1,ystyle=1,xtitle='V2',ytitle='V3',xcharsize=1.3,ycharsize=1.3,xmargin=12,ymargin=5,title=channel
for i=0,nslices-1 do begin
  oplot,[v2_corners[*,i],v2_corners[0,i]],[v3_corners[*,i],v3_corners[0,i]],color=colors[i]
endfor
oplot,[inscr_v2,inscr_v2[0]],[inscr_v3,inscr_v3[0]]
device,/close
set_plot,'x'

close,lun
free_lun,lun

return
end
