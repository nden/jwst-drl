; This program takes the EC-provided MRS dither patterns
; (which are in alpha,beta relative to the Ch1A pointing)
; and transforms them to XIdl, YIdl offsets relative
; to the pointing origin of each channel.
;
; This is because APT works by applying the dithers optimized
; for Ch2,3,4 to the central point of those channels defining
; the pointing origin.

pro mmrs_dithers,rootdir=rootdir,siafdir=siafdir,outdir=outdir

; This is where the dither input file from the EC lives
if (~keyword_set(rootdir)) then $
  rootdir='~/jwst/trunk/mirimrs/dithers/oct2016/'

; This is where the SIAF parameter files created from mmrs_siaf live
if (~keyword_set(siafdir)) then $
  siafdir='/Users/dlaw/STSCI/MIRI/SIAF/MRS/Sep2016/'

; This is where the results will go
if (~keyword_set(outdir)) then $
  outdir='./'

input=yanny_readone(concat_dir(rootdir,'ditherinput.par'))

; Output in .par format (easy to use programmatically) and
; in .txt format for easy import to APT Excel (Import-textfile)
; Output is your local directory
outpar=concat_dir(outdir,'dithers.par')
outapt=concat_dir(outdir,'apt.txt'); Import this to Excel from Import-textfile

; Need to reorganize Alistair's table to swap each position 3 and 4
dithers=input
for i=1,8 do begin
  dithers[i*4-1].alpha1A=input[i*4-2].alpha1A
  dithers[i*4-2].alpha1A=input[i*4-1].alpha1A
  dithers[i*4-1].beta1A=input[i*4-2].beta1A
  dithers[i*4-2].beta1A=input[i*4-1].beta1A
endfor

dithers=jjadd_tag(dithers,'v2',0.)
dithers=jjadd_tag(dithers,'v3',0.)
dithers=jjadd_tag(dithers,'dxidl',0.)
dithers=jjadd_tag(dithers,'dyidl',0.)

ndither=n_elements(dithers)
for i=0,ndither-1 do begin
  ; Convert desired 1A frame alpha,beta position to a v2,v3 position
  mmrs_abtov2v3,dithers[i].alpha1A,dithers[i].beta1A,tempv2,tempv3,'1A',xan=xan,yan=yan
  dithers[i].v2=tempv2
  dithers[i].v3=tempv3

  ; Figure out the zeropoint location of this band
  ; (local alpha=beta=0 for 1A,2A,3A,or 4A) in v2,v3 coordinates
  mmrs_abtov2v3,0.,0.,zpv2,zpv3,strtrim(dithers[i].band,2)

  ; Convert both the dither locations and the zeropoint locations
  ; to the XIdl, YIdl reference frame
  mmrs_v2v3toideal,zpv2,zpv3,zpx,zpy
  mmrs_v2v3toideal,tempv2,tempv3,tempx,tempy
  ; Determine dXIdl, dYIdl offsets.  Note that these dXIdl,dYIdl are OFFSETS
  ; from the base location rather than POSITIONS (as for v2,v3)
  dithers[i].dxidl=tempx-zpx
  dithers[i].dyidl=tempy-zpy
endfor

; Write out new parfile
yanny_write,outpar,ptr_new(dithers)

; Write out text file
openw,lun,outapt,/get_lun
; Make sure dither positions go 1-4
printf,lun,'PosnIndex   dXIdeal(arcsec)   dYIdeal(arcsec)'
for i=0,ndither-1 do $
  printf,lun,dithers[i].dpos-floor(i/4.)*4,dithers[i].dxidl,dithers[i].dyidl
close,lun
free_lun,lun

; Quality control plots
;start=33-1
;stop=40-1
;plot,dithers[start:stop].alpha1A,dithers[start:stop].beta1A,psym=1
;for i=start,stop do xyouts,dithers[i].alpha1A,dithers[i].beta1A,i+1

; Read in SIAF information about the nominal field boundaries
siaf1a=yanny_readone(concat_dir(siafdir,'siaf_1A.par'))
siaf1b=yanny_readone(concat_dir(siafdir,'siaf_1B.par'))
siaf1c=yanny_readone(concat_dir(siafdir,'siaf_1C.par'))
siaf2a=yanny_readone(concat_dir(siafdir,'siaf_2A.par'))
siaf2b=yanny_readone(concat_dir(siafdir,'siaf_2B.par'))
siaf2c=yanny_readone(concat_dir(siafdir,'siaf_2C.par'))
siaf3a=yanny_readone(concat_dir(siafdir,'siaf_3A.par'))
siaf3b=yanny_readone(concat_dir(siafdir,'siaf_3B.par'))
siaf3c=yanny_readone(concat_dir(siafdir,'siaf_3C.par'))
siaf4a=yanny_readone(concat_dir(siafdir,'siaf_4A.par'))
siaf4b=yanny_readone(concat_dir(siafdir,'siaf_4B.par'))
siaf4c=yanny_readone(concat_dir(siafdir,'siaf_4C.par'))

; Define field boundaries
  box1A_v2=[siaf1a[0].v2_ll,siaf1a[0].v2_ul,siaf1a[0].v2_ur,siaf1a[0].v2_lr,siaf1a[0].v2_ll]
  box1A_v3=[siaf1a[0].v3_ll,siaf1a[0].v3_ul,siaf1a[0].v3_ur,siaf1a[0].v3_lr,siaf1a[0].v3_ll]
  box1B_v2=[siaf1b[0].v2_ll,siaf1b[0].v2_ul,siaf1b[0].v2_ur,siaf1b[0].v2_lr,siaf1b[0].v2_ll]
  box1B_v3=[siaf1b[0].v3_ll,siaf1b[0].v3_ul,siaf1b[0].v3_ur,siaf1b[0].v3_lr,siaf1b[0].v3_ll]
  box1C_v2=[siaf1c[0].v2_ll,siaf1c[0].v2_ul,siaf1c[0].v2_ur,siaf1c[0].v2_lr,siaf1c[0].v2_ll]
  box1C_v3=[siaf1c[0].v3_ll,siaf1c[0].v3_ul,siaf1c[0].v3_ur,siaf1c[0].v3_lr,siaf1c[0].v3_ll]
  box2A_v2=[siaf2a[0].v2_ll,siaf2a[0].v2_ul,siaf2a[0].v2_ur,siaf2a[0].v2_lr,siaf2a[0].v2_ll]
  box2A_v3=[siaf2a[0].v3_ll,siaf2a[0].v3_ul,siaf2a[0].v3_ur,siaf2a[0].v3_lr,siaf2a[0].v3_ll]
  box2B_v2=[siaf2b[0].v2_ll,siaf2b[0].v2_ul,siaf2b[0].v2_ur,siaf2b[0].v2_lr,siaf2b[0].v2_ll]
  box2B_v3=[siaf2b[0].v3_ll,siaf2b[0].v3_ul,siaf2b[0].v3_ur,siaf2b[0].v3_lr,siaf2b[0].v3_ll]
  box2C_v2=[siaf2c[0].v2_ll,siaf2c[0].v2_ul,siaf2c[0].v2_ur,siaf2c[0].v2_lr,siaf2c[0].v2_ll]
  box2C_v3=[siaf2c[0].v3_ll,siaf2c[0].v3_ul,siaf2c[0].v3_ur,siaf2c[0].v3_lr,siaf2c[0].v3_ll]
  box3A_v2=[siaf3a[0].v2_ll,siaf3a[0].v2_ul,siaf3a[0].v2_ur,siaf3a[0].v2_lr,siaf3a[0].v2_ll]
  box3A_v3=[siaf3a[0].v3_ll,siaf3a[0].v3_ul,siaf3a[0].v3_ur,siaf3a[0].v3_lr,siaf3a[0].v3_ll]
  box3B_v2=[siaf3b[0].v2_ll,siaf3b[0].v2_ul,siaf3b[0].v2_ur,siaf3b[0].v2_lr,siaf3b[0].v2_ll]
  box3B_v3=[siaf3b[0].v3_ll,siaf3b[0].v3_ul,siaf3b[0].v3_ur,siaf3b[0].v3_lr,siaf3b[0].v3_ll]
  box3C_v2=[siaf3c[0].v2_ll,siaf3c[0].v2_ul,siaf3c[0].v2_ur,siaf3c[0].v2_lr,siaf3c[0].v2_ll]
  box3C_v3=[siaf3c[0].v3_ll,siaf3c[0].v3_ul,siaf3c[0].v3_ur,siaf3c[0].v3_lr,siaf3c[0].v3_ll]
  box4A_v2=[siaf4a[0].v2_ll,siaf4a[0].v2_ul,siaf4a[0].v2_ur,siaf4a[0].v2_lr,siaf4a[0].v2_ll]
  box4A_v3=[siaf4a[0].v3_ll,siaf4a[0].v3_ul,siaf4a[0].v3_ur,siaf4a[0].v3_lr,siaf4a[0].v3_ll]
  box4B_v2=[siaf4b[0].v2_ll,siaf4b[0].v2_ul,siaf4b[0].v2_ur,siaf4b[0].v2_lr,siaf4b[0].v2_ll]
  box4B_v3=[siaf4b[0].v3_ll,siaf4b[0].v3_ul,siaf4b[0].v3_ur,siaf4b[0].v3_lr,siaf4b[0].v3_ll]
  box4C_v2=[siaf4c[0].v2_ll,siaf4c[0].v2_ul,siaf4c[0].v2_ur,siaf4c[0].v2_lr,siaf4c[0].v2_ll]
  box4C_v3=[siaf4c[0].v3_ll,siaf4c[0].v3_ul,siaf4c[0].v3_ur,siaf4c[0].v3_lr,siaf4c[0].v3_ll]

plotname=concat_dir(outdir,'dithers.ps')
set_plot,'ps'
device,filename=plotname,/color,xsize=16,ysize=15
loadct,39

; Plot field bounding boxes
plot,box1A_v2,box1A_v3,xrange=[-8.29,-8.49],yrange=[-5.43,-5.23],/xstyle,/ystyle,xthick=5,ythick=5,thick=4,charsize=1.5,xtitle='V2 (arcmin)',ytitle='V3 (arcmin)',charthick=4,title='MRS Dithers'
oplot,box1A_v2,box1A_v3,color=60,thick=4
oplot,box1B_v2,box1B_v3,color=60,thick=4
oplot,box1C_v2,box1C_v3,color=60,thick=4
oplot,box2A_v2,box2A_v3,color=140,thick=4
oplot,box2B_v2,box2B_v3,color=140,thick=4
oplot,box2C_v2,box2C_v3,color=140,thick=4
oplot,box3A_v2,box3A_v3,color=200,thick=4
oplot,box3B_v2,box3B_v3,color=200,thick=4
oplot,box3C_v2,box3C_v3,color=200,thick=4
oplot,box4A_v2,box4A_v3,color=250,thick=4
oplot,box4B_v2,box4B_v3,color=250,thick=4
oplot,box4C_v2,box4C_v3,color=250,thick=4

; Plot dither points 33-40 simply
oplot,siaf1a[0].v2_ref-dithers[32:39].dxidl/60.,siaf1a[0].v3_ref+dithers[32:39].dyidl/60.,psym=1,thick=4

; Plot dither points 1-32 with circles showing PSF FWHM
theta=findgen(360)
oplot,siaf1a[0].v2_ref-dithers[0:7].dxidl/60.,siaf1a[0].v3_ref+dithers[0:7].dyidl/60.,psym=1,thick=4,color=60
rad=0.3/60.
for i=0,7 do begin
  xtemp=siaf1a[0].v2_ref-dithers[i].dxidl/60.+rad*cos(theta*!PI/180.),thick=4,color=60
  ytemp=siaf1a[0].v3_ref+dithers[i].dyidl/60.+rad*sin(theta*!PI/180.),thick=4,color=60
  ;oplot,xtemp,ytemp
endfor

oplot,siaf2a[0].v2_ref-dithers[8:15].dxidl/60.,siaf2a[0].v3_ref+dithers[8:15].dyidl/60.,psym=1,color=140,thick=4
rad=0.3/60.
for i=8,15 do begin
  xtemp=siaf2a[0].v2_ref-dithers[i].dxidl/60.+rad*cos(theta*!PI/180.),color=140,thick=4
  ytemp=siaf2a[0].v3_ref+dithers[i].dyidl/60.+rad*sin(theta*!PI/180.),color=140,thick=4
  ;oplot,xtemp,ytemp,color=100
endfor

oplot,siaf3a[0].v2_ref-dithers[16:23].dxidl/60.,siaf3a[0].v3_ref+dithers[16:23].dyidl/60.,psym=1,color=200,thick=4
rad=0.3/60.
for i=16,23 do begin
  xtemp=siaf3a[0].v2_ref-dithers[i].dxidl/60.+rad*cos(theta*!PI/180.),color=200,thick=4
  ytemp=siaf3a[0].v3_ref+dithers[i].dyidl/60.+rad*sin(theta*!PI/180.),color=200,thick=4
  ;oplot,xtemp,ytemp,color=150
endfor

oplot,siaf4a[0].v2_ref-dithers[24:31].dxidl/60.,siaf4a[0].v3_ref+dithers[24:31].dyidl/60.,psym=1,color=250,thick=4
rad=0.3/60.
for i=24,31 do begin
  xtemp=siaf4a[0].v2_ref-dithers[i].dxidl/60.+rad*cos(theta*!PI/180.),color=250,thick=4
  ytemp=siaf4a[0].v3_ref+dithers[i].dyidl/60.+rad*sin(theta*!PI/180.),color=250,thick=4
  ;oplot,xtemp,ytemp,color=250
endfor

device,/close
spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))



return
end