; This program takes the EC-provided MRS dither patterns
; (which are in alpha,beta relative to the Ch1A pointing)
; and transforms them to XIdl, YIdl offsets relative
; to the pointing origin of each channel.
;
; This is because APT works by applying the dithers optimized
; for Ch2,3,4 to the central point of those channels defining
; the pointing origin.
;
; As of Feb 2017, it also computes the 'EXTENDED SOURCE'
; patterns based on the EC provided locations.

pro mmrs_dithers,rootdir=rootdir,siafdir=siafdir,outdir=outdir

; This is where the dither input file from the EC lives
if (~keyword_set(rootdir)) then $
  rootdir='~/jwst/trunk/mirimrs/dithers/oct2016/'

; This is where the SIAF parameter files created from mmrs_siaf live
if (~keyword_set(siafdir)) then $
  siafdir='/Users/dlaw/STSCI/MIRI/SIAF/MRS/Nov2016/'

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

; Construct the extended source patterns based on points 1/3 and 5/7
; of the original tables.  Recenter around their midpoint.
dithers=[dithers,dithers[0:15]]
; Loop over channels
for i=0,3 do begin
  ; Dither position indices
  dithers[40+i*4+0].dpos=40+i*4+1
  dithers[40+i*4+1].dpos=40+i*4+2
  dithers[40+i*4+2].dpos=40+i*4+3
  dithers[40+i*4+3].dpos=40+i*4+4
  ; First pair is entries 1/3
  mid_alpha=(dithers[i*8+0].alpha1a+dithers[i*8+2].alpha1a)/2.
  mid_beta=(dithers[i*8+0].beta1a+dithers[i*8+2].beta1a)/2.
  dithers[40+i*4+0].alpha1a=dithers[i*8+0].alpha1a-mid_alpha
  dithers[40+i*4+1].alpha1a=dithers[i*8+2].alpha1a-mid_alpha
  dithers[40+i*4+0].beta1a=dithers[i*8+0].beta1a-mid_beta
  dithers[40+i*4+1].beta1a=dithers[i*8+2].beta1a-mid_beta
  ; Second pair is entries 5/7
  mid_alpha=(dithers[i*8+4].alpha1a+dithers[i*8+6].alpha1a)/2.
  mid_beta=(dithers[i*8+4].beta1a+dithers[i*8+6].beta1a)/2.
  dithers[40+i*4+2].alpha1a=dithers[i*8+4].alpha1a-mid_alpha
  dithers[40+i*4+3].alpha1a=dithers[i*8+6].alpha1a-mid_alpha
  dithers[40+i*4+2].beta1a=dithers[i*8+4].beta1a-mid_beta
  dithers[40+i*4+3].beta1a=dithers[i*8+6].beta1a-mid_beta
  ; Add band information
  if (i eq 0) then dithers[40+i*4:40+(i+1)*4-1].band='1A'
  if (i eq 1) then dithers[40+i*4:40+(i+1)*4-1].band='2A'
  if (i eq 2) then dithers[40+i*4:40+(i+1)*4-1].band='3A'
  if (i eq 3) then dithers[40+i*4:40+(i+1)*4-1].band='4A'
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




thetime=strsplit(systime(),' ',/extract)
thedate=strcompress(thetime[1]+' '+thetime[4])


plotname=concat_dir(outdir,'dithers.ps')
set_plot,'ps'
device,filename=plotname,/color,xsize=16,ysize=15
loadct,39

; Plot field bounding boxes
plot,box1A_v2,box1A_v3,xrange=[-8.29,-8.49]*60,yrange=[-5.43,-5.23]*60,/xstyle,/ystyle,xthick=5,ythick=5,thick=4,charsize=1.5,xtitle='V2 (arcsec)',ytitle='V3 (arcsec)',charthick=4,title=strcompress('MRS Dithers: Pre-flight ('+thedate+')')
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
oplot,siaf1a[0].v2_ref-dithers[32:39].dxidl,siaf1a[0].v3_ref+dithers[32:39].dyidl,psym=1,thick=4

; Plot dither points 1-32 with circles showing PSF FWHM
theta=findgen(360)
oplot,siaf1a[0].v2_ref-dithers[0:7].dxidl,siaf1a[0].v3_ref+dithers[0:7].dyidl,psym=1,thick=4,color=60
rad=0.3/60.
for i=0,7 do begin
  xtemp=siaf1a[0].v2_ref-dithers[i].dxidl+rad*cos(theta*!PI/180.)
  ytemp=siaf1a[0].v3_ref+dithers[i].dyidl+rad*sin(theta*!PI/180.)
  ;oplot,xtemp,ytemp,thick=4,color=60
endfor

oplot,siaf2a[0].v2_ref-dithers[8:15].dxidl,siaf2a[0].v3_ref+dithers[8:15].dyidl,psym=1,color=140,thick=4
rad=0.3/60.
for i=8,15 do begin
  xtemp=siaf2a[0].v2_ref-dithers[i].dxidl+rad*cos(theta*!PI/180.)
  ytemp=siaf2a[0].v3_ref+dithers[i].dyidl+rad*sin(theta*!PI/180.)
  ;oplot,xtemp,ytemp,color=140,thick=4
endfor

oplot,siaf3a[0].v2_ref-dithers[16:23].dxidl,siaf3a[0].v3_ref+dithers[16:23].dyidl,psym=1,color=200,thick=4
rad=0.3/60.
for i=16,23 do begin
  xtemp=siaf3a[0].v2_ref-dithers[i].dxidl+rad*cos(theta*!PI/180.)
  ytemp=siaf3a[0].v3_ref+dithers[i].dyidl+rad*sin(theta*!PI/180.)
  ;oplot,xtemp,ytemp,color=200,thick=4
endfor

oplot,siaf4a[0].v2_ref-dithers[24:31].dxidl,siaf4a[0].v3_ref+dithers[24:31].dyidl,psym=1,color=250,thick=4
rad=0.3/60.
for i=24,31 do begin
  xtemp=siaf4a[0].v2_ref-dithers[i].dxidl+rad*cos(theta*!PI/180.)
  ytemp=siaf4a[0].v3_ref+dithers[i].dyidl+rad*sin(theta*!PI/180.)
  ;oplot,xtemp,ytemp,color=250,thick=4
endfor

device,/close
spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))




; Make a plot illustrating a 4-pt ALL dither
plotname=concat_dir(outdir,'dithers_4ptall.ps')
set_plot,'ps'
device,filename=plotname,/color,xsize=16,ysize=15
loadct,39

; Set up plot
plot,box1A_v2,box1A_v3,xrange=[-8.29,-8.49]*60,yrange=[-5.43,-5.23]*60,/nodata,/xstyle,/ystyle,xthick=5,ythick=5,thick=4,charsize=1.5,xtitle='V2 (arcsec)',ytitle='V3 (arcsec)',charthick=4,title=strcompress('MRS Dithers: Pre-flight ('+thedate+')')
oplot,siaf1a[0].v2_ref-dithers[0:3].dxidl,siaf1a[0].v3_ref+dithers[0:3].dyidl,psym=1,thick=4,color=60
oplot,box1A_v2-dithers[0].dxidl,box1A_v3+dithers[0].dyidl,color=60,thick=4
oplot,box1A_v2-dithers[1].dxidl,box1A_v3+dithers[1].dyidl,color=60,thick=4
oplot,box1A_v2-dithers[2].dxidl,box1A_v3+dithers[2].dyidl,color=60,thick=4
oplot,box1A_v2-dithers[3].dxidl,box1A_v3+dithers[3].dyidl,color=60,thick=4
oplot,box4A_v2-dithers[0].dxidl,box4A_v3+dithers[0].dyidl,color=250,thick=4
oplot,box4A_v2-dithers[1].dxidl,box4A_v3+dithers[1].dyidl,color=250,thick=4
oplot,box4A_v2-dithers[2].dxidl,box4A_v3+dithers[2].dyidl,color=250,thick=4
oplot,box4A_v2-dithers[3].dxidl,box4A_v3+dithers[3].dyidl,color=250,thick=4

device,/close
spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))




; Make a plot illustrating a MOSAIC dither
plotname=concat_dir(outdir,'dithers_mos.ps')
set_plot,'ps'
device,filename=plotname,/color,xsize=16,ysize=15
loadct,39

; Set up plot
plot,box1A_v2,box1A_v3,xrange=[-8.29,-8.49]*60,yrange=[-5.43,-5.23]*60,/nodata,/xstyle,/ystyle,xthick=5,ythick=5,thick=4,charsize=1.5,xtitle='V2 (arcsec)',ytitle='V3 (arcsec)',charthick=4,title=strcompress('MRS Dithers: Pre-flight ('+thedate+')')
oplot,siaf1a[0].v2_ref-dithers[32:33].dxidl,siaf1a[0].v3_ref+dithers[32:33].dyidl,psym=1,thick=4
oplot,box1A_v2-dithers[32].dxidl,box1A_v3+dithers[32].dyidl,color=60,thick=4
oplot,box1A_v2-dithers[33].dxidl,box1A_v3+dithers[33].dyidl,color=60,thick=4
oplot,box4A_v2-dithers[32].dxidl,box4A_v3+dithers[32].dyidl,color=250,thick=4
oplot,box4A_v2-dithers[33].dxidl,box4A_v3+dithers[33].dyidl,color=250,thick=4

device,/close
spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))



; Make a plot illustrating a EXTENDED dither
plotname=concat_dir(outdir,'dithers_ext.ps')
set_plot,'ps'
device,filename=plotname,/color,xsize=16,ysize=15
loadct,39

; Set up plot
plot,box1A_v2,box1A_v3,xrange=[-8.29,-8.49]*60,yrange=[-5.43,-5.23]*60,/nodata,/xstyle,/ystyle,xthick=5,ythick=5,thick=4,charsize=1.5,xtitle='V2 (arcsec)',ytitle='V3 (arcsec)',charthick=4,title=strcompress('MRS Dithers: Pre-flight ('+thedate+')')
oplot,siaf1a[0].v2_ref-dithers[32:33].dxidl,siaf1a[0].v3_ref+dithers[32:33].dyidl,psym=1,thick=4
oplot,box1A_v2-dithers[40].dxidl,box1A_v3+dithers[40].dyidl,color=60,thick=4
oplot,box1A_v2-dithers[41].dxidl,box1A_v3+dithers[41].dyidl,color=60,thick=4
oplot,box4A_v2-dithers[40].dxidl,box4A_v3+dithers[40].dyidl,color=250,thick=4
oplot,box4A_v2-dithers[41].dxidl,box4A_v3+dithers[41].dyidl,color=250,thick=4

device,/close
spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))





return
end
