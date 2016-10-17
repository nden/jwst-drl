pro miri_simcube,directory,band,imonly=imonly

;band='1A'

; Log runtime
stime0 = systime(1)
; Memory tracking
thismem = memory()
maxmem = 0

indir=concat_dir(directory,'preproc/')
outdir=concat_dir(directory,'stack/')
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir
outfile=concat_dir(outdir,'stack.fits')

; Require inputs already background subtracted and hacked
; to include dither keywords
if ((band eq '1A')or(band eq '2A')or(band eq '1B')or(band eq '2B')or(band eq '1C')or(band eq '2C')) then begin
  files=file_search(concat_dir(indir,'*MIRIFUSHORT*'))
endif else begin
  files=file_search(concat_dir(indir,'*MIRIFULONG*'))
endelse

;files=files[0]
nfiles=n_elements(files)

raref=dblarr(nfiles)
decref=dblarr(nfiles)
v2ref=dblarr(nfiles)
v3ref=dblarr(nfiles)

; Loop over inputs reading dither positions
for i=0,nfiles-1 do begin
  hdr=headfits(files[i])
  
  if (i eq 0) then ny=fix(sxpar(hdr,'NAXIS2'))

  raref[i]=sxpar(hdr,'RA_REF')
  decref[i]=sxpar(hdr,'DEC_REF')
  v2ref[i]=sxpar(hdr,'V2_REF')
  v3ref[i]=sxpar(hdr,'V3_REF')
endfor

print,raref,decref

if ((band eq '1A')or(band eq '1B')or(band eq '1C')) then begin
  pwidth=0.196; pixel size along alpha in arcsec
  swidth=0.176; slice width in arcsec
  xmin=8; Minimum x pixel
  xmax=509; Maximum x pixel

  ; Output cube parameters
  rlim_arcsec=0.15; in arcseconds
  ps_x=0.1; arcsec
  ps_y=0.1; arcsec
  ps_z=0.002; micron
endif
if ((band eq '2A')or(band eq '2B')or(band eq '2C')) then begin
  pwidth=0.196; pixel size along alpha in arcsec
  swidth=0.277; slice width in arcsec
  xmin=520; Minimum x pixel
  xmax=1020; Maximum x pixel

  ; Output cube parameters
  rlim_arcsec=0.15; in arcseconds
  ps_x=0.1; arcsec
  ps_y=0.1; arcsec
  ps_z=0.002; micron
endif
if ((band eq '3A')or(band eq '3B')or(band eq '3C')) then begin
  pwidth=0.245; pixel size along alpha in arcsec
  swidth=0.387; slice width in arcsec
  xmin=510; Minimum x pixel
  xmax=1025; Maximum x pixel

  ; Output cube parameters
  rlim_arcsec=0.1; in arcseconds
  ps_x=0.1; arcsec
  ps_y=0.1; arcsec
  ps_z=0.002; micron
endif
if ((band eq '4A')or(band eq '4B')or(band eq '4C')) then begin
  pwidth=0.273; pixel size along alpha in arcsec
  swidth=0.645; slice width in arcsec
  xmin=14; Minimum x pixel
  xmax=480; Maximum x pixel

  ; Output cube parameters
  rlim_arcsec=0.4; in arcseconds
  ps_x=0.2; arcsec
  ps_y=0.2; arcsec
  ps_z=0.002; micron
endif

nx=xmax-xmin+1

; Define base x and y pixel number
basex=rebin(findgen(nx)+xmin,[nx,ny])
basey=transpose(rebin(findgen(ny),[ny,nx]))

; Convert to base alpha,beta,lambda
mmrs_xytoabl,basex,basey,basealpha,basebeta,baselambda,band,slicenum=slicenum

; Crop to only pixels on a real slice
index0=where(slicenum gt 0,nindex0)
slicenum=slicenum[index0]
basex=basex[index0]
basey=basey[index0]
basebeta=basebeta[index0]
basealpha=basealpha[index0]
baselambda=baselambda[index0]

npix=n_elements(index0)

; Convert all alpha,beta base locations to v2,v3 base locations
mmrs_abtov2v3,basealpha,basebeta,basev2,basev3,band,xan=xan,yan=yan

; Create a master vector of fluxes and v2,v3 locations
master_flux=fltarr(nindex0*nfiles)
master_ra=fltarr(nindex0*nfiles)
master_dec=fltarr(nindex0*nfiles)
master_lam=fltarr(nindex0*nfiles)
master_expnum=intarr(nindex0*nfiles)

; Loop over input files reading them into master vectors
for i=0,nfiles-1 do begin
  thisimg=(mrdfits(files[i],0,hdr))[*,*,0]
  ; Crop to correct 1/2 of detector
  thisflux=thisimg[xmin:xmax,*]
  ; Crop to only pixels on real slices
  thisflux=thisflux[index0]
  master_flux[i*nindex0:(i+1)*nindex0-1]=thisflux
  ; Coordinate transform
  jwst_v2v3toradec,basev2,basev3,ra,dec,hdr=hdr,/local

  master_ra[i*nindex0:(i+1)*nindex0-1]=ra
  master_dec[i*nindex0:(i+1)*nindex0-1]=dec
  master_lam[i*nindex0:(i+1)*nindex0-1]=baselambda
  master_expnum[i*nindex0:(i+1)*nindex0-1]=i
endfor

; Safety case; deal with 0-360 range to ensure no problems
; around ra=0 with coordinate wraparound
medra=median(master_ra)
wrapind=where(abs(master_ra - medra) gt 180.,nwrap)
if ((nwrap ne 0)and(medra lt 180.)) then master_ra[wrapind]=master_ra[wrapind]-360.
if ((nwrap ne 0)and(medra ge 180.)) then master_ra[wrapind]=master_ra[wrapind]+360.

; Trim to eliminate any nan fluxes
index1=where(finite(master_flux) eq 1)
master_flux=master_flux[index1]
master_ra=master_ra[index1]
master_dec=master_dec[index1]
master_lam=master_lam[index1]
master_expnum=master_expnum[index1]

; Reference location for output WCS
racen=median(master_ra)
decen=median(master_dec)

cube_x=3600.*(master_ra-min(master_ra))*cos(decen*!PI/180.)/ps_x ; X output cube location in pixels
cube_xsize=fix(max(cube_x)*1.1); 10% oversized in X
cube_x=cube_x+0.05*cube_xsize
xcen=3600.*(racen-min(master_ra))*cos(decen*!PI/180.)/ps_x+0.05*cube_xsize+1; 1-indexed

cube_y=3600.*(master_dec-min(master_dec))/ps_y ; Y output cube location in pixels
cube_ysize=fix(max(cube_y)*1.1); 10% oversized in Y
cube_y=cube_y+0.05*cube_ysize
ycen=3600.*(decen-min(master_dec))/ps_y+0.05*cube_ysize+1; 1-indexed

cube_z=(master_lam-min(master_lam))/ps_z ; Z output cube location in pixels
cube_zsize=fix(max(cube_z))

; squash factors
xpsf_arcsec=0.3
ypsf_arcsec=0.24
zpsf_micron=1; NOT DONE
xpsf=1.;0.3;xpsf_arcsec/ps_x
ypsf=1.0;.24;ypsf_arcsec/ps_y
; Note- arbitrarily squashing x and y too much effectively
; downranks the importance of the z distance; decreasing spectral
; resolution to improve spatial resolution.

; roi
rlimx=rlim_arcsec/ps_x; in pixels
rlimy=rlimx; in pixels
rlimz=2.0
; (Gives about 1-2 spec elements at each spatial element)
rlim=[rlimx,rlimy,rlimz]; in pixels


; Scale correction factor is the ratio of area between an input pixel
; (in arcsec^2) and the output pixel size in arcsec^2
; The result means that the cube will be in calibrated units/pixel
scale=ps_x*ps_y/(pwidth*swidth)

; Just make images at slice 50
slice=125
im=mmrs_cube(cube_x,cube_y,cube_z,master_flux,master_expnum,[cube_xsize,cube_ysize,cube_zsize],rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,slice=slice,wtype=2,expsig=expsig)
;expsig=0.04/ps_x

; Recover gaussian FWHM
imfit=gauss2dfit(im,coeff)
fwhmx=round(coeff[2]*2.355*ps_x*1e3)/1e3
fwhmy=round(coeff[3]*2.355*ps_y*1e3)/1e3
print,'Wavelength (micron): ',slice*ps_z+min(baselambda)
print,'X FWHM (arcsec): ',fwhmx
print,'Y FWHM (arcsec): ',fwhmy

if (~keyword_set(imonly)) then begin
  cube=mmrs_cube(cube_x,cube_y,cube_z,master_flux,master_expnum,[cube_xsize,cube_ysize,cube_zsize],rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,wtype=2,expsig=expsig)

  ; Create header
  mkhdr,cubehdr,cube
  sxaddpar, cubehdr, 'CRPIX1', xcen, 'Reference pixel (1-indexed)'
  sxaddpar, cubehdr, 'CRPIX2', ycen, 'Reference pixel (1-indexed)'
  sxaddpar, cubehdr, 'CRVAL1', racen
  sxaddpar, cubehdr, 'CRVAL2', decen
  sxaddpar, cubehdr, 'CD1_1', -2.77778e-4*ps_x
  sxaddpar, cubehdr, 'CD2_2', 2.77778e-4*ps_y
  sxaddpar, cubehdr, 'CTYPE1', 'RA---TAN'
  sxaddpar, cubehdr, 'CTYPE2', 'DEC--TAN'
  sxaddpar, cubehdr, 'CUNIT1', 'deg'
  sxaddpar, cubehdr, 'CUNIT2', 'deg'

  writefits,outfile,cube,cubehdr
endif

outfile_coll=concat_dir(outdir,'collapse.fits')
collapse=median(cube[*,*,30:cube_zsize-30],dimension=3)
writefits,outfile_coll,collapse

imfit=gauss2dfit(collapse,coeff)
fwhmx=round(coeff[2]*2.355*ps_x*1e3)/1e3
fwhmy=round(coeff[3]*2.355*ps_y*1e3)/1e3
print,'Wavelength (micron): ',slice*ps_z+min(baselambda)
print,'X FWHM (arcsec): ',fwhmx
print,'Y FWHM (arcsec): ',fwhmy



;stop

; Recover gaussian FWHM
;imfit=gauss2dfit(im,coeff)
;fwhmx=round(coeff[2]*2.355*ps_x*1e3)/1e3
;fwhmy=round(coeff[3]*2.355*ps_y*1e3)/1e3
;print,fwhmx,fwhmy
;stop

; Iterate over noised realizations relative to initial
;nrand=30
;tempx=fltarr(nrand)
;tempy=fltarr(nrand)
;test2=ml_meanclip(master_flux,a,b)
;for i=0,nrand-1 do begin
;print,i

 ; tempim=mmrs_cube(cube_x,cube_y,cube_z,master_flux+2*randomn(float(i),n_elements(master_flux))*b,master_expnum,[cube_xsize,cube_ysize,cube_zsize],rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,slice=50,wtype=2,expsig=expsig)

 ; tempimfit=gauss2dfit(tempim,tcoeff)
 ; tempx[i]=tcoeff[2]*2.355*ps_x
 ; tempy[i]=tcoeff[3]*2.355*ps_y
;endfor
;junk=ml_meanclip(tempx,junk1,xsig)
;junk=ml_meanclip(tempy,junk2,ysig)
;xsig=round(xsig*1e3)/1e3
;ysig=round(ysig*1e3)/1e3

;print,'X FWHM (arcsec): ',fwhmx,' +- ',xsig
;print,'Y FWHM (arcsec): ',fwhmy,' +- ',ysig

;stop

return
end

