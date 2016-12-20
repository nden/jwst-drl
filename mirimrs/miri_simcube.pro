; Dec 2016: Update for cube testing with Jane.
; Now works by default on Lvl2b processed data, but
; can work on ramps data too with /rampdata.

; stopx/stopy refer to locations in the output image and will
; optionally stop the code here for debugging

pro miri_simcube,directory,band,imonly=imonly,rampdata=rampdata,slice=slice,stopx=stopx,stopy=stopy

subband=strupcase(strmid(band,1,1));A,B, or C
if (subband eq 'A') then subband_name='SHORT'
if (subband eq 'B') then subband_name='MEDIUM'
if (subband eq 'C') then subband_name='LONG'

if ((subband ne 'A')and(subband ne 'B')and(subband ne 'C')) then begin
  print,'Subband not known!'
  exit
endif

; Log runtime
stime0 = systime(1)
; Memory tracking
thismem = memory()
maxmem = 0

outdir=concat_dir(directory,'stack/')
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir
outcube=concat_dir(outdir,'cube.fits')
outslice=concat_dir(outdir,'slice.fits')
outcollapse=concat_dir(outdir,'collapse.fits')

; If /rampdata is set, assume that the input files are ramps
; otherwise assume they have been reduced by JWST pipeline through the
; end of CALSPEC2
if (keyword_set(rampdata)) then begin
  indir=concat_dir(directory,'preproc/')

  if ((band eq '1A')or(band eq '2A')or(band eq '1B')or(band eq '2B')or(band eq '1C')or(band eq '2C')) then begin
    files=file_search(concat_dir(indir,'*MIRIFUSHORT*'))
  endif else begin
    files=file_search(concat_dir(indir,'*MIRIFULONG*'))
  endelse
endif else begin
  filter=strcompress(directory+'*.fits',/remove_all)
  files = dialog_pickfile( title='Read Files to Process', $
                                     filter=filter, $
                                     get_path=new_path, $
                                     /MUST_EXIST        , $
                                     /MULTIPLE_FILES)
endelse

nfiles=n_elements(files)

raref=dblarr(nfiles)
decref=dblarr(nfiles)
v2ref=dblarr(nfiles)
v3ref=dblarr(nfiles)

; Loop over inputs reading dither positions
for i=0,nfiles-1 do begin
  hdr=headfits(files[i])
  hdr1=headfits(files[i],exten=1)
  
  if (i eq 0) then ny=fix(sxpar(hdr1,'NAXIS2'))

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
  rlimz_mic=0.004;
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

; Ramps data are not pixel area corrected, while Lvl2b
; data are
parea=1.0
if (keyword_set(rampdata)) then parea=pwidth*swidth

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
; Create a master vector of detector pixel locations (debugging only)
master_detx=intarr(nindex0*nfiles)
master_dety=intarr(nindex0*nfiles)

; Loop over input files reading them into master vectors
for i=0,nfiles-1 do begin
  hdr=headfits(files[i])
  raw=mrdfits(files[i],1)

  ; If these were ramps data, super-crude slopes conversion
  if (keyword_set(rampdata)) then begin
    ; Super-crude ramps to slopes
    ndim=size(raw,/dim)
    ; Take difference of first and last frames
    lf=raw[*,*,ndim[2]-1,*]-raw[*,*,0,*]
    ; Median across each integration to reject cosmics
    thisimg=median(lf,dimension=4)
  endif else thisimg=raw

  ; Check that this is correct band
  thisband=strtrim(fxpar(hdr,'BAND'),2)
  if (thisband ne subband_name) then begin
    print,'Band mismatch: input file is not band ',subband_name
    exit
  endif

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
  master_detx[i*nindex0:(i+1)*nindex0-1]=basex
  master_dety[i*nindex0:(i+1)*nindex0-1]=basey
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
master_detx=master_detx[index1]
master_dety=master_dety[index1]

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
cube_wavesol=findgen(cube_zsize)*ps_z+min(master_lam)

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
rlimz=rlimz_mic/ps_z
; (Gives about 1-2 spec elements at each spatial element)
rlim=[rlimx,rlimy,rlimz]; in pixels
print,'rlim=',rlim

; Scale correction factor is the ratio of area between an input pixel
; (in arcsec^2) and the output pixel size in arcsec^2
; The result means that the cube will be in calibrated units/pixel
;scale=ps_x*ps_y/(parea)
; Edit: change to output in units of per arcsec^2
scale=1.0/parea
print,'scale=',scale

if ((keyword_set(slice))and(keyword_set(stopx))and(keyword_set(stopy))) then $
  print,'Will stop at ',stopx,' ',stopy,' ',cube_wavesol[slice],' microns'

; Make images at specified slice
if (~keyword_set(slice)) then slice=50
if (~keyword_set(stopx)) then stopx=-1
if (~keyword_set(stopy)) then stopy=-1
im=mmrs_cube(cube_x,cube_y,cube_z,master_flux,master_expnum,[cube_xsize,cube_ysize,cube_zsize],rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,slice=slice,wtype=2,expsig=expsig,detx=master_detx,dety=master_dety,stopx=stopx,stopy=stopy)

; Recover gaussian FWHM
imfit=gauss2dfit(im,coeff)
fwhmx=round(coeff[2]*2.355*ps_x*1e3)/1e3
fwhmy=round(coeff[3]*2.355*ps_y*1e3)/1e3
print,'Wavelength (micron): ',slice*ps_z+min(baselambda)
print,'X FWHM (arcsec): ',fwhmx
print,'Y FWHM (arcsec): ',fwhmy

; Write file
mkhdr,imhdr,im
cdarr=dblarr(2,2)
cdarr[0,0]=2.77778e-4*ps_x
cdarr[1,1]=2.77778e-4*ps_y
make_astr,astr,cd=cdarr,crpix=[xcen,ycen],crval=[racen,decen]
putast,imhdr,astr
writefits,outslice,im,imhdr

if (~keyword_set(imonly)) then begin
  cube=mmrs_cube(cube_x,cube_y,cube_z,master_flux,master_expnum,[cube_xsize,cube_ysize,cube_zsize],rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,wtype=2,expsig=expsig)
  ; Write file
  mkhdr,cubehdr,cube
  putast,cubehdr,astr
  writefits,outcube,cube,cubehdr

  collapse=median(cube[*,*,30:cube_zsize-30],dimension=3)
  ; Write file
  mkhdr,collhdr,collapse
  putast,collhdr,astr
  writefits,outcollapse,collapse,collhdr
endif


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

