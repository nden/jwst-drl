; Oct 2016: Changed sign on RA WCS, was wrong
; Dec 2016: Update for cube testing with Jane
pro miri_cv3cube,imonly=imonly,rampdata=rampdata,slice=slice,stopx=stopx,stopy=stopy

; Log runtime
stime0 = systime(1)
; Memory tracking
thismem = memory()
maxmem = 0

; Output cube parameters
ps_x=0.1; arcsec
ps_y=0.1; arcsec
ps_z=0.002; micron

outdir='./'
outcube=concat_dir(outdir,'cube.fits')
outslice=concat_dir(outdir,'slice.fits')
outcollapse=concat_dir(outdir,'collapse.fits')

; If /rampdata is set, assume that the input files are ramps
; otherwise assume they have been reduced by JWST pipeline through the
; end of CALSPEC2
if (keyword_set(rampdata)) then begin
  ; This hasn't been tested lately
endif else begin
; This case assumes already pipeline processed to calibrated slopes
  files = dialog_pickfile( title='Read Files to Process', $
                                     filter='*.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST        , $
                                     /MULTIPLE_FILES)
  nfiles=n_elements(files)

  ; Assume inputs have been calibrated to specific intensity
  ; (flux/arcsec^2)
  parea=1.0
endelse

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

; Set up rough channel parameters
; Let's work just with channel 1
xmin=0; Minimum x pixel for ch1
xmax=511; Maximum x pixel for ch1
nx=xmax-xmin+1

; Define base x and y pixel number
basex=rebin(findgen(nx)+xmin,[nx,ny])
basey=transpose(rebin(findgen(ny),[ny,nx]))

; Convert to base alpha,beta,lambda
mmrs_xytoabl,basex,basey,basealpha,basebeta,baselambda,'1A',slicenum=slicenum

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
mmrs_abtov2v3,basealpha,basebeta,basev2,basev3,'1A',xan=xan,yan=yan

; Create a master vector of fluxes and v2,v3 locations
master_flux=fltarr(nindex0*nfiles)
master_ra=dblarr(nindex0*nfiles)
master_dec=dblarr(nindex0*nfiles)
master_lam=fltarr(nindex0*nfiles)
master_expnum=intarr(nindex0*nfiles)
master_dq=replicate(long64(0),nindex0*nfiles)
; Extra vectors for debugging
master_detx=fltarr(nindex0*nfiles)
master_dety=fltarr(nindex0*nfiles)
master_v2=dblarr(nindex0*nfiles)
master_v3=dblarr(nindex0*nfiles)

; Loop over input files reading them into master vectors
for i=0,nfiles-1 do begin
  hdr=headfits(files[i])
  thisimg=(mrdfits(files[i],'SCI',hdr1))[*,*,0]
  thisdq=(mrdfits(files[i],'DQ',hdr2))[*,*,0]
  bzero=long64(fxpar(hdr2,'BZERO'))
  ; Crop to correct 1/2 of detector
  thisflux=thisimg[xmin:xmax,*]
  thisdq=thisdq[xmin:xmax,*]+bzero
  ; Crop to only pixels on real slices
  thisflux=thisflux[index0]
  thisdq=thisdq[index0]
  master_flux[i*nindex0:(i+1)*nindex0-1]=thisflux
  master_dq[i*nindex0:(i+1)*nindex0-1]=thisdq
  ; Coordinate transform
  jwst_v2v3toradec,double(basev2),double(basev3),ra,dec,hdr=hdr,/local

  master_ra[i*nindex0:(i+1)*nindex0-1]=ra
  master_dec[i*nindex0:(i+1)*nindex0-1]=dec
  master_lam[i*nindex0:(i+1)*nindex0-1]=baselambda
  master_expnum[i*nindex0:(i+1)*nindex0-1]=i
  master_detx[i*nindex0:(i+1)*nindex0-1]=basex
  master_dety[i*nindex0:(i+1)*nindex0-1]=basey
  master_v2[i*nindex0:(i+1)*nindex0-1]=basev2
  master_v3[i*nindex0:(i+1)*nindex0-1]=basev3
endfor

; Declare minimum wavelength *before* doing any QA cuts for specific exposures
lmin=min(master_lam)
lmax=max(master_lam)

badnum=long64(2)^0+long64(2)^3+long64(2)^9+long64(2)^10+long64(2)^11+long64(2)^14+long64(2)^16
bad=where((master_dq and badnum)ne 0,complement=good)
master_flux=master_flux[good]
master_ra=master_ra[good]
master_dec=master_dec[good]
master_lam=master_lam[good]
master_expnum=master_expnum[good]
master_dq=master_dq[good]
master_detx=master_detx[good]
master_dety=master_dety[good]
master_v2=master_v2[good]
master_v3=master_v3[good]

; Trim to eliminate any nan fluxes
index1=where(finite(master_flux) eq 1,complement=tnan)
master_flux=master_flux[index1]
master_ra=master_ra[index1]
master_dec=master_dec[index1]
master_lam=master_lam[index1]
master_expnum=master_expnum[index1]
master_dq=master_dq[index1]
master_detx=master_detx[index1]
master_dety=master_dety[index1]
master_v2=master_v2[index1]
master_v3=master_v3[index1]

dec_min=min(master_dec)
dec_max=max(master_dec)
dec_ave=(dec_min+dec_max)/2.
ra_min=min(master_ra)
ra_max=max(master_ra)
ra_ave=(ra_min+ra_max)/2.

ra_range=(ra_max-ra_min)*3600.*cos(dec_ave*!PI/180.)
dec_range=(dec_max-dec_min)*3600.
cube_xsize=ceil(ra_range/ps_x)
cube_ysize=ceil(dec_range/ps_y)

xi_min=3600.*(ra_min-ra_ave)*cos(dec_ave*!PI/180.)
xi_max=3600.*(ra_max-ra_ave)*cos(dec_ave*!PI/180.)
eta_min=3600.*(dec_min-dec_ave)
eta_max=3600.*(dec_max-dec_ave)

xi_min=xi_min - ps_x/2.
eta_min=eta_min - ps_y/2.
xi_max=xi_max + ps_x/2.
eta_max=eta_max + ps_y/2.

n1a=ceil(abs(xi_min)/ps_x)
n1b=ceil(abs(xi_max)/ps_x)
n2a=ceil(abs(eta_min)/ps_y)
n2b=ceil(abs(eta_max)/ps_y)
cube_xsize=n1a+n1b
cube_ysize=n2a+n2b

xi=3600.*(master_ra-ra_ave)*cos(dec_ave*!PI/180.)
eta=3600.*(master_dec-dec_ave)
cube_x=(xi-xi_min-ps_x/2.)/ps_x
cube_y=(eta-eta_min-ps_y/2.)/ps_y

racen=ra_ave
decen=dec_ave
xcen=n1a
ycen=n2a


;ramin=44.999011d
;ramax=45.000927d
;demin=-0.00049387922d
;demax=0.00086723190d
;racen=44.9999550191d
;decen=0.000172787451335d
;xcen=34+1
;ycen=24+1
;cube_xsize=70
;cube_ysize=50
;cube_x=3600.*(master_ra-ramin)*cos(decen*!PI/180.)/ps_x
;cube_y=3600.*(master_dec-demin)/ps_y




; Reference location for output WCS
; DRL: Hard override my reference location to match Jane
;racen=44.99995501914714d;median(master_ra)
;decen=0.000172787451335265d;median(master_dec)

;obound=0.05; Oversizing
;janeshift_x=-1.3
;janeshift_y=-1.3

;cube_x=3600.*(master_ra-min(master_ra))*cos(decen*!PI/180.)/ps_x ; X output cube location in pixels
;cube_xsize=fix(max(cube_x)*(1.+obound)); 10% oversized in X
;cube_x=cube_x+obound/2.*cube_xsize+janeshift_x
;xcen=3600.*(racen-min(master_ra))*cos(decen*!PI/180.)/ps_x+obound/2.*cube_xsize+1+janeshift_x; 1-indexed

;cube_y=3600.*(master_dec-min(master_dec))/ps_y ; Y output cube location in pixels
;cube_ysize=fix(max(cube_y)*(1.+obound)); 10% oversized in Y
;cube_y=cube_y+obound/2.*cube_ysize+janeshift_y
;ycen=3600.*(decen-min(master_dec))/ps_y+obound/2.*cube_ysize+1+janeshift_y; 1-indexed

zrange=lmax-lmin
cube_zsize=ceil(zrange/ps_z)
lamcen=(lmax+lmin)/2.
lamstart=lamcen-(cube_zsize/2.)*ps_z
lamstop=lamstart+cube_zsize*ps_z
cube_z=(master_lam-lamstart)/ps_z ; Z output cube location in pixels
wavevec=findgen(cube_zsize)*ps_z+min(master_lam)

; Cull single exposures
;index1=where(master_expnum eq 8)
;master_flux=master_flux[index1]
;master_ra=master_ra[index1]
;master_dec=master_dec[index1]
;master_lam=master_lam[index1]
;master_expnum=master_expnum[index1]
;cube_x=cube_x[index1]
;cube_y=cube_y[index1]
;cube_z=cube_z[index1]

; squash factors
xpsf_arcsec=0.3
ypsf_arcsec=0.24
zpsf_micron=1; NOT DONE
xpsf=1.;0.3;xpsf_arcsec/ps_x
ypsf=1.;0.24;ypsf_arcsec/ps_y
; Note- arbitrarily squashing x and y too much effectively
; downranks the importance of the z distance; decreasing spectral
; resolution to improve spatial resolution.

; roi
rlim_arcsec=0.15; in arcseconds
rlimx=rlim_arcsec/ps_x; in pixels
rlimy=rlimx; in pixels
rlimz=0.004/ps_z
; (Gives about 1-2 spec elements at each spatial element)
rlim=[rlimx,rlimy,rlimz]; in pixels


; Scale correction factor is the ratio of area between an input pixel
; (in arcsec^2) and the output pixel size in arcsec^2
; The result means that the cube will be in calibrated units/pixel
;scale=ps_x*ps_y/(parea)
scale=1.0 ; Output is in flux/solid angle
print,'scale=',scale

if ((keyword_set(slice))and(keyword_set(stopx))and(keyword_set(stopy))) then $
  print,'Will stop at ',stopx,' ',stopy,' ',wavevec[slice],' microns'

; Make images at specified slice
if (~keyword_set(slice)) then slice=50
if (~keyword_set(stopx)) then stopx=-1
if (~keyword_set(stopy)) then stopy=-1
im=mmrs_cube(cube_x,cube_y,cube_z,master_flux,master_expnum,[cube_xsize,cube_ysize,cube_zsize],rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,slice=slice, wtype=2,expsig=expsig,detx=master_detx,dety=master_dety,detlam=master_lam,stopx=stopx,stopy=stopy)

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
  fxaddpar,cubehdr,'CD3_3',ps_z
  fxaddpar,cubehdr,'CDELT3',ps_z
  fxaddpar,cubehdr,'CRPIX3',1
  fxaddpar,cubehdr,'CTYPE3','WAVE'
  fxaddpar,cubehdr,'CRVAL3',lamstart
  fxaddpar,cubehdr,'CUNIT3','um'
  writefits,outcube,cube,cubehdr

  collapse=median(cube[*,*,30:cube_zsize-30],dimension=3)
  ; Write file
  mkhdr,collhdr,collapse
  putast,collhdr,astr
  writefits,outcollapse,collapse,collhdr
endif



; Iterate over noised realizations relative to initial
;nrand=30
;tempx=fltarr(nrand)
;tempy=fltarr(nrand)
;test2=ml_meanclip(master_flux,a,b)
;for i=0,nrand-1 do begin
;print,i

;  tempim=mmrs_cube(cube_x,cube_y,cube_z,master_flux+2*randomn(float(i),n_elements(master_flux))*b,master_expnum,[cube_xsize,cube_ysize,cube_zsize],rlim,xsquash=xpsf,ysquash=ypsf,scale=scale,slice=50,wtype=2,expsig=expsig)

;  tempimfit=gauss2dfit(tempim,tcoeff)
;  tempx[i]=tcoeff[2]*2.355*ps_x
;  tempy[i]=tcoeff[3]*2.355*ps_y
;endfor
;junk=ml_meanclip(tempx,junk1,xsig)
;junk=ml_meanclip(tempy,junk2,ysig)
;xsig=round(xsig*1e3)/1e3
;ysig=round(ysig*1e3)/1e3
;
;print,'X FWHM (arcsec): ',fwhmx,' +- ',xsig
;print,'Y FWHM (arcsec): ',fwhmy,' +- ',ysig

;stop

return
end

