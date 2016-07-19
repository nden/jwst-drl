; Preprocess Lvl2b CV3 MRS data
; to add dither information.
; Also do background subtraction.

pro mmrs_cv3_preprocess,directory,outdir=outdir

; directory is a local directory path inside which it will look
; for LVL2b files
rootdir=directory

; Default output directory
if (~keyword_set(outdir)) then outdir=concat_dir(indir,'Converted')
; Ensure output directory exists
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir

; Pretend we're pointing at RA=45 degrees, dec=0 degrees
; with local roll 0
razp=45.D
deczp=0.D
ROLLREF=0.D

; Define V2REF, V3REF at middle of 1A field (alpha=beta=0.0)
V2REF = -8.3942412d ; In arcmin
V3REF = -5.3123744d ; In arcmin

; Only care about ch1a, find those files and their background
; exposures
files=file_search(concat_dir(rootdir,'*-SHORT-*'))
; background exposures
bfiles=file_search(concat_dir(rootdir,'*-SHORTB-*'))
nfiles=n_elements(files)

; Sort based on Q number order
qnum=intarr(nfiles)
for i=0,nfiles-1 do begin
  qname=(strsplit(files[i],'-',/extract))[1]
  qnum[i]=reform(fix(strsplit(qname,'Q',/extract)))
endfor
theorder=sort(qnum)
files=files[theorder]
bfiles=bfiles[theorder]

; Read in and write out files
for i=0,nfiles-1 do begin
  outfile=concat_dir(outdir,fileandpath(files[i]))

  image=readfits(files[i],hdr)
  imageb=readfits(bfiles[i],hdrb)

  ; Subtract backround and zero bad pixels
  imtemp=image[*,*,0]-imageb[*,*,0]
  sigtemp=image[*,*,1]
  mask=image[*,*,2]
  index=where(mask ne 0,nindex)
  if (nindex ne 0) then begin
    imtemp[index]=0.
    sigtemp[index]=0.
  endif
  image[*,*,0]=imtemp
  image[*,*,1]=sigtemp

  ; Derive the offsets from the XACTPOS, YACTPOS in the headers
  osim_xact=fxpar(hdr,'XACTPOS')
  osim_yact=fxpar(hdr,'YACTPOS')
  ; Convert to V2,V3 positions of the OSIM point source
  v2=osim_xact
  v3=-(osim_yact+7.8)
  ; Calculate the deltav2 and deltav3 relative to the reference point
  dv2=(v2-V2REF)/60.*cos(V3REF/60.*!PI/180.); Arc offset in degrees
  dv3=(v3-V3REF)/60.; Arc offset in degrees
  dra=dv2*cos(ROLLREF*!PI/180.)+dv3*sin(ROLLREF*!PI/180.)
  ddec=-dv2*sin(ROLLREF*!PI/180.)+dv3*cos(ROLLREF*!PI/180.)
  ; Figure out corresponding RA,DEC location
  ; (subtract because we're pretending that the telescope moved
  ; instead of the point source?)
  RAREF=razp-dra/cos(deczp*!PI/180.)
  DECREF=deczp-ddec

  ; Add dither offset keywords
  fxaddpar,hdr,'V2_REF',V2REF
  fxaddpar,hdr,'V3_REF',V3REF
  fxaddpar,hdr,'ROLL_REF',ROLLREF
  fxaddpar,hdr,'RA_REF',RAREF
  fxaddpar,hdr,'DEC_REF',DECREF

  print,RAREF,DECREF,v2,v3

  ; Write out file
  writefits,outfile,image,hdr
endfor


return
end
