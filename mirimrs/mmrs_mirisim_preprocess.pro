; Take the output of a mirisim simulation and convert it to a format
; useful for building cubes.  This means adding some WCS keywords
; and doing a very rough ramp-slope correction.

function mmrs_refdark,det,band,dark_dir=dark_dir

if (~keyword_set(dark_dir)) then dark_dir='./darks/'

if ((strtrim(det,2) eq 'MIRIFUSHORT')and(strtrim(band,2) eq 'SHORT')) then darkfile='dark_12A.fits'
if ((strtrim(det,2) eq 'MIRIFUSHORT')and(strtrim(band,2) eq 'MEDIUM')) then darkfile='dark_12B.fits'
if ((strtrim(det,2) eq 'MIRIFUSHORT')and(strtrim(band,2) eq 'LONG')) then darkfile='dark_12C.fits'
if ((strtrim(det,2) eq 'MIRIFULONG')and(strtrim(band,2) eq 'SHORT')) then darkfile='dark_34A.fits'
if ((strtrim(det,2) eq 'MIRIFULONG')and(strtrim(band,2) eq 'MEDIUM')) then darkfile='dark_34B.fits'
if ((strtrim(det,2) eq 'MIRIFULONG')and(strtrim(band,2) eq 'LONG')) then darkfile='dark_34C.fits'

darkfile=concat_dir(dark_dir,darkfile)

dark=readfits(darkfile)

return,dark
end

; Convert a nominal dither offset in alpha,beta
; to offset in ra, dec
pro mmrs_mirisim_dloc,aoff,boff,ROLLREF,raoff,decoff

; Convert alpha,beta offsets to v2, v3
mmrs_abtov2v3,0.,0.,v2ref,v3ref,'1A'
mmrs_abtov2v3,aoff,boff,v2off,v3off,'1A'
v2off=-(v2off-v2ref)*60.
v3off=-(v3off-v3ref)*60.

raoff=v2off*cos(ROLLREF*!PI/180.)+v3off*sin(ROLLREF*!PI/180.)
decoff=-v2off*sin(ROLLREF*!PI/180.)+v3off*cos(ROLLREF*!PI/180.)

;print,aoff,boff,raoff,decoff

return
end

pro mmrs_mirisim_preprocess,directory, dark_dir=dark_dir

; directory is a local directory path, inside which it will look
; for det_images/*.fits and convert each of those files
indir=concat_dir(directory,'det_images/')
outdir=str_replace(indir,'det_images','preproc')
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir

infiles=file_search(indir,'*.fits')
nfiles=n_elements(infiles)
outfiles=infiles
for i=0,nfiles-1 do begin
  outfiles[i]=str_replace(infiles[i],'det_image_','preproc_')
  outfiles[i]=str_replace(outfiles[i],'det_images','preproc')
endfor

ditherfile=concat_dir(getenv('MIRISIM_DIR'),'obssim/data/mrs_recommended_dither.dat')
; Table of possible dither positions (2 x N)
dithertable=(read_ascii(ditherfile)).field1

; Read the log file to pull out key information
logfile=concat_dir(directory,'mirisim.log')
openr,lun,logfile,/get_lun
line=''
while not eof(lun) do begin
  ; Read one line
  readf,lun,line
  ; Look for a line of info about the dither pattern start index (1-indexed)
  a=stregex(line,'Starting index of executed dither pattern:')
  if (a ge 0) then begin
    dstart=fix((strsplit(line,'Starting index of executed dither pattern:',/extract,/regex))[1])-1 ;0-indexed
    break
  endif
endwhile
close,lun
free_lun,lun

; Pretend we're pointing at RA=45 degrees, dec=0 degrees
; with local roll 0
razp=45.D
deczp=0.D
ROLLREF=0.D

; Define V2REF, V3REF at middle of 1A field (alpha=beta=0.0)
V2REF = -8.3942412d ; In arcmin
V3REF = -5.3123744d ; In arcmin


for i=0,nfiles-1 do begin
  hdr=headfits(infiles[i])
  raw=mrdfits(infiles[i],1)
  thisdet=fxpar(hdr,'DETECTOR')
  thisband=fxpar(hdr,'BAND')
  ditherno=dstart+fxpar(hdr,'PNTG_SEQ')-1; Add starting point, convert to 0 indexed

  ndim=size(raw,/dim)
  ; Take only the last frame
  lf=raw[*,*,ndim[2]-1,*]
  ; Median across each integration to reject cosmics
  frame=median(lf,dimension=4)

  ; Subtract off a reference dark
  dark=mmrs_refdark(thisdet,thisband,dark_dir=dark_dir)
  frame=frame-dark

  ; Add dither offset keywords.
  aoffset=dithertable[0,ditherno]; In alpha
  boffset=dithertable[1,ditherno]; In beta
  ; Convert alpha,beta offsets to ra,dec
  mmrs_mirisim_dloc,aoffset,boffset,ROLLREF,dra,ddec
  ; Figure out corresponding RA,DEC location
  ; (subtract because we're pretending that the telescope moved
  ; instead of the point source?)
  RAREF=razp-dra/cos(deczp*!PI/180.)/3600.D
  DECREF=deczp-ddec/3600.D
  fxaddpar,hdr,'V2_REF',V2REF
  fxaddpar,hdr,'V3_REF',V3REF
  fxaddpar,hdr,'ROLL_REF',ROLLREF
  fxaddpar,hdr,'RA_REF',RAREF
  fxaddpar,hdr,'DEC_REF',DECREF

  ; Write out
  writefits,outfiles[i],frame,hdr
endfor

return
end
