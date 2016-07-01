; Take the output of a mirisim simulation and convert it to a format
; useful for building cubes.  This means adding some WCS keywords
; and doing a very rough ramp-slope correction.

pro mmrs_mirisim_preprocess,directory, dark_dir=dark_dir

; directory is a local directory path, inside which it will look
; for det_images/*.fits and convert each of those files
indir=concat_dir(directory,'det_images/')
infiles=file_search(indir,'*.fits')
nfiles=n_elements(infiles)
outfiles=infiles
for i=0,nfiles-1 do $
  outfiles[i]=str_replace(infiles[i],'det_image_','preproc_')

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
    dstart=fix((strsplit(line,'Starting index of executed dither pattern:',/extract,/regex))[1])
    break
  endif
endwhile
close,lun
free_lun,lun

stop
for i=0,nfiles-1 do begin

raw=mrdfits(file,1,hdr)
thisdet=fxpar(hdr,'DETECTOR')
thisband=fxpar(hdr,'BAND')
ditherno=fxpar(hdr,'PNTG_SEQ')



ndim=size(raw,/dim)
; Take only the last frame
lf=raw[*,*,ndim[2]-1,*]
; Median across each integration to reject cosmics
frame=median(lf,dimension=4)

; Subtract off a reference dark
endfor

return
end

; ignore integrateion #1
