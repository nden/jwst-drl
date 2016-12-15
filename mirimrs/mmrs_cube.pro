; REQUIRES rlim[0]=rlim[1]

; rlim is 3-element vector
function mmrs_cube, x, y, z, f, expnum, dim_out, rlim, scale=scale, xsquash=xsquash, ysquash=ysquash, zsquash=zsquash, maskcube=maskcube,slice=slice, wtype=wtype, expsig=expsig

; wtype:
; 1: 1/d weighting
; 2: 1/d^2 weighting
; 3: gaussian weighting (define expsig)
if (~keyword_set(wtype)) then wtype=1

if (~keyword_set(xsquash)) then xsquash=1.
if (~keyword_set(ysquash)) then ysquash=1.
if (~keyword_set(zsquash)) then zsquash=1.
if (~keyword_set(scale)) then scale=1.



; If the slice keyword is set, make only a single slice
; and override the nominal output dimensions (ensure we
; don't pass this backwards)
thisdim_out=dim_out
if (keyword_set(slice)) then begin
  if ((slice ge 0)and(slice lt thisdim_out[2])) then begin
    thisdim_out[2]=1
    arr_zcoord=slice
  endif
endif

; Dimensions
ntot=n_elements(f) ; Number of total samples
all_dim=[thisdim_out,ntot] ; Output X,Y,Z size, tot sample size

; Output cubes
fcube=dblarr(thisdim_out)
maskcube=intarr(thisdim_out)
maskcube[*]=1

; XYZ output pixel coordinate arrays
arr_xcoord = dindgen(thisdim_out[0])
arr_ycoord = dindgen(thisdim_out[1])
; Don't overwrite arr_zcoord it if was already set to one slice
if (~keyword_set(arr_zcoord)) then $
  arr_zcoord = dindgen(thisdim_out[2])

;define epsilon
eps=1e-2

; Loop over output cube building the cube slice by slice
for k=0,thisdim_out[2]-1 do begin
  ; Print a message every 5% of the loop
  if (k mod fix(thisdim_out[2]/20) eq 0) then $
    print,'Constructing cube: '+strcompress(string(round(k*100./thisdim_out[2])),/remove_all)+'% complete'

  ; First pass cut: trim to only stuff within rlim of this z location
  indexk=where(abs(z-arr_zcoord[k]) le rlim[2],nindexk)

  tempx=x[indexk]
  tempy=y[indexk]
  tempz=z[indexk]
  tempf=f[indexk]
  tempenum=expnum[indexk]

; QA plot
plot,tempx,tempy,psym=1
circxcen=median(tempx)
circycen=median(tempy)
circphi=findgen(360)/180.*!PI
circx=circxcen+rlim[0]*cos(circphi)
circy=circycen+rlim[1]*sin(circphi)
oplot,circx,circy,color=250
;stop


  ; Loop over output image, building the image row by row
  for j=0,thisdim_out[1]-1 do begin
;print,j
    ; Second pass cut: trim to only stuff within rlim of this y location
    indexj=where(abs(tempy-arr_ycoord[j]) le rlim[1],nindexj)

    ; If nothing makes the cut, do nothing.  Otherwise
    ; build the row
    if (nindexj gt 0) then begin
      tempx2=tempx[indexj]
      tempy2=tempy[indexj]
      tempz2=tempz[indexj]
      tempf2=tempf[indexj]

      ; Now do a 1d build within this slice, looping over input points
      arr_weights=dblarr(thisdim_out[0],nindexj)

      for q=0L,nindexj-1 do begin
        arr_radius=replicate(rlim[0]+1,thisdim_out[0])
        arr_sradius=replicate(rlim[0]+1,thisdim_out[0])

        ; Which output pixels are affected by input points, i.e.
        ; within rlim of this x location?
        ; Don't go outside output array boundaries
        xmin=floor(tempx2[q]-rlim[0]) > 0
        xmax=ceil(tempx2[q]+rlim[0]) < (thisdim_out[0]-1)

        ; Number of points within box
        nbox=[n_elements(arr_xcoord[xmin:xmax])]

        ; Calculate physical spatial radius for ROI determination
        rx=arr_xcoord[xmin:xmax]-tempx2[q]
        ry=arr_ycoord[j]-tempy2[q]
        arr_radius[xmin:xmax]=sqrt(rx^2+replicate(ry^2,nbox))

        ; Determine points within the final circular ROI
        tocalc=where(arr_radius le rlim[0],ncalc)

        ; Squashed radii for weights
        srx=rx/xsquash
        sry=ry/ysquash
        srz=(arr_zcoord[k]-tempz2[q])/zsquash


        ; Combine normalized radii inside ROI
        arr_sradius[xmin:xmax] = sqrt( $
          (srx^2) + $
          replicate(sry^2,nbox) + $
          replicate(srz^2,nbox)  )

;stop

        ; Ensure no divide by zero
        if (ncalc gt 0) then begin

;if ((j eq 20)and(q eq nindexj/2) )then stop
          if (wtype eq 0) then $
            arr_weights[tocalc+q*thisdim_out[0]]=1.
          if (wtype eq 1) then $
            arr_weights[tocalc+q*thisdim_out[0]]=1./(arr_sradius[tocalc] > eps)
          if (wtype eq 2) then $
            arr_weights[tocalc+q*thisdim_out[0]]=1./(arr_sradius[tocalc]^2 > eps)
          if (wtype eq 3) then $
            arr_weights[tocalc+q*thisdim_out[0]]=exp(-0.5/expsig^2*arr_sradius[tocalc]^2)
          if (wtype eq 4) then $
            arr_weights[tocalc+q*thisdim_out[0]]=1./(arr_sradius[tocalc]^4 > eps)
        endif
      endfor

      ; Normalization matrix
      if (nindexj eq 1) then matr_norm=arr_weights $
      else matr_norm = total(arr_weights,2)

      ; Flag where the normalization matrix is zero; there is no good data here
      nodata=where(matr_norm eq 0,complement=gooddata,ncomplement=ngood,nnodata)
      ; Mark good locations in output mask
      if (ngood ne 0) then maskcube[gooddata,j,k]=0

      ; We don't want to divide by zero where there is no data; set the normalization
      ; matrix to 1 in these cases
      if (nnodata gt 0) then matr_norm[nodata]=1.

      ; Apply the weights to calculate the output flux in this row
      frow=dblarr(thisdim_out[0])
      for q=0L,nindexj-1 do begin
        alpha=arr_weights[*,q]/matr_norm
        frow+=tempf2[q]*alpha
      endfor
;if (j eq 29) then stop
      ; Put the row into the final cube
      fcube[*,j,k]=frow*scale
    endif

 endfor

endfor

return,fcube
end
