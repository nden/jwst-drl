function mmrs_cube, x, y, z, f, dim_out, rlim, sigma, scale=scale, xsquash=xsquash, ysquash=ysquash, zsquash=zsquash

if (~keyword_set(xsquash)) then xsquash=1.
if (~keyword_set(ysquash)) then ysquash=1.
if (~keyword_set(zsquash)) then zsquash=1.
if (~keyword_set(scale)) then scale=1.

; Dimensions
ntot=n_elements(f) ; Number of total samples
all_dim=[dim_out,ntot] ; Output X,Y,Z size, tot sample size

; Output cubes
fcube=dblarr(dim_out)

; XYZ output pixel coordinate arrays
arr_xcoord = dindgen(dim_out[0])
arr_ycoord = dindgen(dim_out[1])
arr_zcoord = dindgen(dim_out[2])

;define epsilon
eps=1e-2

; Loop over output cube building the cube slice by slice
; test starting at wave=20
for k=0,dim_out[2]-1 do begin
  ; Print a message every 5% of the loop
  if (k mod fix(dim_out[2]/20) eq 0) then $
    print,'Constructing cube: '+strcompress(string(round(k*100./dim_out[2])),/remove_all)+'% complete'

  ; First pass cut: trim to only stuff within rlim of this z location
  indexk=where(abs(z-arr_zcoord[k])/zsquash le rlim,nindexk)

  tempx=x[indexk]
  tempy=y[indexk]
  tempz=z[indexk]
  tempf=f[indexk]

  ; Loop over output image, building the image row by row
  for j=0,dim_out[1]-1 do begin
    ; Second pass cut: trim to only stuff within rlim of this y location
    indexj=where(abs(tempy-arr_ycoord[j])/ysquash le rlim,nindexj)

    ; If nothing makes the cut, do nothing.  Otherwise
    ; build the row
    if (nindexj gt 0) then begin
      tempx2=tempx[indexj]
      tempy2=tempy[indexj]
      tempz2=tempz[indexj]
      tempf2=tempf[indexj]

      ; Now do a 1d build within this slice, looping over input points
      arr_weights=dblarr(dim_out[0],nindexj)
      for q=0L,nindexj-1 do begin
        arr_radius=replicate(rlim+1,dim_out[0])

        xmin=floor(tempx2[q]-rlim*xsquash) > 0
        xmax=ceil(tempx2[q]+rlim*xsquash) < (dim_out[0]-1)

        dim_c=[n_elements(arr_xcoord[xmin:xmax])]

        rx=(arr_xcoord[xmin:xmax]-tempx2[q])/xsquash
        ry=(arr_ycoord[j]-tempy2[q])/ysquash
        rz=(arr_zcoord[k]-tempz2[q])/zsquash

        arr_radius[xmin:xmax] = sqrt( $
          (rx^2) + $
          replicate(ry^2,dim_c) + $
          replicate(rz^2,dim_c)  )

        ; Calculate weights in the region where radius < rlim
        tocalc = where(arr_radius le rlim,ncalc)
        ; Ensure we don't divide by zero
        if (ncalc gt 0) then $
          arr_weights[tocalc+q*n_elements(arr_radius)]=1./(arr_radius[tocalc] > eps)
      endfor

      ; Normalization matrix
      if (nindexj eq 1) then matr_norm=arr_weights $
      else matr_norm = total(arr_weights,2)

      ; Flag where the normalization matrix is zero; there is no good data here
      nodata=where(matr_norm eq 0,nnodata)
      ; We don't want to divide by zero where there is no data; set the normalization
      ; matrix to 1 in these cases
      if (nnodata gt 0) then matr_norm[nodata]=1.

      ; Apply the weights to calculate the output flux in this row
      frow=dblarr(dim_out[0])
      for q=0L,nindexj-1 do begin
        alpha=arr_weights[*,q]/matr_norm
        frow+=tempf2[q]*alpha
      endfor

      ; Put the row into the final cube
      fcube[*,j,k]=frow*scale
    endif

  endfor

endfor

return,fcube
end
