function mirim_noise,ngroupin,nintin

signalrate=0.4*1*15./17.*1.05; e-/s ; dial this to match what my baseline should be for a given test
ngroup=float(ngroupin)
nint=float(nintin)

det_rn=14
det_dkcurr=0.2
nsample=10; for slow mode

; Each sample is 10 us apart
; MIRI set by nsample (1 or 10, fast or slow), ngroups, nint
; 'frames' and 'groups' interchangeable for miri
tgroup=27.8 ; seconds
tframe=tgroup

; rauscher notation
; m=# of frames/group=1
; n=#groups for integration
mframe=1
bg=(det_dkcurr+signalrate/8.); Signal is divided by 8 to ROUGHLY account for detector footprint



    term1=12.*(ngroup-1)/(mframe*ngroup*(ngroup+1))*det_rn*det_rn
    term2=6.*(ngroup*ngroup+1)/(5*ngroup*(ngroup+1))*(ngroup-1)*tgroup*bg
    term3=2.*(mframe*mframe-1)*(ngroup-1)/(mframe*ngroup*(ngroup+1))*tframe*bg
    noisevec=sqrt(term1+term2-term3)

signal=signalrate*tgroup*ngroup*nint
finalnoise=sqrt(nint)*noisevec
snr=signal/finalnoise
print,signal,finalnoise,snr


return,snr
end
