pro mirim_v2v3toideal,v2,v3,xideal,yideal

; Per Colin Cox:
;V2 = V2Ref + VIdlParity*XIdl*cos(a) + YIdl*sin(a)
;V3 = V3Ref - VIdlParity*XIdl*sin(a) + YIdl*cos(a)
; a being the 4.449705 degrees when rotation is evaluated
; at the middle of the detector

; Parity
P=-1

; Reference point for the Imager
xref=693.5; SIAF convention
yref=512.5; SIAF convention
mirim_xytov2v3,xref-1,yref-1,v2ref,v3ref,'F770W'
print,'v2ref=',v2ref
print,'v3ref=',v3ref

; Derive the rotation between Imager and V2V3
mirim_xytov2v3,516.-1.,512.-1.,v2_0,v3_0,'F770W'
mirim_xytov2v3,516.-1.,512.-1.+1.,v2_1,v3_1,'F770W'
dv2=v2_1-v2_0
dv3=v3_1-v3_0
a=atan(-dv3,dv2)+!PI/2.

xideal=P*((v2-v2ref)-(v3-v3ref)*tan(a))/(cos(a)+sin(a)*tan(a))
yideal=((v2-v2ref)-P*xideal*cos(a))/sin(a)

return
end
