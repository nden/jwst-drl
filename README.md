# jwst

This product contains a bunch of homegrown IDL based routines by D. Law
to test and evaluate coordinate transforms associated with the MIRI
instrument on board JWST.

manga/ contains some convenience functions taken from the IDL code library
of the SDSS MaNGA IFU pipeline.

mirim/ contains routines to apply the MIRI Imager distortion solution to
convert between pixel coordinates and JWST v2,v3 telescope-frame coordinates

mirimrs/ contains routines to apply the MIRI Medium Resolution Spectrometer
distortion solutions to convert between pixel coordinates, local 
alpha/beta coordinates, and JWST v2,v3 telescope-frame coordinates

tel/ contains routines designed to convert between v2,v3 and RA,DEC coordinates
(and vice versa) given suitable attitude matrix information.

In order to use any of these tools, you will need to add some information
to your shell path.  In the case of .cshrc you should define the 
environmental variable pointing to the code location:
setenv JWSTTOOLS_DIR /Users/dlaw/jwst/trunk
And then be sure to add this to your IDL path:
setenv IDL_PATH ${IDL_PATH}:+$JWSTTOOLS_DIR

Additionally, if using this code to build cubes from the mirisim simulator
you should also set the environmental variable:
setenv MIRISIM_DIR /Users/dlaw/anaconda2/envs/mirisim/lib/python2.7/site-packages/mirisim/
pointing to your mirisim install location.  This is necessary in order to find
the dither reference table assumed by mirisim.


Tag 201609 contains code from September 2016, which is pertinent to when
mirisim did not include position information in the header data.