# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Modified by David Law based on code by Nadia Dencheva
"""
Tools to create WCS reference files for MIRI MRS using IDT reference
products delivered with CDP-6:
e.g., MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_06.04.00.fits
MIRI MRS uses 4 reference files of type:
DISTORTION
SPECWCS
REGIONS
WAVELENGTHRANGE
create_cdp6_references() creates all reference files.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os as os
import numpy as np
from asdf import AsdfFile
import pdb as pdb
from astropy.io import fits
from astropy.modeling import models
from astropy import units as u
from jwst import datamodels
from jwst.assign_wcs import miri
from numpy.testing import utils

import mmrs_tools as mmrs_tools
import drltimer as drltimer

from jwst.datamodels import *


# Function to loop over all 6 MIRI MRS distortion files
# making reference files for all of them
# create_cdp6_all('./')
@drltimer.fn_timer
def create_cdp6_all(outdir):
    detbands='12A','12B','12C','34A','34B','34C'
    nbands=len(detbands)
    for i in range(nbands):
        create_cdp6_setfiles(detbands[i],outdir)

# Function to automatically figure out the input/output required to make
# a CDP-6 reference file for a particular detector band (e.g., 12A)
# create_cdp6_setfiles('12A','./')
@drltimer.fn_timer
def create_cdp6_setfiles(detband,outdir):
    def_rootdir=os.path.expandvars('$JWSTTOOLS_DIR')

    if (detband == '12A'):
        fname=def_rootdir+'/cdp/cdp6/MIRI_FM_MIRIFUSHORT_12SHORT_DISTORTION_06.04.00.fits'
    elif (detband == '34A'):
        fname=def_rootdir+'/cdp/cdp6/MIRI_FM_MIRIFULONG_34SHORT_DISTORTION_06.04.00.fits'
    elif (detband == '12B'):
        fname=def_rootdir+'/cdp/cdp6/MIRI_FM_MIRIFUSHORT_12MEDIUM_DISTORTION_06.04.00.fits'
    elif (detband == '34B'):
        fname=def_rootdir+'/cdp/cdp6/MIRI_FM_MIRIFULONG_34MEDIUM_DISTORTION_06.04.00.fits'
    elif (detband == '12C'):
        fname=def_rootdir+'/cdp/cdp6/MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_06.04.00.fits'
    elif (detband == '34C'):
        fname=def_rootdir+'/cdp/cdp6/MIRI_FM_MIRIFULONG_34LONG_DISTORTION_06.04.00.fits'

    distfile=outdir+'jwst_miri_mrs'+detband+'_distortion_cdp6.asdf'
    regfile=outdir+'jwst_miri_mrs'+detband+'_regions_cdp6.asdf'
    specfile=outdir+'jwst_miri_mrs'+detband+'_specwcs_cdp6.asdf'
    wavefile=outdir+'jwst_miri_mrs_wavelengthrange_cdp6.asdf'
    refs={'distortion': distfile, 'regions':regfile, 'specwcs':specfile, 'wavelengthrange':wavefile}
    print('Working on: '+detband)
    create_cdp6_onereference(fname,refs)
    print('Testing: '+detband)
    test_cdp6_onereference(detband,refs)
    print('Done testing: '+detband)

def create_cdp6_onereference(fname, ref):
    """
    Create ASDF WCS reference files for MIRI MRS data from a single CDP6 reference file.
    Parameters
    ----------
    fname : str
        name of reference file
    ref : dict
        A dictionary {reftype: refname}, e.g.
        {'distortion': 'jwst_miri_distortion_0001.asdf',
         'regions': 'jwst_miri_regions_0001.asdf',
         'specwcs': 'jwst_miri_specwcs_0001.asdf',
         'wavelengthrange': 'jwst_miri_wavelengthrange_0001.asdf'}
    Examples
    --------
    >>> create_cdp6_references('MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_06.04.00.fits', ref)
    """
    with fits.open(fname) as f:
        channel = f[0].header['CHANNEL']
        band = f[0].header['BAND']
        detector = f[0].header['DETECTOR']
        ch1 = 'CH{0}'.format(channel[0])
        ch2 = 'CH{0}'.format(channel[1])
        slices = f[1].data
        fov1 = f[2].data
        fov2 = f[3].data
        alpha1 = f[('Alpha_'+ch1, 1)].data
        lam1 = f[('Lambda_'+ch1, 1)].data
        alpha2 = f[('Alpha_'+ch2, 1)].data
        lam2 = f[('Lambda_'+ch2, 1)].data
        x1 = f[('X_'+ch1, 1)].data
        y1 = f[('Y_'+ch1, 1)].data
        x2 = f[('X_'+ch2, 1)].data
        y2 = f[('Y_'+ch2, 1)].data
        ab_v23 = f[('albe_to_XANYAN', 1)].data.copy()
        v23_ab = f[('XANYAN_to_albe', 1)].data.copy()
        b0_ch1 = f[0].header['B_ZERO'+ch1[2]]
        bdel_ch1 = f[0].header['B_DEL'+ch1[2]]
        b0_ch2 = f[0].header['B_ZERO'+ch2[2]]
        bdel_ch2 = f[0].header['B_DEL'+ch2[2]]
    # Get channel names, e.g. 1LONG, 2LONG
    channels = [c + band for c in channel]
    # Note that now 'channel' is (e.g.) 12, while 'channels' is (e.g.) '1SHORT','2SHORT'

    bzero = {}
    bdel = {}
    for c in channel:
        cb = c+band
        bzero[cb] = f[0].header['B_ZERO' + c]
        bdel[cb] = f[0].header['B_DEL' + c]

    # MRS reference files are long enough that keeping tables as inline
    # text is impractical
    outformat='inline'

    coeff_names = build_coeff_names(alpha1.names)
    amodel1 = create_poly_models(alpha1, int(channel[0]), coeff_names, name='det2local')
    lmodel1 = create_poly_models(lam1, int(channel[0]), coeff_names, name='det2local')
    amodel2 = create_poly_models(alpha2, int(channel[1]), coeff_names, name='det2local')
    lmodel2 = create_poly_models(lam2, int(channel[1]), coeff_names, name='det2local')
    # reverse models

    # 'x' in the report corresponds to y in python and 'y' to x,
    # The x/y models take (lam, alpha)
    xmodel1 = create_xy_models(x1, int(channel[0]), coeff_names, name='x')
    ymodel1 = create_xy_models(y1, int(channel[0]), coeff_names, name='y')
    xmodel2 = create_xy_models(x2, int(channel[1]), coeff_names, name='x')
    ymodel2 = create_xy_models(y2, int(channel[1]), coeff_names, name='y')
    amodel1.update(amodel2)
    xmodel1.update(xmodel2)
    ymodel1.update(ymodel2)
    lmodel1.update(lmodel2)

    bmodel1 = create_beta_models(b0_ch1, bdel_ch1, int(channel[0]), len(alpha1))
    bmodel2 = create_beta_models(b0_ch2, bdel_ch2, int(channel[1]), len(alpha2))

    bmodel1.update(bmodel2)
    useafter = "2000-01-01T00:00:00"
    author =  'Adrian M. Glauser, David R. Law'  #Author of the data
    description = 'MIRI MRS CDP6 distortion reference data.'

    create_distortion_file(reftype='distortion', detector=detector, band=band, channel=channel, channels=channels,
                           data=(amodel1, bmodel1, xmodel1, ymodel1, bzero, bdel, ab_v23, v23_ab), name=ref['distortion'],
                           author=author, useafter=useafter, description=description, outformat=outformat)

    create_specwcs_file('specwcs', detector, band, channel, lmodel1, ref['specwcs'], author,
                        useafter, description, outformat)

    create_regions_file(slices, detector, band, channel, ref['regions'], author, useafter,
               description, outformat)

    create_wavelengthrange_file(ref['wavelengthrange'], detector, author, useafter,
                                description, outformat)


def create_regions_file(slices, detector, band, channel, name, author, useafter, description, outformat):
    model = RegionsModel()
    model = create_reffile_header(model, detector, band, channel, author, useafter,
                                  description)
    model.meta.filename = os.path.split(name)[-1]
    model.regions = slices
    #f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    model.save(name)


def create_reffile_header(model, detector, band, channel, author, useafter,
                          description=""):
    model.meta.description = description
    model.meta.author = author
    model.meta.useafter = useafter
    model.meta.pedigree = 'GROUND'
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = detector
    model.meta.instrument.channel = channel
    model.meta.instrument.band = band
    model.meta.exposure.type = "MIR_MRS"

    return model


def create_distortion_file(reftype, detector,  band, channel, channels, data, name, author,
                           useafter, description, outformat):
    dist = DistortionMRSModel()
    description = 'MIRI MRS Distortion Maps'
    dist = create_reffile_header(dist, detector, band, channel, author, useafter,
                                 description)

    dist.meta.filename = os.path.split(name)[-1]
    # Split the provided data vector into its pieces
    adata, bdata, xdata, ydata, bzero, bdel, ab_v23, v23_ab = data

    slices = list(xdata.keys())
    dist.slices = slices
    xd = [xdata[sl] for sl in slices]
    dist.x_model = xd
    yd = [ydata[sl] for sl in slices]
    dist.y_model = yd
    ad = [adata[sl] for sl in slices]
    dist.alpha_model = ad
    bd = [bdata[sl] for sl in slices]
    dist.beta_model = bd
    dist.bzero = {'channel_band': list(bzero.keys()), 'beta_zero': list(bzero.values())}
    dist.bdel = {'channel_band': list(bdel.keys()), 'delta_beta': list(bdel.values())}

    """
    Create the transform from MIRI Local to telescope V2/V3 system for all channels.
    """
    channel = "".join([ch[0] for ch in channels])

    # Read in coefficients from the tables.  Note that we'll flip the coefficient
    # ordering since they were set up for column-major indexing (IDL) but we're working in
    # python (row-major)
    # ab -> v2 transform for first channel
    c0_0, c1_0, c0_1, c1_1 = ab_v23[0][1:]
    ch1_v2 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    # v2,v3 -> a transform for first channel
    c0_0, c1_0, c0_1, c1_1 = v23_ab[0][1:]
    ch1_a = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")
    # ab -> v3 transform for first channel
    c0_0, c1_0, c0_1, c1_1 = ab_v23[1][1:]
    ch1_v3 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    # v2,v3 -> b transform for first channel
    c0_0, c1_0, c0_1, c1_1 = v23_ab[1][1:]
    ch1_b = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")
    # ab -> v2 transform for second channel
    c0_0, c1_0, c0_1, c1_1 = ab_v23[2][1:]
    ch2_v2 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    # v2,v3 -> a transform for second channel
    c0_0, c1_0, c0_1, c1_1 = v23_ab[2][1:]
    ch2_a = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")
    # ab -> v3 transform for second channel
    c0_0, c1_0, c0_1, c1_1 = ab_v23[3][1:]
    ch2_v3 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    # v2,v3 -> b transform for second channel
    c0_0, c1_0, c0_1, c1_1 = v23_ab[3][1:]
    ch2_b = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")

    # Transforms from the CDP only went to XAN,YAN, now need a transform to V2,V3
    # Mapping to transform from XAN,YAN in arcmin to V2,V3 in arcsec
    xanyan_to_v2v3 = models.Identity(1) & (models.Scale(-1) | models.Shift(-7.8)) | models.Scale(60.) & models.Scale(60.)

    # Since the matrix transforms need a 4-element input we need a mapping
    # to go from (0,1) to (0,1,0,1)
    mapping = models.Mapping((0, 1, 0, 1), n_inputs=2)

    # Put the mappings all together
    ch1 = mapping | ch1_v2 & ch1_v3 | xanyan_to_v2v3
    ch2 = mapping | ch2_v2 & ch2_v3 | xanyan_to_v2v3

    # And make the inverse mapping
    ch1.inverse =  xanyan_to_v2v3.inverse | mapping | ch1_a & ch1_b
    ch2.inverse =  xanyan_to_v2v3.inverse | mapping | ch2_a & ch2_b

    #pdb.set_trace()
    # save to file
    dist.abv2v3_model = {'channel_band': channels, 'model': [ch1, ch2]}
    dist.meta.input_units = u.pix
    dist.meta.output_units = u.arcsec
    dist.save(name)


def create_specwcs_file(reftype, detector, band, channel, lmodel, name, author, useafter, description, outformat):
    spec = SpecwcsModel()
    spec = create_reffile_header(spec, detector, band, channel, author, useafter,
                                 description)

    spec.meta.subarray.name = "N/A"
    spec.meta.filename = os.path.split(name)[-1]

    slices = list(lmodel.keys())
    spec.slices = slices
    lam_data = [lmodel[sl] for sl in slices]
    spec.model = lam_data

    #f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    spec.save(name)

# Create the x,y to a,b models
def create_poly_models(data, channel, coeff_names, name):
    """
    Create a 2D polynomial model for the transformation
    detector --> local MIRI frame
    Works for alpha and lambda coordinates.
    """
    nslices = len(data)
    sl = channel * 100 + np.arange(1, nslices+1)

    transforms = {}
    for i in range(nslices):
        sl = channel * 100 + i +1
        al = data[i]
        xs = al[0]
        coeffs = {}
        for c, val in zip(coeff_names, al[1:]):
            coeffs[c] = val

        # First we need to shift from the 0-indexed pixels input by the pipeline to the
        # 1-indexed pixels used by Adrians transforms
        # (both Adrian and the pipeline include reference pixels)
        # Remember that the coordinates here are in order y,x
        thisshift=models.Shift(1) & models.Shift(1)
        # Now do Adrians transforms
        thisxform=models.Identity(1) & models.Shift(-xs) | models.Polynomial2D(8, name=name, **coeffs)
        # Put the models together
        transforms[sl] = thisshift | thisxform

    return transforms

# Create the a,b, to x,y models
def create_xy_models(data, channel, coeff_names, name):
    """
    Create a 2D polynomial model for the transformation
    local_MIRI --> detector frame.
    """
    nslices = len(data)
    sl = channel * 100 + np.arange(1, nslices+1)
    shname = "shift_{0}".format(name)
    pname = "polynomial_{0}".format(name)
    transforms = {}
    for i in range(nslices):
        sl = channel * 100 + i +1
        al = data[i]
        xs = al[0]
        coeffs = {}
        for c, val in zip(coeff_names, al[1:]):
            coeffs[c] = val

        # First we do Adrian's transforms
        thisxform=models.Shift(-xs, name=shname) & models.Identity(1) | models.Polynomial2D(8, name=pname, **coeffs)
        # Then we need to shift from the 1-indexed pixels used by transforms to the
        # 0-indexed pixels used by the pipeline
        # (both Adrian and the pipeline include reference pixels)
        # Only a single output so only a single shift
        thisshift=models.Shift(-1)
        # Put the models together
        transforms[sl] =  thisxform | thisshift

    return transforms


def build_coeff_names(names):
    names = names[1:]
    names = [name.replace('VAR2_', "c") for name in names]
    return names

def create_beta_models(b0, bdel, channel, nslices):
    beta = {}
    for s in range(nslices):
        sl = channel * 100 + s +1
        beta_s = b0 + s * bdel
        m = models.Const1D(beta_s, name='det2local') #xy2beta and xy2lam
        beta[sl] = m
    return beta


def create_wavelengthrange_file(name, detector, author, useafter, description, outformat):
    model = WavelengthrangeModel()

    # Relaxing the range to match the distortion. The table above
    # comes from the report and is "as designed".
    wavelengthrange = {'1SHORT': (4.68, 5.97),
                        '1MEDIUM': (5.24, 6.87),
                        '1LONG': (6.2, 7.90),
                        '2SHORT': (7.27, 9.03),
                        '2MEDIUM': (8.43, 10.39),
                        '2LONG': (9.76, 11.97),
                        '3SHORT': (11.29, 13.75),
                        '3MEDIUM': (13.08, 15.86),
                        '3LONG': (15.14, 18.29),
                        '4SHORT': (17.40, 21.20),
                        '4MEDIUM': (20.31, 24.68),
                        '4LONG': (23.72, 28.75)
                        }
    channels = ['1SHORT', '1MEDIUM', '1LONG', '2SHORT', '2MEDIUM', '2LONG',
                '3SHORT', '3MEDIUM', '3LONG', '4SHORT', '4MEDIUM', '4LONG']

    model = create_reffile_header(model, detector, band="N/A", channel="N/A", author=author,
                                 useafter=useafter, description=description)
    model.meta.filename = os.path.split(name)[-1]
    model.meta.author = ''
    model.meta.instrument.detector = "N/A"
    model.waverange_selector = channels
    wr = [wavelengthrange[ch] for ch in channels]
    model.wavelengthrange = wr
    model.meta.wavelength_units = u.micron
    #f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    model.save(name)

# Function to test the implemented transforms and ASDF files
# Detband is (e.g.) '12A'
def test_cdp6_onereference(detband,refs):
    ch1,ch2=mmrs_tools.channel(detband)# Convert to (e.g.) '1A' and '2A'

    # Do the first channel of the file
    input_model = datamodels.ImageModel()
    # Convert input of type '12A' into the band and channel that pipeline needs
    theband,thechan=mmrs_tools.bandchan(ch1)
    # Set the filter in the data model meta header
    input_model.meta.instrument.band = theband
    input_model.meta.instrument.channel = thechan
    ref_data=mrs_ref_data[ch1]
    ref_x,ref_y=ref_data['x']-1,ref_data['y']-1# Convert to 0-indexed
    ref_a=ref_data['alpha']
    ref_b=ref_data['beta']
    ref_l=ref_data['lam']
    ref_xan=ref_data['xan']
    ref_yan=ref_data['yan']
    ref_v2,ref_v3=mmrs_tools.xanyan_to_v2v3(ref_xan,ref_yan)
    xform1=mmrs_tools.xytoabl(ch1,refs=refs)
    xform2=mmrs_tools.abtov2v3(ch1,refs=refs)
    for i, s in enumerate(ref_data['s']):
      a,b,l=xform1(ref_x[i],ref_y[i])
      v2,v3=xform2(a,b)
      ainv,binv=xform2.inverse(v2,v3)
      xinv,yinv=xform1.inverse(ainv,binv,l)
      # Test a,b agreement to 5e-3
      utils.assert_allclose(ref_a[i], a, atol=5*10**-3)
      utils.assert_allclose(ref_a[i], ainv, atol=5*10**-3)
      utils.assert_allclose(ref_b[i], b, atol=5*10**-3)
      utils.assert_allclose(ref_b[i], binv, atol=5*10**-3)
      # Test v2,v3,l agreement to 5e-3
      utils.assert_allclose(ref_v2[i], v2, atol=5*10**-3)
      utils.assert_allclose(ref_v3[i], v3, atol=5*10**-3)
      utils.assert_allclose(ref_l[i], l, atol=5*10**-3)
      # Test x,y agreement to 5e-2
      utils.assert_allclose(ref_x[i], xinv, atol=5*10**-2)
      utils.assert_allclose(ref_y[i], yinv, atol=5*10**-2)


    # Do the second channel of the file
    input_model = datamodels.ImageModel()
    # Convert input of type '12A' into the band and channel that pipeline needs
    theband,thechan=mmrs_tools.bandchan(ch2)
    # Set the filter in the data model meta header
    input_model.meta.instrument.band = theband
    input_model.meta.instrument.channel = thechan
    ref_data=mrs_ref_data[ch2]
    ref_x,ref_y=ref_data['x']-1,ref_data['y']-1# Convert to 0-indexed
    ref_a=ref_data['alpha']
    ref_b=ref_data['beta']
    ref_l=ref_data['lam']
    ref_xan=ref_data['xan']
    ref_yan=ref_data['yan']
    ref_v2,ref_v3=mmrs_tools.xanyan_to_v2v3(ref_xan,ref_yan)
    xform1=mmrs_tools.xytoabl(ch2,refs=refs)
    xform2=mmrs_tools.abtov2v3(ch2,refs=refs)
    for i, s in enumerate(ref_data['s']):
      a,b,l=xform1(ref_x[i],ref_y[i])
      v2,v3=xform2(a,b)
      ainv,binv=xform2.inverse(v2,v3)
      xinv,yinv=xform1.inverse(ainv,binv,l)
      # Test a,b agreement to 5e-3
      utils.assert_allclose(ref_a[i], a, atol=5*10**-3)
      utils.assert_allclose(ref_a[i], ainv, atol=5*10**-3)
      utils.assert_allclose(ref_b[i], b, atol=5*10**-3)
      utils.assert_allclose(ref_b[i], binv, atol=5*10**-3)
      # Test v2,v3,l agreement to 5e-3
      utils.assert_allclose(ref_v2[i], v2, atol=5*10**-3)
      utils.assert_allclose(ref_v3[i], v3, atol=5*10**-3)
      utils.assert_allclose(ref_l[i], l, atol=5*10**-3)
      # Test x,y agreement to 5e-2
      utils.assert_allclose(ref_x[i], xinv, atol=5*10**-2)
      utils.assert_allclose(ref_y[i], yinv, atol=5*10**-2)

# MRS test reference data
# This stuff is all 1-indexed for x,y from Adrian's report
mrs_ref_data = {
    '1A': {'x': np.array([28.310396, 475.02154, 493.9777, 41.282537, 58.998266]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([11, 1, 1, 21, 21]),
           'alpha': np.array([0, -1.66946, 1.65180, -1.70573, 1.70244]),
           'beta': np.array([0, -1.77210, -1.77210, 1.77210, 1.77210]),
           'lam': np.array([5.34437, 4.86642, 4.95325, 5.65296, 5.74349]),
           'xan': np.array([-8.39424, -8.41746, -8.36306, -8.42653, -8.37026]),
           'yan': np.array([-2.48763, -2.52081, -2.51311, -2.46269, -2.45395]),
           },
    '1B': {'x': np.array([28.648221, 475.07259, 493.98157, 41.559386, 59.738296]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([11, 1, 1, 21, 21]),
           'alpha': np.array([0., -1.70796, 1.60161, -1.70854, 1.78261]),
           'beta': np.array([0., -1.77204, -1.77204, 1.77204, 1.77204]),
           'lam': np.array([6.17572, 5.62345, 5.72380, 6.53231, 6.63698]),
           'xan': np.array([-8.39426, -8.41808, -8.36368, -8.42682, -8.36899]),
           'yan': np.array([-2.48492, -2.51808, -2.51040, -2.46001, -2.45126])
           },
    '1C': {'x': np.array([30.461871, 477.23742, 495.96228, 43.905314, 60.995683]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([11, 1, 1, 21, 21]),
           'alpha': np.array([0., -1.60587, 1.67276, -1.60766, 1.68720]),
           'beta': np.array([0., -1.77202, -1.77202, 1.77202, 1.77202]),
           'lam': np.array([7.04951, 6.42424, 6.53753, 7.45360, 7.57167]),
           'xan': np.array([-8.39357, -8.41570, -8.36165, -8.42457, -8.36996]),
           'yan': np.array([-2.48987, -2.52271, -2.51525, -2.46467, -2.45649])
           },
    '2A': {'x': np.array([992.158, 545.38386, 525.76143, 969.29711, 944.19303]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([9, 1, 1, 17, 17]),
           'alpha': np.array([0., -2.11250, 2.10676, -2.17239, 2.10447]),
           'beta': np.array([0., -2.23775, -2.23775, 2.23775, 2.23775]),
           'lam': np.array([8.20797, 7.52144, 7.64907, 8.68677, 8.83051]),
           'xan': np.array([-8.39393, -8.42259, -8.35355, -8.43583, -8.36499]),
           'yan': np.array([-2.48181, -2.52375, -2.51357, -2.44987, -2.44022])
           },
    '2B': {'x': np.array([988.39977, 541.23447, 521.60207, 964.91753, 940.10325]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([9, 1, 1, 17, 17]),
           'alpha': np.array([0., -2.10593, 2.10015, -2.08817, 2.10422]),
           'beta': np.array([0., -2.23781, -2.23781, 2.23781, 2.23781]),
           'lam': np.array([9.44205, 8.65341, 8.79991, 9.99257, 10.15795]),
           'xan': np.array([-8.39645, -8.42502, -8.35603, -8.43716, -8.36742]),
           'yan': np.array([-2.47773, -2.51972, -2.50938, -2.44554, -2.43626])
           },
    '2C': {'x': np.array([990.89693, 543.82344, 524.34514, 967.98318, 942.77564]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([9, 1, 1, 17, 17]),
           'alpha': np.array([0., -2.07490, 2.11234, -2.14704, 2.14196]),
           'beta': np.array([0., -2.23778, -2.23778, 2.23778, 2.23778]),
           'lam': np.array([10.90225, 9.99162, 10.16079, 11.53780, 11.72887]),
           'xan': np.array([-8.39303, -8.42129, -8.35221, -8.43454, -8.36352]),
           'yan': np.array([-2.47869, -2.52052, -2.51036, -2.44668, -2.43712])
           },
    '3A': {'x': np.array([574.80828, 1001.0602, 984.6387, 547.27479, 518.89992]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([8, 1, 1, 16, 16]),
           'alpha': np.array([0., -2.86745, 3.20982, -3.01230, 2.96643]),
           'beta': np.array([-0.19491, -2.92360, -2.92360, 2.92360, 2.92360]),
           'lam': np.array([12.5335, 13.49968, 13.33846, 11.77148, 11.52350]),
           'xan': np.array([-8.40590, -8.44849, -8.34906, -8.46070, -8.36174]),
           'yan': np.array([-2.48992, -2.54104, -2.52854, -2.44547, -2.43112])
           },
    '3B': {'x': np.array([574.26012, 1001.7349, 985.30166, 548.016, 519.98]),
           'y': np.array([512., 10., 100, 900, 1014]),
           's': np.array([8, 1, 1, 16, 16]),
           'alpha': np.array([0, -3.17728, 2.92434, -3.29402, 2.60797]),
           'beta': np.array([-0.19491, -2.92360, -2.92360, 2.92360, 2.92360]),
           'lam': np.array([14.53997, 15.66039, 15.47355, 13.65622, 13.36833]),
           'xan': np.array([-8.40044, -8.44785, -8.34786, -8.46088, -8.36211]),
           'yan': np.array([-2.48588, -2.53771, -2.52512, -2.44219, -2.42776])
           },
    '3C': {'x': np.array([573.25446, 1000.21721, 983.92918, 546.00285, 518.2782]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([8, 1, 1, 16, 16]),
           'alpha': np.array([0., -2.94573, 3.09057, -3.07810, 2.73161]),
           'beta': np.array([-0.19491, -2.92360, -2.92360, 2.92360, 2.92360]),
           'lam': np.array([16.79017, 18.08441, 17.86845, 15.76948, 15.43724]),
           'xan': np.array([-8.40205, -8.44574, -8.34664, -8.45859, -8.36196]),
           'yan': np.array([-2.48627, -2.53761, -2.52502, -2.44221, -2.42787]),
           },
    '4A': {'x': np.array([80.987181, 434.34987, 461.90855, 26.322503, 53.674656]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([6, 1, 1, 12, 12]),
           'alpha': np.array([0., -3.74625, 3.72621, -3.94261, 3.62762]),
           'beta': np.array([-0.32802, -3.60821, -3.60821, 3.60821, 3.60821]),
           'lam': np.array([19.34914, 20.93078, 20.6464, 18.07975, 17.67221]),
           'xan': np.array([-8.38446, -8.43506, -8.31378, -8.46256, -8.33609]),
           'yan': np.array([-2.48058, -2.5444, -2.52426, -2.42449, -2.40839])
           },
    '4B': {'x': np.array([77.625553, 431.57061, 458.86869, 23.559111, 50.632416]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([6, 1, 1, 12, 12]),
           'alpha': np.array([0., -3.64817, 3.73313, -3.73558, 3.74096]),
           'beta': np.array([-0.32802, -3.60821, -3.60821, 3.60821, 3.60821]),
           'lam': np.array([22.38267, 24.21212, 23.88327, 20.91426, 20.44279]),
           'xan': np.array([-8.38581, -8.43443, -8.3141, -8.46152, -8.33604]),
           'yan': np.array([-2.48185, -2.54526, -2.52568, -2.42513, -2.40959])
           },
    '4C': {'x': np.array([79.662059, 433.73384, 460.75026, 25.820431, 52.412219]),
           'y': np.array([512., 10, 100, 900, 1014]),
           's': np.array([6, 1, 1, 12, 12]),
           'alpha': np.array([0., -3.61682, 3.69713, -3.66259, 3.69888]),
           'beta': np.array([-0.32802, -3.60819, -3.60819, 3.60819, 3.60819]),
           'lam': np.array([26.18343, 28.32354, 27.93894, 24.46574, 23.91417]),
           'xan': np.array([-8.38603, -8.43509, -8.31524, -8.45888, -8.33707]),
           'yan': np.array([-2.48315, -2.54647, -2.52661, -2.42721, -2.41060])
           }
}




