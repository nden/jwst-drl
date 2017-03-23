# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Edited by DRL based on code by Nadia Dencheva
"""
Tools to create WCS reference files for MIRI MRS using IDT reference
products delivered with CDP-4:
MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_04.02.00.fits
MIRI MRS uses 5 reference files of type:
DISTORTION
SPECWCS
REGIONS
WAVELENGTHRANGE
V2V3
create_cdp4_references() creates all reference files.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from asdf import AsdfFile

from astropy.io import fits
from astropy.modeling import models

from matplotlib import image as mplimage
from matplotlib import pyplot as plt


def create_cdp6_references(fname, ref):
    """
    Create ASDF WCS reference files for MIRI MRS data from a CDP6 reference file.
    Parameters
    ----------
    fname : str
        name of reference file
    ref : dict
        A dictionary {reftype: refname}, e.g.
        {'distortion': 'jwst_miri_distortion_0001.asdf',
         'regions': 'jwst_miri_regions_0001.asdf',
         'specwcs': 'jwst_miri_specwcs_0001.asdf',
         'v2v3': 'jwst_miri_v2v3_00001.asdf',
         'wavelengthrange': 'jwst_miri_wavelengthrange_0001.asdf'}
    Examples
    --------
    >>> create_cdp6_references('MIRI_FM_MIRIFUSHORT_12LONG_DISTORTION_04.02.00.fits', ref)
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

    # Mandate that output ASDF files contain numbers inline
    outformat='inline'

    coeff_names = build_coeff_names_cdp5(alpha1.names)
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
    bmodel1, smodel1 = create_beta_models(b0_ch1, bdel_ch1, int(channel[0]), len(alpha1))
    bmodel2, smodel2 = create_beta_models(b0_ch2, bdel_ch2, int(channel[1]), len(alpha2))
    bmodel1.update(bmodel2)
    useafter = "2000-01-01T00:00:00"
    author =  'Adrian M. Glauser, David R. Law'  #Author of the data
    description = 'MIRI MRS CDP6 distortion reference data.'
    create_distortion_file(reftype='DISTORTION', detector=detector, band=band, channel=channel,
                           data=(amodel1, bmodel1, xmodel1, ymodel1, smodel1, smodel2), name=ref['distortion'],
                           author=author, useafter=useafter, description=description, outformat=outformat)

    create_specwcs_file('SPECWCS', detector, band, channel, lmodel1, ref['specwcs'], author,
                        useafter, description, outformat)

    create_v23('V2V3', detector, band, channels, (ab_v23, v23_ab), ref['v2v3'], author, useafter,
               description, outformat)

    create_regions_file(slices, detector, band, channel, ref['regions'], author, useafter,
               description, outformat)

    create_wavelengthrange_file(ref['wavelengthrange'], detector, author, useafter,
                                description, outformat)


def create_regions_file(slices, detector, band, channel, name, author, useafter, description, outformat):
    tree = create_reffile_header("REGIONS", detector, band, channel, author, useafter,
                                 description)
    tree['filename'] = name
    f = AsdfFile()
    tree['regions'] = slices
    f.tree = tree
#    f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    f.write_to(name,all_array_storage=outformat)


def create_reffile_header(reftype, detector, band, channel, author, useafter,
                          description=""):
    tree = {"author": author,
            "description": description,
            "detector": detector,
            "exp_type": "MIR_MRS",
            "filename": "",
            "instrume": "MIRI",
            "pedigree": "GROUND",
            "reftype": reftype,
            "telescope": "JWST",
            "useafter": useafter,
            "band": band,
            "channel": channel,
            "title": "MIRI IFU model - based on CDP-6",
            #"history": "DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;"
            }
    return tree


def create_v23(reftype, detector, band, channels, data, name, author, useafter, description, outformat):
    """
    Create the transform from MIRI Local to telescope V2/V3 system for all channels.
    """
    channel = "".join([ch[0] for ch in channels])
    tree = create_reffile_header(reftype, detector, band, channel, author, useafter, description)
    """
    tree = {"detector": detector,
            "instrument" : "MIRI",
            "band": band,
            "channel": channel,
            "exp_type": "MIR_MRS",
            "pedigree": "GROUND",
            "title": "MIRI IFU model - based on CDP-6",
            "reftype": reftype,
            "author": author,
            "useafter": useafter,
            "description": description
            }
    """
    tree['filename'] = name
    ab_v23 = data[0]
    v23_ab = data[1]
    m = {}
    c0_0, c0_1, c1_0, c1_1 = ab_v23[0][1:]
    ch1_v2 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[0][1:]
    ch1_a = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")

    c0_0, c0_1, c1_0, c1_1 = ab_v23[1][1:]
    ch1_v3 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[1][1:]
    ch1_b = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")
    c0_0, c0_1, c1_0, c1_1 = ab_v23[2][1:]
    ch2_v2 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[2][1:]
    ch2_a = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")

    c0_0, c0_1, c1_0, c1_1 = ab_v23[3][1:]
    ch2_v3 = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                 name="ab_v23")
    c0_0, c0_1, c1_0, c1_1 = v23_ab[3][1:]
    ch2_b = models.Polynomial2D(2, c0_0=c0_0, c1_0=c1_0, c0_1=c0_1, c1_1=c1_1,
                                name="v23_ab")
    ch1_for =  ch1_v2 & ch1_v3
    ch2_for = ch2_v2 & ch2_v3
    ch1_for.inverse =  ch1_a & ch1_b
    ch2_for.inverse =  ch2_a & ch2_b
    m[channels[0]] = ch1_for
    m[channels[1]] = ch2_for
    tree['model'] = m

    f = AsdfFile()
    f.tree = tree
#    f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    f.write_to(name,all_array_storage=outformat)


def create_distortion_file(reftype, detector,  band, channel, data, name, author,
                           useafter, description, outformat):

    description = 'MIRI MRS Distortion Maps'
    tree = create_reffile_header(reftype, detector, band, channel, author, useafter,
                                 description)
    tree['filename'] = name
    adata, bdata, xdata, ydata, sdata1, sdata2 = data
    tree['alpha_model'] = adata
    tree['beta_model'] = bdata
    tree['x_model'] = xdata
    tree['y_model'] = ydata
    tree['slice_model'] = {str(channel[0])+band: sdata1, str(channel[1])+band: sdata2}
    f = AsdfFile()
    f.tree = tree
#    f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    f.write_to(name,all_array_storage=outformat)


def create_specwcs_file(reftype, detector, band, channel, lmodel, name, author, useafter, description, outformat):
    tree = create_reffile_header(reftype, detector, band, channel, author, useafter,
                                 description)
    tree['subarray'] = "N/A"
    tree['filename'] = name
    tree['model'] = lmodel
    f = AsdfFile()
    f.tree = tree
#    f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    f.write_to(name,all_array_storage=outformat)


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
        transforms[sl] = models.Identity(1) & models.Shift(-xs) | \
                  models.Polynomial2D(8, name=name, **coeffs)
    return transforms


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

        transforms[sl] =  models.Shift(-xs, name=shname) & models.Identity(1) | \
                  models.Polynomial2D(8, name=pname, **coeffs)
    return transforms


def build_coeff_names_cdp5(names):
    names = names[1:]
    names = [name.replace('VAR2_', "c") for name in names]
    return names


def build_coeff_names(names):
    names = names[1:]
    names = [name.replace('VAR2(', "c") for name in names]
    names = [name.replace(')', "") for name in names]
    names = [name.replace(',', "_") for name in names]
    return names


def create_beta_models(b0, bdel, channel, nslices):
    beta = {}
    slices = {}
    for s in range(nslices):
        sl = channel * 100 + s +1
        beta_s = b0 + s * bdel
        m = models.Const1D(beta_s, name='det2local') #xy2beta and xy2lam
        beta[sl] = m
        inv = models.Const1D(sl)
        slices[beta_s] = models.Mapping([1,]) | inv
    return beta, slices


def create_wavelengthrange_file(name, detector, author, useafter, description, outformat):
    f = AsdfFile()

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

    tree = create_reffile_header("WAVELENGTHRANGE", detector, band="N/A", channel="N/A", author=author, 
                                 useafter=useafter, description=description)
    tree['filename'] = name
    tree['author'] = 'Nadia Dencheva'
    tree['detector'] = "N/A"
    tree['channels'] = channels

    f.tree = tree
    vr = np.empty((12, 2), dtype=np.float)
    for i, ch in enumerate(channels):
        vr[i] = wavelengthrange[ch]
    f.tree['wavelengthrange'] = vr
#    f.add_history_entry("DOCUMENT: MIRI-TN-00001-ETH; SOFTWARE: polyd2c_CDP5.pro; DATA USED: Data set of: - FM Test Campaign relevant to MRS-OPT-01, MRS-OPT-02, MRS-OPT-04, MRS-OPT-08; - CV1 Test Campaign relevant to MRS-OPT-02; - CV2 Test Campaign relevant to MRS-OPT-02; - Laboratory measurement of SPO; ============ DIFFERENCES: - New file structure: Change of Extention names and Table Column Headers.; - Replaced V2/V3 with XAN/YAN;")
    f.write_to(name,all_array_storage=outformat)
