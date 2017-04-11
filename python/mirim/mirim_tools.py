# Useful python tools for working with the MIRI Imager

import numpy as np
from astropy.modeling import models
from asdf import AsdfFile
from jwst import datamodels
from jwst.assign_wcs import miri

# Convenience function to convert input X,Y pixel coordinates
# to V2V3 given a filter and an input reference file
def imager_detector_to_v2v3_cdp6(filter, refs):
    """
    Return the distortion, followed by the filter offset,
    the Xan_Yan to V2V3 transform and scaling to arcsec.

    """
    # Construct the reference data model in general JWST imager type
    input_model = datamodels.ImageModel()
    # Set the filter in the data model meta header
    input_model.meta.instrument.filter = filter
    # Call the pipeline code to make a distortion object given these inputs
    # This calls code at 
    # /Users/dlaw/anaconda2/envs/jwstb7rc2/lib/python2.7/site-packages/jwst-0.7.7.dev261-py2.7-macosx-10.7-x86_64.egg/jwst/assign_wcs/miri.py
    # Note that I've started modifying this code too
    distortion = miri.imaging_distortion(input_model, refs)
    # Return the distortion object that can then be queried
    return distortion

# Convenience function to convert input X,Y pixel coordinates
# to V2V3 given a filter and an input reference file
def imager_detector_to_v2v3_cdp7b(filter, refs):
    """
    Return the distortion, followed by the filter offset,
    the Xan_Yan to V2V3 transform and scaling to arcsec.

    """
    # Construct the reference data model in general JWST imager type
    input_model = datamodels.ImageModel()
    # Set the filter in the data model meta header
    input_model.meta.instrument.filter = filter
    # Call the pipeline code to make a distortion object given these inputs
    # This calls code at 
    # /Users/dlaw/anaconda2/envs/jwstb7rc2/lib/python2.7/site-packages/jwst-0.7.7.dev261-py2.7-macosx-10.7-x86_64.egg/jwst/assign_wcs/miri.py
    # Note that I've started modifying this code too
    distortion = miri.imaging_distortion(input_model, refs)
    distortion = distortion | models.Scale(3600.) & models.Scale(3600.)
    # Return the distortion object that can then be queried
    return distortion

# Convert v2,v3 in arcsec to xan,yan in arcmin
def v2v3_to_xanyan(v2,v3):
    xan=v2/60.
    yan=-(v3+7.8*60.)/60.
    return xan,yan

# Convert xan,yan in arcmin to v2,v3 in arcsec
def xanyan_to_v2v3(xan,yan):
    v2=xan*60.
    v3=(-yan-7.8)*60.
    return v2,v3

