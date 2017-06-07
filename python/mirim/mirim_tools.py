# Useful python tools for working with the MIRI Imager

import os as os
import numpy as np
from astropy.modeling import models
from asdf import AsdfFile
from jwst import datamodels
from jwst.assign_wcs import miri

########################################################################################
# General convenience tools
########################################################################################

# Convenience function to set default most-recent CRDS reference files from
# jwst-drl/trunk/CRDS/
def setreffiles_default():
    def_rootdir=os.path.expandvars('$JWSTTOOLS_DIR')

    # See https://jwst-crds.stsci.edu//display_result/52cef902-ad77-4792-9964-d26a0a8a96a8
    distfile=def_rootdir+'/CRDS/jwst_miri_distortion_0023.asdf'
    offfile=def_rootdir+'/CRDS/jwst_miri_filteroffset_0003.asdf'

    refs = {"distortion": distfile, "filteroffset": offfile}

    return refs

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

# Function to convert input X,Y pixel coordinates
# to V2V3 given a filter and the most recent reference files
def xytov2v3(filter):
    # Construct the reference data model in general JWST imager type
    input_model = datamodels.ImageModel()
    # Set the filter in the data model meta header
    input_model.meta.instrument.filter = filter
    # Set default reference files
    refs=setreffiles_default()
    # Call the pipeline code to make a distortion object given these inputs
    distortion = miri.imaging_distortion(input_model, refs)
    # Convert to arcsec
    distortion = distortion | models.Scale(3600.) & models.Scale(3600.)
    # Return the distortion object that can then be queried
    return distortion

########################################################################################
# CDP-specific tools
########################################################################################

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

