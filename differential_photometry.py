# Import the required libraries

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils.aperture import SkyCircularAperture
from astropy.visualization.stretch import SinhStretch, LinearStretch
from astropy.visualization import ImageNormalize
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
from astroquery.astrometry_net import AstrometryNet
from photutils.aperture import SkyCircularAperture
from astropy.time import Time
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry
from astropy.visualization import simple_norm
import glob
import os

# User inputs

AstrometryNet_key = 'edqxeasvvonajjpl'

Target_RA = 194.0849
Target_DEC = 21.2909

Ref_RA = 194.1270
Ref_DEC = 21.2909

ap_size = 10 # aperture size, in pixels, 5
an_small = 20 # size of inner ring of annulus, in pixels, 10
an_large = 30 #size of outer ring of annulus, in pixels, 15
imsz = 150 # number of pixels from the middle to the edge of the image, in x and y

plotpath = "/Users/emilywhittaker/AnacondaProjects/CodeAstro/group_project/codeastro_project/" #path where you want to save the lightcurve and csv
fitpath = "/Users/emilywhittaker/AnacondaProjects/CodeAstro/group_project/codeastro_project/Small" #path where the fits files are

# Plate Solve the FIT file

def AstrometryNet_Plate_Solve(AstrometryNet_key = None, FIT_File = None, Target_RA = None, Target_DEC = None, Ref_RA = None, Ref_DEC = None):

    ast = AstrometryNet()
    ast.api_key = AstrometryNet_key

    try_again = True
    submission_id = None

    while try_again:
        try:
            if not submission_id:
                wcs_header = ast.solve_from_image(FIT_File,
                                                submission_id = submission_id)
            else:
                wcs_header = ast.monitor_submission(submission_id,
                                                    solve_timeout = 3000000000000)
        except TimeoutError as e:
            submission_id = e.args[1]
        else:
            # got a result, so terminate
            try_again = False

        w = WCS(wcs_header)

    Target_Pix = w.wcs_world2pix(Target_RA, Target_DEC, 1)

    Ref_Pix = w.wcs_world2pix(Ref_RA, Ref_DEC, 1)

    # Return the locations of the target star and the reference star in pixels

    return(Target_Pix, Ref_Pix)

def star_counts(pixels, data):
    """ Star Counts

    Using aperture photometry to find the number of counts measured from a star.

    Args:
        pixels (array): numpy vector. Location of the target star in pixels.
        data (array): numpy 2D array. Image data of the star from the fits file.

    Returns:
        float: instrument counts from the star
    """
    #Getting the mean background 
    aperture = CircularAperture(pixels, r=ap_size)
    annulus_aperture = CircularAnnulus(pixels, r_in=an_small, r_out=an_large)
    aperstats = ApertureStats(data, annulus_aperture)
    bkg_mean = aperstats.mean

    #Getting the aperture area
    phot_table = aperture_photometry(data, aperture)
    aperture_area = aperture.area_overlap(data)

    #Getting the photometry within the aperture
    total_bkg = bkg_mean * aperture_area
    phot_bkgsub = phot_table['aperture_sum'] - total_bkg
    return phot_bkgsub, aperture, annulus_aperture #total counts within the aperture, aperture object, annulus_object