from codeastro_project.differential_photometry import  AstrometryNet_Plate_Solve, star_counts, do_dif_photometry, plot_lightcurve

# User inputs

# API Key from Astrometry.net
AstrometryNet_key = 'edqxeasvvonajjpl'

# RA and DEC of the target star
Target_RA = 194.0849
Target_DEC = 21.2909

# RA and DEC of the reference star
Ref_RA = 194.1270
Ref_DEC = 21.2909

# Settings for the differential photometry
ap_size = 10 # aperture size, in pixels, 5
an_small = 20 # size of inner ring of annulus, in pixels, 10
an_large = 30 #size of outer ring of annulus, in pixels, 15
imsz = 150 # number of pixels from the middle to the edge of the image, in x and y

# Wheather or not to plot the images with their apertures.
plotting = True

# Set the paths
plotpath = "Plot/" #path where you want to save the lightcurve and csv
fitpath = "Small/" #path where the fits files are

# Name of the output file
output_file = "Differential_Photometry.txt"

data = do_dif_photometry(AstrometryNet_key, Target_RA, Target_DEC, Ref_RA, Ref_DEC, ap_size, an_small, an_large, imsz, fitpath, plotpath, plotting = True, output_file = "Differential_Photometry.txt")

plot_lightcurve(data)