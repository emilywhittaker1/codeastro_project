# differential_photometry

Creates a photometric lightcurve of a variable star.


# Overview:
This package is designed to input a collection of fits files containing 
reduced images containing a target variable star and the coordinates of the 
target star and a non-varying reference star at different points in time, and
output the photometric lightcurve of the target star.  


# Inputs: 
The user will put all of the necessary inputs into the user inputs 
block of the Jupyter notebook.  For example, your inputs might look
like this:

User inputs

an individual account key that the user must have on astrometry.net.
AstrometryNet_key = 'edqxeasvvonajjpl'   

The RA and Dec (in degrees) of the target star
Target_RA = 194.0849 
Target_DEC = 21.2909

The RA and Dec (in degrees) of the reference star
Ref_RA = 194.1270
Ref_DEC = 21.2909

aperture size, in pixels
ap_size = 10 
size of inner ring of annulus, in pixels
an_small = 20 
size of outer ring of annulus, in pixels
an_large = 30 
number of pixels from the middle to the edge of the image, in x and y
imsz = 150 

path where you want to save the lightcurve and csv
plotpath = "/Users/emilywhittaker/AnacondaProjects/CodeAstro/group_project/codeastro_project/" 
path where the fits files are
fitpath = "/Users/emilywhittaker/AnacondaProjects/CodeAstro/group_project/codeastro_project/Small" 


# Outputs:
The package will output a light curve of the variable star, with 
time listed in minutes since the first observation, and relative flux 
of the variable star measured in instrument counts.  This will be a 
.png file.

It will also output a .csv file containing the name of the fits file, the 
time of observations (MJD), the relative flux of the variable star (counts), 
the measured flux of the variable star (counts), the measured flux of the 
reference star (counts), the location of the target star in pixels, and the 
location of the reference star in pixels.

# Documentation:
Documentation for all of our functions is listed in the html file under
docs2/_build/html/index.html (double click to open in web browser).

In an ideal world, the documentation would also be available here: 
https://codeastro-project.readthedocs.io/en/latest/ 

# Contacts:
If you have any questions about the package, contact us at:
emilywhittaker at g.ucla.edu

![logo](https://github.com/emilywhittaker1/codeastro_project/blob/main/differential_photometry_logo.png?raw=true)
