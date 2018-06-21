
# coding: utf-8

# # CIG96 Pseudoring
# 
# ##GENERAL INFORMATION
# 
# Notebook with the analysis of the pseudoring, based on the calibrated R,G,B images by Peris.
# 
# ### AMIGA WIKI PAGES:
# 
# **Optical data**
# - Info from the CAHA campaign can be found <a href="http://amiga.iaa.es:8888/display/science/CIG96+-+Deep+optical+and+HI+images#CIG96-DeepopticalandHIimages-InterestingHIandopticalplots">here</a>.
# 
# - The information about CIG96 paper (2016) can be found <a href="http://amiga.iaa.es:8888/display/science/CIG96+-+Deep+optical+and+HI+images#CIG96-DeepopticalandHIimages-Paperstructureandcomments">here</a>.
# 
# ###a) Import necessary modules

# In[1]:

import aplpy
from astropy.utils.data import get_pkg_data_filename as fitsWCS
from astropy.wcs import WCS as wcs
import astropy.io.fits as fits
from   astropy import units as u
import csv
from   glob import glob as ls
import imexam
#from   kapteyn import maputils
import lemon.astrometry as astrometry
#import lemon.photometry as photometry
import matplotlib as mpl
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import montage_wrapper as montage
import numpy as np
from   numpy.linalg import eig, inv
import os
#import pyds9 as ds9
import pylab
import pyraf
import pyregion as pyreg
#import random
#import repipy.combine as combine
import repipy.find_sky as find_sky
import scipy
from scipy import interpolate as interpolate
from scipy import ndimage
import shutil
#from   spectral_cube import SpectralCube as spc
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")
#import pickle
#import pandas
#import sqlite3
#import photutils as phot
#import seaborn as sns
#import repipy.astroim as astroim
#import repipy.rename as rename
#import repipy.utilities as utilities
#import repipy.combine as combine
#import repipy.create_masks as create_masks
#import repipy.remove_cosmics as remove_cosmics
#import repipy.arith as arith
#import repipy.scale_to_ref as scale_to_ref
#import repipy.extract_mag_airmass_common as extract


# In[2]:

# renders interactive figures in the Notebook
get_ipython().magic(u'matplotlib nbagg')


# #0. Functions
# 
# All final and working functions used throughout the notebook are stored in this cell.

# In[3]:

def go_to(dir):
    """Function to select and switch to optical or HI folders"""
    work = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16/" + dir
    os.chdir(work)
    #print "Work/save directory:", work

##########

# Ellipse fitting function

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def ellipse_angle_of_rotation2( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a < c:
            return 0
        else:
            return np.pi/2
    else: 
        if a < c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2

###############

# Isophotes regions cleaning

def clean_isoph(raw_isoph_file, regions, clean_file):
    '''This function *removes* all the points that fall inside the specified
    rectangular region/s of the image.
    'raw_isoph_file' is the input file with the raw isophote contours described in 2 columns
    'regions' must be a list of lists consisting on four floats ordered as:
    x0 = X initial coordinate, same raw file coordinates
    y0 = Y initial coordinate, same raw file coordinates
    x1 = X final coordinate, same raw file coordinates
    y1 = Y final coordinate, same raw file coordinates
    'clean_file' is the output file with the clean isophote contours described in 2 columns as well
    '''    
    with open(raw_isoph_file) as infile, open(clean_file, "w+") as outfile:
        raw = infile.readlines()
        for line in raw:
            parts = line.split()
            cx = float(parts[0])
            cy = float(parts[1])
            j = True                                     # all coords are good unless they belong to a region (see below)
            for reg in regions:
                x0 = float(reg[0])
                y0 = float(reg[1])
                x1 = float(reg[2])
                y1 = float(reg[3])
                if ((x0 < cx < x1) and (y0 < cy < y1)):  # if the coordinates fall in any region box -> the coords are not good
                    j = False                            # bad coords
            if j:
                outfile.write(line)
                
###############

# Isophotes fitting

def isofit(iso_file):
    '''This function will fit an ellipse over any two-columns x,y (pix) file provided.
    'iso_file' is the input coordinates file
    'surf_bri' is the surface brightness of the isophote'''
    x, y = np.loadtxt(iso_file, unpack=True, usecols=(0,1))
    fit = fitEllipse(x, y)
    center = ellipse_center(fit)
    phi = ellipse_angle_of_rotation(fit)
    axes = ellipse_axis_length(fit)
    return  center[0], center[1], phi, axes[0], axes[1], max(axes)/min(axes)
               
###############

# Emission region aperture sum in a ds9 image

def apersum(image, coordinates, radius, sums, errs):
    '''This function will calculate the sum of the pixels contained in
    certain apertures defined by "coords" file from a given "image".'''
    # read coordinates
    coord = np.genfromtxt(coordinates)

    im = fits.open(image)
    data = im[0].data
    dimx = data.shape[1]
    dimy = data.shape[0]
    hist = []
        
    with open(sums, "w+") as output, open(errs, "w+") as errors:
        for pair in coord:
            cx = pair[0]   # center x coordinate of the aperture (in pixels)
            cy = pair[1]   # center y coordinate of the aperture (in pixels)
            r = float(radius)          # radius of the aperture (in pixels)
            y, x = np.ogrid[-cy:dimy-cy, -cx:dimx-cx]   # circular shape, must be defined like this to fit im and data shapes
            mask = x**2 + y**2 <= r**2    # mask condition: all points within 'r', described as right triangles: x² + y² <= r²
            hist.append(np.array(data[mask]).tolist())
            mean  = np.mean(data[mask])       # mean value
            sd = np.std(data[mask])        # st.dev. value
            output.write(str(mean) + "\n")
            errors.write(str(sd) + "\n")
    return hist, mean, sd
            
###############

# Sky aperture sum in a ds9 image

def apersumsky(image, coordinates, radius, sums):
    '''This function will calculate the rms of the pixels contained in
    certain apertures defined by "coords" file from a given "image".'''
    # read coordinates
    coord = np.genfromtxt(coordinates)

    im = fits.open(image)
    data = im[0].data
    dimx = data.shape[1]
    dimy = data.shape[0]
        
    with open(sums, "w+") as output:
        for pair in coord:
            cx = pair[0]   # center x coordinate of the aperture (in pixels)
            cy = pair[1]   # center y coordinate of the aperture (in pixels)
            r = float(radius)          # radius of the aperture (in pixels)
            y, x = np.ogrid[-cy:dimy-cy, -cx:dimx-cx]   # circular shape, must be defined like this to fit im and data shapes
            mask = x**2 + y**2 <= r**2    # mask condition: all points within 'r', described as right triangles: x² + y² <= r²
            rms  = np.sqrt(np.mean(np.square(data[mask])))       # rms value
            output.write(str(rms) + "\n")
            
###############

# Rotate FITS image X degrees

def rotate(image, degrees):
    '''This function will rotate a FITS image by n degrees,
    either clockwise or anticlockwise, and write the result.'''
    # read image
    image = image
    im = fits.open(image)
    data = im[0].data
    # rotation
    data = ndimage.interpolation.rotate(data, float(degrees), reshape=False)
    # save rotated fits
    fits.writeto(image[:-5]+"_"+str(degrees)+"rot.fits", data, clobber=True)
    
###############

# Discrete colormap
# Source: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
# By Jake VanderPlas
# License: BSD-style

import matplotlib.pyplot as plt
import numpy as np


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

if __name__ == '__main__':
    N = 5

    x = np.random.randn(40)
    y = np.random.randn(40)
    c = np.random.randint(N, size=40)

    # Edit: don't use the default ('jet') because it makes @mwaskom mad...
    #plt.scatter(x, y, c=c, s=50, cmap=discrete_cmap(N, 'cubehelix'))
    #plt.colorbar(ticks=range(N))
    #plt.clim(-0.5, N - 0.5)
    #plt.show()
    
###############

def ping():
    os.system('spd-say "Ping"')


# To see if there is a gradient along the ring, we follow this process:
# 
# - We define regions containing the emission over the ring in the RGB image.
# - They are overlaid in the B image and R image and their medians and rms are extracted, respectively.
# - Color is extracted from subtracting the R region value to the B region value, both in magnitudes.
# 
# **NOTE**: this method assumes same exposure time and calibration of the images.
# 
# #1. Images calculations and deprojections
# 
# - CIG96 morphological type = Sc = 5:      **q (Sc) = 0.175**
# 
# - Axis relation (minor over major axis):    **b/a = 0.698821608** (extracted from HI data)
# 
# - Inclination (sources: Jakoby, Masters et al. 2014):   $cos^{2}(i) = \frac{(b/a)^{2} - q^{2}}{1 - q^{2}}$:       **i = 46.594368772 deg = 46.59 deg**
# 
# ##1.1 Extinction correction
# 
# Explained in <a href="http://amiga.iaa.es:8888/display/science/Extinction+correction+for+CIG+96+colours">Extinction correction for CIG96</a>.
# 
# - SDSS conversions to extinction (g,r bands):
# $$A_x(g) = 3.793 * E_{B-V}$$
# $$A_x(r) = 2.751 * E_{B-V}$$
# 
# - From AMIGA Data release 2012 (B band):
# $$A_{int}(B) = 0.276$$
# $$A_k(B) = 0.009$$
# 
# - From Savage & Mathis 1979:
# $$A_x(B) = 4.1 * E_{B-V}$$
# 
# such as:
# 
# $$A_x(g) = 3.793 * A(B) / 4.1$$
# 
# In the case of the <b>internal extinction</b>:
# $$A_{int}(g) = 3.793 * 0.276 / 4.1 $$
# $$A_{int}(r) = 2.751 * 0.276 / 4.1 $$
# 
# In the case of the <b>k correction</b>:
# $$A_k(g) = 3.793 * 0.009 / 4.1$$
# $$A_k(r) = 2.751 * 0.009 / 4.1$$
# 
# So:
# 
# $$A_{int}(g) = 0.255$$
# $$A_{int}(r) = 0.185$$
# $$A_k(g) = 0.008$$
# $$A_k(r) = 0.006$$
# 
# The conversion to g,r SDSS magnitudes, after extinction and k corrections are:
# 
# $$m_{r_{SDSS}} = 1.01*mag_{R_{phot}} - 9.83 - (0.185 + 0.06) \pm0.15$$
# $$m_{g_{SDSS}} = 0.99*mag_{B_{phot}} - 9.70 - (0.255 + 0.08) \pm0.33$$
# 
# simplified as:
# 
# $$m_{r_{SDSS}} = 1.01*mag_{R_{phot}} - 10.02 \pm0.15$$
# $$m_{g_{SDSS}} = 0.99*mag_{B_{phot}} - 9.96 \pm0.33$$
# 
# 
# ##1.2 g, r magnitude, flux and color index g over r images.
# 
# Conversion of R, B images into r, g images both in flux and magnitudes.

# In[4]:

# Select working directory
go_to('deep_peris')

# convert B,R to g,r in magnitudes
#r = 1.01*(-2.5*np.log10(float(elemr))) - 9.83 - 0.191
#g = 0.99*(-2.5*np.log10(float(elemg))) - 9.70 - 0.263

# image r in mag and flux
image = "CIG96_R_koords_skysub_al_psf.fits"
im = fits.open(image)
im[0].data = 1.01*(-2.5*np.log10(im[0].data)) - 10.02 + 25
im.writeto("cig96_r.fits", clobber=True)
im[0].data = 10**(-0.4*im[0].data)
np.nan_to_num(im[0].data)
im.writeto("cig96_r_flux.fits", clobber=True)

# image g in mag and flux
image = "CIG96_B_koords_skysub_al_psf.fits"
im = fits.open(image)

im[0].data = 0.99*(-2.5*np.log10(im[0].data)) - 9.96 + 25
im.writeto("cig96_g.fits", clobber=True)
im[0].data = 10**(-0.4*im[0].data)
np.nan_to_num(im[0].data)
im.writeto("cig96_g_flux.fits", clobber=True)


                In IRAF via imarith:

g over r in flux

imarith:

operand1=    cig96_g_flux.fits  Operand image or numerical constant
op      =                    /  Operator
operand2=    cig96_r_flux.fits  Operand image or numerical constant
result  = cig96_g_over_r_flux.fits  Resultant image
                
# In[133]:

# Select working directory
go_to('deep_peris')

image = "cig96_g_over_r_flux.fits"
im = fits.open(image)
im[0].data = -2.5*np.log10(im[0].data)
im.writeto("cig96_g_minus_r_mag.fits", clobber=True)


# In[195]:

# Select working directory
go_to('deep_peris')

result = "cig96_g_minus_r_mag.fits"

#print np.amin(fits.open(result)[0].data)
#print np.amax(fits.open(result)[0].data)

fig = plt.figure(figsize=(8, 8))
plt.axis('off')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.08

goverr = aplpy.FITSFigure(result, figure=fig, north=True)
goverr.show_colorscale(stretch='linear', cmap='RdBu_r', vmin=0.35, vmax=1.2)
goverr.add_colorbar()
goverr.recenter(RA, DEC, radius=rad)
##goverr.save("cig96_colorindex_goverr.png")


# In[151]:

plt.close()


# ##1.2 Deprojections

                IRAF commands for all deprojections necessary:

geotran input=CIG96_B_SB.fits output=CIG96_B_SB_dep.fits xin=371.37686 yin=416.86709 xrot=-20.39 yrot=-20.39 xscale=1.430615
geotran input=CIG96_R_SB.fits output=CIG96_R_SB_dep.fits xin=371.37686 yin=416.86709 xrot=-20.39 yrot=-20.39 xscale=1.430615

geotran input=CIG96_B_koords_skysub_al_psf.fits output=CIG96_B_koords_skysub_al_psf_dep.fits xin=371.37686 yin=416.86709 xrot=-20.39 yrot=-20.39 xscale=1.430615
geotran input=CIG96_R_koords_skysub_al_psf.fits output=CIG96_R_koords_skysub_al_psf_dep.fits xin=371.37686 yin=416.86709 xrot=-20.39 yrot=-20.39 xscale=1.430615

geotran input=cig96_g_flux.fits output=cig96_g_flux_dep.fits xin=371.37686 yin=416.86709 xrot=-20.39 yrot=-20.39 xscale=1.430615
geotran input=cig96_r_flux.fits output=cig96_r_flux_dep.fits xin=371.37686 yin=416.86709 xrot=-20.39 yrot=-20.39 xscale=1.430615

geotran input=cig96_g_over_r_flux.fits[0] output=cig96_g_over_r_flux_dep.fits xin=506.456 yin=539.199 xrot=-20.39 yrot=-20.39 xscale=1.430615

DS9 region file:   cig96_pseudo_reg_dep_BandR.reg
Circular areas radius = 12.47 arcsec
                
# #2. Apertures selection

# In[262]:

# Select working directory
go_to('deep_peris')

# figure structure
fig = plt.figure(figsize=(12, 6))
plt.suptitle("De-projected images (g and r) of CIG96 with emission extraction apertures")
plt.axis('off')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.08   # zoom radius so, the larger rad, the larger FoV in the image

# plot window design
vw = 0.80 # vertical width of the subplots in %
hw = 0.38 # horizontal width of the subplots in % 

# deprojected B
Bdep = aplpy.FITSFigure('cig96_g_flux_dep.fits', figure=fig, subplot=[0.11, 0.1, hw, vw])
Bdep.show_colorscale(cmap="Blues", stretch='log', exponent=2, vmin=1e-18, vmax=1e-10, vmid=1e-20)
Bdep.show_regions("cig96_pseudo_reg_dep_BandR_wcs.reg")
# zoom
Bdep.recenter(RA, DEC, radius=rad)

# deprojected R
Rdep = aplpy.FITSFigure('cig96_r_flux_dep.fits', figure=fig, subplot=[0.52, 0.1, hw, vw])
Rdep.show_colorscale(cmap="Reds", stretch='log', exponent=2, vmin=1e-18, vmax=1e-10, vmid=1e-20)
Rdep.show_regions("cig96_pseudo_reg_dep_BandR_wcs.reg")
# zoom
Rdep.recenter(RA, DEC, radius=rad)

# to hide
Rdep.hide_yaxis_label()
Rdep.hide_ytick_labels()


# In[263]:

plt.close()


# #3. Pseudoring g - r color and HI measurements
# 
# ##3.1 Pseudoring apertures
# 
# Here we measure the apertures selected for the pseudoring in both g,r flux images, via 'apersum' function.

# In[7]:

# Select working directory
go_to('deep_peris')

# PSEUDORING, g and r

# images 
img = "cig96_g_flux_dep.fits"
imr = "cig96_r_flux_dep.fits"

# coordinates file
coordsgr = "cig96_pseudo_reg_dep_BandR_xy.reg"

# apertures calculations
aperture_rad = 8.33 # 12.46 arcsec
hist, mean, err = apersum(img, coordsgr, aperture_rad, "g.txt", "errg.txt")
apersum(imr, coordsgr, aperture_rad, "r.txt", "errr.txt")

color_gr = []

with open("g.txt", "r") as gsum, open("r.txt", "r") as rsum:
    for elemg, elemr in zip(gsum, rsum):
        g = -2.5*np.log10(float(elemg))
        r = -2.5*np.log10(float(elemr))
        color_gr.append(float(g) - float(r))
    meancolor_gr = np.mean(color_gr)

print "max g-r =", round(max(color_gr),4)
print "min g-r =", round(min(color_gr),4)

###

# pseudoring error
with open("g.txt", "r") as g, open("r.txt", "r") as r, open("errg.txt", "r") as errg, open("errr.txt", "r") as errr:
    abserr_gr = []
    for elemg, elemr, errorg, errorr in zip(g, r, errg, errr):
        #error = 2.5/np.log(10.)*(float(errorg)/float(elemg) + float(errorr)/float(elemr))
        error = 2.5*np.log10(float(errorg)) - 2.5*np.log10(float(errorr))
        abserr_gr.append(round(error,4))
        
print np.std(abserr_gr)
print "color_gr = ", color_gr


# - Histograms of each aperture (g - r)

# In[10]:

fig, ax = plt.subplots(7,5, figsize=(15, 20))#, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.5, wspace=0.02)
fig.tight_layout()

ax = ax.ravel()
medians = []  # median values of each histogram/aperture

for elem in range(0,len(hist)):
    num = np.random.random()
    ax[elem].hist(hist[elem], color=(num, num, num), bins=20)
    ax[elem].set_title(str(elem), color='blue')
    medians.append(round(np.median(hist[elem]),13))

a = np.arange(0, len(medians), 1)
for elem, med in zip(a, medians):
    if med > 2.2e-11 and med < 2.7e-11:
        print "Typical RED aperture: ", elem, med
    elif med > 4.3e-11 and med < 4.8e-11:
        print "Typical BLUE aperture: ", elem, med


# In[148]:

plt.close()


                flux = [3.7e-11, 5.25e-11, 4.38e-11, 3.84e-11, 6.4e-11, 5.76e-11, 4.24e-11, 4.73e-11, 4.77e-11, 5.16e-11, 9.46e-11, 4.91e-11, 3.7e-11, 3.14e-11, 3.99e-11, 5.28e-11, 4.43e-11, 4.51e-11, 6.11e-11, 5e-11, 4.95e-11, 2.97e-11, 2.35e-11, 2.07e-11, 2.81e-11, 2.65e-11, 4.95e-11, 1.6e-11, 2.2e-11, 1.84e-11, 2.6e-11, 2.77e-11, 3.34e-11]
mags = [-2.5*np.log10(elem) for elem in flux]

# all median values
plt.plot(flux, 'b*')
                
                plt.close()
                
                # all median values in magnitudes
plt.plot(mags, 'C4d')
                
                plt.close()
                
                blue_mean = []
red_mean = []
for elem in flux:
    if elem <= 3e-11:
        red_mean.append(elem)
    else:
        blue_mean.append(elem)

print "blue mean =", np.mean(blue_mean)
print "red mean =",  np.mean(red_mean)
                
# In[149]:

fig = plt.figure(figsize=(10,7))
plt.hist(hist[8], bins=20, color='blue', alpha=0.8)
plt.annotate(s='Mean BLUE apertures = 4.87 * 10e-11', xy=[0.487e-10,52], xytext=(0.65e-10,55), size=22, color='blue', arrowprops=dict(facecolor='blue'))
plt.hist(hist[23], bins=20, color='red', alpha=0.5)
plt.annotate(s='Mean RED apertures = 2.39 * 10e-11', xy=[0.2386e-10,63], xytext=(0.05e-10,67), size=22, color='red', arrowprops=dict(facecolor='red'))
#plt.annotate("North", xy=(150, 349), xytext=(160, 350), size=20, color='white', arrowprops=dict(facecolor='white'))
plt.ylim(0.0, 75)
plt.xlabel("Flux", size=24)
plt.show()


# In[150]:

plt.close()


# ##3.2 Sky apertures
# 
# Here we measure the surrounding sky, to see later if there is any correlation with the pseudoring and the sky nearby, indicating an artificial color pattern in the pseudoring.

# In[151]:

# SKY, g and r

# Select working directory
go_to('deep_peris')

# images 
img = "cig96_g_flux_dep.fits"
imr = "cig96_r_flux_dep.fits"

# coordinates file
skygr = "cig96_pseudo_reg_dep_ring_sky_xy.reg"

# apertures calculations
aperture_rad = 8.33 # 12.46 arcsec
apersumsky(img, skygr, aperture_rad, "skyg.txt")
apersumsky(imr, skygr, aperture_rad, "skyr.txt")

skycolor_gr = []

with open("skyg.txt") as gsum, open("skyr.txt") as rsum:
    for elemg, elemr in zip(gsum, rsum):
        skyg = elemg
        skyr = elemr
        skycolor_gr.append(float(skyg)/float(skyr))                      # g - r in fluxes

###

# sky error
skyerr = np.std(skycolor_gr)
print np.median(skycolor_gr), skyerr


# ####Angles distribution:

                # clockwise
#ang = [0.0, 7.8, 17.6, 26.6, 36.1, 43.97, 51.99, 60.2, 68.2, 75.7, 86.9, 96.63, 105.53, 119.5, 128.91, 136.29, 144.47, 150.31, 158, 171.04, 178.55, 191.22, 199.6, 208.5, 222.9, 240.6, 249.29, 265.25, 273.63, 284.26, 295.33, 313.16, 322.1]

# anticlockwise
ang = [360.0, 352.2, 342.4, 333.4, 323.9, 316.03, 308.01, 299.8, 291.8, 284.3, 273.1, 263.37, 254.47, 240.5, 231.09, 223.71, 215.53, 209.69, 202.0, 188.96, 181.45, 168.78, 160.4, 151.5, 137.1, 119.4, 110.71, 94.79, 86.37, 75.74, 64.67, 46.84, 37.9]
                
# In[152]:

# Pseudoring angles (not equidistant) for both optical B,R and HI pseudoring apertures
# setting 0 deg in the Northern feature, leftmost aperture and increasing ANTIclockwise:
ang = [360.0, 352.2, 342.4, 333.4, 323.9, 316.03, 308.01, 299.8, 291.8, 284.3, 273.1, 263.37, 254.47, 240.5, 231.09, 223.71, 215.53, 209.69, 202.0, 188.96, 181.45, 168.78, 160.4, 151.5, 137.1, 119.4, 110.71, 94.79, 86.37, 75.74, 64.67, 46.84, 37.9]

# Sky angles (equidistant) for 360 coverage apertures
angsk = range(0, len(skycolor_gr))
angsky = []
for elem in angsk:
    angle = round(float(360./float(len(skycolor_gr))*elem),2)
    angsky.append(angle) 

angsky = angsky[::-1]  # invert angle so they increase ANTIclockwise
print angsky


# ##3.3 HI apertures
# 
# From <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96%20-%20HI%20maps.ipynb#">HI maps notebook</a> we take 'HI_sums.txt' file.
# 
# Green valley boundaries from <a href="http://iopscience.iop.org/article/10.1088/0004-637X/775/2/129/pdf">Walker et al. 2013</a>, Figure 3.

# In[29]:

# Select working directory
go_to('cig')

# HI: 0th moment 3sigma
HI_mean_val = []
with open("HI_sums.txt") as HI:
    for elem in HI:
        HI_mean_val.append(float(elem))

#conversion to N(HI) in 10²¹ at/cm²
NHI_apers = [elem*1.411429847 for elem in HI_mean_val]   # WARNING: this is x10²¹ at/cm²

# invert list to start from NorthEast aperture running anticlockwise until Northernmost aperture
#NHI_apers = NHI_apers[::-1]

# NHI rms
NHIrms = 0.06*1.411429847


# In[30]:

plt.figure(figsize=(11, 6))
for nhi, color, angle in zip(NHI_apers, color_gr, ang):
    plt.plot(nhi, color, 'bD')
    plt.annotate(s=int(angle), xy=(nhi, color), xytext=(nhi+0.005,color), fontsize=8)
plt.xlabel("N(HI) 10^21 at/cm-2")
plt.ylabel("g-r")
plt.title("Color vs N(HI) relation")
plt.axhspan(0.55, 0.7, color='g', alpha=0.2, label="Green valley")
plt.axhspan(0.1, 0.55, color='b', alpha=0.05, label="Bluer")
plt.axhspan(0.7, 1.1, color='r', alpha=0.05, label="Redder")
plt.ylim(0.1, 1.1)
plt.legend()

# linear regression
slope, interc, R, P, std_Err = scipy.stats.linregress(NHI_apers , color_gr)
print "Slope =", slope
print "R2 =", R**2
line = [slope*elem + interc for elem in NHI_apers]
plt.plot(NHI_apers, line, 'r-')
plt.annotate(s="R2 = " + str(round(R**2,4)), xy=(0.6, 0.2))


# In[31]:

plt.close()


# ###3.2.1 g - r (pseudoring) + g - r (sky) + N(HI)
# 
# More on y axis label ticks: <a href="http://stackoverflow.com/questions/20350503/remove-first-and-last-ticks-label-of-each-y-axis-subplot">here</a>.

# In[35]:

# plot
fig, (ax0, ax1, ax2) = plt.subplots(3, sharex=True, squeeze=True, figsize=(10,10))

plt.suptitle('g - r and HI distribution along the pseudoring', fontsize=16)
plt.xlim(-10, 380)
plt.xlabel("Angle (deg)")
plt.xticks(np.arange(0, 370, 20), rotation=0)
plt.subplots_adjust(wspace=0.01, hspace=0.05)

# pseudoring error: yerr=abserr_gr
ax0.errorbar(ang, color_gr,  capsize=2, ecolor='red', c='brown', fmt='o-', markersize=6, markeredgecolor='black', label="g - r  (mag, pseudoring regions)")
ax0.axvline(x=0.0, color='b', alpha=0.3, linestyle='--')
ax0.axvline(x=360.0, color='b', alpha=0.3, linestyle='--')
ax0.xaxis.set_tick_params(labeltop='on')
#fillbottom = [fsa - err for fsa,err in zip(color_gr, abserr_gr)]
#filltop    = [fsb + err for fsb,err in zip(color_gr, abserr_gr)]
#ax0.fill_between(ang, fillbottom, filltop, color='red', alpha=0.02)
ax0.axes.set_ylabel("g - r (mag)")
ax0.legend(numpoints=1, loc=0, fontsize=10)
ax0.grid(alpha=0.2)

# sky
ax1.errorbar(angsky, skycolor_gr,  yerr=skyerr, alpha=0.1, fmt='s-', label="g - r  (mag, sky 360deg coverage)", linewidth=2)
ax1.axvline(x=0.0, color='b', alpha=0.3, linestyle='--')
ax1.axvline(x=360.0, color='b', alpha=0.3, linestyle='--')
ax1.axes.set_ylabel("g - r (mag)")
fillSbottom = [fsa - skyerr for fsa in skycolor_gr]
fillStop    = [fsb + skyerr for fsb in skycolor_gr]
ax1.fill_between(angsky, fillSbottom, fillStop, color='blue', alpha=0.02)
ax1.grid(alpha=0.2)
#cut = 1.12
#ax1.axhline(y=cut, color='k', alpha=0.8, linestyle=':')
#for angle,elem in zip(angsky, skycolor_gr):
#    clean = []
#    if float(elem) > cut:
#        ax1.plot(angle, elem, 'yh', markersize=8) #, label='Aperture contains star')
#ax1.plot(-100, 0, 'yh', label='Aperture contains star (g-r > ' + str(cut) + ' mag)') # plots a ghost hexagon se I can use the label and aovid repeating it in the loop
#ax1.axhspan(cut + 0.01, 2.5, xmin=0, xmax=360, alpha=0.1, color='yellow')   # yellow region
ax1.legend(numpoints=1, loc=0, fontsize=10)

# HI
ax2.errorbar(ang, NHI_apers, yerr=NHIrms, alpha=0.9, color='green', fmt='o-', markersize=6, markeredgecolor='black', label="N(HI) (pseudoring regions)", linewidth=1)
ax2.axvline(x=0.0, color='b', alpha=0.3, linestyle='--')
ax2.axvline(x=360.0, color='b', alpha=0.3, linestyle='--')
ax2.axes.set_ylabel("Mean N(HI) 10e21 at/cm2")
ax2.legend(numpoints=1, loc=0, fontsize=10)
ax2.grid(alpha=0.2)


#ax2.errorbar(ang, HInorm, yerr=HIerr, alpha=0.9, color='green', fmt='-', label="HI (pseudoring regions)", linewidth=1)
#ax2.axvline(x=0.0, color='b', alpha=0.3, linestyle='--')
#ax2.axvline(x=360.0, color='b', alpha=0.3, linestyle='--')
#ax2.axes.set_ylabel("Mean subtracted integrated flux")
#ax2.legend(numpoints=1, loc=0, fontsize=10)


# In[36]:

plt.close()


# ###3.2.2 g - r (pseudo) and 360deg sky

                print np.median(color_gr[:18])
print np.std(color_gr[:18])

print np.median(color_gr[18:])
print np.std(color_gr[18:])
                
# To have a nice view of the two-colored sides of the pseudoring, we must shift all data (both from the pseudoring and HI-sky):
# 
# #### PSEUDORING (red points)

# In[153]:

# PA angles
ang = [360.0, 352.2, 342.4, 333.4, 323.9, 316.03, 308.01, 299.8, 291.8, 284.3, 273.1, 263.37, 254.47, 240.5, 231.09, 223.71, 215.53, 209.69, 202.0, 188.96, 181.45, 168.78, 160.4, 151.5, 137.1, 119.4, 110.71, 94.79, 86.37, 75.74, 64.67, 46.84, 37.9]
ang_shifted = [(x-37.9+15.5) for x in ang] # 37.9 to dismantle the previous shift and 15.5 to match the first value to 60 degrees, later on, the x-axis of the final plot will also be shifted 60 degrees as well to produce a circular plot
#print ang_shifted


# #### SKY (blue points)

# In[154]:

# ANGLES
angsky = [354.19, 348.39, 342.58, 336.77, 330.97, 325.16, 319.35, 313.55, 307.74, 301.94, 296.13, 290.32, 284.52, 278.71, 272.9, 267.1, 261.29, 255.48, 249.68, 243.87, 238.06, 232.26, 226.45, 220.65, 214.84, 209.03, 203.23, 197.42, 191.61, 185.81, 180.0, 174.19, 168.39, 162.58, 156.77, 150.97, 145.16, 139.35, 133.55, 127.74, 121.94, 116.13, 110.32, 104.52, 98.71, 92.9, 87.1, 81.29, 75.48, 69.68, 63.87, 58.06, 52.26, 46.45, 40.65, 34.84, 29.03, 23.23, 17.42, 11.61, 5.81, 0.0]
angsky_shifted = [angle-60 for angle in angsky] # to match the x-value of the pseudoring regions, we shift the sky data the same 60 degrees, later on, the x-axis of the final plot will also be shifted 60 degrees as well to produce a circular plot

# converting negative values to positive (+360 deg)
for pos, elem in enumerate(angsky_shifted):
    if elem < 0:
        angsky_shifted[pos] = elem + 360

print angsky_shifted


# In[270]:

print np.median(skycolor_gr)
print np.std(skycolor_gr)
print len(skycolor_gr)


# In[264]:

# Select working directory
go_to('cig')

# figure
plt.figure(figsize=(13,7))

plt.ylim(-0.4, 2.50)
#xaxis_names = np.append(list(np.arange(60, 380, 20)), np.arange(20, 80, 20))
plt.xticks(np.arange(0.0, 380, 20), xaxis_names, size=11)
plt.yticks(np.arange(-0.4,2.5,0.2), size=11)
plt.xlabel("PA (deg)", size=14)
plt.ylabel("$g - r$ (mag)", size=14)

# pseudoring
plt.errorbar(ang_shifted, color_gr, yerr=np.std(abserr_gr), capsize=2, ecolor='brown', fmt='o', c='red', label="g - r (pseudoring)") #yerr=abserr_gr,
#plt.xticks(ang_shifted)
fillbottom = [fsa - err for fsa,err in zip(color_gr, abserr_gr)]
filltop    = [fsb + err for fsb,err in zip(color_gr, abserr_gr)]
#plt.fill_between(ang, fillbottom, filltop, color='red', alpha=0.04)

# sky colors and plot
cskydots    = 'skyblue'
cskyerrbars = 'skyblue'
cskyfill    = 'skyblue'
plt.errorbar(angsky_shifted, skycolor_gr, yerr=skyerr, fmt='s', alpha=0.45, c=cskydots, ecolor=cskyerrbars, label="g - r (sky)", linewidth=2)

# fill between values needs redefinition of the two sides with an extra value as a consequence of the redefinition of the angles (above 360)
leftside_angsky_shifted  = angsky_shifted[0:-11]
rightside_angsky_shifted = angsky_shifted[-11:]
# includes the redundant value (2.75 + 360) where the break happens
a = 300
leftside_angsky_shifted = [a] + leftside_angsky_shifted
# we add the same y-value to preserve the dimensions in y data as well
leftside_skycolor_gr  = skycolor_gr[0:-11]
rightside_skycolor_gr = skycolor_gr[-11:]
b =  0.9711650812627145
leftside_skycolor_gr = [b] + leftside_skycolor_gr
# color filling command
fillSbottom = [fsa - skyerr for fsa in leftside_skycolor_gr]
fillStop    = [fsb + skyerr for fsb in leftside_skycolor_gr]
plt.fill_between(leftside_angsky_shifted, fillSbottom, fillStop, color=cskyfill, alpha=0.15)

fillSbottom = [fsa - skyerr for fsa in rightside_skycolor_gr]
fillStop    = [fsb + skyerr for fsb in rightside_skycolor_gr]
plt.fill_between(rightside_angsky_shifted, fillSbottom, fillStop, color=cskyfill, alpha=0.15)

# SDSS red/blue VS Green Valley boundaries according to Walker et al. 2013
plt.axhspan(0.55, 0.7, color='g', alpha=0.2)

# angle limits and regions
plt.axvline(x=198, ymin=0, ymax=1.0, color='k', alpha=0.8, linestyle=':')
plt.axvline(x=0.0, color='b', alpha=0.3, linestyle='--')
plt.axvline(x=360.0, color='b', alpha=0.3, linestyle='--')

# grid
plt.grid(alpha=0.1)

# legend
plt.legend(numpoints=1, loc=(0.30,0.01), fontsize=12)
plt.show()

# factor for the next histograms
plt.annotate(s='x 10$^{-10}$', xy=(150, 1.325), fontsize=10)
plt.annotate(s='x 10$^{-10}$', xy=(331, 1.325), fontsize=10)

# Typical BLUE and RED apertures histograms
# we eliminate the 10**10 factor from the histogram
hred = []
for elem in hist[30]:
    a = elem*10**10
    hred.append(a)
# and now print the histogram
red = plt.axes([0.22, 0.60, 0.22, 0.21])
n, bins, patches = plt.hist(hred, bins=20, color='red', alpha=0.7)
plt.title('Flux histogram of a typical \n aperture in PA = 70$\degree$ - 258$\degree$ arc', size=10)

# we eliminate the 10**10 factor from the histogram
hblue = []
for elem in hist[16]:
    b = elem*10**10
    hblue.append(b)
# and now print the histogram
blue = plt.axes([0.59, 0.60, 0.22, 0.21])
n, bins, patches = plt.hist(hblue, bins=20, color='darkblue', alpha=0.8)
plt.title('Flux histogram of a typical \n aperture in PA = 258$\degree$ - 38$\degree$ arc', size=10)

# save png
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_pseudo_apertures.png', bbox_inches='tight')


# In[271]:

plt.close()


# ###3.2.3 Scaled g - r and HI (both in pseudoring)
# 
# Scale = mean subtracted and sigma weighted

# In[278]:

# figure
#plt.figure(figsize=(12,7))
#plt.title("Scaled g - r and HI along pseudoring | Scale = mean subtracted and sigma weighted")
#plt.xlim(-10, 380)

# scaled pseudoring
sigma_gr = np.std(color_gr)
scalecolor_gr = [(elem - meancolor_gr)/sigma_gr for elem in color_gr]  # scale = data - mean divided by sigma

# scaled HI
HIsums  = [elem for elem in NHI_apers]
meanHI  = np.mean(HIsums)                                     # mean of the HI data
sigmaHI = np.std(HIsums)                                      # standard deviation of the HI data
HIsumscale = [(elem - meanHI)/sigmaHI for elem in HIsums]

# errors
#scalecolor_gr_err = 
HIsumscale_err = np.abs((NHIrms - meanHI)/sigmaHI)
print HIsumscale_err
print meanHI
# HI - log and downsized to fit the Y axis (irrelevant units)
#HIsumslog  = [np.log10(elem) for elem in NHI_apers]                    # HI in logarithmic scale
#meanHIlog  = np.mean(HIsumslog)                                     # mean of the logarithmic HI data
#sigmaHIlog = np.std(HIsumslog)                                      # standard deviation of the log. HI data
#HIsumscale = [(elem - meanHI)/sigmaHI for elem in HIsums]           # arithmetic scale = data - mean divided by sigma
#HIsumscale = [(elem - meanHIlog)/sigmaHIlog for elem in HIsumslog]  # logarithmic scale = data - mean divided by sigma

# plot properties
#plt.errorbar(ang, scalecolor_gr, capsize=2, ecolor='brown', marker='o', c='red', linewidth=1, label="Scaled g - r")
#plt.plot(ang, HIsumscale, color='green', marker='o', alpha=0.9, label="Scaled N(HI)", linewidth=1)
#plt.xlabel("Angle (deg)")
#plt.ylabel("Scaled g - r (mag) / Scaled N(HI) (10e21 at/cm2)")
#plt.legend(numpoints=1, loc=0, fontsize=11)

#plt.axvspan(240,360, color='b', alpha=0.05) # most blue region of the pseudoring
#plt.axvspan(0,239, color='r', alpha=0.05) # rest of the ring, mainly redder

#plt.legend()

#plt.savefig("cig96_pseudo_HI.png")


                plt.close()
                
# ###Combined plot: HI vs color

# In[33]:

# Select working directory
go_to('cig')

# HI: 0th moment 3sigma
HI_mean_val = []
with open("HI_sums.txt") as HI:
    for elem in HI:
        HI_mean_val.append(float(elem))

#conversion to N(HI) in 10^20 at/cm²
NHI_apers = [elem*14.11429847 for elem in HI_mean_val]   # WARNING: this is x10²¹ at/cm²

# invert list to start from NorthEast aperture running anticlockwise until Northernmost aperture
#NHI_apers = NHI_apers[::-1]

# NHI rms
NHIrms = 0.06*14.11429847


# In[281]:

# PA angles
ang = [38, 30, 20, 11, 1, 353, 345, 337, 329, 322, 311, 301, 292, 278, 268, 261, 253, 247, 239, 226, 219, 206, 198, 189, 175, 157, 148, 132, 124, 113, 102, 84, 75]

# g over r figure
fig, [ax0, ax1] = plt.subplots(nrows=2, ncols=1, figsize=(9,8.5))
plt.subplots_adjust(left=0.125, bottom=0.125, right=0.90, top=0.9, hspace=0.3)

####################

# g-r vs PA

# plot
ax0.errorbar(ang, scalecolor_gr, color='brown', marker='o', markeredgecolor='k', markersize=10, linewidth=0, label="Scaled $g - r$")
ax0.plot(ang, HIsumscale, color='g', marker='X', markersize=10, markeredgecolor='k', linewidth=0, alpha=0.9, label="Scaled N$_{HI}$")
ax0.set_xlabel("PA (deg)", size=18)
xti = np.arange(0,370,30)
ax0.set_xticks(xti)
ax0.set_xticklabels(xti, fontsize=14)
ax0.set_ylabel("Scaled $g-r$ \n Scaled N$_{HI}$", size=18)
yti = np.arange(-2.6,3.0,0.6)
ax0.set_yticks(yti)
ax0.set_yticklabels(yti, fontsize=14)
ax0.grid(alpha=0.15)

# marks
ax0.axhline(0, color='k', linestyle='--', alpha=0.2)
ax0.set_xlim(-5, 375)

# anticorrelation broken in 70-180 deg
ax0.axvspan(90, 180, color='black', alpha=0.10)

####################

# g-r vs NHI

# discrete scaled colormap
cmap_discrete = discrete_cmap(18, base_cmap='hsv')

# plot
ax = ax1.scatter(NHI_apers, color_gr, s=100, c=ang, cmap=cmap_discrete, edgecolor='k')
ax1.set_xlabel("N$_{HI}$ (10$^{20}$ at cm$^{-2})$", size=18)
xt = np.arange(5, 14, 1)
ax1.set_xticks(xt)
ax1.set_xticklabels(xt, fontsize=14)
ax1.set_ylabel("$g-r$ (mag)", size=18)
yt = np.arange(0.0, 1.2, 0.1)
ax1.set_yticks(yt)
ax1.set_yticklabels(yt, fontsize=14)
ax1.grid(alpha=0.15)
ax1.axhspan(0.55, 0.7, color='g', alpha=0.2, label="Green valley")
ax1.axhspan(0.0, 0.55, color='b', alpha=0.02, label="Bluer")
ax1.axhspan(0.7, 1.2,  color='r', alpha=0.02, label="Redder")
ax1.set_ylim(0.05, 1.15)

# colorbar
cbar = fig.add_axes([0.91, 0.13, 0.01, 0.30])
#cbar.set_yscale(np.arange(0, 360, 10))
cbar.set_label(np.arange(50, 480, 10))
#cbar.set_yticklabels()
cbar = fig.colorbar(ax, cax=cbar, ticks=np.arange(0, 360, 60))
cbar.ax.tick_params(labelsize=14)
ax1.annotate(s="  PA (deg)", xy=(13.5, 1.05), fontsize=16)



# linear regression
#slope, interc, R, P, std_Err = scipy.stats.linregress(NHI_apers , color_gr)
#print "Slope =", slope
#print "R =", R
#print "P =", P
#print "std_Err = ", std_Err
#line = [slope*elem + interc for elem in NHI_apers]
#ax1.plot(NHI_apers, line, 'r')
#ax1.annotate(s="R2 = " + str(round(R**2,4)), xy=(6, 0.2))

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_HIvsColor.png', bbox_inches='tight')


# In[282]:

plt.close()


# ####Correlation between scaled HI and g - r data:

# In[145]:

# linear regression
slope, interc, R, P, std_Err = scipy.stats.linregress(HIsumscale, scalecolor_gr)
print "Slope =", slope
print "R2 =", R**2

# plot
plt.plot(HIsumscale, scalecolor_gr, 'co')
line = [slope*elem + interc for elem in HIsumscale]
plt.plot(HIsumscale, line, 'r-')
plt.annotate(s="R2 = "+str(round(R**2,4)), xy=(0.0, -0.6))
plt.xlabel("Scaled HI")
plt.ylabel("Scaled g - r (mag)")

plt.savefig("cig96_pseudo_HI_corr.png")


# In[146]:

plt.close()


# ###3.2.4 Interpolation of scaled pseudoring and full (360 degrees) HI regions

# In[147]:

# apertures sums
# HI: 0th moment 3sigma

# image
HI = "cig96_HI_mom0_3.0rms_dep.fits"

#coordinates file
coordsHI = "cig96_HI360_dep_xy.reg"   # 360deg regions on deprojected HI image
aperture_rad = 2.4465 # = 14 arcsec

# apertures calculations
apersum(HI, coordsHI, aperture_rad, "HI360sums.txt", "HI360err.txt")

HI360sums     = []
HI360sumslog  = []
HI360scale    = []
HI360scalelog = []

with open("HI360sums.txt") as HI:
    for sumap in HI:
        HI360sums.append(float(sumap))
        HI360sumslog.append(np.log10(float(sumap)))
    meanHI = np.mean(HI360sums)
    sigmaHI = np.std(HI360sums)
    meanHIlog = np.mean(HI360sumslog)
    sigmaHIlog = np.std(HI360sumslog)
    for elem in HI360sums:
        HI360scale.append((elem - meanHI)/sigmaHI)
    for elemlog in HI360sumslog:
        HI360scalelog.append((elemlog - meanHIlog)/sigmaHIlog)

# HI angles (ONLY FOR circular coverage)
HI360ang = range(0, len(HI360sums))
HI360angles = []
for elem in HI360ang:
    HI360angles.append(float(360./float(len(HI360sums)))*elem)        # equidistant angle distribution


# In[148]:

# 1D interpolation function
HIf = interpolate.interp1d(HI360angles + [360.0], HI360scale + [HI360scale[0]])
HIflog = interpolate.interp1d(HI360angles + [360.0], HI360scalelog + [HI360scalelog[0]])

# run interpolation
interpHI = HIf(ang)
interpHIlog = HIflog(ang)


# In[149]:

fig = plt.figure(figsize=(12,8))
plt.xlim(-10, 370)
plt.xlabel("")
plt.ylabel("")
plt.plot(ang, scalecolor_gr, 'ro-', label="g - r (mag)")                              # pseudoring B-R measurements
plt.plot(HI360angles, HI360scalelog, 'g.-', label="Original HI measurements (42 apertures)")  # pseudoring HI measurements (42)
plt.plot(ang, interpHIlog, 'bo-', label="Interpolated HI (33 apertures)")                     # interpolated HI measurements
plt.xlabel("Angle (deg)")
plt.ylabel("Scaled g - r (mag) and HI")
plt.legend(numpoints=1, loc=0, fontsize=10)

plt.savefig("cig96_pseudo_HIinterp.png")


# In[150]:

plt.close()


# ####Correlation between scaled and interpolated HI and g - r data

# In[151]:

# linear regression
slope, interc, R, P, std_Err = scipy.stats.linregress(interpHIlog, scalecolor_gr)
print "Slope =", slope
print "R2 =", R**2

# plot
plt.plot(interpHIlog, scalecolor_gr, 'co')
line = [slope*elem + interc for elem in interpHIlog]
plt.plot(interpHIlog, line, 'r-')
plt.annotate(s="R2 = " + str(round(R**2,4)), xy=(0.0, -0.6))
plt.xlabel("Scaled HI")
plt.ylabel("Scaled g - r")

plt.savefig("cig96_pseudo_HIinterp_corr.png")


# In[152]:

plt.close()


# ##3.3 Optical g - r pseudoring evolution

                # reproject the FITS image - perform this only once
montage.reproject(image, "cig96_RdRGB_SB_oriented.fits", north_aligned=True)
                
# In[153]:

# Select working directory
go_to('deep_peris')

# read text file:
regfile = open("cig96_pseudo_reg_BandR_xyorient.reg")

# define x,y
x = []
y = []

for col in (raw.strip().split() for raw in regfile):
    x.append(float(col[0]))
    y.append(float(col[1]))

# select and open image
image = "cig96_RdRGB_SB_oriented.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.07

# figure
fig = plt.figure(figsize=(12, 12))

# image plot
pseudo = aplpy.FITSFigure(im, auto_refresh=True, figure=fig, north=True)
pseudo.set_title("CIG96 pseudoring's g - r color evolution")
pseudo.show_colorscale(stretch='log', cmap='bone', vmin=22, vmax=28)   # for RGB-image visualization
#pseudo.show_regions("cig96_pseudo_reg_BandR_wcs.reg", cmap="seismic")

# plot   
plt.scatter(x, y, s=550, c=color_gr, cmap="seismic", edgecolors='black')

# zoom
pseudo.recenter(RA, DEC, rad)

# colorbar
cbar = plt.colorbar(fraction=0.04, pad=0.02)
cbar.set_label("g - r (mag)")

plt.savefig("cig96_opt_regs.png")


# In[154]:

plt.close()


# ##3.4 Optical g - r pseudoring evolution over HI

# In[581]:

# Select working directory
go_to('deep_peris')

# read text file:
regfile = open("cig96_pseudo_reg_HI_xy.reg")

# define x,y
x = []
y = []

for col in (raw.strip().split() for raw in regfile):
    x.append(float(col[0]))
    y.append(float(col[1]))
    
# select and open image
image = "cig96_HI_mom0_3rms_coldens.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.3

# figure
fig = plt.figure(figsize=(12, 12))

# image plot
pseudoHI = aplpy.FITSFigure(im, auto_refresh=True, figure=fig, north=True, zorder=1)
pseudoHI.set_title("CIG96 pseudoring's g - r color evolution over HI")
pseudoHI.show_colorscale(stretch='linear', cmap='RdYlGn_r', vmin=0, vmax=11)   # for RGB-image visualization

# beam definition
pseudoHI.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, edgecolor='k', facecolor=(0.9,0.9,0.2,0.99))

# N(HI) contours levels in 10²¹ cm⁻²
levs = [1.0]
pseudoHI.show_contour(levels=levs, colors='red', linewidth=6, zorder=2)
levs = [0.25, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 3.0, 5.0, 7.0, 10.0]
pseudoHI.show_contour(levels=levs, cmap="bone", linewidth=6, zorder=2)

# plot   
plt.scatter(x, y, s=100, c=color_gr, cmap="seismic", edgecolors='black', zorder=3)  # s = 550
# colorbar
cbar = plt.colorbar(fraction=0.04, pad=0.02)
cbar.set_label("g - r (mag)")

# zoom
pseudoHI.recenter(RA, DEC, rad)

plt.savefig("cig96_NHI_regs.png")


# In[582]:

pseudoHI.close()


# ##3.5 Optical image with HI contours

# In[157]:

# Select working directory
go_to('deep_peris')

# read text file:
regfile = open("cig96_pseudo_reg_BandR_xyorient.reg")

# define x,y
x = []
y = []

for col in (raw.strip().split() for raw in regfile):
    x.append(float(col[0]))
    y.append(float(col[1]))

# select images
image  = "cig96_RdRGB_SB_oriented.fits"
imcont = "cig96_HI_mom0_3.0rms.fits"
im   = fits.open(image)
cont = fits.open(imcont)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.07

# figure
fig = plt.figure(figsize=(12, 12))

# image plot
map_im   = "bone"
map_cont = "plasma"
opt = aplpy.FITSFigure(im, auto_refresh=True, figure=fig, north=True)
opt.set_title("CIG96 pseudoring's g - r color evolution over HI")
opt.show_colorscale(stretch='linear', cmap=map_im , vmin=20, vmax=28)
opt.show_contour(cont, levels=np.arange(1250, 4000, 250), cmap=map_cont, smooth=1, kernel="gauss")

# plot   
plt.scatter(x, y, s=550, c=color_gr, cmap="seismic", edgecolors='black')

# zoom
opt.recenter(RA, DEC, rad)

# colorbar
cbar = plt.colorbar(fraction=0.04, pad=0.04)
cbar.set_label("g - r (mag)")

plt.savefig("cig96_optHIcont_regs.png")


# In[158]:

plt.close()


# #4. Radial profile in specific directions

# ##4.1 Individual g,r bands
# 
# a) Rotation

# In[94]:

# Select working directory
go_to('deep_peris')

ims = ["cig96_r_flux_dep.fits", "cig96_g_flux_dep.fits"]

for im in ims:
    rotate(im, -8)
    rotate(im, -23)
    rotate(im, -100)
    rotate(im, -153)
    rotate(im, -174)


# b) Selection:
# 
# - rotation(NORTH-oriented image) = -100 deg
# - rotation(NORTH cut) = -8 deg
# - rotation(SOUTHWEST cut) = -153 deg
# - rotation(SOUTHWEST cut, external feature) = -174 deg

# In[8]:

# Select working directory
go_to('deep_peris')

PA1 = 6
PA2 = 25
PA3 = 52
PA4 = 19

# CIG96 coordinates in deg
RA = 262
DEC = 600
rad = 200

oriented = "cig96_g_flux_dep_-100rot.fits"
image = fits.open(oriented)

# image plot
rotim = aplpy.FITSFigure(image)
#rotim.set_title("CIG96 radial profiles")
rotim.show_colorscale(stretch='linear', cmap='bone_r', vmin=1e-11, vmax=9e-11)   # for RGB-image visualization

# profile cuts
plt.plot(RA, DEC, 'wx', zorder=4)
plt.plot([RA,RA+35],[DEC,DEC+145], linewidth=4, color='g',  label="PA = 6$\degree$")
plt.plot([RA,RA+2],[DEC,DEC+150],  linewidth=4, color='C6', label="PA = 25$\degree$")
plt.plot([RA,RA+50],[DEC,DEC-180], linewidth=4, color='y', label="PA = 199$\degree$")
plt.plot([RA,RA+95],[DEC,DEC-120], linewidth=4, color='C0', label="PA = 232$\degree$")
plt.legend()

# zoom
rotim.recenter(RA, DEC, rad)

# save figure
#plt.savefig('/home/prm/Desktop/cig96_images/cig96_radprofiles_cuts.png', bbox_inches='tight')


# In[15]:

rotim.close()


# Angluar-physical distance scales:

# In[86]:

akpc = 0.098417 # kpc/arcsec
pkpc = 0.103015 # kpc/pix


# In[93]:

# Select working directory
go_to('deep_peris')

# PA = 6 deg
images = ["cig96_g_flux_dep_-23rot.fits", "cig96_r_flux_dep_-23rot.fits"]

# box size
width = 10  # pixels
length = 140

# distance
dist = [pkpc*elem for elem in np.arange(0,length)]

for image in images:
    radprof = np.zeros(length, dtype=list)  # empty list for the radial profile
    # open image
    im = fits.open(image)
    data = im[0].data
    # collapse box
    for i in range(0,width):
        sel = data[598+i,433:433+length]
        radprof = [float(x + y) for x,y in zip(radprof, sel.tolist())]
    radprof = [-2.5*np.log10(elem/width) for elem in radprof]
    plt.plot(dist, radprof, label=image[6:7])   # plots the profile

# lines
plt.axvline(0, linestyle=':', color='k', alpha=0.3)
plt.axvline(101*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
plt.axvline(111*pkpc, linestyle='-.', color='g', alpha=0.5)
plt.axvline(131*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')
plt.axhline(27.5, linestyle='-', color='k', alpha=0.2, label="SB limit")
    
plt.ylim(max(radprof)+1, min(radprof)-1)
plt.legend(fontsize=8)
plt.xlabel("Distance (kpc)")
plt.ylabel("$\mu$ (mag/arcsec$^{2}$)")
plt.title("PA = 6 deg")


# In[94]:

plt.close()


# In[95]:

# PA = 25 deg
images = ["cig96_g_flux_dep_-8rot.fits", "cig96_r_flux_dep_-8rot.fits"]

# box size
width = 10    # pixels
length = 140  # pixels

# distance
dist = [pkpc*elem for elem in np.arange(0,length)]

for image in images:
    radprof = np.zeros(length, dtype=list)  # empty list for the radial profile
    # open image
    im = fits.open(image)
    data = im[0].data
    # collapse box
    for i in range(0,width):
        sel = data[574+i,457:457+length]
        radprof = [float(x + y) for x,y in zip(radprof, sel.tolist())]
    radprof = [-2.5*np.log10(elem/width) for elem in radprof]
    plt.plot(dist, radprof, label=image[6:7])   # plots the profile

# lines
plt.axvline(0, linestyle=':', color='k', alpha=0.3)
plt.axvline(92*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
plt.axvline(103*pkpc, linestyle='-.', color='g', alpha=0.5)
plt.axvline(128*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')
plt.axhline(27.5, linestyle='-', color='k', alpha=0.2, label="SB limit")

plt.ylim(max(radprof)+1, min(radprof)-1)
plt.legend(fontsize=8)
plt.title("PA = 25 deg")


# In[96]:

plt.close()


# In[97]:

# PA = 199 deg
images = ["cig96_g_flux_dep_-174rot.fits", "cig96_r_flux_dep_-174rot.fits"]

# box size
width = 10 # pixels
length = 200

# distance
dist = [pkpc*elem for elem in np.arange(0,length)]

for image in images:
    radprof = np.zeros(length, dtype=list)  # empty list for the radial profile
    # open image
    im = fits.open(image)
    data = im[0].data
    # collapse box
    for i in range(0,width):
        sel = data[485+i,213:213+length]
        radprof = [float(x + y) for x,y in zip(radprof, sel.tolist())]
    radprof = [-2.5*np.log10(elem/width) for elem in radprof]
    plt.plot(dist, radprof, label=image[6:7])   # plots the profile

# lines
plt.axvline(0, linestyle=':', color='k', alpha=0.3)
plt.axvline(106*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
plt.axvline(112*pkpc, linestyle='-.', color='g', alpha=0.5)
plt.axvline(147*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')
plt.axhline(27.5, linestyle='-', color='k', alpha=0.2, label='SB limit')
    
plt.ylim(max(radprof)+1, min(radprof)-1)
plt.legend(fontsize=8)
plt.xlabel("Distance (kpc)")
plt.ylabel("$\mu$ (mag/arcsec$^{2}$)")
plt.title("PA = 199 deg")


# In[98]:

plt.close()


# In[99]:

# PA = 232 deg
images = ["cig96_g_flux_dep_-153rot.fits", "cig96_r_flux_dep_-153rot.fits"]

# box size
width = 10 # pixels
length = 160

# distance
dist = [pkpc*elem for elem in np.arange(0,length)]

for image in images:
    radprof = np.zeros(length, dtype=list)  # empty list for the radial profile
    # open image
    im = fits.open(image)
    data = im[0].data
    # collapse box
    for i in range(0,width):
        sel = data[485+i,213:213+length]
        radprof = [float(x + y) for x,y in zip(radprof, sel.tolist())]
    radprof = [-2.5*np.log10(elem/width) for elem in radprof]
    plt.plot(dist, radprof, label=image[6:7])   # plots the profile

# lines
plt.axvline(0, linestyle=':', color='k', alpha=0.3)
plt.axvline(106*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
plt.axvline(112*pkpc, linestyle='-.', color='g', alpha=0.5)
plt.axvline(147*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')
plt.axhline(27.5, linestyle='-', color='k', alpha=0.2, label='SB limit')
    
plt.ylim(max(radprof)+1, min(radprof)-1)
plt.legend(fontsize=8)
plt.xlabel("Distance (kpc)")
plt.ylabel("$\mu$ (mag/arcsec$^{2}$)")
plt.title("PA = 232 deg")


# In[100]:

plt.close()


# ##4.2 Color g-r images
# 

# In[101]:

# Select working directory
go_to('deep_peris')

# g over r figure
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=4, sharey=True, squeeze=True, figsize=(20,4))
plt.subplots_adjust(wspace=0.02, hspace=0.05)
#fig.tight_layout()

##################

# box size
width = 15 # pixels
length = 140

# distance
dist = [pkpc*elem for elem in np.arange(0,length)]

##################

# PA = 6 deg
images = ["cig96_g_flux_dep_-23rot.fits", "cig96_r_flux_dep_-23rot.fits"]

# coordinates
col_x = 433
lin_y = 598

radprof_g = np.zeros(length, dtype=list)  # empty list for the radial profile
radprof_r = np.zeros(length, dtype=list)  # empty list for the radial profile

# open images
img = fits.open(images[0])
imr = fits.open(images[1])
datag = img[0].data
datar = imr[0].data

# collapse box
for i in range(0,width):
    selg = datag[lin_y+i,col_x:col_x+length]
    radprof_g = [float(x + y) for x,y in zip(radprof_g, selg.tolist())]
    selr = datar[lin_y+i,col_x:col_x+length]
    radprof_r = [float(x + y) for x,y in zip(radprof_r, selr.tolist())]
    g_r6 = np.divide(radprof_g, radprof_r)
    g_r6 = [-2.5*np.log10(elem) for elem in g_r6]

ax1.plot(dist, g_r6, label="PA = 6 deg", color="g")   # plots the profile

# lines
ax0.axvline(0, linestyle=':', color='k', alpha=0.3)
ax0.axvline(101*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
ax0.axvline(111*pkpc, linestyle='-.', color='g', alpha=0.5)
ax0.axvline(131*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')

# SDSS red/blue VS green valley boundaries according to Walker et al. 2013
ax0.axhspan(0.55, 0.7, color='y', alpha=0.5, label="Green valley")

# hide Y ticks
ax0.tick_params(axis='y', length=0)

##################

# PA = 25 deg
images = ["cig96_g_flux_dep_-8rot.fits", "cig96_r_flux_dep_-8rot.fits"]

# coordinates
col_x = 457
lin_y = 574

radprof_g = np.zeros(length, dtype=list)  # empty list for the radial profile
radprof_r = np.zeros(length, dtype=list)  # empty list for the radial profile

# open images
img = fits.open(images[0])
imr = fits.open(images[1])
datag = img[0].data
datar = imr[0].data

# collapse box
for i in range(0,width):
    selg = datag[lin_y+i,col_x:col_x+length]
    radprof_g = [float(x + y) for x,y in zip(radprof_g, selg.tolist())]
    selr = datar[lin_y+i,col_x:col_x+length]
    radprof_r = [float(x + y) for x,y in zip(radprof_r, selr.tolist())]
    g_r25 = np.divide(radprof_g, radprof_r)
    g_r25 = [-2.5*np.log10(elem) for elem in g_r25]

ax0.plot(dist, g_r25, label="PA = 25 deg", color="C6")   # plots the profile

# lines
ax1.axvline(0, linestyle=':', color='k', alpha=0.3)
ax1.axvline(92*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
ax1.axvline(103*pkpc, linestyle='-.', color='g', alpha=0.5)
ax1.axvline(128*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')

# SDSS red/blue VS green valley boundaries according to Walker et al. 2013
ax1.axhspan(0.55, 0.7, color='y', alpha=0.5, label="Green valley")

##################

# PA = 232 deg
images = ["cig96_g_flux_dep_-153rot.fits", "cig96_r_flux_dep_-153rot.fits"]

# coordinates
col_x = 213
lin_y = 485

radprof_g = np.zeros(length, dtype=list)  # empty list for the radial profile
radprof_r = np.zeros(length, dtype=list)  # empty list for the radial profile

# open images
im_g = fits.open(images[0])
im_r = fits.open(images[1])
datag = im_g[0].data
datar = im_r[0].data

# collapse box
for i in range(0,width):
    selg = datag[lin_y+i,col_x:col_x+length]
    radprof_g = [float(x + y) for x,y in zip(radprof_g, selg.tolist())]
    selr = datar[lin_y+i,col_x:col_x+length]
    radprof_r = [float(x + y) for x,y in zip(radprof_r, selr.tolist())]
    g_r52 = np.divide(radprof_g, radprof_r)
    g_r52 = [-2.5*np.log10(elem) for elem in g_r52]
    
ax2.plot(dist, g_r52, label="Southwestern cut", color='C0')   # plots the profile

# lines
ax2.axvline(0, linestyle=':', color='k', alpha=0.3)
ax2.axvline(106*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
ax2.axvline(112*pkpc, linestyle='-.', color='g', alpha=0.5)
ax2.axvline(147*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')

# SDSS red/blue VS green valley boundaries according to Walker et al. 2013
ax2.axhspan(0.55, 0.7, color='y', alpha=0.5, label="Green valley")

# hide Y ticks
ax2.tick_params(axis='y', length=0)

####################

# Overlay of both the profiles
ax3.plot(dist, g_r6,   label="PA = 6 deg",   color="g")   # plots the profile
ax3.plot(dist, g_r25,  label="PA = 25 deg",  color="C6")   # plots the profile
ax3.plot(dist, g_r52, label="PA = 232 deg", color='C0')   # plots the profile

# lines: cuts selected over the mean value of the 
ax3.axvline(0, linestyle=':', color='k', alpha=0.3)
ax3.axvline(100*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
ax3.axvline(109*pkpc,linestyle='-.', color='g', alpha=0.5)
ax3.axvline(135*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')

# SDSS red/blue VS green valley boundaries according to Walker et al. 2013
ax3.axhspan(0.55, 0.7, color='y', alpha=0.5, label="Green valley")

# hide Y ticks
ax3.tick_params(axis='y', length=0)

####################

# Titles and legend

ax0.set_ylabel("g-r (mag/arcsec$^{2}$)")
ax0.set_xlabel("Distance (kpc)")
ax0.set_title("PA = 6 deg")
ax1.set_xlabel("Distance (kpc)")
ax1.set_title("PA = 25 deg")
ax2.set_xlabel("Distance (kpc)")
ax2.set_title("PA = 232 deg")
ax3.set_xlabel("Distance (kpc)")
ax3.set_title("Overlay")

# grid
for graph in [ax0, ax1, ax2, ax3]:
    graph.grid(alpha=0.1)
    
# legend
plt.legend(fontsize=8)


# In[102]:

plt.close()


# In[5]:

# Select working directory
go_to('deep_peris')

##################

# radial profile box size
width = 15      # pixels
length = 180    # pixels

# distance
dist = [pkpc*elem for elem in np.arange(0,length)]

##################

# PA = 6 deg
images = ["cig96_g_flux_dep_-23rot.fits", "cig96_r_flux_dep_-23rot.fits"]

# coordinates
col_x = 433
lin_y = 598

radprof_g = np.zeros(length, dtype=list)  # empty list for the radial profile
radprof_r = np.zeros(length, dtype=list)  # empty list for the radial profile

# open images
img = fits.open(images[0])
imr = fits.open(images[1])
datag = img[0].data
datar = imr[0].data

# collapse box
for i in range(0,width):
    selg = datag[lin_y+i,col_x:col_x+length]
    radprof_g = [float(x + y) for x,y in zip(radprof_g, selg.tolist())]
    selr = datar[lin_y+i,col_x:col_x+length]
    radprof_r = [float(x + y) for x,y in zip(radprof_r, selr.tolist())]
    g_r6 = np.divide(radprof_g, radprof_r)
    g_r6 = [-2.5*np.log10(elem) for elem in g_r6]

##################

# PA = 25 deg
images = ["cig96_g_flux_dep_-8rot.fits", "cig96_r_flux_dep_-8rot.fits"]

# coordinates
col_x = 457
lin_y = 574

radprof_g = np.zeros(length, dtype=list)  # empty list for the radial profile
radprof_r = np.zeros(length, dtype=list)  # empty list for the radial profile

# open images
img = fits.open(images[0])
imr = fits.open(images[1])
datag = img[0].data
datar = imr[0].data

# collapse box
for i in range(0,width):
    selg = datag[lin_y+i,col_x:col_x+length]
    radprof_g = [float(x + y) for x,y in zip(radprof_g, selg.tolist())]
    selr = datar[lin_y+i,col_x:col_x+length]
    radprof_r = [float(x + y) for x,y in zip(radprof_r, selr.tolist())]
    g_r25 = np.divide(radprof_g, radprof_r)
    g_r25 = [-2.5*np.log10(elem) for elem in g_r25]

##################

# PA = 199 deg
images = ["cig96_g_flux_dep_-174rot.fits", "cig96_r_flux_dep_-174rot.fits"]

# coordinates
col_x = 225
lin_y = 441

# new box size
length = 190    # pixels
# new distance
dist = [pkpc*elem for elem in np.arange(0,length)]

radprof_g = np.zeros(length, dtype=list)  # empty list for the radial profile
radprof_r = np.zeros(length, dtype=list)  # empty list for the radial profile

# open images
im_g = fits.open(images[0])
im_r = fits.open(images[1])
datag = im_g[0].data
datar = im_r[0].data

# collapse box
for i in range(0,width):
    selg = datag[lin_y+i,col_x:col_x+length]
    radprof_g = [float(x + y) for x,y in zip(radprof_g, selg.tolist())]
    selr = datar[lin_y+i,col_x:col_x+length]
    radprof_r = [float(x + y) for x,y in zip(radprof_r, selr.tolist())]
    g_r19 = np.divide(radprof_g, radprof_r)
    g_r19 = [-2.5*np.log10(elem) for elem in g_r19]
    
##################

# PA = 232 deg
images = ["cig96_g_flux_dep_-153rot.fits", "cig96_r_flux_dep_-153rot.fits"]

# coordinates
col_x = 213
lin_y = 485

radprof_g = np.zeros(length, dtype=list)  # empty list for the radial profile
radprof_r = np.zeros(length, dtype=list)  # empty list for the radial profile

# open images
im_g = fits.open(images[0])
im_r = fits.open(images[1])
datag = im_g[0].data
datar = im_r[0].data

# collapse box
for i in range(0,width):
    selg = datag[lin_y+i,col_x:col_x+length]
    radprof_g = [float(x + y) for x,y in zip(radprof_g, selg.tolist())]
    selr = datar[lin_y+i,col_x:col_x+length]
    radprof_r = [float(x + y) for x,y in zip(radprof_r, selr.tolist())]
    g_r52 = np.divide(radprof_g, radprof_r)
    g_r52 = [-2.5*np.log10(elem) for elem in g_r52]


# In[6]:

print len(dist)
print len(g_r6)
print len(g_r25)
print len(g_r52)
print len(g_r19)


# In[7]:

akpc = 0.098417 # kpc/arcsec
pkpc = 0.103015 # kpc/pix


# In[ ]:

# g over r figure
fig, [ax0, ax1, ax2, ax3] = plt.subplots(nrows=4, ncols=1, sharex=True, squeeze=True, figsize=(8,12))
plt.subplots_adjust(left=0.125, bottom=0.125, right=0.95, top=0.9)

# plots
length = 180    # pixels
dist = [pkpc*elem for elem in np.arange(0,length)]
ax0.plot(dist, g_r6, linewidth=3,  label="PA = 6 deg",  color="g")
ax1.plot(dist, g_r25, linewidth=3, label="PA = 25 deg", color="C6")
plt.subplots_adjust(wspace=0.0, hspace=0.0)

length = 190    # pixels
dist = [pkpc*elem for elem in np.arange(0,length)]
ax2.plot(dist, g_r19, linewidth=3, label="PA = 199 deg", color='y')
ax3.plot(dist, g_r52, linewidth=3, label="PA = 232 deg", color='C0')


# disc lines
ax0.axvline(0, linestyle=':', color='k', alpha=0.3)
ax0.axvline(101*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
ax1.axvline(0, linestyle=':', color='k', alpha=0.3)
ax1.axvline(92*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
ax2.axvline(0, linestyle=':', color='k', alpha=0.3)
ax2.axvline(100*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')
ax3.axvline(0, linestyle=':', color='k', alpha=0.3)
ax3.axvline(106*pkpc, linestyle=':', color='k', alpha=0.3, label='Disc')

# pseudoring lines
ax0.axvline(111*pkpc, linestyle='-.', color='g', alpha=0.5)
ax0.axvline(131*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')
ax1.axvline(103*pkpc, linestyle='-.', color='g', alpha=0.5)
ax1.axvline(128*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')
ax2.axvline(115*pkpc, linestyle='-.', color='g', alpha=0.5)
ax2.axvline(137*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')
ax3.axvline(112*pkpc, linestyle='-.', color='g', alpha=0.5)
ax3.axvline(145*pkpc, linestyle='-.', color='g', alpha=0.5, label='Pseudoring')

# bulge lines
ax0.axvline(24.27*pkpc, linestyle='--', color='b', alpha=0.3, label='Bulge')
ax1.axvline(24.27*pkpc, linestyle='--', color='b', alpha=0.3, label='Bulge')
ax2.axvline(24.27*pkpc, linestyle='--', color='b', alpha=0.3, label='Bulge')
ax3.axvline(24.27*pkpc, linestyle='--', color='b', alpha=0.3, label='Bulge')

# hide X axis
ax0.tick_params(axis='x', length=0)
ax1.tick_params(axis='x', length=0)
ax2.tick_params(axis='x', length=0)
plt.xticks(np.arange(0,21,1))

# X and Y axis limits
ax0.set_xlim(0, 20, 1)
ax1.set_xlim(0, 20, 1)
ax2.set_xlim(0, 20, 1)
ax3.set_xlim(0, 20, 1)

#ax0.set_ylim(-0.8,1.5) # original values
#ax1.set_ylim(-0.8,1.5) # original values
#ax2.set_ylim(-0.8,2.2) # original values
#ax3.set_ylim(-0.1,1.5) # original values

ax0.set_ylim(-0.8,2.2)
ax1.set_ylim(-0.8,2.2)
ax2.set_ylim(-0.8,2.2)
ax3.set_ylim(-0.8,2.2)


####################

# SDSS red/blue VS green valley boundaries according to Walker et al. 2013
ax0.axhspan(0.55, 0.7, color='g', alpha=0.2, label="Green valley")
ax1.axhspan(0.55, 0.7, color='g', alpha=0.2, label="Green valley")
ax2.axhspan(0.55, 0.7, color='g', alpha=0.2, label="Green valley")
ax3.axhspan(0.55, 0.7, color='g', alpha=0.2, label="Green valley")

####################

# Labels, ticks and legend
ax0.set_ylabel("$g-r$ (mag)", fontsize=16)
ax0.set_yticklabels(np.arange(-1, 2.5, 0.5), fontsize=14)
ax0.annotate(s="PA = 6$\degree$", xy=(5.3, -0.6), fontsize=15)


ax1.set_ylabel("$g-r$ (mag)", fontsize=16)
ax1.set_yticklabels(np.arange(-1, 2.5, 0.5), fontsize=14)
ax1.annotate(s="PA = 16$\degree$", xy=(5.3, -0.6), fontsize=15)

ax2.set_ylabel("$g-r$ (mag)", fontsize=16)
ax2.set_yticklabels(np.arange(-1, 2.5, 0.5), fontsize=14)
ax2.annotate(s="PA = 30$\degree$", xy=(5.3, -0.5), fontsize=15)
ax2.axvspan(17.5, 18.4, color='red', alpha=0.2, label="")

ax3.set_ylabel("$g-r$ (mag)", fontsize=16)
ax3.set_yticklabels(np.arange(-1, 2.5, 0.5), fontsize=14)
ax3.set_xticklabels(np.arange(0, 21, 1), fontsize=14)
ax3.set_xlabel("Distance (kpc)", fontsize=16)
ax3.annotate(s="PA = 55$\degree$", xy=(5.3, 0.0), fontsize=15)

# grid
for graph in [ax0, ax1, ax2, ax3]:
    graph.grid(alpha=0.1)
    
# legend
#plt.legend(fontsize=8)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_radprofiles.png', bbox_inches='tight')


# In[274]:

plt.close()


# #5. Disc and pseudoring mean color
# 
# Displaying two areas over the disc (removing bright features and the inner **30 arcsec** of the bulge, according to **Espada+ 2011a**) and the pseudoring, respectively, the mean and median values of each one are:

# In[10]:

#from image:  cig96_g_over_r_flux_dep.fits
mean_disc = 0.50012729
median_disc = 0.49885345

mean_full_pseudoring = 1.4642415
median_full_pseudoring = 0.64481032

# bluer vs redder regions of pseudoring
mean_blue = 1.4634161
median_blue = 0.75250494

mean_red = 0.99061446
median_red = 0.53974068


# #6. Color index g-r, regions color and cuts

# In[5]:

# Select working directory
go_to('deep_peris')

# images
goverr = "cig96_g_minus_r_mag.fits"
regions = "cig96_RdRGB_SB_oriented.fits"
oriented = "cig96_g_flux_dep_-100rot.fits"

# reminder: color_gr variable
color_gr =  [0.6259716006291072, 0.36300770477647504, 0.307734966852955, 0.2355621588626633, 0.29228524203696793, 0.3526424828164778,             0.2790106947835227, 0.297460811984255, 0.3542155036889838, 0.2792155690809395, 0.19178883178920358, 0.30237066554180103,              0.4545786013884836, 0.3211585015780365, 0.45062802025530146, 0.23381907560446535, 0.5021268976871376, 0.5416916825805913,              0.6379768719413406, 0.73355361149655, 0.7439927192003424, 0.7276399374276785, 0.771708591029256, 0.9043116435588878,              0.5636875793307716, 0.7267345496621722, 0.49615698168135225, 0.47082067425720453, 1.00273286274248, 0.918469795711097,             0.6094727700934222, 0.4969637363756334, 0.7452621620966937]

# regions: read files and define x,y for both files
regfile = open("cig96_pseudo_reg_BandR_xyorient.reg")
x = []
y = []
for col in (raw.strip().split() for raw in regfile):
    x.append(float(col[0]))
    y.append(float(col[1]))

armsfile = open("cig96_pseudo_reg_arms_BandR_xyorient.reg")
a_x = []
a_y = []
for col in (raw.strip().split() for raw in armsfile):
    a_x.append(float(col[0]))
    a_y.append(float(col[1]))


# In[62]:

# figure
fig = plt.figure(figsize=(21, 7))
plt.axis('off')

# color index plot
gr = aplpy.FITSFigure(goverr, figure=fig, north=True, subplot=[0.06,  0.15, 0.29, 0.77])
gr.show_colorscale(stretch='linear', cmap='RdBu_r', vmin=0.25, vmax=1.0)
gr.axis_labels.set_font(size=16)

# grid and ticks
#gr.add_grid()
#gr.grid.set_color('w')
#gr.grid.set_alpha(0.6)
#gr.grid.show()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
gr.tick_labels.set_xformat('hh:mm:ss')
gr.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=2, color='w')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=2, color='w')

# arrows
plt.annotate('', xy=(310, 480), xytext=(260, 430), arrowprops=dict(facecolor='orange', shrink=0.0))
plt.annotate('', xy=(580, 275), xytext=(663, 275), arrowprops=dict(facecolor='cyan', shrink=0.0))

# profile cuts
gr_x = 471
gr_y = 504
plt.plot(gr_x, gr_y, 'wo', zorder=4)
plt.plot([gr_x,gr_x-10],  [gr_y,gr_y+195],  linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 6$\degree$")
plt.plot([gr_x,gr_x-38],  [gr_y,gr_y+190],  linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 16$\degree$")
plt.plot([gr_x,gr_x+98],  [gr_y,gr_y-245],  linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 210$\degree$")
plt.plot([gr_x,gr_x+135], [gr_y,gr_y-165],  linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 235$\degree$")

# colorbar
gr.add_colorbar(location='right', width=0.22, pad=0.0, axis_label_text="g - r (mag)")
gr.colorbar.set_axis_label_font(size=15)
gr.colorbar.set_font(size=14)

# g-r image: CIG96 coordinates in pix
RAgr = 500
DECgr = 500
rgr = 200

#--------------------------------------

# regions plot
im = fits.open(regions)
reg = aplpy.FITSFigure(im, auto_refresh=True, figure=fig, north=True, subplot=[0.375,  0.15, 0.29, 0.77], zorder=1)
reg.show_colorscale(stretch='log', cmap='bone', vmin=22, vmax=28)   # for RGB-image visualization
reg.axis_labels.set_font(size=16)
reg.axis_labels.hide_y()
reg.tick_labels.hide_y()

# grid and ticks
#reg.add_grid()
#reg.grid.set_color('k')
#reg.grid.set_alpha(0.3)
#reg.grid.show()
plt.xticks(fontsize=15)
reg.tick_labels.set_xformat('hh:mm:ss')
reg.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='k')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='k')


# regions overlay   
regs = plt.scatter(x, y, s=150, c=color_gr, cmap="RdBu_r", edgecolors='black', zorder=4)

# arms regions overlay
plt.scatter(a_x, a_y, s=75, marker='X', color='y', edgecolors='black', zorder=4)

# arrow and plot colorbar (do nor alter order with the next "arms regions overlay" or the plot will crash)
plt.annotate('', xy=(582, 305), xytext=(665, 305), arrowprops=dict(facecolor='cyan', shrink=0.0))
cbar = plt.colorbar(mappable=regs, fraction=0.07, pad=0.0)
cbar.set_label("g - r (mag)", fontsize=15)
cbar.ax.tick_params(labelsize=14)

# profile cuts
reg_x = 478
reg_y = 534
plt.plot(reg_x, reg_y, 'wo', zorder=4)
plt.plot([reg_x,reg_x-10],  [reg_y,reg_y+195],  linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 6$\degree$")
plt.plot([reg_x,reg_x-38],  [reg_y,reg_y+190],  linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 16$\degree$")
plt.plot([reg_x,reg_x+98],  [reg_y,reg_y-245],  linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 210$\degree$")
plt.plot([reg_x,reg_x+135], [reg_y,reg_y-165],  linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 235$\degree$")

#--------------------------------------

# deprojected plot
imor = fits.open(oriented)
dep = aplpy.FITSFigure(imor, auto_refresh=True, figure=fig, subplot=[0.695, 0.15, 0.29, 0.77])
dep.show_colorscale(stretch='linear', cmap='bone_r', vmin=1e-11, vmax=9e-11)   # for RGB-image visualization
dep.axis_labels.hide_x()
dep.axis_labels.hide_y()
dep.tick_labels.hide_x()
dep.tick_labels.hide_y()

# cuts position angles
PA1 = 25
PA2 = 6
PA3 = 227

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.09

# deprojected image: CIG96 coordinates in pix
RAx = 262
DECy = 600
r = 200

# profile cuts
plt.plot(RAx, DECy, 'wo', zorder=4)
plt.plot([RAx,RAx+35],[DECy,DECy+145], linewidth=2, color='forestgreen',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 6$\degree$")
plt.plot([RAx,RAx+2], [DECy,DECy+150], linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 16$\degree$")
plt.plot([RAx,RAx+50],[DECy,DECy-180], linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 30$\degree$")
plt.plot([RAx,RAx+95],[DECy,DECy-120], linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 55$\degree$")

# legend
plt.legend(fontsize=16)

# zooms
gr.recenter(RA, DEC, radius=rad-0.01)
reg.recenter(RA, DEC, radius=rad-0.01)
dep.recenter(RAx, DECy, radius=r)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_gr_pseudo_cuts_all.png', bbox_inches='tight')


# In[63]:

plt.close()


# In[58]:

# figure
fig = plt.figure(figsize=(14, 7))
plt.axis('off')

# color index plot
gr = aplpy.FITSFigure(goverr, figure=fig, north=True, subplot=[0.0,  0.02, 0.50, 0.95])
gr.show_colorscale(stretch='linear', cmap='RdBu_r', vmin=0.25, vmax=1.0)
gr.axis_labels.set_font(size=16)

# grid and ticks
#gr.add_grid()
#gr.grid.set_color('w')
#gr.grid.set_alpha(0.6)
#gr.grid.show()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
gr.tick_labels.set_xformat('hh:mm:ss')
gr.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=2, color='w')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=2, color='w')

# arrows
plt.annotate('', xy=(310, 480), xytext=(260, 430), arrowprops=dict(facecolor='orange', shrink=0.0))
plt.annotate('', xy=(580, 275), xytext=(663, 275), arrowprops=dict(facecolor='cyan', shrink=0.0))

# profile cuts
gr_x = 471
gr_y = 504
plt.plot(gr_x, gr_y, 'wo', zorder=4)
plt.plot([gr_x,gr_x-10],  [gr_y,gr_y+195],  linewidth=2, color='forestgreen',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 6$\degree$")
plt.plot([gr_x,gr_x-38],  [gr_y,gr_y+190],  linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 16$\degree$")
plt.plot([gr_x,gr_x+98],  [gr_y,gr_y-245],  linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 210$\degree$")
plt.plot([gr_x,gr_x+135], [gr_y,gr_y-165],  linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 235$\degree$")

# colorbar
gr.add_colorbar(location='right', width=0.22, pad=0.0, axis_label_text="g - r (mag)")
gr.colorbar.set_axis_label_font(size=15)
gr.colorbar.set_font(size=14)

# g-r image: CIG96 coordinates in pix
RAgr = 500
DECgr = 500
rgr = 200

#--------------------------------------

# regions plot
im = fits.open(regions)
reg = aplpy.FITSFigure(im, auto_refresh=True, figure=fig, north=True, subplot=[0.55,  0.02, 0.45, 0.95], zorder=1)
reg.show_colorscale(stretch='log', cmap='bone', vmin=22, vmax=28)   # for RGB-image visualization
reg.axis_labels.set_font(size=16)
reg.axis_labels.hide_y()
reg.tick_labels.hide_y()

# grid and ticks
#reg.add_grid()
#reg.grid.set_color('k')
#reg.grid.set_alpha(0.3)
#reg.grid.show()
plt.xticks(fontsize=15)
reg.tick_labels.set_xformat('hh:mm:ss')
reg.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='k')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='k')


# regions overlay   
regs = plt.scatter(x, y, s=150, c=color_gr, cmap="RdBu_r", edgecolors='black', zorder=4)

# arms regions overlay
plt.scatter(a_x, a_y, s=75, marker='X', color='y', edgecolors='black', zorder=4)

# arrow and plot colorbar (do nor alter order with the next "arms regions overlay" or the plot will crash)
plt.annotate('', xy=(582, 305), xytext=(665, 305), arrowprops=dict(facecolor='cyan', shrink=0.0))
cbar = plt.colorbar(mappable=regs, fraction=0.07, pad=0.0)
cbar.set_label("g - r (mag)", fontsize=15)
cbar.ax.tick_params(labelsize=14)

# profile cuts
reg_x = 478
reg_y = 534
plt.plot(reg_x, reg_y, 'wo', zorder=4)
plt.plot([reg_x,reg_x-10],  [reg_y,reg_y+195],  linewidth=2, color='forestgreen',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 6$\degree$")
plt.plot([reg_x,reg_x-38],  [reg_y,reg_y+190],  linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 16$\degree$")
plt.plot([reg_x,reg_x+98],  [reg_y,reg_y-245],  linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 210$\degree$")
plt.plot([reg_x,reg_x+135], [reg_y,reg_y-165],  linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 235$\degree$")

#--------------------------------------

# deprojected plot
#imor = fits.open(oriented)
#dep = aplpy.FITSFigure(imor, auto_refresh=True, figure=fig, subplot=[0.695, 0.15, 0.29, 0.77])
#dep.show_colorscale(stretch='linear', cmap='bone_r', vmin=1e-11, vmax=9e-11)   # for RGB-image visualization
#dep.axis_labels.hide_x()
#dep.axis_labels.hide_y()
#dep.tick_labels.hide_x()
#dep.tick_labels.hide_y()

# cuts position angles
PA1 = 25
PA2 = 6
PA3 = 227

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.09

# deprojected image: CIG96 coordinates in pix
#RAx = 262
#DECy = 600
#r = 200

# profile cuts
#plt.plot(RAx, DECy, 'wo', zorder=4)
#plt.plot([RAx,RAx+35],[DECy,DECy+145], linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 6$\degree$")
#plt.plot([RAx,RAx+2], [DECy,DECy+150], linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 16$\degree$")
#plt.plot([RAx,RAx+50],[DECy,DECy-180], linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 30$\degree$")
#plt.plot([RAx,RAx+95],[DECy,DECy-120], linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 55$\degree$")

# legend
plt.legend(fontsize=16)

# zooms
gr.recenter(RA, DEC, radius=rad-0.01)
reg.recenter(RA, DEC, radius=rad-0.01)
#dep.recenter(RAx, DECy, radius=r)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_gr_pseudo_cuts_all_two.png', bbox_inches='tight')
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_gr_pseudo_cuts_all_two.eps', bbox_inches='tight')


# In[59]:

plt.close()


# In[60]:

# figure
fig = plt.figure(figsize=(10, 10))
plt.axis('off')

# deprojected plot alone
imor = fits.open(oriented)
dep = aplpy.FITSFigure(imor, auto_refresh=True, figure=fig)
dep.show_colorscale(stretch='linear', cmap='bone_r', vmin=1e-11, vmax=9e-11)   # for RGB-image visualization
dep.axis_labels.hide_x()
dep.axis_labels.hide_y()
dep.tick_labels.hide_x()
dep.tick_labels.hide_y()

# deprojected image: CIG96 coordinates in pix
RAx = 262
DECy = 600
r = 200

# profile cuts
plt.plot(RAx, DECy, 'wo', zorder=4)
plt.plot([RAx,RAx+35],[DECy,DECy+145], linewidth=2, color='forestgreen',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 6$\degree$")
plt.plot([RAx,RAx+2], [DECy,DECy+150], linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 16$\degree$")
plt.plot([RAx,RAx+50],[DECy,DECy-180], linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 30$\degree$")
plt.plot([RAx,RAx+95],[DECy,DECy-120], linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 55$\degree$")

# legend
plt.legend(fontsize=16)

# zooms

dep.recenter(RAx, DECy, radius=r)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_gr_pseudo_cuts_all_deproj.png', bbox_inches='tight')
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_gr_pseudo_cuts_all_deproj.eps', bbox_inches='tight')


# In[56]:

plt.close()


# In[39]:

# figure
fig = plt.figure(figsize=(7, 21))
plt.axis('off')

# color index plot
gr = aplpy.FITSFigure(goverr, figure=fig, north=True, subplot=[0.0,  0.66, 1.0, 0.33]) # subplot=[x_original_position, y_original_position, width, height]
gr.show_colorscale(stretch='linear', cmap='afmhot', vmin=0.6, vmax=1.5)
gr.axis_labels.set_font(size=16)

# grid and ticks
#gr.add_grid()
#gr.grid.set_color('w')
#gr.grid.set_alpha(0.6)
#gr.grid.show()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
gr.tick_labels.set_xformat('hh:mm:ss')
gr.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=2, color='w')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=2, color='w')

# arrows
plt.annotate('', xy=(310, 480), xytext=(260, 430), arrowprops=dict(facecolor='orange', shrink=0.0))
plt.annotate('', xy=(580, 275), xytext=(663, 275), arrowprops=dict(facecolor='cyan', shrink=0.0))

# profile cuts
gr_x = 471
gr_y = 504
plt.plot(gr_x, gr_y, 'wo', zorder=4)
plt.plot([gr_x,gr_x-10],  [gr_y,gr_y+195],  linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 6$\degree$")
plt.plot([gr_x,gr_x-38],  [gr_y,gr_y+190],  linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 16$\degree$")
plt.plot([gr_x,gr_x+98],  [gr_y,gr_y-245],  linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 210$\degree$")
plt.plot([gr_x,gr_x+135], [gr_y,gr_y-165],  linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 235$\degree$")

# colorbar
#gr.add_colorbar(location='right', width=0.22, pad=0.0, axis_label_text="g - r (mag)")
#gr.colorbar.set_axis_label_font(size=15)
#gr.colorbar.set_font(size=14)

# g-r image: CIG96 coordinates in pix
RAgr = 500
DECgr = 500
rgr = 200

#--------------------------------------

# regions plot
im = fits.open(regions)
reg = aplpy.FITSFigure(im, auto_refresh=True, figure=fig, north=True, subplot=[0.0,  0.33, 1.0, 0.33], zorder=1)
reg.show_colorscale(stretch='log', cmap='bone', vmin=22, vmax=28)   # for RGB-image visualization
reg.axis_labels.set_font(size=16)
reg.axis_labels.hide_y()
reg.tick_labels.hide_y()

# grid and ticks
#reg.add_grid()
#reg.grid.set_color('k')
#reg.grid.set_alpha(0.3)
#reg.grid.show()
plt.xticks(fontsize=15)
reg.tick_labels.set_xformat('hh:mm:ss')
reg.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='k')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='k')


# regions overlay   
regs = plt.scatter(x, y, s=150, c=color_gr, cmap="RdBu_r", edgecolors='black', zorder=4)

# arms regions overlay
plt.scatter(a_x, a_y, s=75, marker='X', color='y', edgecolors='black', zorder=4)

# arrow and plot colorbar (do nor alter order with the next "arms regions overlay" or the plot will crash)
plt.annotate('', xy=(582, 305), xytext=(665, 305), arrowprops=dict(facecolor='cyan', shrink=0.0))
cbar = plt.colorbar(mappable=regs, fraction=0.07, pad=0.0)
cbar.set_label("g - r (mag)", fontsize=15)
cbar.ax.tick_params(labelsize=14)

# profile cuts
reg_x = 478
reg_y = 534
plt.plot(reg_x, reg_y, 'wo', zorder=4)
plt.plot([reg_x,reg_x-10],  [reg_y,reg_y+195],  linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 6$\degree$")
plt.plot([reg_x,reg_x-38],  [reg_y,reg_y+190],  linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 16$\degree$")
plt.plot([reg_x,reg_x+98],  [reg_y,reg_y-245],  linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 210$\degree$")
plt.plot([reg_x,reg_x+135], [reg_y,reg_y-165],  linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 235$\degree$")

#--------------------------------------

# deprojected plot
imor = fits.open(oriented)
dep = aplpy.FITSFigure(imor, auto_refresh=True, figure=fig, subplot=[0.0, 0.0, 1.0, 0.33])
dep.show_colorscale(stretch='linear', cmap='bone_r', vmin=1e-11, vmax=9e-11)   # for RGB-image visualization
dep.axis_labels.hide_x()
dep.axis_labels.hide_y()
dep.tick_labels.hide_x()
dep.tick_labels.hide_y()

# cuts position angles
PA1 = 25
PA2 = 6
PA3 = 227

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.09

# deprojected image: CIG96 coordinates in pix
RAx = 262
DECy = 600
r = 200

# profile cuts
plt.plot(RAx, DECy, 'wo', zorder=4)
plt.plot([RAx,RAx+35],[DECy,DECy+145], linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 6$\degree$")
plt.plot([RAx,RAx+2], [DECy,DECy+150], linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 16$\degree$")
plt.plot([RAx,RAx+50],[DECy,DECy-180], linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 30$\degree$")
plt.plot([RAx,RAx+95],[DECy,DECy-120], linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 55$\degree$")

# legend
plt.legend(fontsize=16)

# zooms
gr.recenter(RA, DEC, radius=rad-0.01)
reg.recenter(RA, DEC, radius=rad-0.01)
dep.recenter(RAx, DECy, radius=r)

# save figure
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_gr_pseudo_cuts_all_vert.png', bbox_inches='tight')


# In[40]:

plt.close()


# In[30]:

# figure
fig = plt.figure(figsize=(7, 7))
plt.axis('off')

# color index plot
gr = aplpy.FITSFigure(goverr, figure=fig, north=True, subplot=[0.0,  0.0, 1.0, 1.0]) # subplot=[x_original_position, y_original_position, width, height]
gr.show_colorscale(stretch='linear', cmap='RdBu', vmin=0.6, vmax=1.5)
gr.axis_labels.set_font(size=16)


# grid and ticks
#gr.add_grid()
#gr.grid.set_color('w')
#gr.grid.set_alpha(0.6)
#gr.grid.show()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
gr.tick_labels.set_xformat('hh:mm:ss')
gr.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=2, color='w')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=2, color='w')

# arrows
plt.annotate('', xy=(310, 480), xytext=(260, 430), arrowprops=dict(facecolor='orange', shrink=0.0))
plt.annotate('', xy=(580, 275), xytext=(663, 275), arrowprops=dict(facecolor='cyan', shrink=0.0))

# profile cuts
gr_x = 471
gr_y = 504
plt.plot(gr_x, gr_y, 'wo', zorder=4)
plt.plot([gr_x,gr_x-10],  [gr_y,gr_y+195],  linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 6$\degree$")
plt.plot([gr_x,gr_x-38],  [gr_y,gr_y+190],  linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 16$\degree$")
plt.plot([gr_x,gr_x+98],  [gr_y,gr_y-245],  linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)   #, label="PA = 210$\degree$")
plt.plot([gr_x,gr_x+135], [gr_y,gr_y-165],  linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], zorder=3)  #, label="PA = 235$\degree$")

# colorbar
gr.add_colorbar(location='right', width=0.22, pad=0.0, axis_label_text="g - r (mag)")
gr.colorbar.set_axis_label_font(size=15)
gr.colorbar.set_font(size=14)

# g-r image: CIG96 coordinates in pix
RAgr = 500
DECgr = 500
rgr = 200


# cuts position angles
PA1 = 25
PA2 = 6
PA3 = 227

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.09

# deprojected image: CIG96 coordinates in pix
RAx = 262
DECy = 600
r = 200

# profile cuts
plt.plot(RAx, DECy, 'wo', zorder=4)
plt.plot([RAx,RAx+35],[DECy,DECy+145], linewidth=2, color='g',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 6$\degree$")
plt.plot([RAx,RAx+2], [DECy,DECy+150], linewidth=2, color='C6', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 16$\degree$")
plt.plot([RAx,RAx+50],[DECy,DECy-180], linewidth=2, color='y',  path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 30$\degree$")
plt.plot([RAx,RAx+95],[DECy,DECy-120], linewidth=2, color='C0', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()], label="PA = 55$\degree$")

# legend
plt.legend(fontsize=16)

# zooms
gr.recenter(RA, DEC, radius=rad-0.01)
reg.recenter(RA, DEC, radius=rad-0.01)
dep.recenter(RAx, DECy, radius=r)

# save figure
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_gr_pseudo_cuts_all_vert.png', bbox_inches='tight')


# In[31]:

plt.close()


# In[ ]:



