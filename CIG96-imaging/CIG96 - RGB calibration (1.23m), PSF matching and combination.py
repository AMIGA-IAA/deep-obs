
# coding: utf-8

# # CIG96 - RGB calibration, PSF matching and combination with R-deep
# 
# ##GENERAL INFORMATION
# 
# Notebook with the calibration of the R,G,B images by Peris with the *g* and *r* images from SDSS.
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
from   astropy import units as u
from   astropy import wcs as wcs
import astropy.io.fits as fits
import csv
from   glob import glob as ls
import imexam
import lemon.astrometry as astrometry
import lemon.photometry as photometry
import matplotlib as mpl
import matplotlib.colors
import matplotlib.pyplot as plt
import montage_wrapper as montage
import numpy as np
from   numpy.linalg import eig, inv
import os
import pylab
import pyraf
import pyregion as pyreg
import repipy.combine as combine
import repipy.find_sky as find_sky
import scipy
from   scipy import stats
import shutil
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")
#from   kapteyn import maputils
#import pandas
#import pickle
#import pyds9 as ds9
#import random
#from   spectral_cube import SpectralCube as spc
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


# ###b) "Go to" function

# In[3]:

def go_to(dir):
    """Function to select and switch to optical or HI folders"""
    # options
    opt = ["opt","Opt","OPT"]
    peris = ["per", "peris", "Peris"]
    hi = ["HI", "Hi", "hi"]
    
    if dir in opt:
        # raw data and work directories:
        raw = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16/cig"
        work = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16/cig"
    
    elif dir in peris:
        # raw data and work directories:
        raw = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16/deep_peris/"
        work = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16/deep_peris/"
        
    elif dir in hi:
        # raw data and work directories:
        raw = "/home/prm/Desktop/evla-vla-2015"
        work = "/home/prm/Desktop/evla-vla-2015"
        
    else:
        return "Invalid option, please select 'opt' or 'HI' or 'peris'."

    # switch to selected working directory
    os.chdir(os.path.expandvars(work))
    print "Work/save directory: ", work


# In[4]:

go_to('peris')


# #0. Functions

# In[5]:

# function definitions

##########

def imex_coord_mag(coordfile, magsfile, sky_rms, outputf):
    '''This function will extract the RA, DEC (deg) coordinates from a standard file,
    the mangitude with error from a typical imexam file and save the 4 parameters in an output file.'''
    
    with open(coordfile) as coordsf, open(magsfile) as magsf, open(outputf, 'w') as outfile:
        ra = []    # RA coordinate
        dec = []   # DEC coordinate
        f = []     # flux
        m = []     # empty lists to append later to the output file
        e = []     # error
        
        data_coords = coordsf.readlines()   
        data_mags = magsf.readlines()[2:]   # read from line 3 on (2 in pyton count), skipping column names
        
        for line in data_coords:
            parts = line.split()
            ra.append(parts[0])   # ra coordinate
            dec.append(parts[1])  # dec coordinate
                
        for line in data_mags:       # extract magnitudes only from even lines (they have 79 characters, we take the magnitudes from these 
            if len(line) > 34:
                parts = line.split()
                m.append(parts[1])   # magnitude
                f.append(parts[2])   # flux
                err = np.abs(-2.5/np.log(10)*float(sky_rms)/float(parts[2]))
                e.append(round(err,5))
                
        writer = csv.writer(outfile, delimiter='\t')  
        writer.writerows(zip(ra,dec,m,e))   # write file to specified output name
        
##########

def clean_indef(files_with_indefs, files_to_clean, output_files):
    '''This function will point to the line number where any INDEF value is found
       in the magnitude column, delete the line in all files involved and 
       save the outputs in a new file.'''
    
    indef = 'INDEF'  # condition to identify lines
    del_lines = []   # list of lines to be deleted
    
    for file in files_with_indefs:                   # 1st, identify lines with INDEF values and list the lines only once
        with open(file, 'r') as infile:
            for num, line in enumerate(infile, 1):
                if indef in line:
                    if num not in del_lines:
                        del_lines.append(num)         
    
    # order and print lines to delete
    del_lines.sort()
    print ">> The following lines contain 'INDEF' values:", del_lines
    print ">> Deleting", len(del_lines), "lines..."
    
    for infile, outfile in zip(files_to_clean, output_files):                      # 2n, delete the lines for all files
        with open(infile, 'r') as readfile, open(outfile, "w+") as out:
            for num, line in enumerate(readfile,1):
                if num not in del_lines:  # si no esta en mi lista negra...
                    out.write(line)
    
    print ">> DONE:", len(del_lines), "line/s have been deleted from the specified files."

##########

# rms rough calculation - must be refined because as it is now, it takes all pixels into account galaxy and stars included
def rms_RGB(image):
    im = fits.open(image)
    data = im[0].data
    return np.sqrt(np.mean(np.square(data)))

##########

# extract magnitudes from SDSS catalog, generated by Aladin
def mag_from_cat(inputf, outputf):
    '''Extracts RA, DEC, r, g, err_r, err_g from a SDSS catalog generated by Aladin
    and saves them in an output file'''
    
    with open(inputf) as infile, open(outputf, 'w') as outfile:
        ra = []      # RA coordinate (deg)
        dec = []     # DEC coordinate (deg)
        rm = []      # r magnitude (mag)
        gm = []      # g magnitude (mag)
        rme = []     # r error (mag)
        gme = []     # g error (mag)
        
        data = infile.readlines()[2:]   # read from line 3 on (2 in pyton count), skipping column names
        
        for line in data:
                parts = line.split()
                ra.append(parts[7])   # ra coordinate
                dec.append(parts[8])  # dec coordinate
                gm.append(parts[10])
                rm.append(parts[11])
                gme.append(parts[15])
                rme.append(parts[16])
                
        writer = csv.writer(outfile, delimiter='\t')  
        writer.writerows(zip(ra,dec,rm,rme,gm,gme))   # write file to specified output name
        
##########

# matching SDSS catalog with RGB data

def match_SDSS_RGB(SDSSfile, RGBfile, outfile):
    '''This function will match the stars from the two sources, the SDSS catalog and the R,G,B files.'''
    
    with open(SDSSfile) as sdssf, open(RGBfile) as rgbf, open(outfile, "w+") as finalfile:
        
        datasdss = sdssf.readlines()
        datargb  = rgbf.readlines()
        rgb_ra  = []
        rgb_dec = []
        rgb_mag = []
        rgb_err = []
 
        for line in datargb:            # opens R/G/B data file and stores coordinated in two lists
            parts = line.split()
            rgb_ra.append(parts[0])
            rgb_dec.append(parts[1])
            rgb_mag.append(parts[2])
            rgb_err.append(parts[3])
            
        for line in datasdss:           # opens SDSS catalog file
            parts = line.split()
            ra  = parts[0]
            dec = parts[1]
            for rgbra, rgbdec, rgbmag, rgberr in zip(rgb_ra, rgb_dec, rgb_mag, rgb_err):      # compares the distance between the R/G/B coordinates and
                dist = np.sqrt((float(ra) - float(rgbra))**2 + (float(dec) - float(rgbdec))**2)
                #if dist < 2:
                if dist < 0.0005:
                    finalfile.write(str(rgbra) + "   " + str(rgbdec) + "   " + str(rgbmag) + "   " + str(rgberr) + "   " + line)
                    


# To combine the R,G,B images with the R deep image, we need to degrade the images with best seeing or smaller PSF (the R,G,B images) to match the PSF of the worst seeing image or larger PSF (R deep).
# 
# For this we need to:
# 
# - extract and subtract sky from all images
# - align R,G,B and Rdeep
# - calibrate R,G,B
# - degrade PSF
# - combine R,G,B and Rdeep
# 
# #1. R,G,B and Rdeep sky subtraction
# 
# ##1.1 Cropping
# 
# First, we crop the deep image to remove the FoV-off regions and have a more accurate sky estimation.

# In[14]:

# cropping of cig96 (2.2mCAHA) image to remove the too large region outside FoV for a better sky estimation

# Select working directory
go_to('opt')

# read the image
image = fits.open('cig96_def_crop_wcs.fits', mode='update')
#image.info()

# select layer (in this case, it is unique)
im = image[0]
im.shape

# crop image from x0 to x1 and from y0 to y1
im.data = im.data[150:850,150:850]

im.writeto('cig96_def_crop_skycrop.fits', clobber=True)


# In[89]:

get_ipython().run_cell_magic(u'bash', u'', u'# find sky via repipy routine over cropped deep image\nfind_sky.py cig96_def_crop_skycrop.fits')


# ##1.2 Sky value estimation
# 
# Second, we extract the sky value of the double-cropped image and store it for later use.
# 
# We make a copy of the deep image ("cig96_def_crop_wcs.fits" version) to the new work folder.

# In[90]:

# Select working directory
go_to('opt')

# sky value extraction from double-cropped image
im = fits.open("cig96_def_crop_skycrop.fits")
skycrop = float(im[0].header["SKY"])

# copy to new folder
shutil.copy2("cig96_def_crop_skycrop.fits", "../deep_peris/")
shutil.copy2("cig96_def_crop_wcs.fits", "../deep_peris/")


# And we find the sky for Peris images:

# In[91]:

# Select working directory
go_to('peris')


# In[92]:

get_ipython().run_cell_magic(u'bash', u'', u'# find sky via repipy routine over Peris images\nfind_sky.py *_koords.fits')


# ##1.3 Sky subtraction
# 
# Third, subtract sky from all images.

# In[54]:

# Select working directory
go_to('peris')

# makes a list containing all fits images
image_list = ["cig96_def_crop_wcs.fits", "CIG96_R_koords.fits", "CIG96_G_koords.fits", "CIG96_B_koords.fits"]

for image in image_list:
    if image  == "cig96_def_crop_wcs.fits":
        im = fits.open(image)
        im[0].data = im[0].data - skycrop     # for this image, biased by the big FoV, we subtract double-cropped image sky value
        im.writeto(image[0:-5] + "_skysub.fits", clobber=False)
        print image, " - sky = ", skycrop
    else:
        im = fits.open(image)
        sky = float(im[0].header["SKY"])
        im[0].data = im[0].data - sky
        im.writeto(image[0:-5] + "_skysub.fits", clobber=False)
        print image, " - sky = ", sky


                ###skipped Normalization

Fourth, we normalize images to their medians.
                
                # Select working directory
go_to('peris')

# makes a list containing all fits images
image_list = ls(os.path.join("*_skysub.fits"))

# extract median value of double-cropped and skysubtracted image
im = fits.open("cig96_def_crop_skycrop_skysub.fits")
median_crop = float(np.median(im[0].data))
print median_crop

# median normalization
for image in image_list:
    if image == "cig96_def_crop_wcs_skysub.fits":
        im = fits.open(image)
        print image, median_crop      # in this case, we divide by the median value of the double-cropped image
        im[0].data = im[0].data / median_crop
        im.writeto(image[0:-5] + "_norm.fits", clobber=True)
    else:
        im = fits.open(image)
        median = float(np.median(im[0].data))
        print image, median
        im[0].data = im[0].data / median
        im.writeto(image[0:-5] + "_norm.fits", clobber=True)
                
# #2. Align R,G,B and Rdeep images
# 
# R,G,B Peris' images are already aligned. We have to center deep image with the rest.
# 
# Options:
# - Via IRAF
# 
# - Via <a href="http://obswww.unige.ch/~tewes/alipy/">Alipy</a> module in Pique's notebook: <a href="http://amiga.iaa.es:9999/ipython/notebooks/pique/AlignFITS/AlignFITS-PabloApril.ipynb">Align FITS</a>.
# 
# - <a href="http://montage-wrapper.readthedocs.io/en/v0.9.5/_generated/montage_wrapper.wrappers.reproject.html#montage_wrapper.wrappers.reproject">Montage: reproject</a> task, powered by Montage + AstroPy.
# 
# *** Warning***: montage does not align correctly.

                ### MONTAGE

# Select working directory
go_to('peris')

R = "CIG96_R_SB.fits"
G = "CIG96_G_SB.fits"
B = "CIG96_B_SB.fits"
ref = "cig96_def_SB.fits"

reference = fits.open(ref)
header = reference[0].header

#montage.wrappers.reproject([R,G,B,ref], [R[:-5]+'_al.fits', G[:-5]+'_al.fits', B[:-5]+'_al.fits', ref[:-5]+'_al.fits'])
montage.wrappers.reproject([R,ref], [R[:-5]+'_al.fits', ref[:-5]+'_al.fits'])

#### IRAF

#### ALIPY
WARNING: to download the resulting fits, you must do it in a terminal, logged as 'amiga' user.
                
# #3. Calibration of R,G,B images
# 
# ##3.1 Stars selection (via IRAF) and photometry
# 
# Via *IRAF + ds9*,  we select the unsaturated stars in R,G,B images:
# 
# - imexam + r: to see and select the unsaturated stars with enough counts and a good Moffat fit,
# - imexam + a: over the selected stars to extract the magnitude (-2.5log10(flux) + 25 mag), repeated for the 3 images (R,G,B)
# 
# From the manual selection, we create a **RA, DEC (in degrees) image coordinates file** as well as an **x,y image coordinates file**  for the R,G,B images:
# 
# > *calib_stars_RGB_wcs_deg.reg*
# 
# > *calib_stars_RGB_xy.reg*
# 
# We use the x,y coordinates file to extract the magnitudes of the R,G,B sky subtracted images via *imexam* photometry task:
# 
# > *imexam_RGB_skysub.cl*
# 
# It runs a display task over each image, followed by an imexam taken over the 'imagecu' coordinates file.
# 
# **NOTE**: *imexam_RGB_skysub.cl* may take a few minutes to finish if it finds some stars that do not converge.

                display CIG96_R_koords_skysub.fits 1
imexam image=CIG96_R_koords_skysub.fits imagecu=calib_stars_RGB_xy.dat defkey=a use_dis=no > calib_R.dat

display CIG96_G_koords_skysub.fits 1
imexam image=CIG96_G_koords_skysub.fits imagecu=calib_stars_RGB_xy.dat defkey=a use_dis=no > calib_G.dat

display CIG96_B_koords_skysubb.fits 1
imexam image=CIG96_B_koords_skysub.fits imagecu=calib_stars_RGB_xy.dat defkey=a use_dis=no > calib_B.dat
                
# There are several versions of this .cl script, depending on the number of stars used for the calibration:
# 
# > imexam_RGB_skysub_10.cl
# 
# > imexam_RGB_skysub_24.cl
# 
# > imexam_RGB_skysub_62.cl
# 
# > imexam_RGB_skysub_67.cl
# 
# ##3.2. Extraction of coordinates and magnitudes
# 
# Since the imexam output is complex, we take only x,y and mangitudes and save them in easily manageable files via '***imex_coord_mag***' function.
# 
# We get rid of INDEF values from unsuccessful Moffat fitting attempts using '***clean_indef***' function.
# 
# ###3.2.1 Imexam data extraction to simple files

                rmsR = rms_RGB("CIG96_R_koords_skysub.fits")
rmsG = rms_RGB("CIG96_G_koords_skysub.fits")
rmsB = rms_RGB("CIG96_B_koords_skysub.fits")
print rmsR, rmsG, rmsB

### these result as imprecise measurements because it does take into account all pixels - must be refined before using it
                
# In[93]:

# median sky value from NON-sky-subtracted images (obtained via ds9)
rms_R = 0.0015519347
rms_G = 0.00081506655
rms_B = 0.00047444035

# variables: number of stars, filters and rms
starcount = [10, 24, 62, 67]
filterRGB = ["R", "G", "B"]
rmslist = [rms_R, rms_G, rms_B]

# run imex_coord_mag extraction function over all permutations
for stars in starcount:
    for fil,rms in zip(filterRGB,rmslist):
        # files/variables to insert
        coord = "calib_stars_RGB_wcs_" + str(stars) + ".reg"
        mags = "calib_" + fil + "_" + str(stars) + ".dat"
        rms = rms
        output = "calib_coords_mag_" + fil + "_" + str(stars) + ".dat"
        # run function
        imex_coord_mag(coord, mags, rms, output)


# ###3.2.2  Removal of INDEF values
# 
# Similar process to what is done above but in this case, every line with an INDEF mangitude will be deleted from all files to perform a precise calibration.

# In[94]:

# files who have INDEF values
stars = [10, 24, 62, 67]
for number in stars:
    filesindef = ["calib_coords_mag_R_" + str(number) + ".dat", "calib_coords_mag_G_" + str(number) + ".dat", "calib_coords_mag_B_" + str(number) + ".dat"]
    # files to clean from INDEF: all of them
    filestoclean = filesindef
    # output names
    outfiles = ["calib_coords_mag_R_clean_" + str(number) + ".dat", "calib_coords_mag_G_clean_" + str(number) + ".dat", "calib_coords_mag_B_clean_" + str(number) + ".dat"]
    # Clean files to avoid any undefined values
    clean_indef(filesindef, filestoclean, outfiles)


# ##3.3. SDSS r and g data
# 
# From Aladin, we extract a catalog with all the stars located in the 15-20 arcmin around the galaxy. The next function will extract the RA, DEC, r, g, err_r, err_g columns we need.

# In[95]:

# Aladin catalog file and output file
inf = "cig96_r_g_sources.txt"
outf = "cig96_SDSS_coords_r_g_errs.txt"

# extract coordinates, magnitudes and errors
mag_from_cat(inf, outf)

print outf


# ##3.4. Filters correspondence and linear fitting
# 
# According to Peris:
# 
# - R: 5800-7000 ---> corresponds to r_SDSS
# - G: 4900-5800 ---> corresponds to g_SDSS and r_SDSS
# - B: 3900-5100 ---> corresponds to g_SDSS
# 
# We match the stars selected in R,G,B images with the SDSS catalog.

# There are four different reference stars files with 10, 24, 62 and 67 stars each. To select any of these, we simply specify the number of stars in the next variable:

# In[96]:

stars = 67

match_SDSS_RGB("cig96_SDSS_coords_r_g_errs.txt", "calib_coords_mag_R_clean_" + str(stars) + ".dat", "calibrated_R_SDSS.dat")
match_SDSS_RGB("cig96_SDSS_coords_r_g_errs.txt", "calib_coords_mag_G_clean_" + str(stars) + ".dat", "calibrated_G_SDSS.dat")
match_SDSS_RGB("cig96_SDSS_coords_r_g_errs.txt", "calib_coords_mag_B_clean_" + str(stars) + ".dat", "calibrated_B_SDSS.dat")


# In[20]:

# Peris mag(R,G,B) and errors
Rmag = np.genfromtxt("calibrated_R_SDSS.dat")[:,2]
Rerr = np.genfromtxt("calibrated_R_SDSS.dat")[:,3]

Gmag = np.genfromtxt("calibrated_G_SDSS.dat")[:,2]
Gerr = np.genfromtxt("calibrated_G_SDSS.dat")[:,3]

Bmag = np.genfromtxt("calibrated_B_SDSS.dat")[:,2]
Berr = np.genfromtxt("calibrated_B_SDSS.dat")[:,3]

# SDSS mags and errors
rmag = np.genfromtxt("calibrated_R_SDSS.dat")[:,6]
rerr = np.genfromtxt("calibrated_R_SDSS.dat")[:,7]
gmag = np.genfromtxt("calibrated_R_SDSS.dat")[:,8]
gerr = np.genfromtxt("calibrated_R_SDSS.dat")[:,9]

# r,g combination for G filter
rgmag02 = 0.2*rmag + 0.8*gmag
rgerr02 = 0.2*rerr + 0.8*gerr

rgmag05 = 0.5*rmag + 0.5*gmag
rgerr05 = 0.5*rerr + 0.5*gerr

rgmag08 = 0.8*rmag + 0.2*gmag
rgerr08 = 0.8*rerr + 0.2*gerr

############

# figure and plots
fig = plt.figure(1, figsize=(21,6))
plt.suptitle("R, G and B transformation to r_SDSS and g_SDSS", fontsize=16)
plt.plot(bbox_inches="tight")
ax1 = plt.subplot(1,3,1)
ax2 = plt.subplot(1,3,2)
ax3 = plt.subplot(1,3,3)

############

# filtering of too shallow magnitudes and too large errorbars
sigmaerr  = 0.5
topmagRGB = 31
topmagrg  = 22.4

# R
fil_R  = Rmag[(Rerr < sigmaerr) & (Rmag < topmagRGB) & (rmag < topmagrg)]
fil_r  = rmag[(Rerr < sigmaerr) & (Rmag < topmagRGB) & (rmag < topmagrg)]
fil_Re = Rerr[(Rerr < sigmaerr) & (Rmag < topmagRGB) & (rmag < topmagrg)]
fil_re = rerr[(Rerr < sigmaerr) & (Rmag < topmagRGB) & (rmag < topmagrg)]
# G
fil_G     = Gmag[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag05 < topmagrg)]
fil_rg02  = rgmag02[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag02 < topmagrg)]
fil_rg05  = rgmag05[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag05 < topmagrg)]
fil_rg08  = rgmag08[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag08 < topmagrg)]
fil_Ge    = Gerr[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag05 < topmagrg)]
fil_rg02e = rgerr02[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag02 < topmagrg)]
fil_rg05e = rgerr05[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag05 < topmagrg)]
fil_rg08e = rgerr08[(Gerr < sigmaerr) & (Gmag < topmagRGB) & (rgmag08 < topmagrg)]
# B
fil_B  = Bmag[(Berr < sigmaerr) & (Bmag < topmagRGB) & (gmag < topmagrg)]
fil_g  = gmag[(Berr < sigmaerr) & (Bmag < topmagRGB) & (gmag < topmagrg)]
fil_Be = Berr[(Berr < sigmaerr) & (Bmag < topmagRGB) & (gmag < topmagrg)]
fil_ge = gerr[(Berr < sigmaerr) & (Bmag < topmagRGB) & (gmag < topmagrg)]

############

# filtered data: weighting after converting error to sigma (1/err)
fil_Rew = [1/elem for elem in fil_Re]
fil_Gew = [1/elem for elem in fil_Ge]
fil_Bew = [1/elem for elem in fil_Be]

# linear fitting
# R 
fitRw = np.polyfit(fil_r, fil_R, 1, w = fil_Rew)
# G 
#fitG02w = np.polyfit(fil_rg02, fil_G, 1, w = fil_Gew)
fitG05w = np.polyfit(fil_rg05, fil_G, 1, w = fil_Gew)
#fitG08w = np.polyfit(fil_rg08, fil_G, 1, w = fil_Gew)
# B 
fitBw = np.polyfit(fil_g, fil_B, 1, w = fil_Bew)

############

# fitting data
Rslope, Rintercept, Rr_value, Rp_value, Rstd_err = scipy.stats.linregress(fil_r, fil_R)
Gslope05, Gintercept05, Gr_value05, Gp_value05, Gstd_err05 = scipy.stats.linregress(fil_rg05, fil_G)
Bslope, Bintercept, Br_value, Bp_value, Bstd_err = scipy.stats.linregress(fil_g, fil_B)
Rr_value = Rr_value**2
Gr_value05 = Gr_value05**2
Br_value = Br_value**2

print "R slope = ", fitRw[0], "- R intercept = ", fitRw[1]
print "G slope = ", fitG05w[0], "- R intercept = ", fitG05w[1]
print "B slope = ", fitBw[0], "- R intercept = ", fitBw[1]

print "R slope error =", Rstd_err, "- R2 =",   Rr_value
print "G slope error =", Gstd_err05, "- R2 =", Gr_value05
print "B slope error =", Bstd_err, "- R2 =",   Br_value

# errors
Rmagerr = [Rstd_err*elem for elem in fil_Rew]
Gmagerr = [Gstd_err05*elem for elem in fil_Gew]
Bmagerr = [Bstd_err*elem for elem in fil_Bew]
print "err_mag(R) =", np.median(Rmagerr)
print "err_mag(G) =", np.median(Gmagerr)
print "err_mag(B) =", np.median(Bmagerr)

############

# R,G,B data plot
rg02color = (0.1, 0.9, 0.1)
rg05color = (0.2, 0.7, 0.2)
rg08color = (0.4, 0.6, 0.2)

ax1.errorbar(fil_r, fil_R, fmt='o', markersize=5, markeredgecolor='brown', color='r', xerr=fil_re, yerr=fil_Re, ecolor='brown', label='')
ax2.errorbar(fil_rg05, fil_G, fmt='o', markersize=5, color=rg05color, markeredgecolor=rg05color, xerr=fil_rg05e, yerr=fil_Ge, ecolor=rg05color,              label='50% r$_{SDSS}$ + 50% g$_{SDSS}$')
ax3.errorbar(fil_g, fil_B, fmt='o', markersize=5, markeredgecolor=(0.1, 0.4, 0.9), color='blue', xerr=fil_ge, yerr=fil_Be, ecolor=(0.1, 0.4, 0.7), label='')

############

# R,G,B linear fitting plot
ax1.plot(fil_r, fitRw[0]*fil_r + fitRw[1], 'orange', label='mag$_{R}$ = ' + str(round(fitRw[0],2)) + '*m$_{r}$ + ' + str(round(fitRw[1],2)))
ax1.annotate("R$^{2}$ =", xy=(19.6, 28))
ax1.annotate(round(Rr_value,4), xy=(20.2, 28))
ax2.plot(fil_rg05, fitG05w[0]*fil_rg05 + fitG05w[1], 'r:', linewidth=1,          label='mag$_{G}$ = ' + str(round(fitG05w[0],2)) + '*m$_{r,g}$ + ' + str(round(fitG05w[1],2)) + ' (50% r_SDSS + 50% g_SDSS)')
ax2.annotate("R$^{2}$ =", xy=(19.6, 28.5))
ax2.annotate(round(Gr_value05,4), xy=(20.2, 28.5))
ax3.plot(fil_g, fitBw[0]*fil_g + fitBw[1], 'cyan', label='mag$_{B}$ = ' + str(round(fitBw[0],2)) + '*m$_{g}$ + ' + str(round(fitBw[1],2)))
ax3.annotate("R$^{2}$ =", xy=(19.2, 28.5))
ax3.annotate(round(Br_value,4), xy=(19.6, 28.5))

############

# axes
#ax1.set_title("R")
#ax2.set_title("G")
#ax3.set_title("B")
ax1.set_xlabel("r_SDSS (mag)")
ax1.set_ylabel("R (mag)")
ax2.set_xlabel("r_SDSS and g_SDSS (mag)")
ax2.set_ylabel("G (mag)")
ax3.set_xlabel("g_SDSS (mag)")
ax3.set_ylabel("B (mag)")

# legend
ax1.legend(loc=2, numpoints=1, fontsize=12)
ax2.legend(loc=2, numpoints=1, fontsize=10)
ax3.legend(loc=2, numpoints=1, fontsize=12)


# In[21]:

plt.close()


# The transformation equations between R,G,B and r,g SDSS filters are:
# 
# $$mag_{R_{phot}} = 0.99*m_{r_{SDSS}} + 9.74 \pm0.15$$
# 
# $$mag_{G_{phot}} = 0.99*m_{0.5(r+g)_{SDSS}} + 9.96 \pm0.19$$
# 
# $$mag_{B_{phot}} = 1.01*m_{g_{SDSS}} + 9.80 \pm0.33$$
# 
# The opposite being:
# 
# $$m_{r_{SDSS}} = 1.01*mag_{R_{phot}} - 9.83 \pm0.15$$
# 
# $$m_{0.5(r+g)_{SDSS}} = 1.01*mag_{G_{phot}} - 10.06 \pm0.19$$
# 
# $$m_{g_{SDSS}} = 0.99*mag_{B_{phot}} - 9.70 \pm0.33$$

# ###3.4.1 Testing the linear regression with the ODR (orthogonal distance regression)

# In[22]:

import scipy.odr

# Data !!!! These must be changed to produce the desired plot
x = fil_g
y = fil_B

# Fit the data using numpy.polyfit
fit_np = np.polyfit(x, y, 1, w = fil_Bew) # !!!! This must be changed to produce the desired plot

# Define model function for scipy.odr
def fit_func(pol, t):
  return pol[0]*t + pol[1]

#Fit the data using scipy.odr
Model = scipy.odr.Model(fit_func)
Data = scipy.odr.RealData(x, y)
Odr = scipy.odr.ODR(Data, Model, [1, 10], maxit = 10000)
output = Odr.run()

#output.pprint()
beta = output.beta
betastd = output.sd_beta

print "poly [slope, indep.term] =", fit_np
print "ODR [slope, indep.term] =", beta
print "% differences in slope and indep.term = ", round(((max(fit_np[0],beta[0])/min(fit_np[0],beta[0])-1)*100),2), "%,", round(((max(fit_np[1],beta[1])/min(fit_np[1],beta[1])-1)*100),2), "%"

plt.plot(x, y, "ro")
plt.plot(x, np.polyval(fit_np, x), "y--", lw = 4, zorder=1)
plt.plot(x, fit_func(beta, x), "g--", lw = 1, zorder=2)

plt.tight_layout()

plt.show()


# In[23]:

plt.close()


# Setting **no restriction** in error or r, g, R, G, B magnitudes:
# 
# - **R** parameters show a       **0.1%** difference between weighted linear regression and ODR methods.
# - **G** parameters show a **1.0 - 2.0%** difference between weighted linear regression and ODR methods.
# - **B** parameters show a **2.0 - 4.0%** difference between weighted linear regression and ODR methods.
# 
# Setting a **top of 0.5 mag for R,G,B mag. errors** and **no restriction** r, g, R, G, B magnitudes:
# 
# - **R** parameters show a       **0.1%** difference between weighted linear regression and ODR methods.
# - **G** parameters show a **2.6 - 4.6%** difference between weighted linear regression and ODR methods.
# - **B** parameters show a **0.8 - 1.8%** difference between weighted linear regression and ODR methods.
# 
# Setting a **top of 0.5 mag for R,G,B mag. errors**, a **top of 22.4 mag for r, g magnitudes** and **31 mag for R, G, B magnitudes**:
# 
# - **R** parameters show a **0.2 - 0.4%** difference between weighted linear regression and ODR methods.
# - **G** parameters show a **2.0 - 3.6%** difference between weighted linear regression and ODR methods.
# - **B** parameters show a **2.4 - 4.6%** difference between weighted linear regression and ODR methods.

# #4. Extinction correction
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
# $$A(B) = 4.1 * E_{B-V}$$
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
# $$m_{r_{SDSS}} = 1.01*mag_{R_{phot}} - 10.021 \pm0.15$$
# $$m_{g_{SDSS}} = 0.99*mag_{B_{phot}} - 9.963 \pm0.33$$

# #5. PSF matching: degrading R,G,B PSFs to match Rdeep PSF

# The images with better seeing (smaller, PSFs, R,G,B) are degraded to the one with the worst seeing (R deep): **Peris' R, G and B** images will be convolved with a function calculated by the psfmatch task in IRAF.
# 
# ***WARNING***: the images MUST have the same dimensions so the parameter 'psfdata' contains the same X,Y coordinates of the same star in the images to match.

                In IRAF: psfmatch task

PACKAGE = immatch
   TASK = psfmatch

input   =             @in.list  Input images
referenc= cig96_def_crop_wcs_al.fits  Reference images or reference psfs
psfdata =     611.12 416.53  Objects lists or input psfs 
kernel  =                       Input/output convolution kernels
(output =            @out.list) Output convolved images
(convolu=                image) Kernel computation method
(dnx    =                   12) X width of data region to extract
(dny    =                   12) Y width of data region to extract
(pnx    =                    9) X width of convolution kernel
(pny    =                    9) Y width of convolution kernel

Where:

in.list  contains: R, G, B images aligned:
CIG96_B_koords_skysub_al.fits
CIG96_G_koords_skysub_al.fits
CIG96_R_koords_skysub_al.fits

out.list contains: R, G, B images aligned AND corrected:
CIG96_B_koords_skysub_al_psf.fits
CIG96_G_koords_skysub_al_psf.fits
CIG96_R_koords_skysub_al_psf.fits
                
# #6. Creating calibrated R,G,B images
# 
# The calibration has been done with the original skysubtracted images.
# 
# The SB images will be created with the _**aligned and PSF corrected** images as input_:
# 
# - CIG96_R_koords_skysub_al.fits
# - CIG96_G_koords_skysub_al.fits
# - CIG96_B_koords_skysub_al.fits
# 
# ***WARNING***: there are some negative values in the backgorund of these images, np.abs() must be applied to the image data to avoid fake bright background results.

# In[125]:

# Select working directory
go_to('peris')

# R to SB
image = fits.open("CIG96_R_koords_skysub_al_psf.fits", mode='update')
# transform data to SB in magnitudes/arcsec2
im = image[0]
im.data = (25 -9.74 - 2.5*np.log10(np.abs(im.data)))/Rslope
# write image in SB
im.writeto("CIG96_R_SB.fits", clobber=True)
image.close

# G to SB
image = fits.open("CIG96_G_koords_skysub_al_psf.fits", mode='update')
# transform data to SB in magnitudes/arcsec2
im = image[0]
im.data = (25 -9.96 - 2.5*np.log10(np.abs(im.data)))/Gslope05
# write image in SB
im.writeto("CIG96_G_SB.fits", clobber=True)
image.close

# B to SB
image = fits.open("CIG96_B_koords_skysub_al_psf.fits", mode='update')
# transform data to SB in magnitudes/arcsec2
im = image[0]
im.data = (25 -9.80 - 2.5*np.log10(np.abs(im.data)))/Bslope
# write image in SB
im.writeto("CIG96_B_SB.fits", clobber=True)
image.close


# #7. Combine R,G,B and Rdeep images
# 
# #WARNING: this combined image is NOT VALID for science, they are not in the same band/filter!
# 
# R,G,B and R-deep now have a matching PSF, same orientation and size.
# 
# R-deep has already been converted to SB in <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96%20-%20Analysis.ipynb">Analysis notebook</a>, Section 3.1 so we use its result here: "*cig96_def_SB.fits*".

# In[112]:

images = ["CIG96_R_SB.fits", "CIG96_G_SB.fits", "CIG96_B_SB.fits", "cig96_def_SB.fits"]


                ### IRAF

imcombine task:

combine = median, scale = none, zero = none, weight = none, blank = 0.
                Images       N
        CIG96_R_SB.fits      1
        CIG96_G_SB.fits      1
        CIG96_B_SB.fits      1
      cig96_def_SB.fits     71

  Output image = cig96_RdRGB_SB.fits, ncombine = 74
                
# In[55]:

# median Rdeep & RGB combined image
image = "cig96_RdRGB_SB.fits"


                Only the R Cous (in r SDSS) and the R Peris (in r SDSS):

combine = median, scale = none, zero = none, weight = none, blank = 0.
                Images       N
        CIG96_R_SB.fits      1
      cig96_def_SB.fits     71

  Output image = cig96_RdR_rSDSS.fits, ncombine = 72
                
# In[ ]:

# median (or mean) Rdeep & R combined image
image = "cig96_RdR_rSDSS_SB.fits"


# ##7.1 Cropping median Rdeep & RGB combined image

# In[56]:

# WARNING: cropping DOES NOT respect the WCS, it must be checked!

# Select working directory
go_to('peris')

# read the image
image = fits.open("test.fits", mode='update')
#image.info()

# select layer (in this case, it is unique)
im = image[0]
im.shape
dimx = im.shape[1]
dimy = im.shape[0]

# crop image
x = 130
im.data = im.data[x:dimx-x, x:dimy-x]

# update image
im.writeto('test_cropped.fits', clobber=True)


# ##7.2 Trim method, repipy

# In[151]:

get_ipython().run_cell_magic(u'bash', u'', u'trim_images.py --suffix " _cropped" --region 100 980 80 900 cig96_RdRGB_SB_oriented.fits\ntrim_images.py --suffix " _cropped_small" --region 111 960 51 900 cig96_RdRGB_SB_oriented.fits')


# Visualization of trimmed image.

# In[37]:

# Select working directory
go_to('peris')

# data
image = "cig96_RdRGB_SB_oriented_cropped.fits"
data = fits.open(image)[0].data

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.08

# plot
cig96 = aplpy.FITSFigure(image, north=True)
cig96.axis_labels.set_font(size=18)
cig96.set_tick_labels_font(size=12)
stretch='linear'
vmin=22
vmax=28
cig96.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone')
cig96.add_colorbar()# show_colorbar()
cig96.colorbar.set_font(size=16)
#cig96.add_colorbar()
#cig96.colorbar("$\mu_{r\ SDSS}$")
cig96.show_contour(levels=[19, 20, 21, 22, 23, 24, 25, 26, 27.0, 27.5], smooth=1, cmap="spectral_r", returnlevels=True, inline=True)

# zoom
cig96.recenter(RA, DEC, rad)

# save figure
plt.savefig("/home/prm/Desktop/cig96_images/cig96_"+str(stretch)+str(vmin)+"_"+str(vmax)+".png", bbox_inches='tight')


# In[20]:

cig96.close()


# In[21]:

# Select working directory
go_to('peris')

# data
image = "cig96_RdRGB_SB_oriented_cropped.fits"

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.07

# plot
cig96 = aplpy.FITSFigure(image, north=True)
plt.tight_layout()
cig96.axis_labels.set_font(size=12)
cig96.show_colorscale(stretch='linear', vmin=20.5, vmax=27.7, cmap="bone")
cig96.show_colorbar()
cig96.show_contour(levels=[25, 26.5, 27.0, 27.5], cmap="RdYlGn")

# zoom
cig96.recenter(RA, DEC, rad)

# save figure


# In[22]:

cig96.close()


# In[14]:

# Select working directory
go_to('cig')

# data
#image = "cig96_RdRGB_SB_oriented_cropped.fits"
#image = "cig96_def_crop_wcs.fits"
image = "cig96_def_SB.fits"

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.07

# plot
cig96 = aplpy.FITSFigure(image, north=True)
cig96.axis_labels.set_font(size=18)
cig96.set_tick_labels_font(size=12)
stretch='log'
vmin=25
vmax=28
cig96.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone')

# colorbar
cig96.add_colorbar()# show_colorbar()
cig96.colorbar.set_font(size=16)
cig96.colorbar.set_axis_label_text("mag arcsec$^{-2}$")
cig96.colorbar.set_axis_label_font(size=14)
#cig96.add_colorbar()
#cig96.colorbar("$\mu_{r\ SDSS}$")

# contours
lev = [19, 20, 21, 22, 23, 24, 25.0, 25.5, 26.0, 26.5, 27.0, 27.5]
cig96.show_contour(levels=lev, smooth=3, cmap="spectral_r", linewidths=2, returnlevels=True)

# marks
plt.plot(300, 300, 'yo', markersize=13, markeredgecolor='k', zorder=4)

# zoom
cig96.recenter(RA, DEC, rad)

# save figure
plt.savefig("/home/prm/Desktop/cig96_images/cig96_2_2inv_"+str(stretch)+str(vmin)+"_"+str(vmax)+".png", bbox_inches='tight')


# In[88]:

cig96.close()

