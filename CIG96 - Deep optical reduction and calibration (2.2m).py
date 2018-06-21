
# coding: utf-8

# # CIG96 - Deep optical reduction and calibration
# 
# ### Observation date: 11sep12
# ### Telescope: CAHA2.2m
# 
# #I. GENERAL INFORMATION
# 
# Notebook with the full data reduction and calibration plus external IRAF steps.
# 
# - Further info from this CAHA campaign and log can be found <a href="http://amiga.iaa.es:8888/display/science/CAHA+2.2m+Observations+11sep12">here</a>.
# 
# - Notes on the data reduction and calibration can be found <a href="http://amiga.iaa.es:8888/display/science/Data+reduction+CIG96+CAHA2.2">here</a>.
# 
# - Information about CIG96 paper (2016) can be found <a href="http://amiga.iaa.es:8888/display/science/CIG96+-+Deep+optical+and+HI+images#CIG96-DeepopticalandHIimages-Paperstructureandcomments">here</a>.
# 
# #II. PRE-REDUCTION SETUP
# 
# **NOTE**: when a box is previously labeled as ***IRAF***, it means the notebook needs an external computation in IRAF before proceeding with further cells.
# 
# ##a) Import necessary modules

# In[1]:

import matplotlib as mpl
import repipy.astroim as astroim
import repipy.rename as rename
import repipy.utilities as utilities
import repipy.combine as combine
import repipy.create_masks as create_masks
import repipy.remove_cosmics as remove_cosmics
import repipy.arith as arith
import repipy.find_sky as find_sky
import repipy.scale_to_ref as scale_to_ref
import astropy.io.fits as fits
import lemon.astrometry as astrometry
import lemon.photometry as photometry
import numpy as np
import os
import sys
import shutil
import random
import subprocess
import warnings
import pandas
import sqlite3
import repipy.extract_mag_airmass_common as extract
import glob as glob
import matplotlib.pyplot as plt
from glob import glob as ls
from astropy.io import fits
import pyregion
import scipy
from scipy import stats
from stsci import convolve
#import seaborn as sns


# In[2]:

# renders interactive figures in the Notebook
get_ipython().magic(u'matplotlib nbagg')

# renders figures as static PNG
#%matplotlib inline


# In[3]:

mpl.pylab.rcParams['figure.figsize'] = 14, 7  # that's default image size for this interactive session


# ##b) Raw data and work directories definition

# In[4]:

# Directory where data are:
raw_data_dir = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16"
# Working directory:
work_dir = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16"
print "Raw data path: ", raw_data_dir
print "Work path:     ", work_dir

# Switch to working directory
os.chdir(os.path.expandvars(work_dir))
print "\nYou will be working on the following path: ", work_dir


# ## IMPORTANT NOTE: if previously reduced:
# 
# - AVOID repeating steps 1 to 9 (included) from Section III.
# 
# - **MUST REPEAT ONLY** steps **12 and 16** from Section III.

# #III. DATA REDUCTION and CALIBRATION
# 
# ##1. Headers fix
# Unfortunately, many observatories include WCS keywords in the headers, **but not the proper WCS info**. So, not only you have useless, and sometimes old-fashion keywords in the header, but they confuse programs trying to read the headers AFTER you have used astrometry.net. Best solution is to remove them first. Also, some **keywords** might be doubled or show ambiguous names (e.g. skyflats show up as "flat" instead of "skyflats").
# 
# We make corrections for all these.
# 
# ###1.1 WCS removal

# In[5]:

get_ipython().run_cell_magic(u'bash', u'', u'pwd')


# In[6]:

get_ipython().run_cell_magic(u'bash', u'', u'# removing WCS from all fits files\nwipe_wcs.py cig/cig96_def_crop-noWCS.fits')


# ###1.2 Proper keywords reading

# In[56]:

# makes a list containing all fits images
image_list = ls(os.path.join(work_dir, "*.fits"))
im = astroim.Astroim(image_list[0])


# Now we read from the object astroim the keywords we need. 

# In[55]:

# assign the specified keywords to proper variables
filterk = im.header.filterk
exptimek = im.header.exptimek
objectk = im.header.objectk
airmassk = im.header.airmassk
datek = im.header.datek
telescope = im.header.telescope
gaink = im.header.gaink
read_noisek = im.header.ccdronk
timek = im.header.timek
filtersysk = im.header.filtersysk
print exptimek
print filterk


# ###1.3 'OBJECT' keyword renaming in skyflats
# All the flats are actually sky flats this night, so we will rebrand them to make sure we have all the information. 

# In[8]:

for im_name in image_list:
    im = fits.open(im_name, 'update')
    object_im = im[0].header["OBJECT"]
    if 'flat' in object_im.lower() and 'sky' not in object_im.lower():
        im[0].header["OBJECT"] = "skyflat"  
    im.flush()
    im.close()


# ##2. WCS restoration
# Proper astrometry is added again to all the images using astrometry.net. This will take some minutes, so take it easy ;). 

# In[ ]:

get_ipython().run_cell_magic(u'bash', u'', u'astrom.py --overwrite --radius 1 *.fits')


# ##3. Filter specification
# 'CousR' filter is not recognized by passband.py subroutine so the **filter keyword is changed to 'Cousins R'**, something that can be read.
# Also, a **new keyword named "FILTER" is added** since "INSFLNAM" is not usually recognized as the filter keyword by *any* normal routine.

# In[10]:

for fitsim in image_list:
    data, header = fits.getdata(fitsim, header=True)
    header["INSFLNAM"] = "Cousins R"  # change filter keyword to 'Cousins R' so rename.py can read it correctly
    header.set("FILTER","Cousins R")  # add a redundant filter keyword ("FILTER") just in case INSFLANM is not identified
    fits.writeto(fitsim, data, header, clobber=True)


# ##4. Renaming all the files
# Now we rename the files with a recognizable filter.

# In[11]:

get_ipython().run_cell_magic(u'bash', u'', u'rename.py --copy --overwrite .')


# ###4.1 Manual check and deletion of useless files
# Among the dataset there might be bad images or images from another targets observed during the night. All are manually removed:

# In[12]:

get_ipython().run_cell_magic(u'bash', u'', u"# CIG1019 images in 'cig' folder\nrm -Rf cig/*cig1019*")


# ## 5. Mask  images
# Mask out pixels with values above a maximum number of counts

# In[13]:

get_ipython().run_cell_magic(u'bash', u'', u'pwd\ncreate_masks.py --circular --max_val 55000 skyflat/*.fits \ncreate_masks.py --circular --max_val 55000 cig/*.fits \n#create_masks.py --circular --max_val 55000 standards/*.fits\n#create_masks.py --circular --max_val 10000 blanks/*.fits # for these, the max counts number can be lowered drastically')


# ## 6. Bias
# We will median combine all the bias images to get a master bias. Since it is a bias, no scaling is necessary, and we have chosen the median for the combination because even at zero exposure time you get some cosmic rays in each image. The "--all_together" keyword will tell the program combine.py to ignore the filters of the images, and combine all the images together no matter what filters they have. 
# 
# ###6.1 Bias combination

# In[14]:

get_ipython().run_cell_magic(u'bash', u'', u'cd bias\ncombine.py --output "masterbias.fits" --average "median" --scale "none" --all_together bias*.fits')


# ###6.2 Bias subtraction
# Given that the bias has some small amount of structure, we will remove the whole image instead of an average. We will remove it from all other images, which include cig images, flats, standards, bias and even unknown objects, in case they are needed later. We will add a suffix "-b" before the fits extension. The message "masterbias image removed" will be included in the header and the subtraction will not make use of the keywords "--median" or "--mean", which would subtract only the median or mean of the second image. 

# In[15]:

get_ipython().run_cell_magic(u'bash', u'', u'arith.py --suffix " -b" --message "Masterbias image removed" cig/cig*.fits - bias/masterbias.fits\narith.py --suffix " -b" --message "Masterbias image removed" skyflat/skyflat_*.fits - bias/masterbias.fits\narith.py --suffix " -b" --message "Masterbias image removed" bias/bias*.fits - bias/masterbias.fits\n#arith.py --suffix " -b" --message "Masterbias image removed" standards/*.fits - bias/masterbias.fits\n#arith.py --suffix " -b" --message "Masterbias image removed" blanks/blank*.fits - bias/masterbias.fits')


# ##7. Overscan subtraction
# The difference between the median level of the overscan and the rest of the image is usually a constant. Despite the mean/median amount of counts of a bias image may vary throughout the night, such increase or decrease is usually the same for both the bias itself and the overscan region. This way, after subtracting the master-bias, we can subtract the overscan value to avoid any over/infrasubtraction of the real bias values due to those variations throughout the night.
# 
# For CAFOS at CAHA2.2m, the overscan region is located on the left side of the image. Its physical extension (in X, Y pixels) is a little larger than the values specified below: [1028:1080,6:1020], but these are good enough for a full subtraction.

# In[16]:

get_ipython().run_cell_magic(u'bash', u'', u'subtract_overscan.py --region  1028 1080 6 1020  --suffix " -o" cig/*-b.fits\nsubtract_overscan.py --region  1028 1080 6 1020  --suffix " -o" skyflat/*-b.fits\nsubtract_overscan.py --region  1028 1080 6 1020  --suffix " -o" bias/*-b.fits\n#subtract_overscan.py --region  1028 1080 6 1020  --suffix " -o" standards/*-b.fits\n#subtract_overscan.py --region  1028 1080 6 1020  --suffix " -o" blanks/*-b.fits')


# ## 8. Flatfields
# ### 8.1 Flatfields combination
# All flats are combined. We do not need to separate them by filters, the routine combine will take care of that. We just need to give all the bias-subtracted flats (xxx-b.fits) as input and it will separate them by filter and rename them accordingly. 

# In[17]:

get_ipython().run_cell_magic(u'bash', u'', u'pwd\ncd skyflat/\ncombine.py --output "masterflats.fits" --average "median" --scale "median" skyflat*-o.fits ')


# ### 8.2 Flatfield correction
# One by one, we will flat-field correct all the cig images, standards and skyflats. This would be easiest to do using a program, but we prefer to explicitly show the operations. A lot of "RuntimeWarning: divide by zero [...]" happens in this operation, because it's true, there are a lot of zeroes in the outter parts of the flats images. This will cause no trouble at all, though. They're warnings, not errors.

# In[ ]:

get_ipython().run_cell_magic(u'bash', u'', u'arith.py --suffix " -f" --message "Divided by masterflat image" cig/*-o.fits / skyflat/masterflats_CousinsR.fits\narith.py --suffix " -f" --message "Divided by masterflat image" skyflat/*-o.fits / skyflat/masterflats_CousinsR.fits\n#arith.py --suffix " -f" --message "Divided by masterflat image" blanks/*-o.fits / skyflat/masterflats_SDSSr.fits\n#arith.py --suffix " -f" --message "Divided by masterflat image" standards/*-o.fits / skyflat/masterflats_SDSSr.fits')


# ## 9. Remove cosmics
# We wil use cosmics.py, routine that performs for python what Van Dokkum's L.A.cosmic algorithm (http://arxiv.org/abs/astro-ph/0108003) did in IRAF. This is, by far, the most time-comsuming part of the data reduction. Doing 3 iterations per object, it can take easily ~40 minutes, so this might be a good time to go for a coffee. 
# 
# NOTE: if available, it is extremely unlikely that a cosmic ray falls upon a standard star in any of the 3 or 4 seconds of integration time. Should it happen, we would definitely see this in the flux calibration so, then, we would correct it.

# In[21]:

get_ipython().run_cell_magic(u'bash', u'', u'remove_cosmics.py --suffix " -c" --sigclip 5 --maxiter 2 cig/*-f.fits')


# ### NOTE: no need to repeat anything before this point (from 1 to 9 included).

# ##10. Alignment and PSF/FWHM estimation
# First, we run IRAF to algin the images. Then, we compute the FWHM and seeing measurements for each image and the night mean. The scripts and outputs can be found <a href="http://amiga.iaa.es:8888/display/science/Data+reduction+CIG96+CAHA2.2">here</a>.
# 
# More info on the IRAF commands used in this section:
# - <a href="http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?imexamine">imexam</a>
# - <a href="http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?rimcursor">rimcursor</a>
# - <a href="http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?psfmeasure.hlp">psfmeasure</a>
# - <a href="http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?imalign.hlp">imalign</a>
# 
# ###10.1 Images shifts
# In IRAF, selection of the **physical coordinates of one reference star** over the reference image "*cig96_20120911_CousinsR_001-b-o-f-c.fits*":
# 
# ##### *IRAF*

                iraf_prompt> !ds9 &
iraf_prompt> display cig96_20120911_CousinsR_001-b-o-f-c.fits 1
iraf_prompt> imexam 

# Ref. star coordinates (physical) on image 1:
# 263.947    615.66
                
# Manual selection of the reference star in the rest of images (plus coordinates saving):
# 
# ##### *IRAF*

                iraf_prompt> cl < refstar-coords.cl

# This script consists on a display+rimcursor loop over the 71 images with these commands:
# | iraf_prompt> display image(i) 1
# | iraf_prompt> rimcursor >> refstar-coords.dat      # SPACE BAR to mark the position of the star; CTRL+D to skip to next image.
# The output file may have to be edited: only the 2 first columns are useful and indicate x and y position of the reference star.
                
# Shift computation from all the images according to the reference star (in this case, x=263.947,y=615.66) in the reference image:

# In[245]:

# reference star physical x,y coordinates
refx = 263.947
refy = 615.660

# load data
coords = np.genfromtxt("cig/refstar-coords.dat", delimiter=' ', dtype=float)

# compute shifts
shiftx = refx - coords[:,0]
shifty = refy - coords[:,1]

# create file with shifts
shifts = np.column_stack((shiftx,shifty))

#print shifts   # type: array; shape: 70,2
np.savetxt("cig/refstar-shifts.dat", shifts)


# ### 10.2 Images alignment
# We select a number of stars that **must not be saturated** (refkey=a in *imexam* task provides PEAK value) and **must not have a close stellar companion** (pressing 'r' over the star provides a radial profile).
# 
# They **must have a successful IRAF fitting** (*imexam* step).
# 
# 
# ##### *IRAF*

                iraf_prompt> display cig96_20120911_CousinsR_001-b-o-f-c.fits 1          # open the image
iraf_prompt> imexam                # 'r' key returns fittings (must NOT be saturated, DO have successful fitting and NO companions)
iraf_prompt> rimcursor > alignstars-coords.dat            # record the positions of in the .dat file 
                
# With the shifts file (*refstar-shifts.dat*) and the stars coordinates just selected for the photometry (*alignstars-coords.dat*) we 
# can align all the images.
# 
# Task "imalign" parameters used are:
# 
# ##### *IRAF*

                iraf_prompt> epar imalign

input   =         @images.list                      Input images
referenc= cig96_20120911_CousinsR_001-b-o-f-c.fits  Reference image
coords  = alignstars-coords.dat                     Reference coordinates file
output  =      @imagesout.list                      Output images
(shifts =   refstar-shifts.dat)                     Initial shifts file
(boxsize=                    7)                     Size of the small centering box
(bigbox =                   11)                     Size of the big centering box
...

iraf_prompt> imalign
                
# This task will return:
# 
#    > Trim_Section = [23:1080,15:1004]
#    
#    > Shifting images:
#    
#    > Trimming images:  corrected section = [23:1080,15:1004]
# 
# And the output images (according to imagesout.list file) will have an "**-a**" (after 'aligned') added to the end of the filename.

# ### 10.3 M-FWHM of the selected stars
# 
# Now we compute the Moffat-FWHM of each image.
# 
# To select the **appropriate stars for the FWHM**, *psfmeasure* and or *imexam (refkey=r)* should be tested previously in each candidate star individually. 
# 
# **Example**: a candidate star should have been checked with imexam first, to see if it is not saturated, companion-less and present in all images. If so, printing the output and running the following lines may help see the peak counts in each image.

                # Example to check an imexam output over a certain star:

star = open("cig/star_imexam_output.dat", "r")
peaks = open("cig/star4_just_peaks.dat", "wa")

for index,line in enumerate(star):
    if index in range(3,284,4):
        peaks.write(line[32:38] + "\n")
peaks.close()

peaks = np.genfromtxt("cig/star_just_peaks.dat")
s_peaks = np.sort(peaks)
for i in s_peaks:
    if i > 60000:
        print i
                
# #### ¬ M-FWHM of the night:
# 
# We compute an **all-night average M-FWHM** ("M-" stands for Moffat fit) from all aligned images ("*@images_a.list*") based on the positions of the selected stars "*photstars-coords.list*" (same as ds9 region file: "*photstars-psfmeasure-8stars.reg*") by running a short script "*mfwhm-allnight.cl*" that will set up the following parameters in *psfmeasure* command:
# 
# ##### *IRAF*

                iraf_prompt> epar psfmeasure

images  =       @images_a.list    List of images         # list of all aligned images "-a.fits".
(coords =                mark1)   Object coordinates
(wcs    =              logical)   Coordinate system
(display=                  yes)   Display images?
(frame  =                    1)   Display frame to use
(level  =                  0.5)   Measurement level (fraction or percent)
(size   =                 MFWHM)  Size to display        # MFWHM = Moffat fit; GFWHM = Gaussian fit;
...
(imagecu= photstars-coords.list)  Image cursor input

iraf_prompt> cl < mfwhm-allnight.cl                      # we measure the PSF and MFWHM over the selected stars
                
# The output should be something like this:

                [...]
cig96_20120912_  460.55  787.49    2.93   2.225  2.3    0.53      73
                 373.27  749.23    2.55   2.594  4.1    0.81       5
                 291.08  739.81    0.54   1.697  4.7    0.11      80
                 241.82  602.04   -0.07   1.641  3.6    0.03     -78
                 280.14  464.85    2.93   1.967  3.8    0.99      76
                 417.12  205.36    2.53   2.190  2.7    0.61       8
                 461.52  249.73    2.75   2.061  4.6    0.00       0
                 535.95  280.05    3.34   2.136  3.9    0.97     -27
                 635.77  262.64    3.71   2.112  2.7    0.58      16

  Average star @ (460.45, 787.45): MFWHM=2.37, m=3.17
  Average star @ (373.14, 749.42): MFWHM=2.70, m=2.60
  Average star @ (290.98, 739.81): MFWHM=1.80, m=0.53
  Average star @ (241.86, 601.96): MFWHM=1.82, m=0.00
  Average star @ (280.12, 464.74): MFWHM=2.14, m=2.94
  Average star @ (417.15, 205.41): MFWHM=2.37, m=2.53
  Average star @ (461.56, 249.72): MFWHM=2.19, m=2.74
  Average star @ (536.00, 280.13): MFWHM=2.32, m=3.31
  Average star @ (635.88, 262.80): MFWHM=2.28, m=3.79

  Average full width at half maximum (MFWHM) of 1.9348
                
# #### · Illustrative results:
# 
# The mean-averaged all-night M-FWHM shows up in the last line:     **Average full width at half maximum (MFWHM) of 1.9348** but this is not a value we trust since it is the mean of the means.
# 
# The median, min and max values of the mean-averaged all-night M-FWHM are:

# In[33]:

allmfwhm = np.genfromtxt("cig/mfwhm-per-image-values.dat", dtype=float)

pixscale = 1.0461051

print "Minimum M-FWHM of the night = ", min(allmfwhm), "pix"
print "Maximum M-FWHM of the night = ", max(allmfwhm), "pix"
print "Median all-night M-FWHM = ", round(np.median(allmfwhm),4), "+/-", round(np.std(allmfwhm),4), "pix \n"

print "Minimum M-FWHM of the night = ", min(allmfwhm), "pix = ", min(allmfwhm)/pixscale, "arcsec"
print "Maximum M-FWHM of the night = ", max(allmfwhm), "pix = ", max(allmfwhm)/pixscale, "arcsec"
print "Median all-night M-FWHM = ", round(np.median(allmfwhm)/pixscale,4), "+/-", round(np.std(allmfwhm)/pixscale,4), "arcsec"


# #### ¬ M-FWHM of each image:
# 
# We compute now an **average M-FWHM per image** ("M-" stands for Moffat fit again). This means that a **M-FWHM is computed for each image**.
# 
# The coordinates of the stars used to measure the M-FWHM are in the file *"photstars-coords.list"*. They are extracted from any aligned image via "*imexam*" task and are the same as in ds9 region file "*photstars-psfmeasure-8stars.reg*":
# 
# ##### *IRAF*

                iraf_prompt> display cig96_20120911_CousinsR_001-b-o-f-c.fits 1 
iraf_prompt> imexam > photstars-coords.list                           # click 'x' on each star and end with 'q'
                
# In IRAF, we run a file named "*mfwhm-per-image.cl*" that measures the M-FWHM for each image with the "*psfmeasure*" task, with the following structure:
# 
# ##### *IRAF*

                # psfmeasure images=cig96_20120911_CousinsR_001-b-o-f-c-a.fits coords=markall size=MFWHM imagecu=photstars-coords.list > mfwhm-per-image.dat
[...]
# psfmeasure images=cig96_20120911_CousinsR_070-b-o-f-c-a.fits coords=markall size=MFWHM imagecu=photstars-coords.list > mfwhm-per-image.dat

iraf_prompt> cl < mfwhm-per-image.cl
                
# The **average M-FWHM per image** is computed for each image and can be found at the end of each loop in a specific sentence as: "*Average full width at half maximum (MFWHM) of 1.3750*" (or any other value).

# ### 10.4 M-FWHM distribution
# 
# We plot the evolution throughout the night and a histogram to see the trend of the seeing:

# In[5]:

imdates = ['00:00:00', '00:01:19', '00:05:21', '00:09:24', '00:13:27', '00:17:30', '00:21:39', '00:25:42', '00:29:45', '00:33:48',          '00:37:51', '00:41:59', '00:46:01', '00:50:04', '00:54:07', '00:58:10', '01:02:17', '01:06:19', '01:10:22', '01:14:25',          '01:18:27', '01:22:36', '01:26:38', '01:30:41', '01:34:44', '01:38:47', '01:43:03', '01:47:05', '01:51:08', '01:55:11',          '01:59:14', '02:03:21', '02:07:23', '02:11:26', '02:15:29', '02:19:31', '02:23:40', '02:27:43', '02:31:46', '02:35:49',          '02:39:52', '02:44:10', '02:48:13', '02:52:16', '02:56:19', '03:00:21', '03:04:29', '03:08:32', '03:12:35', '03:16:38',          '03:20:41', '03:24:49', '03:28:51', '03:32:54', '03:36:57', '03:40:59', '03:45:06', '03:49:09', '03:53:12', '03:57:14',          '04:01:17', '04:05:26', '04:09:28', '04:13:31', '04:17:34', '04:21:37', '04:25:48', '04:29:51', '04:33:53', '04:37:57',          '04:41:59']
dates = range(71)

mfwhm = np.genfromtxt("cig/mfwhm-per-image-values.dat", dtype=float)

plt.xticks(dates, imdates, rotation=70, size=12)
plt.plot(dates, mfwhm, 'gd-')
plt.title('Moffat-FWHM evolution')
plt.legend(['M-FWHM'], loc='upper right', fontsize='x-large')
plt.show()

plt.hist(mfwhm, range=[1.3,1.9], bins=12, color="g")
plt.title('Moffat FWHM distribution')
plt.legend(['M-FWHM'], loc='upper right', fontsize='x-large')
plt.show()


# In[6]:

plt.close()


# ## 11. Photometry
# ###11.1 From IRAF output (M-FWHM file) to python-readable table
# 
# To work with the IRAF output as a python-readable table, we need to trim some lines from that file so a table can be created and read easily.

# In[21]:

# number of stars used to measure FWHM and to be used for photometry, and total number of images
numstars = int(sum(1 for line in open("cig/photstars-coords.list")))
numimages = int(sum (1 for line in open("cig/images_a.list")))

# python-readable table file creation
fileout = open("cig/mfwhm-allnight-pyready.dat", "w")
fileout.write('Image\tColumn\tLine\tMag\tMFWHM\tBeta\tEllip\tPA\n')
fileout.close()

# number of heading lines to be cut (unuseful information, out of the table-to-be)
heading = 5

# input/output files definitions
filein = np.loadtxt("cig/mfwhm-allnight.dat", skiprows=heading+1, dtype=np.str, delimiter='\t')
numlines = sum(1 for line in filein)
fileout = open("cig/mfwhm-allnight-pyready.dat", "a")

# replace blanks in table
newfile = [filein[i].replace("                 ", "cig   ") for i in range(0, len(filein) - numstars - 1)]

# write the trimmed file 'newfile' in 'fileout'
for i in range(0,len(newfile)):
    fileout.write(newfile[i]+'\n')

fileout.close()


# ### 11.2 Median M-FWHM per image
# Now we may work with the file with the following columns: Image, Column, Line, Mag, MFWHM, Beta, Ellip and PA.

# In[22]:

# read the MFWHM file as a table
mfwhm_tab = pandas.read_csv("cig/mfwhm-allnight-pyready.dat", delim_whitespace=True, dtype=np.str)

# read the proper MFWHM values
mfwhm_values = np.array(mfwhm_tab['MFWHM'].values)

# put the MFWHM values as floats in a list
mfwhm_list = [np.float(mfwhm_values[i]) for i in range(0, len(mfwhm_values))]

# make the median with a regular expression (do X for 'variable' in 'range')
#### median_MFWHM = [np.median(FWHM[numstars*n:numstars*(n+1):1]) for n in range(0,71)]

# or, equivalent, do a normal "for" loop:
filemedian = open("cig/mfwhm-per-image-medians.dat", "wa")
medians = np.zeros(numimages)

# calculate and write the medians for each image in a new file ('filemedian' variable)
for n in range(0,numimages):
    per_image_MFWHM = mfwhm_list[numstars*n:numstars*(n+1):1]     # takes the X measurements of each image (X = numstars)
    medians = np.median(per_image_MFWHM)                          # calculates median value of the X measurements of MFWHM 
    mstdev = np.std(per_image_MFWHM)                              # calculates the st.dev. of the X measurements
    filemedian.write(str(medians) + "\t" + str(mstdev) + "\n")    # writes medians and st.dev. in file 

filemedian.close()


# ###11.3 Magnitudes calculation (apertures definition + qphot + pdump)
# Definition of apertures and annulus (for each image, according to each median M-FWHM) and creation of IRAF scripts to perform *qphot* task and magnitude and errors extraction from the qphot-output databases via *pdump* command.

# In[52]:

# lists of medians and images 
medians = np.genfromtxt("cig/mfwhm-per-image-medians.dat", dtype=float)
images_a = np.genfromtxt("cig/images_a.list", dtype=str)
numimages = int(sum (1 for line in open("cig/images_a.list")))

# aperture and annulus radii definition (in pix)
rad_ap = medians*2
rad_in = medians*4
rad_out = medians*6
width = medians*2

# open IRAF script files
qphotscript = open("cig/qphot_per_image.cl", "w")
pdumpscript = open("cig/pdump_mag.cl","wa")

# qphot and pdump commands written on qphot- and pdump- files
for i in range(0,numimages):
    qphotscript.write("qphot image=" + str(images_a[i])                       + " coords=photstars-coords.list cbox=5. zmag=0.0 exposure=EXPTIME airmass=AIRMASS filter=FILTER obstime=DATE interactive=no aperture="                       + str(rad_ap[i][0]) + " annulus=" + str(rad_in[i][0]) + " dannulus=" + str(width[i][0]) + "\n")
    pdumpscript.write("pdump infiles=" + str(images_a[i]) + ".mag.1 fields=MAG,MERR expr=yes > mag_merr_" + str(images_a[i]) + ".dat" + "\n")

# last line in pdump_mag.cl file: creates a file with all mag_merr databases files listed
pdumpscript.write("ls mag_merr_cig96*.dat > mag_tables.list")

qphotscript.close()
pdumpscript.close()


# #####IRAF
# 
# + Running ***qphot_per_image.cl*** script in IRAF results in the creation of a database per image ("*...mag.1*" files), containing the *qphot* output for each of the 8 stars of *photstars-coords.dat*.
# 
# + By running ***pdump_mag.cl*** script, we extract the necessary variables: magnitude (**MAG**) and its error (**MERR**) and store them in individual files (table shaped).
# 
# ### WARNING: IRAF introduces a 25 mag shift via *zmag* parameter if left as default. Care to reset.
# 
# ###11.4 Magnitudes per star
# Now we create files for each star, containing the star's magnitude in each image.

# In[59]:

# reading of number of stars (redundant, previously defined) and magnitude files
numstars = int(sum(1 for line in open("cig/photstars-coords.list")))  # number of stars selected for the photometry
mag_files = np.genfromtxt("cig/mag_tables.list", dtype=str)           # list with magnitude files "mag_merr_cig96...dat"

# mag and merr are separated by star and added to separate files (one per star)
for star in range(0,numstars):
    mag_star = open("cig/mag_star_00" + str(star+1) + ".dat", "wa")   # new file for each star
    for table in mag_files:
        tab = open("cig/"+table, "r")                                 # open each of the 71 tables
        value = float(tab.readlines()[star][0:6])                     # read 
        #val_err = float(tab.readlines()[star][8:11])
        mag_star.write(str(value) + "\n") #+ "\t" + str(val_err) + "\n")
        tab.close()
    mag_star.close()


# ## 12. Extinction coefficient and ZP of the night
# 
# ### 12.1 Extinction coefficient
# 
# First, we take from SDSS the ***r*** and ***i*** magnitudes of each star and convert them to R using Lupton (2005) conversion formula:
# 
# $$R_{ri} = r - 0.2936*(r-i) - 0.1439$$
# 
# Alternatively, it can be also done by selecting ***r*** and ***g*** magnitudes and converting to R with this formula (also from Lupton 2005):
# 
# $$R_{gr} = r - 0.1837*(g - r) - 0.0971$$

# In[5]:

# SDSS calibrated r,i,g magnitudes for the 8 stars
i_sdss = [19.30, 18.86, 16.94, 19.06, 18.68, 18.79, 19.48, 20.14]
r_sdss = [19.78, 19.20, 17.14, 19.66, 19.27, 19.35, 20.08, 20.66]
g_sdss = [21.03, 20.00, 17.63, 21.26, 20.77, 20.88, 21.84, 22.14]
numstars = int(sum(1 for line in open("cig/photstars-coords.list")))  # number of stars selected for the photometry

# r,i conversion -------------------------------------------------------

# list of R magnitudes from the conversion using r_sdss, i_sdss
Rsdss_ri = []

# conversion to R from r,i
for i in range(0,numstars):
    R = r_sdss[i] - 0.2936*(r_sdss[i] - i_sdss[i]) - 0.1439
    Rsdss_ri.append(round(R,4))
    
Rsdss = Rsdss_ri
print "R (from r,i_SDSS) = ", Rsdss


## r,g conversion --------------------------------------------------------

## list of R magnitudes from the conversion using r_sdss, g_sdss
#Rsdss_rg = []

# conversion to R from r,g
#for n in range(0,numstars):
#    R = r_sdss[n] - 0.1837*(g_sdss[n] - r_sdss[n]) - 0.0971
#    Rsdss_rg.append(round(R,4))
    
#print "R (from rg_SDSS) = ", Rsdss_rg, "\n"
#Rsdss = Rsdss_rg


# Once we have the magnitudes of each star in each image in independent files, we may extract the **airmasses**, plot them against such **magnitudes** and make a Bouguer fit to find the extinction coefficient and Zero Point.

# In[6]:

# list of images
images_a = np.genfromtxt("cig/images_a.list", dtype=str)

# extraction of airmasses
airmass = []
for im in images_a:
    image = fits.open("cig/backup_aligned_images/" + im)
    airmass_item = image[0].header['AIRMASS']
    airmass.append(airmass_item)
#print airmass
#print len(airmass)

# magnitude files and Bouguer fit parameters (extinction coefficient and ZP)
kext = []
zp = []

# 'airmass' list must be converted to an array to plot it with magnitudes
airmass = np.array(airmass)

# number of stars selected for the photometry
numstars = int(sum(1 for line in open("cig/photstars-coords.list")))

# magnitude per date
imdates = ['00:00:00', '00:01:19', '00:05:21', '00:09:24', '00:13:27', '00:17:30', '00:21:39', '00:25:42', '00:29:45', '00:33:48',          '00:37:51', '00:41:59', '00:46:01', '00:50:04', '00:54:07', '00:58:10', '01:02:17', '01:06:19', '01:10:22', '01:14:25',          '01:18:27', '01:22:36', '01:26:38', '01:30:41', '01:34:44', '01:38:47', '01:43:03', '01:47:05', '01:51:08', '01:55:11',          '01:59:14', '02:03:21', '02:07:23', '02:11:26', '02:15:29', '02:19:31', '02:23:40', '02:27:43', '02:31:46', '02:35:49',          '02:39:52', '02:44:10', '02:48:13', '02:52:16', '02:56:19', '03:00:21', '03:04:29', '03:08:32', '03:12:35', '03:16:38',          '03:20:41', '03:24:49', '03:28:51', '03:32:54', '03:36:57', '03:40:59', '03:45:06', '03:49:09', '03:53:12', '03:57:14',          '04:01:17', '04:05:26', '04:09:28', '04:13:31', '04:17:34', '04:21:37', '04:25:48', '04:29:51', '04:33:53', '04:37:57',          '04:41:59']
dates = range(71)

# plot airmass vs magnitude - ext.coeff. and ZP fit
for star in range(0,numstars):
    mag = np.genfromtxt("cig/mag_star_00" + str(star+1) + ".dat")
    
    # linear fit per star
    k , z = np.polyfit(airmass,mag,1)
    print "k_ext = ", k
    print "ZP    =", z
    print "Max. variation in magnitude = ", str(round(max(mag)-min(mag),2)), "mag"
    # save k_ext and zp
    kext.append(k)
    zp.append(z)

    # plot mag vs airmass and linear fit
    plt.plot(airmass, mag, 'o', airmass, z + k*airmass, 'r-')
    plt.title('Magnitude vs airmass - Star ' + str(star+1))
    plt.xlabel("AIRMASS")
    plt.ylabel("MAG (Cousins R)")
    plt.legend(['MAGNITUDE','Linear fit'],loc=0, fontsize='x-large')
    plt.show()
    
    # plot magnitude evolution with time
    plt.xticks(dates, imdates, rotation=70, size=12)
    plt.plot(dates, mag, 'gd-')
    plt.title('Magnitude evolution of star ' + str(star+1))
    plt.xlabel("DATE")
    plt.ylabel("MAG (Cousins R)")
    plt.legend(['MAGNITUDE'],loc='upper right', fontsize='x-large')
    plt.show()


# In[ ]:

plt.close()


# In[7]:

kext_list = []
zp_list = []

for star in range(0,numstars):
    mag = np.genfromtxt("cig/mag_star_00" + str(star+1) + ".dat")

    # linear fit
    k , z = np.polyfit(airmass,mag,1)
    
    # save k_ext and zp
    kext_list.append(k)
    zp_list.append(z)
    
    # plot data and linear fit
    if star == 1:
        plt.plot(airmass, mag, 'go', airmass, z + k*airmass, 'r-')
    else:
        plt.plot(airmass, mag, 'o', airmass, z + k*airmass, 'r-')
        plt.title('Magnitude vs airmass - all stars ')
        plt.xlabel("AIRMASS")
        plt.ylabel("MAG (Cousins R)")
        #plt.ylim([20.4,21.0])
        plt.legend(['MAGNITUDE','Linear fit'],loc=0, fontsize='x-large')
        
        # printing summaries
        print "Star ", star+1
        print "k_ext =", round(k,3)
        print "ZP    =", round(z,3)
        print "Max. magnitude variation =", str(round(max(mag)-min(mag),3)), "mag \n"    
#plt.show()


# In[8]:

print "Extinction coefficients: ",
print kext_list, "\n"
print "k_ext median = ", np.median(kext_list), "+/-(std.dev.)", np.std(kext_list)


# ### 12.2. Zero Point of the night
# By subtracting the converted R magnitude from SDSS to the list of zero points of each fit ("zp"), we can obtain the Zero Point of the night.

# In[9]:

# Zero Point of the Night = SDSS magnitudes 'Rsdss' subtracted to the individual zero points
ZPN = [(ZP - zp) for zp, ZP in zip(zp_list, Rsdss)]

print "Zero point: "
print ZPN, "\n"
print "ZP median = ", np.median(ZPN), "+/-(std.dev.)", np.std(ZPN)


# ### 12.3. Summary: Extinction Coefficient and Zero Point of the night

# In[10]:

# extinction coefficient of the night rounded to 6 decimals and ZP of the night
KEXTN = round(np.median(kext_list),4)
ZPN_N = round(np.median(ZPN),4)

print "k_ext median =", KEXTN
print "ZP median =", ZPN_N


# ##13. Extinction correction
# First, the **extinction coefficient** and **zero point of the night** are **added to the headers** of all images.

# In[12]:

# list of images
images_list = np.genfromtxt("cig/images_a.list", dtype=str)

# ext.coeff. addition to headers
for im in images_list:
    data, header = fits.getdata("cig/"+im, header=True)
    header.set("K_EXT_N", KEXTN, "Extinction coefficient of the night - Aug16 calibration")  # adds KEXTN as a new header line
    header.set("ZP_N", ZPN_N, "Zero Point of the night - Aug16 calibration")  # adds ZPN as a new header line
    fits.writeto("cig/"+im, data, header, clobber=True)


# Finally, now that the extinction coefficient of the night is known and written in each header, all images must be **corrected from extinction** before combination.
# 
# The equation for the extinction correction is:  $$m_{corr} = m_{obs} - \kappa_{ext} * X $$  where $\kappa_{ext}$ is the extinction coefficient, $m_{corr}$ is the magnitude as measured without extinction and $m_{obs}$ the extinguished version of it, we can get:  $$ m_{corr} - m_{obs} = -2.5 * \log \left( \frac{flux_{corr}}{flux_{obs}} \right) \Rightarrow flux_{corr} = flux_{obs} * 10^\left( \frac{\kappa_{ext} X}{2.5}\right)$$
# 
# **To achieve this, we create an IRAF script that will run *imarith* over all images and do the previous operation over all pixels.**

# In[68]:

# list of images
imlist = np.genfromtxt("cig/images_a.list", dtype=str)

# IRAF script file
extcorr = open("cig/extinction_correction.cl", "wa")

# IRAF commands to apply imarith
for fitsim in imlist:
    data, header = fits.getdata("cig/" + fitsim, header=True)   # reads the header and data of the fits files
    airm = header["AIRMASS"]
    kextn = header["K_EXT_N"]
    #print airm, kextn, airm*kextn, 10**(airm*kextn/2.5)
    operand2 = 10**(kextn*airm/2.5)                             # flux correction factor: 10^(k*X/2.5)
    output = fitsim[:-5] + "-e.fits"                            # output name for flux corrected image
    extcorr.write("imarith operand1=" + fitsim +" op=* operand2=" + str(operand2) + " result=" + output + "\n")   # IRAF script


# ##### *IRAF*
# The script created will perform operations on each image, like the following:

                imarith operand1=cig96_20120911_CousinsR_001-b-o-f-c-a.fits op=* operand2=1.24686313448 result=cig96_20120911_CousinsR_001-b-o-f-c-a-e.fits
imarith operand1=cig96_20120912_CousinsR_001-b-o-f-c-a.fits op=* operand2=1.23647091537 result=cig96_20120912_CousinsR_001-b-o-f-c-a-e.fits
imarith operand1=cig96_20120912_CousinsR_002-b-o-f-c-a.fits op=* operand2=1.23244038749 result=cig96_20120912_CousinsR_002-b-o-f-c-a-e.fits
imarith operand1=cig96_20120912_CousinsR_003-b-o-f-c-a.fits op=* operand2=1.22859909598 result=cig96_20120912_CousinsR_003-b-o-f-c-a-e.fits
[...]

-----------------------------

iraf_prompt> cl < extinction_correction.cl
                
# ## 14. Combination by median
# ###NOTE: IF EXPOSURE TIMES ARE DIFFERENT, STEP 15 'DIVISION BY TIME' ***MUST*** BE DONE BEFORE THIS STEP.
# In ***IRAF***, we make the following operation:

                imcombine input=*-e.fits output=cig96-b-o-f-c-a-e.fits combine=median scale=median
                
# ## 15. Division by exposure time
# In ***IRAF***, we make the following operation:

                imarith operand1=cig96-b-o-f-c-a-e.fits op=/ operand2=200 result=cig96-b-o-f-c-a-e-t-def.fits
                
# This operation will generate a fully reduced and calibrated image of CIG96.
# 
# ### 15.1 Header keywords redefinitions (exptime and airmass) for final image
# 
# The header of the resultant image must have now the correct EXPTIME (=1s) and AIRMASS (=1)

# In[14]:

# EXPTIME and AIRMASS modification

# specify resultant image
im = "cig96-b-o-f-c-a-e-t-def.fits"

data, header = fits.getdata("cig/"+im, header=True)

header.set("AIRMASS", 1, "Image corrected from extinction - Aug16 calibration")  # redefines AIRMASS keyword
header.set("EXPTIME", 1, "Image divided by original exposure time of 200s - Aug16 calibration")  # redefines EXPTIME keyword
header.set("FILENAME", im, "Calibrated image") # redefines FILENAME keyword

fits.writeto("cig/"+im, data, header, clobber=True)


# ## 16. Surface brightness limit
# 
# In order to know the depth of the combined image, we have to find out the residuals rms, below which we cannot detect anything.
# 
# Since we have **NOT** subtracted the sky we apply the *physical rms* defined as:
# 
# $$phys.rms =  \sqrt{ \frac{ \sum_{i}^{n} \left( x_{i} - x_{median.sky} \right) ^{2} }{n} }$$
# 
# where $x_{median.sky}$ is the median value of the sky of the combined image  and $x_{i}$ is the total sky in $counts/s$.
# 
# The pixel scale of the instrument is:

# In[11]:

#pixel scale in arcsec²:
pix = 1.0466293
pixarea = pix**2


# We perform the calculations over 3 different areas of the image, since they clearly present different intensity (likely artifacts caused by CAFOS).
# 
# ### 16.1 Region 1: southern region (dark)
# 
# We take a the median values of a large number boxes placed in the dark southern area, compute the average and calculate the *phys.rms* and then the SB of this area:

# In[12]:

# insert the median values of the boxes and divide by pixel area
medians = [50.40, 50.34, 50.32, 50.25, 50.21, 50.29, 50.38, 50.35, 50.25, 50.30, 50.24, 50.34, 50.30, 50.37, 50.32, 50.38]
meas1 = [value / pixarea for value in medians]

# compute average of the median values per arcsecond²
mean1 = np.average(meas1)
print "Average of the median values (region 1) =", mean1, "counts/s/arcsec²"

#calculate the physical rms
meansq1 = []
for elem in meas1:
    sq = (elem - mean1)**2
    meansq1.append(sq)

rms1 = np.sqrt(np.average(meansq1))
print "rms(region 1) =", rms1, "counts/s/arcsec²"

#calculate the SB depth
term1 = -2.5*np.log10(rms1)
print "Depth (region 1)=", term1, "mag/arcsec² \n"
SB1 = ZPN_N + term1
print "SB(region 1) =", SB1, "mag/arcsec²"


# ### 16.2 Region 2: eastern region (bright)
# 
# We take a the median values of a large number of boxes placed in the *bright eastern area*, compute the average and calculate the *phys.rms* and then the SB of this area:

# In[13]:

# insert the median values of the boxes
medians = [50.61, 50.62, 50.64, 50.66, 50.65, 50.66, 50.72, 50.70, 50.69, 50.72, 50.71, 50.72, 50.73, 50.75, 50.76, 50.74, 50.74, 50.74, 50.72, 50.73, 50.96, 50.73]
meas2 = [value / pixarea for value in medians]

# compute average of the median values
mean2 = np.average(meas2)
print "Average of the median values (region 2) =", mean2, "counts/s/arcsec²"

#calculate the physical rms
meansq2 = []
for elem in meas2:
    sq = (elem - mean2)**2
    meansq2.append(sq)

rms2 = np.sqrt(np.average(meansq2))
print "rms(region 2) =", rms2, "counts/s/arcsec²"

#calculate the SB depth
term2 = -2.5*np.log10(rms2)
print "Depth (region 2) =", term2, "mag/arcsec² \n"
SB2 = ZPN_N + term2
print "SB(region 2) =", SB2, "mag/arcsec²"


# ### 16.3 Region 3: northern region (darkish)
# 
# We take a the median values of a large number of boxes placed in the *darkish northern area*, compute the average and calculate the *phys.rms* and then the SB of this area:

# In[14]:

# insert the median values of the boxes
medians = [50.60, 50.58, 50.47, 50.46, 50.45, 50.55, 50.60, 50.61, 50.59, 50.54, 50.45, 50.43, 50.37, 50.40, 50.48, 50.52, 50.50, 50.43]
meas3 = [value / pixarea for value in medians]

# compute average of the median values
mean3 = np.average(meas3)
print "Average of the median values (region 3) =", mean3, "counts/s/arcsec²"

#calculate the physical rms
meansq3 = []
for elem in meas3:
    sq = (elem - mean3)**2
    meansq3.append(sq)

rms3 = np.sqrt(np.average(meansq3))
print "rms(region 3) =", rms3, "counts/s/arcsec²"

#calculate the SB depth
term3 = -2.5*np.log10(rms3)
print "Depth (region 3) =", term3, "mag/arcsec² \n"
SB3 = ZPN_N + term3
print "SB(region 1) =", round(ZPN_N + term1, 3), "mag/arcsec²"
print "SB(region 2) =", round(ZPN_N + term2, 3), "mag/arcsec²"
print "SB(region 3) =", round(ZPN_N + term3, 3), "mag/arcsec²"
SB = [SB1, SB2, SB3]
print "Standard deviation from the 3 values = +/-", round(np.std(SB), 3)


# #END
# 
# ------
# 
# #Appendix A: SB isophotes comparison with SDSS
# 
# SDSS images are in r (or any other SDSS) filter so we have to compute the difference to compare the isophotes between CAHA R Cousins image and SDSS r band image.
# 
# SDSS magnitudes are computed by equation: $m = 22.5 - 2.5*log_{10} (F)$ , measured in $mag$.
# 
# ***IRAF***
# 
# Via qphot and pdump over the final image, we extract the magnitudes of the selected stars and correct them with the SP of the night (ZPN_N):

# In[15]:

#CAHA stars magnitudes
R_CAHA = [-6.34, -5.569, -7.663, -5.194, -5.541, -5.394, -4.845, -4.25]
R_CAHA = [elem + 24.28 for elem in R_CAHA]

print R_CAHA


# In[17]:

#we compute the difference between R Cousins and r SDSS filters so we can compare correctly the two images
print "r_sdss =", r_sdss
print "R_CAHA =", R_CAHA, "\n"
diff_Rr = np.median([a - b for a,b in zip(r_sdss, R_CAHA)])
print "Median difference between r and R:", diff_Rr, "mag"

#median sky for CAHA image
#NOTE: the value 1.0466293**2 is the quadratic pixel scale so we convert the flux from counts/s/arcsec2 to counts/s
skies = [mean1, mean2, mean3]
medsky = np.median(skies)*(1.0466293**2)
print "The median sky value is", round(medsky,4), "counts/arcsec²\n"

#desired depth
SB = 23.5

#values for the contours
#NOTE: the value 0.39609712**2 is the quadratic pixel scale so we convert the flux from counts/s/arcsec2 to counts/s
fSDSS = (10**((SB - 22.5)/(-2.5)))*(0.39609712**2)
fCAHA = medsky + 10**((SB - ZPN_N - diff_Rr)/(-2.5))
print "In order to get a SB of", SB, "mag/arcsec², these are the necessary fluxes for CAHA and SDSS images:\n"
print "flux_CAHA =", fCAHA, "counts/s"
print "flux_SDSS =", fSDSS, "counts/s"


# **NOTE**: in CAHA R-Deep image:
# 
# - SB(**Features in the west**) between 26.2 and 26.6 mag/arcsec2
# - SB(**Pseudoring**) between 25.4 and 25.7 mag/arcsec2

# ------
# 
# #Appendix B: cropping FITS file
# 
# The image has the overscan region present in the right side (non astrometrical orientation). We may crop it to avoid weird visualizations as follows:

# In[23]:

# read the image
image = fits.open('cig/cig96_def.fits', mode='update')
image.info()

# select layer (in this case, it is unique)
im = image[0]
im.shape

# crop image from x0 to x1 and from y0 to y1
im.data = im.data[0:990,0:1000]

im.writeto('cig/cig96_def_crop.fits')#, clobber=True)


# In[ ]:

def cropper(image, x0, x1, y0, y1, outfile):
    '''
    Crops the selected region of any 
    FITS image provided and saves it.
    '''
    img = fits.open(image, mode='update')
    im = img[0].data
    im.data = im.data[x0:x1,y0:y1]
    im.writeto(outfile, clobber=True)


# -------------------
# 
# # Appendix C: dimension keywords cleaning, image coordinates resetting and backup
# 
# When using APLpy, if the dimension keywords "A_0_0" (and alike) are repeated, APLpy can not just choose one so it crashes.
# 
# ***IMPORTANT***: each fits can only contain **one and only one** set of A_0_0 and related dimension keywords.
# 
# **Solution**: remove all the dimension keywords with **missfits** task and reassign coordinates as follows.

                %%bash
missfits -REMOVE_KEYWORD A_?,B_?,AP_?,BP_? cig96_def.fits
wipe_wcs.py cig96_def.fits
astrom.py --overwrite --radius 1 cig96_def.fits
                
# **Backup** of the final image.

                %%bash
cd cig
#mkdir final_backup
cp cig96_def_crop.fits ./final_backup
                
# -------------------
# 
# # Appendix D: image deprojection
# 
# Using IRAF, the image can be deprojected. The necessary parameters are the morphological type, the axis relation and inclination.
# 
# 
# - CIG96 morphological type = Sc = 5:      **q (Sc) = 0.175**
# 
# - Axis relation (minor over major axis):    **b/a = 0.698821608** (extracted from HI data)
# 
# - Inclination (sources: Jakoby, Masters et al. 2014):   cos² i = ( (b/a)² - q² ) / (1 - q²):       **i = 46,594368772 deg = 46,59 deg**

                IRAF command:

geotran  input=cig96_nodep.fits  output=cig96_dep.fits  xin=507  yin=538  xrot=69.61  yrot=69.61  xscale=0.699
                
# ------
# 
# #Appendix E: contrast images
# 
# Contrast images created to enhance structures by removing a smoothed version of the image itself.
# 
# ###E1. Gaussian and Boxcar smoothing
# 
# Define a new smoothed image by setting up a number of kernels (in pixels) and the new images names.
# 
# - <a href="http://stackoverflow.com/questions/17595912/gaussian-smoothing-an-image-in-python">Gaussian smoothing an image in python</a>, where: FWHM = kernel/(2*np.sqrt(2*np.log(2)))
# - <a href="https://matthew-brett.github.io/teaching/smoothing_intro.html">Gaussian kernel and formula</a>
# - <a href="https://docs.scipy.org/doc/scipy-0.7.x/reference/generated/scipy.stsci.convolve.boxcar.html">Boxcar kernel definition in SCIPY</a>

# In[76]:

# gaussian kernels (size in pixels)
kernels = [1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100]

# image selection
image = "cig/cig96_def_crop.fits"

# gaussian smoothing
im = fits.open(image)

for kernel in kernels:
    im[0].data = scipy.ndimage.gaussian_filter(im[0].data, kernel/(2*np.sqrt(2*np.log(2))))
    # smoothed image
    im.writeto(image[0:-5] + "_gauss" + str(kernel) + ".fits", clobber=True)

im.close()

# gaussian kernels (size in pixels)
kernelsbox = [1,2,3,4,5,6,7,8,9,10,15,20]

# image selection
image = "cig/cig96_def_crop.fits"

# boxcar smoothing
im = fits.open(image)

for kernel in kernelsbox:
    boxshape = [kernel,kernel]
    im[0].data = convolve.boxcar(im[0].data, boxshape, output=None, mode='wrap', cval=0.0)
    # smoothed image
    im.writeto(image[0:-5] + "_boxcar" + str(kernel) + ".fits", clobber=True)

im.close()


# In[28]:

# VPeris images gaussian filter

# gaussian kernels (size in pixels)
kernels = [1,2,3,4,5,10,15,20,25,30]

# image selection
image = "v_peris/cig96_suma_peris.fit"

# gaussian smoothing
im = fits.open(image)

for kernel in kernels:
    im[0].data = scipy.ndimage.gaussian_filter(im[0].data, kernel/(2*np.sqrt(2*np.log(2))))
    # smoothed image
    im.writeto(image[0:-5] + "_gauss" + str(kernel) + ".fits", clobber=True)

im.close()


# ###E2. Subtraction of smoothed images
# 
# Each smoothed image is subtracted to the original image.

# In[77]:

# open original image selection
image = "cig/cig96_def_crop.fits"
img = fits.open(image)
im = img[0].data

for kernel in kernels:
    # open smoothed image
    smooth = fits.open(image[0:-5] + "_gauss" + str(kernel) + ".fits")
    sm = smooth[0].data
    # subtraction
    sub = im - sm
    # writing the resulting array as a PrimaryHDU element in the new image
    newhdu = fits.PrimaryHDU(data=sub, header=img[0].header)
    newhdu.writeto(image[0:-5] + "_smooth_gauss" + str(kernel) +".fits", clobber=True)
# source to bring array to HDU:
# http://stackoverflow.com/questions/30858393/fits-image-array-pixel-out-of-bounds?rq=1

for kernel in kernelsbox:
    # open smoothed image
    smooth = fits.open(image[0:-5] + "_boxcar" + str(kernel) + ".fits")
    sm = smooth[0].data
    # subtraction
    sub = im - sm
    # writing the resulting array as a PrimaryHDU element in the new image
    newhdu = fits.PrimaryHDU(data=sub, header=img[0].header)
    newhdu.writeto(image[0:-5] + "_smooth_boxcar" + str(kernel) +".fits", clobber=True)


# In[ ]:



