
# coding: utf-8

# #Visualization of galaxies observed with VST
# 
# - Notebook for cropping VST images <a href="http://localhost:8888/notebooks/DeSIGN/Utilities/Crop%20VST%20fits%20files%20and%20SB%20calculation.ipynb">here</a>.

# In[10]:

import aplpy
import astropy.io.fits as fits
from   astropy import units as u
from   astropy import stats as sts
import copy
from   glob import glob as ls
from   kapteyn import maputils
import lemon.astrometry as astrometry
import lemon.photometry as photometry
import math
import matplotlib as mpl
import matplotlib.colors
import matplotlib.pyplot as plt
import montage_wrapper as montage
import numpy as np
from   numpy.linalg import eig, inv
import os
import pyfits
import pylab
import pyraf as pyraf
import pyregion as pyreg
import random
import repipy.combine as combine
import repipy.find_sky as find_sky
import scipy
from   scipy import stats
from   scipy.misc import imread
import shutil
from   spectral_cube import SpectralCube as spc
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")
#import pyds9 as ds9
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

# renders interactive figures in the Notebook
get_ipython().magic(u'matplotlib nbagg')

# renders figures as static PNG
#%matplotlib inline

# manually defines plots sizes
#mpl.pylab.rcParams['figure.figsize'] = 12, 12  # that's default image size for this interactive session

# FUNCTIONS
###############

def go_to(dir):
    """Function to select and switch to optical or HI folders"""
    # options
    opt = ["opt","Opt","OPT"]
    peris = ["per", "peris", "Peris"]
    hi = ["HI", "Hi", "hi"]
    rgb = ["rgb", "RGB"]
    
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
        raw = "/home/prm/Desktop/cig96-evla-vla-2015"
        work = "/home/prm/Desktop/cig96-evla-vla-2015"
        
    elif dir in rgb:
        # raw data and work directories:
        raw = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16/v_peris/"
        work = "/home/prm/Desktop/optical/optical/CAHA/cig96_jun16/v_peris/"
        
    else:
        return "Invalid option, please select 'opt', 'peris' or 'HI'."

    # switch to selected working directory
    os.chdir(os.path.expandvars(work))
    #print "Work/save directory: ", work
    
go_to('HI')
go_to('opt')
go_to('peris')
go_to('rgb')

###############

def ping():
    os.system('spd-say "Ping"')


# In[11]:

ping()


# ##SB calculation

# In[86]:

# predefined variables
w = 100        # w x w square_width_of_the_box
threshold = 1  # sigma clipping low and high values
iters = 10     # number of iterations for the sigma clipping

# function to calculate the SB limit (according to the rms of the sky)

def SBlim(path, target, bins, col):
    frames_x0 = np.arange(8000,24000,250)   # initial X position of all boxes, scattered from X1 to X2 increments of Z1
    frames_y0 = np.arange(5000,22000,250)   # initial Y position of all boxes, scattered from Y1 to Y2 increments of Z2
    rms = []
    # Select working directory
    os.chdir(path)
    # target image
    image = fits.open(target, mode='update')
    im = image[0] # loads array of data

    # rms calculation
    for y in frames_y0:
        for x in frames_x0:   # loop over the X * Y boxes
            box = im.data[x:x+w,y:y+w]
            sc_box = sts.sigma_clip(box, iters=iters, sigma_lower=threshold, sigma_upper=threshold)
            rms_val = np.sqrt(np.mean(np.square(sc_box)))   # n sigma clipped rms calculation of each box
            #rms_val = np.mean(im.data[x:x+w,y:y+w])   # median calculation of each box
            rms.append(rms_val)
            
    # SB numerical calculation
    rms_med = np.median(rms)
    
    
    # SB calculation by histogram
    SB = [-2.5*np.log10(elem) + 2.5*np.log10(0.21**2) for elem in rms]
    n, b, patches = plt.hist(SB, bins=bins, range=(25,31), color=col, alpha=0.3)
    plt.show()

    # Finding the X mag/arcsec2
    for val in range(0,len(n)):
        elem = n[val]
        if elem == n.max():
         break
    else:   # ideally this should never be tripped
        val = none
    # useful prints
    print len(rms)
    print "rms_clipped_median value =", rms_med
    print "SB_num from rms_clipped_median =", -2.5*np.log10(rms_med) + 2.5*np.log10(0.21**2)
    print "SB limit of "+str(target)+" =", b[val], "mag/arcsec2 - in", str(col)
    # close image
    image.close()


# In[87]:

bins = 60
SBlim("/home/prm/Desktop/VST/cig11", "cig11_vdef.fits", bins, 'blue')
SBlim("/home/prm/Desktop/VST/cig33", "cig33_vdef.fits", bins, 'red')
SBlim("/home/prm/Desktop/VST/cig59", "cig59_vdef.fits", bins, 'magenta')
SBlim("/home/prm/Desktop/VST/cig96", "cig96_VST_coadd.fits", bins, 'forestgreen')
SBlim("/home/prm/Desktop/VST/cig152", "cig152_vdef.fits", bins, 'gold')
SBlim("/home/prm/Desktop/VST/cig154", "cig154_vdef.fits", bins, 'lime')
SBlim("/home/prm/Desktop/VST/cig279", "cig279_vdef.fits", bins, 'maroon')
SBlim("/home/prm/Desktop/VST/cig568", "cig568_vdef.fits", bins, 'black')
SBlim("/home/prm/Desktop/VST/cig1002", "cig1002_vdef.fits", bins, 'orange')
SBlim("/home/prm/Desktop/VST/cig1027", "cig1027_vdef.fits", bins, 'cyan')
SBlim("/home/prm/Desktop/VST/cig1047", "cig1047_vdef.fits", bins, 'pink')

ping()


# In[88]:

plt.close()


# ##Visualization of the targets

# In[90]:

# common colour map
colmap = "rainbow"


# ###CIG11

# In[91]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig11"
os.chdir(path)

# select and open images
image  = "cig11_vdef_7arcmin_cropped.fits"
imagez = "cig11_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 3.6332314
DEC = -0.73685996
rad = 0.14

###########################

# variables definition
stretch='log'
v_min=25
v_max=29
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)
# contours
levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.3, 27.7]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[92]:

cig.close()


# ###CIG33

# In[93]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig33"
os.chdir(path)

# select and open images
image  = "cig33_vdef_7arcmin_cropped.fits"
imagez = "cig33_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 10.866111
DEC = -0.12483331
rad = 0.16

###########################

# variables definition
stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)
# contours
levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5, 28.0, 28.5, 28.9]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[94]:

cig.close()


# ###CIG59

# In[95]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig59"
os.chdir(path)

# select and open images
image  = "cig59_vdef_7arcmin_cropped.fits"
imagez = "cig59_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 24.587031
DEC = 7.5346537
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5, 28.0, 28.4, 28.8]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[96]:

cig.close()


# ###CIG152

# In[120]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig152"
os.chdir(path)

# select and open images
image  = "cig152_vdef_7arcmin_cropped.fits"
imagez = "cig152_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 67.748514
DEC = -2.0033558
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5, 28.0, 28.5]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[121]:

cig.close()


# ###CIG154

# In[99]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig154"
os.chdir(path)

# select and open images
image  = "cig154_vdef_7arcmin_cropped.fits"
imagez = "cig154_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 71.9385
DEC = 1.8191
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5, 28.0, 28.2]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[100]:

cig.close()


# ###CIG279

# In[111]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig279"
os.chdir(path)

# select and open images
image  = "cig279_vdef_7arcmin_cropped.fits"
#imagez = "cig279_vdef_3arcmin_cropped.fits"
imagez = "cig279_vdef_4.5arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 130.38419
DEC = 4.981089
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5, 28.0, 28.2]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[112]:

cig.close()


# ###CIG568

# In[123]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig568"
os.chdir(path)

# select and open images
image  = "cig568_vdef_7arcmin_cropped.fits"
imagez = "cig568_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 196.06271
DEC = 9.223596
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.4, 26.6, 26.8, 27.0, 27.2, 27.3, 27.4]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[124]:

cig.close()


# ###CIG1002

# In[105]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig1002"
os.chdir(path)

# select and open images
image  = "cig1002_vdef_7arcmin_cropped.fits"
imagez = "cig1002_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 345.17068
DEC = 8.4676962
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[106]:

cig.close()


# ###CIG1027

# In[ ]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig1027"
os.chdir(path)

# select and open images
image  = "cig1027_vdef_7arcmin_cropped.fits"
imagez = "cig1027_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 353.859
DEC = 7.322
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[108]:

cig.close()


# ###CIG1047

# In[109]:

# Select working directory
path = "/home/prm/Desktop/optical/optical/VST/cig1047"
os.chdir(path)

# select and open images
image  = "cig1047_vdef_7arcmin_cropped.fits"
imagez = "cig1047_vdef_3arcmin_cropped.fits"
im  = fits.open(image)
imz = fits.open(imagez)

# coordinates in deg
RA = 359.198
DEC = 1.3550 
rad = 0.16

###########################

stretch='log'
v_min=24
v_max=30
vmin = 10**(-0.4*v_max - 1.355561411)
vmax = 10**(-0.4*v_min - 1.355561411)

levs = [22, 23, 24, 25, 26.0, 26.5, 27.0, 27.5, 27.8]
lev = [10**(-0.4*val - 1.355561411) for val in reversed(sorted(levs))]
sm = 3

###########################

# figure structure
fig = plt.figure(figsize=(16, 8))

# plot window design
vw = 1.0 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

###########################

# LEFT plot
cig = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True, subplot=[0.00, 0.00, hw, vw])

#ticks and labels
# ticks
cig.ticks.show()
cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=16)
cig.axis_labels.set_xposition('bottom')
cig.tick_labels.set_xposition('bottom')
cig.tick_labels.set_xformat('hh:mm:ss')
cig.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cig.axis_labels.set_font(size=18)
cig.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')

# colorscale
cig.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cig.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

###########################

# RIGHT plot
cigz = aplpy.FITSFigure(imz, auto_refresh=True, north=True, figure=fig, subplot=[0.50, 0.00, hw, vw])

#ticks and labels
# ticks
cigz.ticks.show()
cigz.axis_labels.set_xposition('bottom')
cigz.tick_labels.set_xposition('bottom')
cigz.axis_labels.set_yposition('right')
cigz.tick_labels.set_yposition('right')
cigz.tick_labels.set_xformat('hh:mm:ss')
cigz.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')

cigz.axis_labels.set_font(size=18)
cigz.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')

# colorscale
cigz.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone_r')

# zoom
#cig.recenter(RA, DEC, rad)

# contours
cigz.show_contour(levels=lev, smooth=sm, cmap=colmap, linewidths=1, returnlevels=True)

# save figure
plt.savefig(path+'/'+image[:-5]+'_contours_sm'+str(sm)+'pix_inner_'+str(levs[-1])+'.png', bbox_inches='tight')
ping()


# In[110]:

cig.close()


# In[ ]:



