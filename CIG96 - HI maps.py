
# coding: utf-8

# # CIG96 - HI maps
# 
# ##GENERAL INFORMATION
# 
# Notebook with some plots of the HI maps of CIG96.
# 
# ### AMIGA WIKI PAGES:
# 
# **HI data**
# - Info from the EVLA campaign can be found <a href="http://amiga.iaa.es:8888/display/science/EVLA+Aug12+CIG96+C+D">here</a>.
# 
# - The data cube calibration notes and scripts can be found <a href="http://amiga.iaa.es:8888/display/science/Data+calibration+Full+%28VLA07+-+EVLA13%29+CIG96">here</a>.
# 
# - Information about CIG96 paper (2016) can be found <a href="http://amiga.iaa.es:8888/display/science/CIG96+-+Deep+optical+and+HI+images#CIG96-DeepopticalandHIimages-Paperstructureandcomments">here</a>.
# 
# ###PYTHON REQUIREMENTS:
# 
# Some VERY interesting (and in most cases necessary) packages are:
# 
# - astropy (calculations and FITS handling)
# 
# - aplpy (plots)
# 
# - photutils (photometric calculations)
# 
# - spectralcube (3D and 4D cubes handling)
# 
# - kapteyn (channels map)
# 
# ###a) Import necessary modules

# In[1]:

import aplpy
import astropy.io.fits as fits
from   astropy import units as u
from   astropy import wcs as wcs
from   glob import glob as ls
from   kapteyn import maputils
import lemon.astrometry as astrometry
import lemon.photometry as photometry
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
from   scipy import stats, optimize
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
    
###########

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

###########

# Emission of a sound when task finished

def ping():
    os.system('spd-say "Ping"')


# In[4]:

ping()


# # 1. HI moments
# 
# After converting a CASA cube to FITS, we may read it, plot it and compute its moments here.
# 
# <a href="http://www.cv.nrao.edu/~aleroy/pytut/topic2/intro_fits_files.py">Useful guide</a>.
# 
# We select CIG96's wavelet filtered datacube with circular beam: *cig96.w.circ.HI.Filt.image.fits*
# 
# ##1.1 0th moment
# 
# ###Calculation via Numpy

                # Select working directory
go_to('cig')

# select and open datacube image
image = "cig96.w.circ.HI.Filt.image.fits"
im = fits.open(image)

# define the cube
cube = im[0].data
#print cube.shape                 # shape = [stokes, channels, RA, DEC]
#pylab.imshow(cube[0,24,:,:])     # will show channel 25

# collapse moment 0
mom0 = np.sum(cube, axis=1)
mom0 = np.sum(mom0, axis=0)       # "sums" stokes axis: necessary to avoid dimension misunderstandings
mom0 = np.transpose(mom0, axes=[-0,1])   # transpose axes to reorient image

# collapse mom0 with masking variable sigma
rms = 0.000126
mask = (cube > 1*rms)
mom0_masked = np.sum(cube*mask, axis=1)
mom0_masked = np.sum(mom0_masked, axis=0)
mom0_masked = np.transpose(mom0_masked, axes=[-0,1])

# mom0 rms
rms0 = np.std(mom0_masked[0:][0:150])
print "rms (mom0) =", rms0 #= 0.000816681

# plot
plt.imshow(mom0_masked, origin='lower', cmap='magma')    # force RA,DEC start bottom left
#plt.contour(mom0_masked, levels=[0.002, 0.003, 0.004, 0.008, 0.012, 0.016, 0.02, 0.03, 0.04, 0.05], cmap='RdYlBu_r')
plt.contour(mom0_masked, cmap='RdYlBu_r',levels=[3*rms0, 3.5*rms0, 5*rms0, 7*rms0])#, 4*rms0, 5*rms0, 6*rms0, 7*rms0, 8*rms0, 9*rms0, 12*rms0, 15*rms0])
#plt.xlim(80,570)
#plt.ylim(80,500)

#######

plt.close()
                
# ###Calculation via "Spectral_cube" module

# In[33]:

# Select working directory
go_to('cig')

# select datacube
image = "cig96.w.circ.HI.Filt.image.fits"
# reads cube
cube = spc.read(image, fill_value=0, wcs="wcs")   
# changes integration axis from frequency to velocity
cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=1.420405751770 * u.GHz)
print cube

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.21   # zoom radius so, the larger rad, the larger FoV in the image

# extract the subcube between 0.000126 and inf Jy
#slab = cube.spectral_slab(0.000126 * u.Jy, 1000 * u.Jy)

# masking, ignore elements fainter than factor*rms (in Jy)
rms = 0.000126  # known from HI data calibration
factor = 3.5
cube_mask = cube.with_mask(cube > factor*rms * u.Jy)  

# compute 0th moment
mom0 = cube_mask.moment(order=0, how='cube')
# NO MASK:
#mom0 = cube.moment(order=0, how='cube')   # NO mask
mom0_save = 'cig96_HI_mom0_' + str(factor) + 'rms.fits'
mom0.write(mom0_save, overwrite=True)

# masked mom0 rms (via numpy computation of 0th moment)
#im = fits.open(image)
#cube = im[0].data
#mom0 = np.sum(cube, axis=1)
#mom0 = np.sum(mom0, axis=0)              # "sums" stokes axis: necessary to avoid dimension misunderstandings
#mom0 = np.transpose(mom0, axes=[-0,1])   # transposes axes to reorient image
#rms = 0.000126                           # known from HI data calibration
#mask = (cube > 1*rms)
#mom0_masked = np.sum(cube > mask, axis=1)
#mom0_masked = np.sum(mom0_masked, axis=0)
#mom0_masked = np.transpose(mom0_masked, axes=[-0,1])
im = fits.open(mom0_save)
data = im[0].data
rms0 = np.sqrt(np.mean(np.square(data[10:590][10:150]))) ### rms computation differs from CASA's value, which is: ~0.00827
print "rms (mom0) =", rms0, "Jy"
im.close()

# plot
m0 = aplpy.FITSFigure(mom0_save, north=True, auto_refresh=True)
m0.show_colorscale(stretch='linear', cmap='jet')#, vmin=0, vmax=3000)
m0.set_nan_color('black')
m0.show_colorbar()
m0.recenter(RA,DEC,rad)

# beam definition
m0.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.99))

# contour levels definition
rms0 = 0.0078
#m0.show_contour(levels=[1.5*rms0, 3*rms0, 5*rms0, 7*rms0, 9*rms0, 15*rms0], cmap='Greys')
#m0.show_contour(levels=[0.034], cmap='Greys')

# save image
#m0.save("cig96-mom0.png")


# In[32]:

m0.close()


# ###Mom0 over optical image

# In[5]:

# Select working directory
go_to('cig')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.19

# optical and mom0 images
optim = "cig96_VST_coadd_12.5arcmin_cropped_SB.fits"
# mom0 = "cig96_HI_mom0_3.5rms_coldens.fits"      # original
#optim = "cig96_RdRGB_SB_oriented.fits"
mom0  = "cig96.w.circ.HI.Filt.image.3.5sigma.mom0.fits"

# plot
fig = plt.figure(figsize=(9,9))
m0opt = aplpy.FITSFigure(optim, figure=fig, north=True, auto_refresh=True)
#m0 = aplpy.FITSFigure(mom0, north=True, auto_refresh=True)
m0opt.axis_labels.set_font(size=20)
m0opt.set_tick_labels_font(size=17)
m0opt.tick_labels.set_xformat('hh:mm:ss')
m0opt.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=8, width=1.5, color='black')

# color scale
m0opt.show_colorscale(vmin=26, vmax=28, stretch='linear', cmap='bone')

# contours
levels = [0.004, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.57, 0.75, 0.9, 1.0] # flux dens.
m0opt.show_contour(data=mom0, levels=levels, cmap='rainbow', linewidths=2)
levels_border = list(np.arange(0.0040, 0.0150, 0.0001, dtype='float'))
m0opt.show_contour(data=mom0, levels=levels_border, colors='purple', linewidths=0.25) # border contours

#m1.show_regions("arms_ALL_def_wcs.reg")
m0opt.recenter(RA, DEC, radius=rad)
#m1.show_axis_labels()

# beam
m0opt.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor='y', edgecolor='k')

# save image
m0opt.save("/home/prm/Desktop/paper_cig96/images/cig96_mom0_overlay_optical.png")


# In[6]:

m0opt.close()


# In[95]:

# conversion to NHI
levels = [0.004, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.57, 0.75, 0.9, 1.0] # flux dens.
lev = []
for elem in levels:
    nhinum = np.round(elem*1.411429847*10*10, 4) ### 10^20 cm-2
    lev.append(nhinum)
print lev


# ##1.2 1st moment

# In[8]:

# Select working directory
go_to('cig')

# select datacube
image = "cig96.w.circ.HI.Filt.image.fits"
cube = spc.read(image, fill_value=0, wcs="wcs")

# masking, ignore elements fainter than rms (in Jy)
rms = 0.000126
factor = 5
cube_mask = cube.with_mask(cube > factor*rms * u.Jy)
#cube_save = 'cig96_HI_mom1.fits'
#cube_mask.write(cube_save, overwrite=True)

# change frequency index to radio velocity
cube_mask_vel = cube_mask.with_spectral_unit(u.km / u.s, velocity_convention='radio')
#print cube_vel

# compute 1st moment
mom1 = cube_mask_vel.moment(order=1, how='cube')
mom1_save = 'cig96_HI_mom1_' + str(factor) + 'rms.fits'
mom1.write(mom1_save, overwrite=True)


# In[9]:

# Select working directory
go_to('cig')

# CIG96 coordinates in deg
RA = 33.935231
DEC = 6.0020581
rad = 0.21

# plot
m1 = aplpy.FITSFigure(mom1_save, north=True, auto_refresh=True)
#m1.show_colorscale(vmin=1440, vmax=1650, stretch='linear', cmap='jet')

# contours
velmin = 1380
velmax = 1730
lev = np.linspace(velmin, velmax, num=(velmax-velmin)/10)
m1.show_contour(levels=lev, cmaps='jet')
#m1.show_regions("arms_ALL_def_wcs.reg")
m1.recenter(RA, DEC, radius=rad)
#m1.show_axis_labels()

# beam
m1.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4))

# save image
#m1.save("cig96-mom1.png")


# In[10]:

m1.close()


                print np.linspace(1440, 1700, 14)
                
                np.arange(1440, 1680, 20)
                
                velmin = 1380
velmax = 1730
lev = np.linspace(velmin, velmax, num=(velmax-velmin)/10+1)
print lev
                
# In[75]:

# Select working directory
go_to('cig')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.14

# figure creation
fig = plt.figure(figsize=(10,10))

# image
img = "cig96_HI_mom1_5rms.fits"
#img = "cig96.w.circ.HI.Filt.image.3.5sigma.mom1.fits"

# plot
m1p = aplpy.FITSFigure(img, figure=fig, north=True, auto_refresh=True)
#m1.show_colorscale(vmin=1440, vmax=1650, stretch='linear', cmap='jet')

# all
velmin = 1380
velmax = 1730
lev = np.linspace(velmin, velmax, num=(velmax-velmin)/10+1)
print lev
m1p.show_contour(levels=lev, cmap='YlGn', linewidth=1, zorder=1)
#m1.show_axis_labels()

# contours
# pairs
chn = 1
pairs = [1550 - 10*chn, 1550 + 10*chn]
m1p.show_contour(levels=pairs, colors=['blue', 'red'], linewidths=3, zorder=2)
#m1.show_regions("arms_ALL_def_wcs.reg")

# zoom
m1p.recenter(RA, DEC, radius=rad)

# beam
m1p.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4))

# save image
#m1p.save('cig96_mom1_ch_'+str(pairs[0])+'_'+str(pairs[1])+'.png')


# In[76]:

m1p.close()


# In[145]:

# Select working directory
go_to('cig')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.16

# figure creation 5 in top row + 5 in bottom row
fig = plt.figure(figsize=(20,8))

# % widths
hw = 0.20
vw = 0.50

# image
img = "cig96_HI_mom1_5rms.fits"
#img = "cig96.w.circ.HI.Filt.image.3.5sigma.mom1.fits"

#### top row

# plot: background + contours + zoom
velmin = 1384
velmax = 1734
lev = np.linspace(velmin, velmax, num=(velmax-velmin)/10+1)

v_x = 0.81
v_up = 0.94
v_down = 0.89

length_tx = 8
length_ty = 4
color_blue = 'royalblue'
color_red  = 'crimson'
labelsize = 16
backcolor = 'YlGn'

# pairs
v_appr  = [1544, 1534, 1524, 1514, 1504, 1494, 1484, 1474, 1464, 1454, 1444]
v_reced = [1544, 1554, 1564, 1574, 1584, 1594, 1604, 1614, 1624, 1634, 1644]

m1p1 = aplpy.FITSFigure(img, figure=fig, subplot=[0.00, 0.50, hw, vw], north=True, auto_refresh=True)
m1p1.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p1.recenter(RA, DEC, radius=rad)
m1p1.set_tick_size(0)
m1p1.add_label(v_x, v_down, str(v_appr[1])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p1.add_label(v_x, v_up, str(v_reced[1])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')


m1p2 = aplpy.FITSFigure(img, figure=fig, subplot=[0.20, 0.50, hw, vw], north=True, auto_refresh=True)
m1p2.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p2.recenter(RA, DEC, radius=rad)
m1p2.set_tick_size(0)
m1p2.add_label(v_x, v_down, str(v_appr[2])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p2.add_label(v_x, v_up, str(v_reced[2])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p3 = aplpy.FITSFigure(img, figure=fig, subplot=[0.40, 0.50, hw, vw], north=True, auto_refresh=True)
m1p3.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p3.recenter(RA, DEC, radius=rad)
m1p3.set_tick_size(0)
m1p3.add_label(v_x, v_down, str(v_appr[3])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p3.add_label(v_x, v_up, str(v_reced[3])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p4 = aplpy.FITSFigure(img, figure=fig, subplot=[0.60, 0.50, hw, vw], north=True, auto_refresh=True)
m1p4.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p4.recenter(RA, DEC, radius=rad)
m1p4.set_tick_size(0)
m1p4.add_label(v_x, v_down, str(v_appr[4])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p4.add_label(v_x, v_up, str(v_reced[4])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p5 = aplpy.FITSFigure(img, figure=fig, subplot=[0.80, 0.50, hw, vw], north=True, auto_refresh=True)
m1p5.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p5.recenter(RA, DEC, radius=rad)
m1p5.set_tick_size(0)
m1p5.add_label(v_x, v_down, str(v_appr[5])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p5.add_label(v_x, v_up, str(v_reced[5])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

#### bottomw row

# plot
m1p6 = aplpy.FITSFigure(img, figure=fig, subplot=[0.00, 0.00, hw, vw], north=True, auto_refresh=True)
m1p6.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p6.recenter(RA, DEC, radius=rad)
m1p6.set_tick_size(0)
m1p6.add_label(v_x, v_down, str(v_appr[6])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p6.add_label(v_x, v_up, str(v_reced[6])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p7 = aplpy.FITSFigure(img, figure=fig, subplot=[0.20, 0.00, hw, vw], north=True, auto_refresh=True)
m1p7.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p7.recenter(RA, DEC, radius=rad)
m1p7.set_tick_size(0)
m1p7.add_label(v_x, v_down, str(v_appr[7])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p7.add_label(v_x, v_up, str(v_reced[7])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p8 = aplpy.FITSFigure(img, figure=fig, subplot=[0.40, 0.00, hw, vw], north=True, auto_refresh=True)
m1p8.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p8.recenter(RA, DEC, radius=rad)
m1p8.set_tick_size(0)
m1p8.add_label(v_x, v_down, str(v_appr[8])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p8.add_label(v_x, v_up, str(v_reced[8])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p9 = aplpy.FITSFigure(img, figure=fig, subplot=[0.60, 0.00, hw, vw], north=True, auto_refresh=True)
m1p9.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p9.recenter(RA, DEC, radius=rad)
m1p9.set_tick_size(0)
m1p9.add_label(v_x, v_down, str(v_appr[9])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p9.add_label(v_x, v_up, str(v_reced[9])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p0 = aplpy.FITSFigure(img, figure=fig, subplot=[0.80, 0.00, hw, vw], north=True, auto_refresh=True)
m1p0.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p0.recenter(RA, DEC, radius=rad)
m1p0.set_tick_size(0)
m1p0.add_label(v_x, v_down, str(v_appr[10])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p0.add_label(v_x, v_up, str(v_reced[10])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

# contours
m1p1.show_contour(levels=[v_appr[1],v_appr[0],v_reced[1]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p2.show_contour(levels=[v_appr[2],v_appr[0],v_reced[2]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p3.show_contour(levels=[v_appr[3],v_appr[0],v_reced[3]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p4.show_contour(levels=[v_appr[4],v_appr[0],v_reced[4]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p5.show_contour(levels=[v_appr[5],v_appr[0],v_reced[5]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p6.show_contour(levels=[v_appr[6],v_appr[0],v_reced[6]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p7.show_contour(levels=[v_appr[7],v_appr[0],v_reced[7]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p8.show_contour(levels=[v_appr[8],v_appr[0],v_reced[8]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p9.show_contour(levels=[v_appr[9],v_appr[0],v_reced[9]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p0.show_contour(levels=[v_appr[10],v_appr[0],v_reced[10]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)

# beam
m1p1.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p2.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p3.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p4.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p5.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p6.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p7.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p8.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p9.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p0.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')


# format of the labels and ticks 
m1p1.set_tick_labels_font(size=14)
m1p1.axis_labels.set_font(size=16)
#m1p1.set_tick_size(0)
m1p1.tick_labels.set_xformat('hh:mm:ss')
m1p1.tick_labels.set_yformat('dd:mm')

m1p6.set_tick_labels_font(size=14)
m1p6.axis_labels.set_font(size=16)
#m1p6.set_tick_size(0)
m1p6.tick_labels.set_xformat('hh:mm:ss')
m1p6.tick_labels.set_yformat('dd:mm')

m1p7.set_tick_labels_font(size=14)
m1p7.axis_labels.set_font(size=16)
#m1p7.set_tick_size(0)
m1p7.tick_labels.set_xformat('hh:mm:ss')
m1p7.tick_labels.set_yformat('dd:mm')

m1p8.set_tick_labels_font(size=14)
m1p8.axis_labels.set_font(size=16)
#m1p8.set_tick_size(0)
m1p8.tick_labels.set_xformat('hh:mm:ss')
m1p8.tick_labels.set_yformat('dd:mm')

m1p9.set_tick_labels_font(size=14)
m1p9.axis_labels.set_font(size=16)
#m1p9.set_tick_size(0)
m1p9.tick_labels.set_xformat('hh:mm:ss')
m1p9.tick_labels.set_yformat('dd:mm')

m1p0.set_tick_labels_font(size=14)
m1p0.axis_labels.set_font(size=16)
#m1p0.set_tick_size(0)
m1p0.tick_labels.set_xformat('hh:mm:ss')
m1p0.tick_labels.set_yformat('dd:mm')

# hide some axes and ticks
m1p1.hide_xaxis_label()
m1p1.tick_labels.hide_x()
m1p2.hide_yaxis_label()
m1p2.hide_xaxis_label()
m1p2.tick_labels.hide_x()
m1p2.tick_labels.hide_y()
m1p3.hide_yaxis_label()
m1p3.hide_xaxis_label()
m1p3.tick_labels.hide_x()
m1p3.tick_labels.hide_y()
m1p4.hide_yaxis_label()
m1p4.hide_xaxis_label()
m1p4.tick_labels.hide_x()
m1p4.tick_labels.hide_y()
m1p5.hide_yaxis_label()
m1p5.hide_xaxis_label()
m1p5.tick_labels.hide_x()
m1p5.tick_labels.hide_y()
m1p7.hide_yaxis_label()
m1p7.tick_labels.hide_y()
m1p8.hide_yaxis_label()
m1p8.tick_labels.hide_y()
m1p9.hide_yaxis_label()
m1p9.tick_labels.hide_y()
m1p0.hide_yaxis_label()
m1p0.tick_labels.hide_y()

# save figure
plt.savefig('/home/prm/Desktop/work/papers/paper_cig96/images/cig96_pairs_velocities.png', bbox_inches='tight')

ping()


# In[146]:

plt.close()


# In[143]:

# Select working directory
go_to('cig')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.16

# figure creation 5 in top row + 5 in bottom row
fig = plt.figure(figsize=(8,20))

# % widths
hw = 0.45
vw = 0.20

# image
img = "cig96_HI_mom1_5rms.fits"
#img = "cig96.w.circ.HI.Filt.image.3.5sigma.mom1.fits"

#### top row

# plot: background + contours + zoom
velmin = 1384
velmax = 1734
lev = np.linspace(velmin, velmax, num=(velmax-velmin)/10+1)

v_x = 0.81
v_up = 0.94
v_down = 0.89

length_tx = 8
length_ty = 4
color_blue = 'royalblue'
color_red  = 'crimson'
color_vcen = 'k'
labelsize = 16
backcolor = 'YlGn'
vcen = 1544 # km/s

# pairs
v_appr  = [1544, 1534, 1524, 1514, 1504, 1494, 1484, 1474, 1464, 1454, 1444]
v_reced = [1544, 1554, 1564, 1574, 1584, 1594, 1604, 1614, 1624, 1634, 1644]

m1p1 = aplpy.FITSFigure(img, figure=fig, subplot=[0.00, 0.80, hw, vw], north=True, auto_refresh=True)
m1p1.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p1.recenter(RA, DEC, radius=rad)
m1p1.set_tick_size(0)
m1p1.add_label(v_x, v_down, str(v_appr[1])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p1.add_label(v_x, v_up, str(v_reced[1])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')


m1p2 = aplpy.FITSFigure(img, figure=fig, subplot=[0.00, 0.60, hw, vw], north=True, auto_refresh=True)
m1p2.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p2.recenter(RA, DEC, radius=rad)
m1p2.set_tick_size(0)
m1p2.add_label(v_x, v_down, str(v_appr[2])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p2.add_label(v_x, v_up, str(v_reced[2])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p3 = aplpy.FITSFigure(img, figure=fig, subplot=[0.00, 0.40, hw, vw], north=True, auto_refresh=True)
m1p3.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p3.recenter(RA, DEC, radius=rad)
m1p3.set_tick_size(0)
m1p3.add_label(v_x, v_down, str(v_appr[3])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p3.add_label(v_x, v_up, str(v_reced[3])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p4 = aplpy.FITSFigure(img, figure=fig, subplot=[0.00, 0.20, hw, vw], north=True, auto_refresh=True)
m1p4.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p4.recenter(RA, DEC, radius=rad)
m1p4.set_tick_size(0)
m1p4.add_label(v_x, v_down, str(v_appr[4])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p4.add_label(v_x, v_up, str(v_reced[4])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p5 = aplpy.FITSFigure(img, figure=fig, subplot=[0.00, 0.00, hw, vw], north=True, auto_refresh=True)
m1p5.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p5.recenter(RA, DEC, radius=rad)
m1p5.set_tick_size(0)
m1p5.add_label(v_x, v_down, str(v_appr[5])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p5.add_label(v_x, v_up, str(v_reced[5])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

#### bottomw row

# plot
m1p6 = aplpy.FITSFigure(img, figure=fig, subplot=[0.50, 0.80, hw, vw], north=True, auto_refresh=True)
m1p6.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p6.recenter(RA, DEC, radius=rad)
m1p6.set_tick_size(0)
m1p6.add_label(v_x, v_down, str(v_appr[6])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p6.add_label(v_x, v_up, str(v_reced[6])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p7 = aplpy.FITSFigure(img, figure=fig, subplot=[0.50, 0.60, hw, vw], north=True, auto_refresh=True)
m1p7.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p7.recenter(RA, DEC, radius=rad)
m1p7.set_tick_size(0)
m1p7.add_label(v_x, v_down, str(v_appr[7])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p7.add_label(v_x, v_up, str(v_reced[7])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p8 = aplpy.FITSFigure(img, figure=fig, subplot=[0.50, 0.40, hw, vw], north=True, auto_refresh=True)
m1p8.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p8.recenter(RA, DEC, radius=rad)
m1p8.set_tick_size(0)
m1p8.add_label(v_x, v_down, str(v_appr[8])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p8.add_label(v_x, v_up, str(v_reced[8])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p9 = aplpy.FITSFigure(img, figure=fig, subplot=[0.50, 0.20, hw, vw], north=True, auto_refresh=True)
m1p9.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p9.recenter(RA, DEC, radius=rad)
m1p9.set_tick_size(0)
m1p9.add_label(v_x, v_down, str(v_appr[9])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p9.add_label(v_x, v_up, str(v_reced[9])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

m1p0 = aplpy.FITSFigure(img, figure=fig, subplot=[0.50, 0.00, hw, vw], north=True, auto_refresh=True)
m1p0.show_contour(levels=lev, cmap=backcolor, linewidth=1, zorder=1)
m1p0.recenter(RA, DEC, radius=rad)
m1p0.set_tick_size(0)
m1p0.add_label(v_x, v_down, str(v_appr[10])+' km/s', relative=True, color=color_blue, size=labelsize)
m1p0.add_label(v_x, v_up, str(v_reced[10])+' km/s', relative=True, color=color_red, size=labelsize)
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=length_tx, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=length_ty, width=1.5, color='black')

# contours
m1p1.show_contour(levels=[v_appr[1],v_appr[0],v_reced[1]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p2.show_contour(levels=[v_appr[2],v_appr[0],v_reced[2]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p3.show_contour(levels=[v_appr[3],v_appr[0],v_reced[3]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p4.show_contour(levels=[v_appr[4],v_appr[0],v_reced[4]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p5.show_contour(levels=[v_appr[5],v_appr[0],v_reced[5]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p6.show_contour(levels=[v_appr[6],v_appr[0],v_reced[6]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p7.show_contour(levels=[v_appr[7],v_appr[0],v_reced[7]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p8.show_contour(levels=[v_appr[8],v_appr[0],v_reced[8]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p9.show_contour(levels=[v_appr[9],v_appr[0],v_reced[9]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)
m1p0.show_contour(levels=[v_appr[10],v_appr[0],v_reced[10]], colors=[color_blue, color_vcen, color_red], linewidths=3, zorder=2)

# beam
m1p1.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p2.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p3.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p4.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p5.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p6.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p7.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p8.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p9.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')
m1p0.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor=(0.2,0.5,0.5,0.4), edgecolor='k')

# format of the labels and ticks 
m1p1.set_tick_labels_font(size=14)
m1p1.axis_labels.set_font(size=16)
m1p1.tick_labels.set_xformat('hh:mm:ss')
m1p1.tick_labels.set_yformat('dd:mm')
m1p2.set_tick_labels_font(size=14)
m1p2.axis_labels.set_font(size=16)
m1p2.tick_labels.set_xformat('hh:mm:ss')
m1p2.tick_labels.set_yformat('dd:mm')
m1p3.set_tick_labels_font(size=14)
m1p3.axis_labels.set_font(size=16)
m1p3.tick_labels.set_xformat('hh:mm:ss')
m1p3.tick_labels.set_yformat('dd:mm')
m1p4.set_tick_labels_font(size=14)
m1p4.axis_labels.set_font(size=16)
m1p4.tick_labels.set_xformat('hh:mm:ss')
m1p4.tick_labels.set_yformat('dd:mm')
m1p5.set_tick_labels_font(size=14)
m1p5.axis_labels.set_font(size=16)
m1p5.tick_labels.set_xformat('hh:mm:ss')
m1p5.tick_labels.set_yformat('dd:mm')
m1p6.set_tick_labels_font(size=14)
m1p6.axis_labels.set_font(size=16)
m1p6.axis_labels.set_yposition('right')
m1p6.tick_labels.set_yposition('right')
m1p6.tick_labels.set_yformat('dd:mm')
m1p7.set_tick_labels_font(size=14)
m1p7.axis_labels.set_font(size=16)
m1p7.axis_labels.set_yposition('right')
m1p7.tick_labels.set_yposition('right')
m1p7.tick_labels.set_yformat('dd:mm')
m1p8.set_tick_labels_font(size=14)
m1p8.axis_labels.set_font(size=16)
m1p8.axis_labels.set_yposition('right')
m1p8.tick_labels.set_yposition('right')
m1p8.tick_labels.set_yformat('dd:mm')
m1p9.set_tick_labels_font(size=14)
m1p9.axis_labels.set_font(size=16)
m1p9.axis_labels.set_yposition('right')
m1p9.tick_labels.set_yposition('right')
m1p9.tick_labels.set_yformat('dd:mm')
m1p0.set_tick_labels_font(size=14)
m1p0.axis_labels.set_font(size=16)
m1p0.axis_labels.set_yposition('right')
m1p0.tick_labels.set_yposition('right')
m1p0.tick_labels.set_xformat('hh:mm:ss')
m1p0.tick_labels.set_yformat('dd:mm')

# hide some axes and ticks
m1p1.hide_xaxis_label()
m1p1.tick_labels.hide_x()
m1p2.hide_xaxis_label()
m1p2.tick_labels.hide_x()
m1p3.hide_xaxis_label()
m1p3.tick_labels.hide_x()
m1p4.hide_xaxis_label()
m1p4.tick_labels.hide_x()
m1p6.hide_xaxis_label()
m1p6.tick_labels.hide_x()
m1p7.hide_xaxis_label()
m1p7.tick_labels.hide_x()
m1p8.hide_xaxis_label()
m1p8.tick_labels.hide_x()
m1p9.hide_xaxis_label()
m1p9.tick_labels.hide_x()

# save figure
plt.savefig('/home/prm/Desktop/work/papers/paper_cig96/images/cig96_pairs_velocities_vert.png', bbox_inches='tight')

ping()


# In[144]:

plt.close()


# ##1.3 Karma position-velocity cuts
# 
# First we change the frequency (Hz) to radio velocity (km/s).

# In[259]:

# Select working directory
go_to('cig')

# select images
images = ["cig96_PAmajor20.39.fits", "cig96_PAminor290.39.fits"]

# change Y axis values to velocity: radvel = c*(fHI - datamajor)/fHI
c = 299792.458  #km/s
fHI = 1420405751.786 # Hz

for image in images:
    data = fits.open(image)[0].data
    head = fits.open(image)[0].header
    head['CRVAL2'] = 1410.0
    head['CTYPE2'] = "VELOCITY"
    head['CRPIX2'] = 10.0
    head['CDELT2'] = 10.0
    head['CRVAL1'] = 0.0
    if image[0]:
        head['CDELT1'] = 0.06666493
    if image[1]:
        head['CDELT1'] = 0.066661538
    fits.writeto(image, data, header=head, clobber=True)


# Second, now we plot the p/V cuts.

# In[4]:

# Select working directory
go_to('cig')

# select images
major = "cig96_PAmajor20.39.fits"
minor = "cig96_PAminor290.39.fits"

# figure structure
fig = plt.figure(figsize=(9, 8))

# plot window design
vw = 0.40 # vertical width of the subplots in %
hw = 0.85 # horizontal width of the subplots in % 

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.20   # zoom radius so, the larger rad, the larger FoV in the image

# upper panel
cutmaj = aplpy.FITSFigure(major, figure=fig, subplot=[0.1, 0.5, hw, vw])
cutmaj.set_axis_labels(xlabel="Offset (arcmin)", ylabel="Velocity (km/s)")
cutmaj.set_axis_labels_size(18)
cutmaj.set_yaxis_coord_type('scalar')
cutmaj.ticks.set_length(length=0)
cutmaj.ticks.set_xspacing(spacing=1)
cutmaj.ticks.set_yspacing(spacing=50)
cutmaj.set_tick_labels_font(size=15)
plt.tick_params(axis='both', which='both', direction='in', right='on', top='on')
plt.grid(which=u'major', alpha=0.4)
# beam
#plt.plot(14, 5, 'yo', markersize=14, markeredgecolor='k', zorder=4)
# SW warp
plt.annotate('', xy=(30, 20), xytext=(30, 25), arrowprops=dict(facecolor='cyan'), zorder=4)
# SW and NE
plt.annotate(s='SW', color='white', fontsize=16, xy=(10, 15), xytext=(10, 20), annotation_clip=False)
plt.annotate(s='NE', color='white', fontsize=16, xy=(230, 30), xytext=(230, 33), annotation_clip=False)

# lower panel
cutmin = aplpy.FITSFigure(minor, figure=fig, subplot=[0.1, 0.1, hw, vw])
cutmin.set_axis_labels(xlabel="Offset (arcmin)", ylabel="Velocity (km/s)")
cutmin.set_axis_labels_size(18)
cutmin.set_yaxis_coord_type('scalar')
cutmin.ticks.set_length(length=0)
cutmin.set_axis_labels_size(18)
cutmin.ticks.set_xspacing(spacing=1)
cutmin.ticks.set_yspacing(spacing=50)
cutmin.set_tick_labels_font(size=15)
plt.tick_params(axis='both', which='both', direction='in', right='on', top='on')
plt.grid(which=u'major', alpha=0.4)

# beam
#plt.plot(14, 5, 'yo', markersize=14, markeredgecolor='k', zorder=4)

# colors
backgr = "PRGn"
cutmaj.show_colorscale(stretch='linear', cmap=backgr, aspect="auto")
cutmaj.set_nan_color('black')
cutmin.show_colorscale(stretch='linear', cmap=backgr, aspect="auto")
cutmin.set_nan_color('black')

# contours
rms = 0.000126
cont = "YlGnBu"
cutmaj.show_contour(levels=[3.5*rms, 5*rms], cmap=cont)
cutmin.show_contour(levels=[3.5*rms, 5*rms], cmap=cont)

# PA
posx, posy = [195, 43]
plt.annotate(s='PA = 20$\degree$',  color='white', fontsize=16, xy=(0.08, 45), xytext=(posx+25, posy+48), annotation_clip=False)
plt.annotate(s='PA = 110$\degree$', color='white', fontsize=16, xy=(0.08, 45), xytext=(posx+20, posy), annotation_clip=False)

# HI features
plt.annotate('', xy=(220, 25), xytext=(220, 30), arrowprops=dict(facecolor='red'), zorder=4)

#legend
#plt.legend(loc=0, facecolor='white')

# save figure
plt.savefig('/home/prm/Desktop/work/papers/paper_cig96/images/cig96_pvcuts.png', bbox_inches='tight')


# In[5]:

plt.close()


# #2. Channels map
# 
# The datacube ingested ***has*** to have only 3 axes: RA, DEC and Frequency. If Stokes axis is present, the function will not work. Stokes axis can be removed via Gipsy. Consult <a href="http://amiga.iaa.es:9999/ipython/notebooks/pique/Cubes/GipsyRecipe/RemoveAxis.ipynb">"Remove Axis" notebook by Pique</a> here for more information.
# 
# <a hreaf="https://www.astro.rug.nl/software/kapteyn/maputils.html">Kapteyn package</a> was used in the function to build the channels map.
# 
# ###a) we define the cube and images:

# In[12]:

# renders interactive figures in the Notebook
get_ipython().magic(u'matplotlib inline')


# In[13]:

# Select working directory
go_to('cig')

# define fits files
#optim_deep = "cig96_def_crop.fits"
optim_deep = "cig96_VST_coadd_12.5arcmin_cropped_SB.fits"
optim_RdRGB = "cig96_RdRGB_SB_oriented_cropped.fits"
cube = "cig96.3D.fits"
mom0 = "cig96_HI_mom0_3.5rms.fits"
mom1 = "cig96_HI_mom1_3.5rms.fits"

# load them
optim_deep_object = maputils.FITSimage(optim_deep)
optim_RdRGB_object = maputils.FITSimage(optim_RdRGB)
cube_object = maputils.FITSimage(cube)
m0_object = maputils.FITSimage(mom0)
m1_object = maputils.FITSimage(mom1)

# Get axes from Cube headers
specaxnum = cube_object.proj.specaxnum
lonaxnum = cube_object.proj.lonaxnum
lataxnum = cube_object.proj.lataxnum


# ###b) we define the channels map function

# In[9]:

# Function to build channel map with planes as background 
'''
    This function makes a channels map and is slightly different from Pique's original function.
    The new arguments introduced are:
        x0, x1, y0, y1 = integers defining the are of the image to be zoomed, in pixels
'''
def channelMap_Bck3D(cols, x0, x1, y0, y1, channels, cube_object, ReprImgobject, contourlevels, mapacol, mapcont, cmin, cmax, Bck3D=True):  

    # Mosaic plot settings
    rows = len(channels) / cols
    if rows*cols < len(channels): rows += 1
    w = cols*5
    h = rows*5
    
    # Zoom in the interesting regions (in pixels) defined in the function
    xlimit = (x0, x1)
    ylimit = (y0, y1)

    # Make Channel Map
    fig = plt.figure(figsize=(w,h))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    
    # Beam position (WCS in deg)
    pos = '33.98 deg, 5.872 deg'
    pos_optical = '34.03 deg, 5.83 deg'

    for i, ch in enumerate(channels):
        frame = fig.add_subplot(rows, cols, i+1)
    
        # Slice channels    
        cube_object.set_imageaxes(lonaxnum, lataxnum, slicepos=ch)
        #cube_object.set_limits(pxlim=xlimit, pylim=ylimit)
    
        if Bck3D:
            # Calculate flux range
            fluxmin, fluxmax = cube_object.get_dataminmax()
            slicedata = cube_object.dat[ch,:,:]
    
            background = []
            background.append(np.std(slicedata[0:30,0:30]))
            background.append(np.std(slicedata[0:30,570:600]))   
            background.append(np.std(slicedata[570:600,0:30]))   
            background.append(np.std(slicedata[570:600,570:600]))
            rms = np.mean(background)
            fluxmin=2.5*rms
            fluxmax=0.95*fluxmax
    
            # Generate images
            bckimage = cube_object.Annotatedimage(frame, cmap=mapacol, clipmin=fluxmin, clipmax=fluxmax, blankcolor='None')
            if (ch == 22):                              # sets the RA, DEC in the bottom-left 22nd channel map
                grat = bckimage.Graticule()             # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_gratline(visible=False)       # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_axislabel(fontsize=14)        # sets the RA, DEC in the bottom-left 22nd channel map
            bckimage.Image()
            bckimage.Beam(28, 28, units='arcsec', pa=0.0, pos=pos_optical, fc='y', linestyle="-", color='k', linewidth=1, fill=True)
                        
            # Contour levels overlay
            overlay = ReprImgobject.Annotatedimage(frame)
            overlay.Contours(levels=contourlevels, colors='red', linewidth=0.5, alpha=0.3)
            
        else:
            bckimage = ReprImgobject.Annotatedimage(frame, cmap=mapacol, clipmin=cmin, clipmax=cmax, blankcolor='None')
            corner_channel = 22       # corner channels are 22 and 34
            if (ch == 6) or (ch == 10) or (ch == 14) or (ch == 18) or (ch == 22) or (ch == 26) or (ch == 30) or (ch == 34):
                # sets the RA, DEC in the bottom-left 22nd channel map
                grat = bckimage.Graticule()             # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_gratline(wcsaxis=(0,1), color='silver', visible=False)       # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_axislabel(fontsize=30)        # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_ticklabel(fontsize=17)        # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_ticklabel(plotaxis="left", rotation=90)
                grat.setp_tickmark(markersize=12, linewidth=2)
            if (ch == 23) or (ch == 24) or (ch == 25) or (ch == 35) or (ch == 36) or (ch == 37):
                grat = bckimage.Graticule()             # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_gratline(visible=False)       # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_axislabel(visible=False) 
                grat.setp_axislabel(fontsize=30)        # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_ticklabel(fontsize=17)        # sets the RA, DEC in the bottom-left 22nd channel map
                grat.setp_tickmark(markersize=12, linewidth=2)
            if (ch != 6) and (ch != 10) and (ch != 14) and (ch != 18) and (ch != 22) and (ch != 26) and (ch != 30) and (ch != 34) and (ch != 35) and (ch != 36) and (ch != 37):
                grat = bckimage.Graticule()             # sets the grid in all channels
                grat.setp_ticklabel(visible=False)      # hides axis ticks
                grat.setp_axislabel(visible=False)      # hides axis labels
                grat.setp_gratline(wcsaxis=(0,1), color='silver', visible=False)
                grat.setp_tickmark(markersize=12, linewidth=2)
            bckimage.Image()            
            bckimage.Beam(28, 28, units='arcsec', pa=0.0, pos=pos_optical, fc='y', linestyle="-", color='k', linewidth=1, fill=True)
                        
            # Contour levels overlay
            overlay = cube_object.Annotatedimage(frame)
            overlay.Contours(levels=contourlevels, cmap=mapcont, linewidth=0.5, alpha=0.7, transparent=True)
    
        # Label velocities upper right corner
        spec = cube_object.proj.sub(specaxnum).spectra("VRAD-???")
        vel = spec.toworld1d(ch) # Calculates velocities
        velinfo = "Ch%d = %.1f km/s" % (ch, vel/1000)    
        frame.text(0.98, 0.98, velinfo, # upper corner
               horizontalalignment='right',
               verticalalignment='top',
               transform = frame.transAxes,
               fontsize=24, color='white',
               bbox=dict(facecolor='black', alpha=0.8), zorder=4)   
        
        # Display
        bckimage.plot()
        overlay.plot()
        
        # Zoom in the designed area
        plt.xlim(xlimit)
        plt.ylim(ylimit)
# DOES NOT WORK: the transparent frame is gone: plt.savefig("/home/prm/Desktop/work/papers/paper_cig96/images/cig96_chmap"+str(start)+"_"+str(end)+"_marks.png", transparent=True)


# ###c) Moment 0 (contours) over cube:

                # Reprojection of 2D WCS into 3D WCS
ReprImgobject = m0_object.reproject_to(cube_object.hdr)

# Max min data values
ReprImgobject.get_dataminmax()
                
                # Get a range of channels in the data cube
start = 21; end = 23; step = 1 # plot-OK channels: start = 6; end = 38; step = 1
channels = range(start, end, step)

# Contour levels
contourlevels = np.linspace(0.015, 1.0, num=5)

# Channels map plot
channelMap_Bck3D(4, 135, 460, 135, 460, channels, cube_object, ReprImgobject, contourlevels, 'spectral_r', 'bone', 1,1)
                
                plt.close()
                
# ###d) optical image (countours) over cube:

# In[14]:

# Reprojection of 2D WCS into 3D WCS  
#ReprImgobject = optim_RdRGB_object.reproject_to(cube_object.hdr)   # CAHA image
ReprImgobject = optim_deep_object.reproject_to(cube_object.hdr)    # VST image

# Display Optical Image
#mplim = ReprImgobject.Annotatedimage(cmap="jet")
#mplim.Beam(28, 28, units='arcsec', pa=0.0, pos="34.1 deg, 5.762 deg", fc='yellow', fill=True)
#mplim.Image()
#mplim.Graticule()
#mplim.plot()

# Max min data values
ReprImgobject.get_dataminmax()


                # HI cube over optical contour

# Get a range of channels in the data cube
start = 26; end = 28; step = 1 # plot-OK channels: start = 6; end = 38; step = 1
channels = range(start, end, step)

# Contour levels
contourlevels = np.linspace(24,26, num=3) # optimal for pseudoring retrieval

# Channels map plot
#channelMap_Bck3D(4, 170, 410, 170, 410, channels, cube_object, ReprImgobject, contourlevels, 'spectral_r', 'bone_r', 1, 1)
channelMap_Bck3D(4, 0, 600, 0, 600, channels, cube_object, ReprImgobject, contourlevels, 'hsv_r', 'bone', 1, 1)
                
# plt.close()

# cdict = {
#   'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
#   'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
#   'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
# }
# 
# cm = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

# In[39]:

# Set a range of channels in the data cube (optimal channels: start = 6; end = 38; step = 1)
start = 6; end = 25; step = 1 # first part of the plot
#start = 26; end = 37; step = 1 # second part of the plot
channels = range(start, end+1, step)

# Contour levels
#in Jy/beam/km*s
rms = 0.000126    
#contourlevels = [3.94*rms, 16.895*rms, 35*rms]#, 15*rms, 30*rms]
#contourlevels = [0.000441, 0.000630, 0.001134, 0.002, 0.005, 0.01]

# in N(HI)
nhi = [6e+18, 7e+18, 8e+18, 9e+18, 1e+19, 5e+19, 1e+20, 2e+20, 4e+20]
contourlevels = [elem/(1.106561*10**24*10/28/28) for elem in nhi]
rmslev = [level/rms for level in contourlevels]

print contourlevels
print nhi
print rmslev

# zoom
down = 130
up = 600 - down

# Channels map plot
channelMap_Bck3D(4, down, up, down, up, channels, cube_object, ReprImgobject, contourlevels, 'Blues_r', 'Reds_r', 26, 28.4, Bck3D=False)


# In[38]:

# Set a range of channels in the data cube (optimal channels: start = 6; end = 38; step = 1)
#start = 6; end = 25; step = 1 # first part of the plot
start = 26; end = 37; step = 1 # second part of the plot
channels = range(start, end+1, step)

# Contour levels
#in Jy/beam/km*s
rms = 0.000126    
#contourlevels = [3.94*rms, 16.895*rms, 35*rms]#, 15*rms, 30*rms]
#contourlevels = [0.000441, 0.000630, 0.001134, 0.002, 0.005, 0.01]

# in N(HI)
nhi = [6e+18, 7e+18, 8e+18, 9e+18, 1e+19, 5e+19, 1e+20, 2e+20, 4e+20]
contourlevels = [elem/(1.106561*10**24*10/28/28) for elem in nhi]
rmslev = [level/rms for level in contourlevels]

print contourlevels
print nhi
print rmslev

# zoom
down = 130
up = 600 - down

# Channels map plot
channelMap_Bck3D(4, down, up, down, up, channels, cube_object, ReprImgobject, contourlevels, 'Blues_r', 'Reds_r', 26, 28.4, Bck3D=False)


# In[15]:

# Function to build channel map with planes as background 
'''
    This function makes a channels map and is slightly different from Pique's original function.
    The new arguments introduced are:
        x0, x1, y0, y1 = integers defining the are of the image to be zoomed, in pixels
'''
def channelMap_Bck3D(cols, x0, x1, y0, y1, start, end, step, channels, cube_object, ReprImgobject, contourlevels, mapacol, mapcont, mapcont2, cmin, cmax, Bck3D=True):  

    # Mosaic plot settings
    rows = len(channels) / cols
    if rows*cols < len(channels): rows += 1
    w = cols*5
    h = rows*5
    
    # Zoom in the interesting regions (in pixels) defined in the function
    xlimit = (x0, x1)
    ylimit = (y0, y1)

    # Make Channel Map
    fig = plt.figure(figsize=(w,h))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    
    # Beam position (WCS in deg)
    pos = '33.98 deg, 5.872 deg'
    pos_optical = '34.03 deg, 5.83 deg'

    for i, ch in enumerate(channels):
        frame  = fig.add_subplot(rows, cols, i+1)
        frame2 = fig.add_subplot(rows, cols, end-start+1-i)
        # Slice channels    
        cube_object.set_imageaxes(lonaxnum, lataxnum, slicepos=ch)
        #cube_object.set_limits(pxlim=xlimit, pylim=ylimit)
    
        if Bck3D:
            # Calculate flux range
            fluxmin, fluxmax = cube_object.get_dataminmax()
            slicedata = cube_object.dat[ch,:,:]
    
            background = []
            background.append(np.std(slicedata[0:30,0:30]))
            background.append(np.std(slicedata[0:30,570:600]))   
            background.append(np.std(slicedata[570:600,0:30]))   
            background.append(np.std(slicedata[570:600,570:600]))
            rms = np.mean(background)
            fluxmin=2.5*rms
            fluxmax=0.95*fluxmax
    
            # Generate images
            bckimage = cube_object.Annotatedimage(frame, cmap=mapacol, clipmin=fluxmin, clipmax=fluxmax, blankcolor='None')
            #if (ch == 22):                              # sets the RA, DEC in the bottom-left 22nd channel map
            #    grat = bckimage.Graticule()             # sets the RA, DEC in the bottom-left 22nd channel map
            #    grat.setp_gratline(visible=False)       # sets the RA, DEC in the bottom-left 22nd channel map
            #    grat.setp_axislabel(fontsize=14)        # sets the RA, DEC in the bottom-left 22nd channel map
            bckimage.Image()
            bckimage.Beam(28, 28, units='arcsec', pa=0.0, pos=pos_optical, fc='y', linestyle="-", color='k', linewidth=1, fill=True)
                        
            # Contour levels overlay
            overlay = ReprImgobject.Annotatedimage(frame)
            overlay.Contours(levels=contourlevels, colors='red', linewidth=0.5, alpha=0.3)
            
        else:
            bckimage = ReprImgobject.Annotatedimage(frame, cmap=mapacol, clipmin=cmin, clipmax=cmax, blankcolor='None')
            corner_channel = 35       # corner channel is frame 34
            RAframes  = [36, 37, 38]  # frames that show RA
            DECframes = [6, 11, 15, 19, 23, 27, 31, 384]  # frames that sow DEC
            if ch == corner_channel:
                grat = bckimage.Graticule()
                grat.setp_gratline(visible=False)
                grat.setp_axislabel("left", fontsize=30)
                grat.setp_ticklabel(plotaxis="left", fontsize=17, rotation=90)
                grat.setp_axislabel(plotaxis="bottom", fontsize=30)
                grat.setp_ticklabel(plotaxis="bottom", fontsize=17, rotation=0)    
            if ch in RAframes: # RA visible
                grat = bckimage.Graticule()
                grat.setp_gratline(visible=False)
                grat.setp_ticklabel(plotaxis="bottom", fontsize=17, rotation=0)
                grat.setp_axislabel(plotaxis="bottom", fontsize=30)
                grat.setp_axislabel(plotaxis="left", visible=False) # removed
                grat.setp_ticklabel(plotaxis="left", visible=False) # removed
                grat.setp_tickmark(plotaxis="left", visible=False)  # removed
            if ch in DECframes: # DEC visible
                grat = bckimage.Graticule()
                grat.setp_gratline(visible=False)
                grat.setp_tickmark(markersize=12, linewidth=2)
                grat.setp_axislabel("left", fontsize=30)
                grat.setp_ticklabel(plotaxis="left", fontsize=17, rotation=90)
                grat.setp_axislabel(plotaxis="bottom", visible=False) # removed
                grat.setp_ticklabel(plotaxis="bottom", visible=False) # removed
                grat.setp_tickmark(plotaxis="bottom", visible=False)  # removed
            if ch not in RAframes or DECframes:
                grat = bckimage.Graticule()
                grat.setp_gratline(visible=False)
                grat.setp_axislabel(visible=False)
                grat.setp_ticklabel(visible=False)
            bckimage.Image()            
            bckimage.Beam(28, 28, units='arcsec', pa=0.0, pos=pos_optical, fc='cyan', linestyle="-", color='k', linewidth=1, fill=True)
                        
            # Contour levels overlay
            overlay = cube_object.Annotatedimage(frame)
            overlay.Contours(levels=contourlevels, cmap=mapcont, linewidth=0.5, alpha=0.7, transparent=True)
            overlay2 = cube_object.Annotatedimage(frame2)
            overlay2.Contours(levels=contourlevels, cmap=mapcont2, linewidth=0.5, alpha=0.7, transparent=True)
            
        # Label velocities upper right corner
        spec = cube_object.proj.sub(specaxnum).spectra("VRAD-???")
        vel = spec.toworld1d(ch) # Calculates velocities
        vel2 = spec.toworld1d(start+end-ch) # Calculates velocities
        velinfo  = "Ch%d = %.0f km/s" % (ch, vel/1000)    
        velinfo2 = "Ch%d = %.0f km/s" % (start+end-ch, vel2/1000)  
        frame.text(0.48, 0.98, velinfo, # upper corner
               horizontalalignment='right',
               verticalalignment='top',
               transform = frame.transAxes,
               fontsize=18, color='red',
               bbox=dict(facecolor='whitesmoke', alpha=0.9), zorder=4)
        frame.text(0.98, 0.98, velinfo2, # upper corner
               horizontalalignment='right',
               verticalalignment='top',
               transform = frame.transAxes,
               fontsize=18, color='blue',
               bbox=dict(facecolor='whitesmoke', alpha=0.9), zorder=4)
        
        # Display and zoom in the designed area
        bckimage.plot()
        overlay.plot()
        plt.xlim(xlimit)
        plt.ylim(ylimit)
        overlay2.plot()
# DOES NOT WORK: the transparent frame is gone: plt.savefig("/home/prm/Desktop/work/papers/paper_cig96/images/cig96_chmap"+str(start)+"_"+str(end)+"_marks.png", transparent=True)


# In[16]:

# Set a range of channels in the data cube (optimal channels: start = 6; end = 38; step = 1)
#start = 6; end = 25; step = 1 # first part of the plot
start = 7; end = 38; step = 1 # second part of the plot
channels = range(start, end+1, step)

# Contour levels
#in Jy/beam/km*s
rms = 0.000126    
#contourlevels = [3.94*rms, 16.895*rms, 35*rms]#, 15*rms, 30*rms]
#contourlevels = [0.000441, 0.000630, 0.001134, 0.002, 0.005, 0.01]

# in N(HI)
nhi = [6e+18, 7e+18, 8e+18, 9e+18, 1e+19, 5e+19, 1e+20, 2e+20, 4e+20]
contourlevels = [elem/(1.106561*10**24*10/28/28) for elem in nhi]
rmslev = [level/rms for level in contourlevels]

print contourlevels
print nhi
print rmslev

# zoom
down = 130
up = 600 - down

# Channels map plot
channelMap_Bck3D(4, down, up, down, up, start, end, step, channels, cube_object, ReprImgobject, contourlevels, 'YlGn_r', 'Reds_r', 'Blues_r', 26, 28.4, Bck3D=False)


# In[17]:

plt.close()


# In[40]:

# renders interactive figures in the Notebook
get_ipython().magic(u'matplotlib nbagg')


# ## 3. Deprojection of 0th moment
# 
# Image selected: "cig96_HI_mom0_3.0rms.fits".

                IN IRAF:
    
geotran input=cig96_HI_mom0_3.0rms.fits output=cig96_HI_mom0_3.0rms_dep.fits xin=300.686 yin=300.887 xrot=69.61 yrot=69.61
                
# In[25]:

# Select working directory
go_to('cig')

# image
image = "cig96_HI_mom0_3.0rms_dep.fits"

# CIG96 coordinates in deg
RAdep  = 33.787073
DECdep = 6.1315624
raddep = 0.09   # zoom radius so, the larger rad, the larger FoV in the image

# deprojected 0th moment (3s)
HIdep = aplpy.FITSFigure(image)
HIdep.show_colorscale(cmap="seismic_r", stretch='linear', exponent=2, vmin=0, vmax=7000, vmid=0)
HIdep.show_regions("cig96_HIpseudo_dep_wcs.reg")
plt.annotate("North", xy=(150, 349), xytext=(160, 350), size=20, color='white', arrowprops=dict(facecolor='white'))
# zoom
HIdep.recenter(RAdep, DECdep, radius=raddep)


# In[58]:

HIdep.close()


# #4. HI 0th moment to column density
# 
# Further info in wiki page <a href="http://amiga.iaa.es:8888/display/science/Instrumental+characteristics+of+HI+observations">'Instrumental characteristics of HI observations'</a>.
# 
# - From flux to **Surface Brightness Temperature** $T_B$ in $K$:
# 
# $$T_B = 6.07 × 10^5 × S / (Maj × Min)$$
# 
# where $S$ is expressed in $Jy/beam$ and $Maj$, $Min$ are, respectively, the major and minor axis of the beam (FWHM) in $arcsec$.
# 
# - The **HI mass** $M_{HI}$ in solar masses $M_{\odot}$ is:
# 
# $$M_{HI} = 2.36 × 10^5 × D^2 \int Sdv$$
# 
# where the integral represents the flux for all channels (each channel with width $dv$), in $Jy · km · s^{−1}$ and $D$ is the distance in $Mpc$.
# 
# - The **column density** $N_{HI}$ is:
# 
# $$N_{HI} = 1.823 × 10^{18} \int T_B  dv$$
# 
# where the integral represents the brightness temperature for all channels (each channel with width $dv$), in $K · km · s^{−1}$.
# 
# 
# - **0th moment** is computed by integrating the flux along all channels:
# 
# $$M_0 = \int S dv$$
# 
# where $S$ is in $Jy/beam$ and $dv$ in $km/s$.
# 
# So, ***to convert a 0th moment map to column density*** (and considering that our beam is $28×28$ $arcsec^2$):
# 
# $$N_{HI} = 1.823 × 10^{18} \int T_B  dv = $$ 
# $$= 1.823 × 10^{18} \int (6.07 × 10^5 × S/(Maj × Min)) dv = $$
# $$ = (1.823 × 10^{18} × 6.07 × 10^5)/(Maj × Min) \int S dv = $$
# $$ = 1.106561 × 10^{24}/(Maj×Min) \int S dv  =$$
# 
# ##4.1 CIG96 EVLA data
# 
# In the case of the EVLA CIG96 data, $Maj×Min = 28 × 28 arcsec$ so:
# $$N_{HI} (cm^{-2})= 1.411429847 × 10^{21} × \int S dv$$

# In[67]:

# Select working directory
go_to('cig')

rms_sel = 3.5
moment0 = "cig96_HI_mom0_"+str(rms_sel)+"rms.fits"
im = fits.open(moment0, mode='update')
im[0].data = 1.411429847*10*(im[0].data)    # conversion to column density in 10e21 at/cm2
im[0].data[np.isnan(im[0].data)] = 0
rms0 = np.sqrt(np.mean(np.square(im[0].data[100:209,125:219])))
im.writeto("cig96_HI_mom0_"+str(rms_sel)+"rms_coldens.fits", clobber=True)

print "rms (mom0) =", rms0, "10e21 at cm-2"


# and:
# $$M_{HI} (M_{\odot}) = 97.25324 × 10^{6} × \int S dv$$

# In[133]:

rms_sel = 5.0
moment0 = "cig96_HI_mom0_"+str(rms_sel)+"rms.fits"
im = fits.open(moment0, mode='update')
im[0].data = 97.25324*10*(im[0].data)    # x 10e6 Msol
im[0].data[np.isnan(im[0].data)] = 0
im.writeto("cig96_HI_mom0_"+str(rms_sel)+"rms_Msol.fits", clobber=True)


# In[102]:

# plot 
# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.20   # zoom radius so, the larger rad, the larger FoV in the image

# plot image
mom0col = aplpy.FITSFigure("cig96_HI_mom0_"+str(rms_sel)+"rms_coldens.fits", north=True, auto_refresh=True)
#mom0col = aplpy.FITSFigure("cig96_HI_mom0_"+str(rms_sel)+"rms.fits", north=True, auto_refresh=True)
mom0col.show_colorscale(stretch='log', cmap='jet', vmin=0.1, vmax=15)
mom0col.set_nan_color('black')
mom0col.recenter(RA,DEC,rad)

# beam definition
mom0col.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, edgecolor='k', facecolor=(0.9,0.9,0.2,0.99))

# Contours

# Espada+05 VLA data
#espada = [0.5, 1.2, 2.4, 3.7, 4.9, 6.1, 7.3, 8.6, 9.8, 11.0, 12.2]    # Espada+ 05 contours

# This work
#mom0col.show_contour(levels=[3*rms0, 4*rms0, 5*rms0], cmap='jet_r')  # *10²¹ at/cm2
levs = [0.15]
# show contours
#mom0col.show_contour(levels=espada, cmap='Blues_r')     # *10²¹ at/cm2
mom0col.show_contour(levels=levs, cmap='hsv')     # *10²¹ at/cm2

#legend
#plt.plot([0,1], color='C0', label="Espada+ 2005: 10 x 10e20 at/cm-2")
#plt.plot([0,1], color='r', label="This work: 2.0, 3.4 x 10e20  at/cm-2")
#plt.legend(loc=0, facecolor='white')


# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/'+'cig96_HImom0_'+str(rms_sel)+'rms_coldens.png', bbox_inches='tight')


# In[103]:

mom0col.close()


# ##4.2 CIG96 VLA data
# 
# In the case of the VLA CIG96 data, $Maj×Min = 49.8 × 46.2 arcsec$ so:
# 
# $$N_{HI} = 1.106561 × 10^{24}/(49.8 × 46.2) \int S dv  =$$
# $$ = 4.809545541 × 10^{20} × \int S dv$$

                im = "CIG96.CBL1-SMO.MOM0.fits"
data, header = fits.getdata(im, header=True)
fits.getheader(im, 0)
header['CUNIT1'] = 'deg'
header['CUNIT2'] = 'deg'
header['CUNIT4'] = 'm/s'
header['RESTFRQ'] = 1420405751.77
header['CRPIX1'] = 301.0
header['CRPIX2'] = 301.0
header['CRVAL1'] = 33.865
header['CRVAL2'] = 6.0025
header['LONPOLE'] = 180.0
header['LATPOLE'] = 6.0025
header['RADESYS'] = 'FK5'
header['EQUINOX'] = 2000.0
header['WCSAXES'] = 2
header['CDELT1'] = -0.002777777845000
header['CDELT2'] = 0.002777777845000
header['EPOCH'] = ""
fits.writeto(im[:-5]+"_new.fits", data, header, clobber=True)
                
                # Select working directory
go_to('cig')

# My image
evla = "cig96_HI_mom0_3rms_coldens.fits"

# VLA image
#vla = "CIG96.CBL1-SMO.MOM0.m_s_new.fits"
vla = "CIG96.CBL1-SMO.MOM0_new.fits"
vlaim = fits.open(vla)
vlaim[0].data = 5.001927363*10.4*vlaim[0].data  # *10²¹ at/cm2
fits.writeto(vla[:-5]+"_coldens.fits", vlaim[0].data, vlaim[0].header, clobber=True)
vla_coldens = vla[:-5]+"_coldens.fits"

# plot 
# CIG96 coordinates in deg
RA = 33.865
DEC = 6.0025
rad = 0.005   # zoom radius so, the larger rad, the larger FoV in the image

#mom0vlaevla = aplpy.FITSFigure(vlaim[0].data)   # VLA data in column density
mom0vlaevla = aplpy.FITSFigure(evla)   # EVLA data in column density
mom0vlaevla.show_colorscale(stretch='linear', cmap='Greys')#, vmin=20, vmax=6000)  # *10²¹ at/cm2
#mom0vlaevla.set_nan_color('black')
#mom0vlaevla.recenter(RA,DEC,rad)

# beam definition
mom0vlaevla.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, edgecolor='k', facecolor=(0.9,0.9,0.2,0.99))

# contour levels definition
#lev = np.arange(0.015,2.0,0.05)
mom0vlaevla.show_contour(levels=[3], cmap='seismic', smooth=1)    # *10²¹ at/cm2
mom0vlaevla.show_contour(levels=[0.5, 1.2, 2.4, 3.7, 4.9, 6.1, 7.3, 8.6, 9.8, 11.0, 12.2], cmap='hsv', smooth=1)    # *10²¹ at/cm2
mom0vlaevla.show_contour(data=vlaim[0].data, levels=[0.5, 1.2, 2.4, 3.7, 4.9, 6.1, 7.3, 8.6, 9.8, 11.0, 12.2], cmap="jet")  # *10²¹ at/cm2
                
                mom0vlaevla.close()
                
# ##4.3 Deprojection of moment 0 (EVLA data)

                In IRAF:
    
geotran input=cig96_HI_mom0_3rms.fits output=cig96_HI_mom0_3rms_dep.fits xin=300.686 yin=300.887 xrot=69.61 yrot=69.61 xmag=1.430615
                
                # Select working directory
go_to('cig')

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

###

image = "cig96_HI_mom0_3rms_dep.fits"
coordinates = "cig96_pseudoring_regions_HI_dep_xy.reg"
radius = 2.5     # in this image 10 arcsec = 2.5 pix

apersum(image, coordinates, radius, "HI_sums.txt", "HI_errs.txt")
                
# These shall be used in the <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96%20-%20Pseudoring%20colors.ipynb#">CIG96 - Pseudoring colors</a> notebook.
