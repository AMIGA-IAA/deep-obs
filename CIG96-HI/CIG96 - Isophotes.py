
# coding: utf-8

# # CIG96 - Optical isophotes
# 
# ##GENERAL INFORMATION
# 
# Notebook with some plots of the isophotes fittings of CIG96.
# 
# ### AMIGA WIKI PAGES:
# 
# **Optical data**
# - Info from the CAHA campaign can be found <a href="http://amiga.iaa.es:8888/display/science/CAHA+2.2m+Observations+11sep12">here</a>.
# 
# - The data reduction and calibration notes and scripts can be found <a href="http://amiga.iaa.es:8888/display/science/Data+reduction+CIG96+CAHA2.2">here</a>.
# 
# - The information about CIG96 paper (2016) can be found <a href="http://amiga.iaa.es:8888/display/science/CIG96+-+Deep+optical+and+HI+images#CIG96-DeepopticalandHIimages-Paperstructureandcomments">here</a>.
# 
# 
# ###a) Import necessary modules

# In[1]:

import aplpy
import astropy.io.fits as fits
from   astropy import units as u
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
import pylab
import pyraf as pyraf
import pyregion as pyreg
import random
import repipy.combine as combine
import repipy.find_sky as find_sky
import scipy
from   scipy import stats
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


# In[3]:

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


# ##0. Functions
# 
# All final and working functions used throughout the notebook are stored in this cell.

# In[4]:

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
    '''
    x, y = np.loadtxt(iso_file, unpack=True, usecols=(0,1))
    fit = fitEllipse(x, y)
    center = ellipse_center(fit)
    phi = ellipse_angle_of_rotation(fit)
    axes = ellipse_axis_length(fit)
    return  center[0], center[1], phi, axes[0], axes[1], max(axes)/min(axes)
               
###############

def ping():
    os.system('spd-say "Ping"')


# In[5]:

ping()


# ##1. Disc isophote fit and comparison with ring elliptical fit in R-Deep image
# 
# ###1.1 Reference isophote selection

# In[8]:

# SB conversion to counts (see reduction notebook, Appendix A)
SB = [23, 23.5, 24, 24.5, 25]
counts = [round(50.5017 + 10**((elem - 24.2823 - 0.5525)/(-2.5)),4) for elem in SB]
SB = [round(elem,2) for elem in SB]

print SB
print counts

# Select working directory
go_to('peris')

# select and open image
image = "cig96_def_SB.fits"
im = fits.open(image)
data = im[0].data

# data smoothing for a better contour selection
data = scipy.ndimage.filters.gaussian_filter(data, 2, mode='nearest')

# plot
plt.figure(figsize=(10,10))
#plt.imshow(data.T, cmap='bone_r', vmin=49.5, vmax=58, origin='lower')  # in counts
plt.imshow(data.T, cmap='bone', vmin=23, vmax=27, origin='lower')  # in magnitudes
plt.colorbar(fraction=0.046, pad=0.04)
plt.gca().invert_xaxis()
cs = plt.contour(data.T, SB, cmap="rainbow")
plt.clabel(cs, inline=True, inline_spacing=-20, fontsize=16)
plt.xlim(700,400), plt.ylim(350,670)
#plt.xlim(800,200), plt.ylim(200,800)
plt.show()


# In[9]:

plt.close()


# Isohpotes **23.0** and **23.5 mag/arcsec2** are selected.
# 
# ###1.2 Elliptical fit of the disc isophotes
# 
# First we extract the 23.5 mag/"2 contour via ds9 and select the disc points. Second, we sample it with PSF resolution. Third, we fit the remaining points.

# In[11]:

# Select working directory
go_to('peris')

# Data selected
# image:            cig96_def_SB.fits
# positions files:   isophote_23.5.dat

#X,Y = np.loadtxt('isophote_23.5.dat', unpack=True, usecols=(0,1))  # all apoints
x,y = np.loadtxt('isophote_23.5_disc_sm3_xy.con', unpack=True, usecols=(0,1))  # no South points
X,Y = np.loadtxt('isophote_23.5_disc_sm3_clean_xy.con', unpack=True, usecols=(0,1))  # clean points selection

a = fitEllipse(X,Y)
center = ellipse_center(a)
phi = ellipse_angle_of_rotation(a)
axes = ellipse_axis_length(a)
    
print "center = ",  center
print "angle of rotation = ",  phi
print "axes = ", axes
    
a, b = axes
arc = np.pi
R = np.arange(0,arc*np.pi, 0.01)
cx = center[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
cy = center[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)


# In[12]:

# Select working directory
go_to('peris')

# CIG96 coordinates in deg
RA = 507
DEC = 538

# plot
fig = plt.figure(figsize=(10,10))
plt.title("Elliptical fit of SB = 23.0 and SB = 23.5 mag/arcsec$^{2}$ disc isophotes", fontsize=16)
plt.xlim(360,660)
plt.ylim(380,680)

# full isophote in BLUE 23.5 mag/"2
plt.plot(x,y,'co', label="Full 23.5 mag/arcsec$^{2}$ isophote")
#plt.plot(xx,yy, color = 'blue', label="23.5 mag/arcsec$^{2}$ isophote")

# selected points of isophote in RED 23.5 mag/"2
plt.plot(X,Y,'ro', zorder=2)
plt.plot(cx,cy, color = 'red', label="Disc fit at 23.5 mag/arcsec$^{2}$ (fit selection points in red)")
plt.plot(RA, DEC, 'y*', label="CIG96 center")
plt.plot(503.64810409, 537.9778365, 'rD', label="Fit center")
plt.show()

plt.legend(loc=0, fontsize=10)

# text in case 23.0 and 23.5 are plotted
#plt.text(380, 400, 'Fit centers (x,y pix): 500-540$_{(23.0mag/arcsec2)}$ vs 502-542$_{(23.5mag/arcsec2)}$\n\
#Axes relation (adim): 1.302$_{(23.0mag/arcsec2)}$ vs 1.298$_{(23.5mag/arcsec2)}$\n\
#Rotation angle (deg): 31.15$_{(23.0mag/arcsec2)}$ vs 28.41$_{(23.5mag/arcsec2)}$', fontsize=14)

#text in case ONLY 23.5 is plotted
plt.text(380, 390, #'Fit centers (x,y pix): 502-542$_{(23.5mag/arcsec2)}$\n\
#Axes relation (adim): 1.298$_{(23.5mag/arcsec2)}$\n\
'Rotation angle (deg): 35.29$_{(23.5mag/arcsec2)}$', fontsize=14)


# In[13]:

plt.close()


# In[14]:

# Select working directory
go_to('peris')

# figure structure
fig = plt.figure(figsize=(12,12))

# CIG96 coordinates in deg
RA  = 33.865231
DEC = 6.0020581
rad = 0.04  # zoom radius so, the larger rad, the larger FoV in the image

fRd = aplpy.FITSFigure('cig96_def_SB.fits', figure=fig, north=True)         # SB  version
fRd.set_tick_labels_font(size='small')
fRd.set_axis_labels_font(size='small')
fRd.set_title('R Deep + Disc isophote 23.5 mag/"2', size=12)
#fRd.show_colorscale(stretch='log', exponent=5, vmin=49.5, vmax=65, cmap="bone_r")    # counts version
fRd.show_colorscale(stretch='log', exponent=2, vmin=18, vmax=29, cmap="bone")         # SB version

# zoom
fRd.recenter(RA, DEC, radius=rad)

#######################

# Disc isophote fit, clean points
d2_cen =  [33.865116, 6.0014161]          # center position in degrees
d2_ang = 35.29                            # 0.507721352092 rad = 29,090296 deg
d2_axes =  [0.027648573, 0.020586471]       # width, height in degrees

#######################

# colors
fit = (0.7, 0.1, 0.2)
cigcen = "yellow"

# 23.5
#fRd.show_regions("cig96_23.5_fit.reg", zorder=1)         # all points fit
fRd.show_regions("isophote_23.5_disc_sm3_clean_wcs.reg", zorder=2) # clean points fit
fRd.show_ellipses(d2_cen[0], d2_cen[1], 2.08*d2_axes[0], 2.08*d2_axes[1], angle=122.29, linewidth=2, edgecolor=fit, facecolor="None", zorder=3)
fRd.show_ellipses(d2_cen[0], d2_cen[1], 0.0015, 0.0015, angle=0.0, linewidth=1, edgecolor=fit, facecolor=fit)
#fRd.show_contour(level=53.920831)                                        # counts version
fRd.show_contour(levels=[23.5], colors="cyan", zorder=1)                # SB version

# CIG96 center
fRd.show_markers(RA, DEC, marker="+", s=250, facecolor=cigcen, edgecolor=cigcen, zorder=3)

#legend
mpl.rcParams['legend.numpoints'] = 1
#blue = mpl.lines.Line2D([], [], color="blue", linestyle="-", linewidth=4, label='23.0 mag/"2 isophote')
red = mpl.lines.Line2D([], [], color=fit, linestyle="-", linewidth=4, label='23.5 mag/"2 isophote')
cen = mpl.lines.Line2D([], [], marker="+", linewidth=0, markeredgewidth=2, markersize=14, markeredgecolor=cigcen, markerfacecolor=cigcen, label='CIG96 center')
plt.legend(handles=[red, cen], fontsize=12)


# In[15]:

fRd.close()
plt.close()


# ###1.3 Pseudoring fit scaled to the disc 23.5mag/"2 isophote fit

# In[165]:

# Select working directory
go_to('peris')

#####################

# figure structure
fig = plt.figure(figsize=(18, 9))

# CIG96 coordinates in deg
RA  = 33.865231
DEC = 6.0020581
rad = 0.05  # zoom radius so, the larger rad, the larger FoV in the image

# plot window design
vw = 0.80 # vertical width of the subplots in %
hw = 0.38 # horizontal width of the subplots in % 


######################

# LEFT: R deep + pseudoring fit + disc isohpote fit
# R deep + pseudoring fit + disc isohpote fit
#fRd = aplpy.FITSFigure('cig96_def_crop_wcs.fits', figure=fig, north=True, subplot=[0.1, 0.1, hw, vw])  # counts version
fRd = aplpy.FITSFigure("cig96_def_SB.fits", figure=fig, north=True, subplot=[0.1, 0.1, hw, vw])         # SB  version
fRd.set_tick_labels_font(size='small')
fRd.set_axis_labels_font(size='small')
fRd.set_title("R Deep + Disc-fit scaled pseudoring fit + Disc isophote", size=12)
#fRd.show_colorscale(stretch='log', exponent=5, vmin=49.5, vmax=65, cmap="bone_r")    # counts version
fRd.show_colorscale(stretch='log', exponent=2, vmin=22, vmax=28, cmap="Greys_r")         # SB version

# zoom
fRd.recenter(RA, DEC, radius=rad)

# RIGHT: R deep - RGB Peris combined + ellipse + disc isohpote
fRdRGB = aplpy.FITSFigure('cig96_RdRGB.not_cal.fits', figure=fig, subplot=[0.5, 0.1, hw, vw], north=True)
fRdRGB.set_tick_labels_font(size='x-small')
fRdRGB.set_axis_labels_font(size='small')
fRdRGB.set_title("R deep-RGB Peris + Disc-fit scaled pseudoring fit + Disc isophote", size=12)
fRdRGB.show_colorscale(stretch='linear', vmin=-5, vmax=200, cmap="Greys")
# to hide
fRdRGB.hide_yaxis_label()
fRdRGB.hide_ytick_labels()
# zoom
fRdRGB.recenter(RA, DEC, radius=rad)

#######################

# Disc isophote fit, clean points
d2_cen =  [33.865116, 6.0014161]          # center position in degrees
d2_ang = 35.29                            # 0.507721352092 rad = 29,090296 deg
d2_axes =  [0.027648573, 0.020586471]       # width, height in degrees

#######################

# colors
fit = (0.7, 0.1, 0.2)    # red
isoph = "yellow"
cigcen = "orange"
pseudo = (0.1, 0.9, 0.3) # green

#######################

# LEFT PLOT

# Rd pseudoring fit data
center1 = [33.867654, 5.9998111]      # center position in degrees
dim1    = [0.032752181, 0.040832538]  # width, height in degrees
rot1    = 20.9557                     # angle of rotation: 0.3657465 rad = 20,955735181 deg

# isophote contour
#fRd.show_regions("cig96_23.5_fit.reg", zorder=1)    # all points
fRd.show_regions("isophote_23.5_disc_sm3_clean_wcs.reg", zorder=2) # clean points fit

# isophote fit
fRd.show_contour(levels=[23.5], colors=isoph)                # SB version
#fRd.show_contour(level=53.920831)                           # counts version

# disc fit
fRd.show_ellipses(d2_cen[0], d2_cen[1], 2.1*d2_axes[0], 2.1*d2_axes[1], angle=d2_ang+87, linewidth=3, edgecolor=fit, facecolor="None")    
fRd.show_ellipses(d2_cen[0], d2_cen[1], 0.0015, 0.0015, angle=0.0, linewidth=1, edgecolor=fit, facecolor=fit)    

# Rd pseudoring fit (major axis = maj.axis of disc isophote; minor axis = maj.axis * axis relation of pseudoring fit)
fRd.show_ellipses(center1[0], center1[1], 1.67*d2_axes[0], 1.67*d2_axes[0]*dim1[1]/dim1[0], angle=rot1, linestyle=":",                   linewidth=3, edgecolor=pseudo, facecolor="None", zorder=2)
fRd.show_ellipses(center1[0], center1[1], 0.001, 0.001, angle=0.0, linewidth=1, edgecolor=pseudo, facecolor=pseudo)

# CIG96 center
fRd.show_markers(RA, DEC, marker="+", s=250, facecolor=cigcen, edgecolor=cigcen)

######################

# RIGHT PLOT

# RdRGB pseudoring fit data
center2 = [33.8669, 6.0001361]       # center position in degrees
dim2    = [0.039267303, 0.048368206] # width, height in degrees
rot2    = 21.94479                   #  angle of rotation: 0.3657465 rad = 21.94479 deg

# Disc isophote
#fRdRGB.show_ellipses(d_cen[0], d_cen[1], 0.0015, 0.0015, angle=0.0, linewidth=1, edgecolor="red", facecolor="red")  # all points
#fRdRGB.show_regions("cig96_23.5_fit.reg", zorder=1)   # all points
fRdRGB.show_ellipses(d2_cen[0], d2_cen[1], 0.0015, 0.0015, angle=0.0, linewidth=1, edgecolor=fit, facecolor=fit) # no south
fRdRGB.show_regions("cig96_23.5_noSouth_fit.reg", zorder=1) # no south

# RdRGB pseudoring fit
fRdRGB.show_ellipses(center2[0], center2[1], 1.65*d2_axes[0], 1.65*d2_axes[0]*dim2[1]/dim2[0], angle=rot2, linestyle=":",                      linewidth=4, edgecolor=pseudo, facecolor="None", zorder=2)
fRdRGB.show_ellipses(center2[0], center2[1], 0.001, 0.001, angle=0.0, linewidth=1, edgecolor=pseudo, facecolor=pseudo)

# CIG96 center
fRdRGB.show_markers(RA, DEC, marker="+", s=250, facecolor="orange", edgecolor="orange")

#######################

#legend
mpl.rcParams['legend.numpoints'] = 1
disc = mpl.lines.Line2D([], [], color=fit, linestyle="-", linewidth=4, label='Disc isophote fit (23.5 mag/arcsec$^{2}$)')
iso = mpl.lines.Line2D([], [], color=isoph, linestyle="-", linewidth=4, label='Disc isophote (23.5 mag/arcsec$^{2}$)')
pse = mpl.lines.Line2D([], [], color=pseudo, linestyle=':', linewidth=4, label='Disc-fit scaled pseudoring fit')
cen = mpl.lines.Line2D([], [], marker="+", linewidth=0, markeredgewidth=2, markersize=14, markeredgecolor=cigcen, markerfacecolor=cigcen, label='CIG96 center')
plt.legend(handles=[iso, disc, pse, cen], fontsize=12)


# In[166]:

fRd.close()
fRdRGB.close()
plt.close()


# Just the galaxy and pseudoring fittings.

# In[94]:

# Select working directory
go_to('peris')

#####################

# figure structure
fig = plt.figure(figsize=(10, 10))

# CIG96 coordinates in deg
RA  = 33.865167
DEC = 6.002611
rad = 0.07  # zoom radius so, the larger rad, the larger FoV in the image

# plot window design
vw = 0.90 # vertical width of the subplots in %
hw = 0.90 # horizontal width of the subplots in % 


######################

# LEFT: R deep + pseudoring fit + disc isohpote fit
# R deep + pseudoring fit + disc isohpote fit
#fRd = aplpy.FITSFigure('cig96_def_crop_wcs.fits', figure=fig, north=True, subplot=[0.1, 0.1, hw, vw])  # counts version
fRd = aplpy.FITSFigure("cig96_def_SB.fits", figure=fig, north=True, subplot=[0.05, 0.05, hw, vw])         # SB  version
fRd.set_tick_labels_font(size='small')
fRd.set_axis_labels_font(size='small')
#fRd.set_title("R Deep + Disc-fit scaled pseudoring fit + Disc isophote", size=12)
#fRd.show_colorscale(stretch='log', exponent=5, vmin=49.5, vmax=65, cmap="bone_r")    # counts version
#fRd.show_colorscale(stretch='log', exponent=5, vmin=20, vmax=28, cmap="Greys_r")         # SB version
fRd.show_colorscale(stretch='log', exponent=10, vmin=22, vmax=28, cmap="bone_r")         # SB version
fRd.set_nan_color('k')

# inner contours
fRd.show_contour(levels=[19.4, 19.6, 19.8, 20.0, 20.2, 20.6, 21.0, 21.4, 21.8, 22.2, 22.6, 23.0], cmap='Pastel2', zorder=1)

# labels
fRd.set_tick_labels_font(size=16)
fRd.axis_labels.set_font(size=16)
fRd.set_tick_size(0)
# ticks
fRd.tick_labels.set_xformat('hh:mm:ss')
fRd.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='off', left='on', right='off', length=12, width=1.5, color='white')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='off', left='on', right='off', length=6, width=1.5, color='white')

# zoom
fRd.recenter(RA, DEC-0.002, radius=rad)

#######################

# Disc isophote fit, clean points
d2_cen =  [33.865116, 6.0014161]          # center position in degrees
d2_ang = 35.29                            # 0.507721352092 rad = 29,090296 deg
d2_axes =  [0.027648573, 0.020586471]       # width, height in degrees

#######################

# colors
isoph_fit = "magenta"
isoph = "maroon"
cigcen = "cyan"
pseudo = "yellow"
pseudo_fit = "orange"

#######################

# Rd pseudoring fit data
center1 = [33.867654, 5.9998111]      # center position in degrees
dim1    = [0.032752181, 0.040832538]  # width, height in degrees
rot1    = 20.9557                     # angle of rotation: 0.3657465 rad = 20,955735181 deg

# isophote contour
#fRd.show_regions("cig96_23.5_fit.reg", zorder=1)    # all points
fRd.show_regions("isophote_23.5_disc_sm3_clean_wcs_black.reg", zorder=2) # clean points fit
# pseudoring regions
fRd.show_regions("pseudoring-WCS-RdeepRGB-yellow.reg", zorder=2)

# isophote plot
fRd.show_contour(levels=[23.5], colors=isoph, linewidths=2)                # SB version
#fRd.show_contour(level=53.920831)                           # counts version

# isophote fit
fRd.show_ellipses(d2_cen[0], d2_cen[1], 2.1*d2_axes[0], 2.1*d2_axes[1], angle=d2_ang+87, linewidth=4, edgecolor=isoph_fit, facecolor='None',zorder=2)    
# black border
fRd.show_ellipses(d2_cen[0], d2_cen[1], 2.1*d2_axes[0], 2.1*d2_axes[1], angle=d2_ang+87, linewidth=7, edgecolor='k', facecolor='None',zorder=1)    

# isophote center
fRd.show_markers(d2_cen[0], d2_cen[1], marker='o', s=75, linewidth=2, edgecolor='k', facecolor=isoph_fit, zorder=1) 

# Rd pseudoring fit (major axis = maj.axis of disc isophote; minor axis = maj.axis * axis relation of pseudoring fit)
factor_inner = 1.67 # factor to multiply the pseudoring fit so it is scaled to the 23.5 isophot
factor_outer = 2.95 # factor to multiply the pseudoring fit so it fits the pseudoring
fRd.show_ellipses(center1[0], center1[1], factor_inner*d2_axes[0], factor_inner*d2_axes[0]*dim1[1]/dim1[0], angle=rot1, linestyle="-",                   linewidth=4, edgecolor=pseudo_fit, facecolor="None", zorder=2)
# black border
fRd.show_ellipses(center1[0], center1[1], factor_inner*d2_axes[0], factor_inner*d2_axes[0]*dim1[1]/dim1[0], angle=rot1, linestyle="-",                   linewidth=7, edgecolor='k', facecolor="None", zorder=1)

fRd.show_ellipses(center1[0], center1[1], factor_outer*d2_axes[0], factor_outer*d2_axes[0]*dim1[1]/dim1[0], angle=rot1, linestyle="-",                   linewidth=4, edgecolor=pseudo_fit, facecolor="None", zorder=2)
#black border
fRd.show_ellipses(center1[0], center1[1], factor_outer*d2_axes[0], factor_outer*d2_axes[0]*dim1[1]/dim1[0], angle=rot1, linestyle="-",                   linewidth=7, edgecolor='k', facecolor="None", zorder=1)

# isophote fitting center
fRd.show_markers(center1[0], center1[1], marker='o', s=75, linewidth=2, edgecolor='k', facecolor=pseudo_fit, zorder=1)
#fRd.show_ellipses(center1[0], center1[1], 0.001, 0.001, angle=0.0, linewidth=6, edgecolor=pseudo_fit, facecolor=pseudo_fit)

# CIG96 center
fRd.show_markers(RA, DEC, marker="X", s=75, linewidth=1, facecolor=cigcen, edgecolor='k', zorder=1)

#######################

# save image
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_pseudoring_fitting_dark.png', transparent=True, bbox_inches='tight')
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_pseudoring_fitting_dark.eps', format='eps', dpi=30, transparent=True, bbox_inches='tight')

ping()


# In[41]:

fRd.close()
plt.close()


# ##2. Multiple isohpotes fitting
# 
# ###2.1 Cleaning the isophotes
# 
# Raw isophotes are full of stars that alter the true galaxy isophote shape. We clean the points of undesiderd regions in **image coordinates**:

# In[5]:

go_to("opt")

# cleaning the following isophote files
raw_files = ["isoph_202.con", "isoph_210.con", "isoph_215.con", "isoph_220.con", "isoph_225.con", "isoph_230.con",              "isoph_235.con", "isoph_240.con", "isoph_245.con", "isoph_250.con", "isoph_255.con", "isoph_260.con", "isoph_264.con"]

reg202 = []                                            # nothing to edit out
reg210 = [[470, 556, 507, 610], [506, 566, 531, 584], [480, 494, 510, 513]]
reg215 = [[456, 556, 507, 606], [433, 500, 462, 531], [435, 482, 493, 520], [515, 552, 560, 589], [542, 544, 559, 553]]
reg220 = [[536, 494, 561, 521], [427, 543, 512, 623], [559, 552, 573, 567], [533, 570, 564, 597], [516, 594, 532, 605]]
reg225 = [[577, 560, 585, 569], [435, 484, 453, 498], [419, 509, 429, 520], [420, 531, 430, 538], [426, 557, 510, 623],           [542, 594, 548, 599]]
reg230 = [[586, 575, 596, 586], [517, 464, 530, 472], [485, 454, 492, 462], [464, 457, 469, 464], [425, 565, 525, 626]]
reg235 = [[594, 521, 614, 605], [570, 601, 588, 619], [388, 473, 426, 550], [414, 558, 539, 664]]
reg240 = [[608, 582, 622, 595], [613, 538, 620, 557], [592, 503, 599, 512], [490, 442, 499, 448], [407, 456, 420, 467],           [385, 501, 393, 528], [393, 547, 408, 570], [404, 578, 528, 674]]
reg245 = [[612, 566, 630, 602], [614, 516, 633, 542], [388, 432, 442, 480], [377, 587, 546, 690], [587, 620, 598, 628]]
reg250 = [[603, 594, 640, 644], [624, 520, 632, 527], [590, 465, 618, 504], [519, 415, 586, 465], [365, 397, 430, 481],           [364, 542, 378, 552], [306, 563, 414, 658], [423, 594, 646, 703]]
reg255 = [[514, 588, 662, 786], [651, 478, 682, 587], [583, 417, 661, 508], [523, 404, 558, 416], [368, 353, 516, 411],           [341, 407, 384, 479], [295, 523, 393, 676], [396, 682, 412, 693]]
reg260 = [[412, 662, 657, 827], [663, 588, 711, 657], [683, 546, 693, 565], [561, 354, 627, 433], [330, 328, 516, 490],           [274, 513, 370, 710], [391, 699, 407, 719]]
reg264 = [[243, 515, 636, 861], [617, 687, 675, 742], [657, 630, 704, 679], [695, 532, 725, 614], [648, 426, 684, 473],           [408, 292, 647, 398], [312, 467, 326, 483]]

regions = [reg202, reg210, reg215, reg220, reg225, reg230, reg235, reg240, reg245, reg250, reg255, reg260, reg264]


# In[6]:

for isoph, regs in zip(raw_files, regions):
    clean_isoph(isoph, regs, isoph[:-4] + "_clean.con")


# ###2.2 Fitting the clean isophotes

# In[7]:

clean_files = ["isoph_202_clean.con", "isoph_210_clean.con", "isoph_215_clean.con", "isoph_220_clean.con",                "isoph_225_clean.con", "isoph_230_clean.con", "isoph_235_clean.con", "isoph_240_clean.con",                "isoph_245_clean.con", "isoph_250_clean.con", "isoph_255_clean.con", "isoph_260_clean.con",                "isoph_264_clean.con"]

for isofile in clean_files:
    cenx, cey, ang, majax, minax, relax = isofit(isofile)
    print cenx, cey, ang, majax, minax, relax


# In[54]:

raw_files = ["isoph_202.con", "isoph_210.con", "isoph_215.con", "isoph_220.con", "isoph_225.con", "isoph_230.con",              "isoph_235.con", "isoph_240.con", "isoph_245.con", "isoph_250.con", "isoph_255.con", "isoph_260.con", "isoph_264.con"]

# clean isophotes
clean_files = ["isoph_202_clean.con", "isoph_210_clean.con", "isoph_215_clean.con", "isoph_220_clean.con",                "isoph_225_clean.con", "isoph_230_clean.con", "isoph_235_clean.con", "isoph_240_clean.con",                "isoph_245_clean.con", "isoph_250_clean.con", "isoph_255_clean.con", "isoph_260_clean.con", "isoph_264_clean.con"]

# magnitudes
magnitudes = ["20.2", "21.0", "21.5", "22.0", "22.5", "23.0", "23.5", "24.0", "24.5", "25.0", "25.5", "26.0", "26.4"]
#fitcolors  = ['darksalmon', 'salmon', 'orangered', 'firebrick', 'crimson', 'mediumvioletred', 'darkviolet', 'slateblue', 'darkblue', 'dodgerblue', 'ightseagreen', 'darkgreen', 'goldenrod']

# plotting isophotes, fittings
#files = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
files = [0,1,2,3,5,8,10,12]
#start = 0

print raw_files[files[5]]


# In[50]:

# Select working directory
go_to('opt')

# CIG96 coordinates in deg
RA = 507
DEC = 538

# figure definition
fig = plt.figure(figsize=(12,13))
#plt.title("Elliptical fittings of 21.0-26.4 mag/arcsec$^{2}$ dics isophotes", fontsize=16)
plt.xlim(250, 770)
plt.ylim(280, 850)


# background image
#image = "cig96_def_SB.fits"
#im = fits.open(image)
#data = im[0].data
#plt.imshow(data.T, cmap='bone', vmin=23, vmax=27, origin='lower', zorder=1)
#plt.gca().invert_xaxis()


# files
# raw isophotes
raw_files = ["isoph_202.con", "isoph_210.con", "isoph_215.con", "isoph_220.con", "isoph_225.con", "isoph_230.con",              "isoph_235.con", "isoph_240.con", "isoph_245.con", "isoph_250.con", "isoph_255.con", "isoph_260.con", "isoph_264.con"]

# clean isophotes
clean_files = ["isoph_202_clean.con", "isoph_210_clean.con", "isoph_215_clean.con", "isoph_220_clean.con",                "isoph_225_clean.con", "isoph_230_clean.con", "isoph_235_clean.con", "isoph_240_clean.con",                "isoph_245_clean.con", "isoph_250_clean.con", "isoph_255_clean.con", "isoph_260_clean.con", "isoph_264_clean.con"]

# magnitudes
magnitudes = ["20.2", "21.0", "21.5", "22.0", "22.5", "23.0", "23.5", "24.0", "24.5", "25.0", "25.5", "26.0", "26.4"]
#fitcolors  = ['darksalmon', 'salmon', 'orangered', 'firebrick', 'crimson', 'mediumvioletred', 'darkviolet', 'slateblue', 'darkblue', 'dodgerblue', 'ightseagreen', 'darkgreen', 'goldenrod']

# plotting isophotes, fittings
files = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
start = 0
end   = 13

#for raw, clean, mag in zip(raw_files[start:end], clean_files[start:end], magnitudes_a[start:end]):
for raw, clean, mag in zip(raw_files[files[0]:files[13]:1], clean_files[files[0]:files[13]:1], magnitudes[files[0]:files[13]:1]):
    ran_col = (round(np.random.random(),2), round(np.random.random(),2), round(np.random.random(),2))  # color selection
    
    # full isophotes points
    rx, ry = np.loadtxt(raw, unpack=True, usecols=(0,1)) # full isophotes
    plt.scatter(rx, ry, color=ran_col, s=2, marker='o', alpha=0.9, zorder=4)#, label=str(mag) + ' mag/"2  full isophote')
    #plt.scatter(rx, ry, color=fitcolor, s=3, marker='o', alpha=0.4, zorder=3, label=str(mag) + ' mag/"2  full isophote')
    
    # selection of points for fitting
    ix, iy   = np.loadtxt(clean, unpack=True, usecols=(0,1)) # clean isopohotes
    plt.scatter(ix, iy, color=ran_col, s=50, marker='o', alpha=0.6, zorder=5, facecolor='k', edgecolor=ran_col)#, label="Fitting points of " + str(mag) + ' mag/"2 isophote')
    #plt.scatter(ix, iy, color=fitcolor, s=15, marker='d', alpha=0.9, zorder=1, facecolor='None', edgecolor=ran_col, label="Fitting points of " + str(mag) + ' mag/"2 isophote')
    
    # elliptical fitting
    cenx, ceny, ang, a, b, rel_ax = isofit(clean)
    arc = np.pi
    R = np.arange(0,arc*np.pi, 0.01)
    cx = cenx + a*np.cos(R)*np.cos(ang) - b*np.sin(R)*np.sin(ang)
    cy = ceny + a*np.cos(R)*np.sin(ang) + b*np.sin(R)*np.cos(ang)
    
    # plot of elliptical fitting
    plt.plot(cx, cy, color=ran_col, linestyle='-', linewidth=3, zorder=5, label=str(mag) + ' mag/"2 fitting')
    plt.plot(cx, cy, color='k', linestyle='-', linewidth=5, zorder=2)
    #plt.plot(cx, cy, color=fitcolor, linestyle='-', linewidth=2, zorder=2, label=str(mag) + ' mag/"2 elliptical fit')
    plt.show()

# legend
plt.legend(loc=0, fontsize=14, bbox_to_anchor=(1,1))
# legend outside the plot, usage sample: ax.legend(bbox_to_anchor=(1.1, 1.05))

# save image
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_isophotes_fittings_'+str(round(np.random.random(),3))+'.png', transparent=True, bbox_inches='tight')

ping()


# In[55]:

plt.close()


# In[92]:

# Select working directory
go_to('opt')

# CIG96 coordinates in deg
RA = 507
DEC = 538

# figure definition
fig = plt.figure(figsize=(12,13))
#plt.title("Elliptical fittings of 21.0-26.4 mag/arcsec$^{2}$ dics isophotes", fontsize=16)
plt.xlim(250, 770)
plt.ylim(280, 850)


# background image
#image = "cig96_def_SB.fits"
#im = fits.open(image)
#data = im[0].data
#plt.imshow(data.T, cmap='bone', vmin=23, vmax=27, origin='lower', zorder=1)
#plt.gca().invert_xaxis()


# files
# raw isophotes
raw_files = ["isoph_202.con", "isoph_210.con", "isoph_215.con", "isoph_220.con", "isoph_225.con", "isoph_230.con",              "isoph_235.con", "isoph_240.con", "isoph_245.con", "isoph_250.con", "isoph_255.con", "isoph_260.con", "isoph_264.con"]

# clean isophotes
clean_files = ["isoph_202_clean.con", "isoph_210_clean.con", "isoph_215_clean.con", "isoph_220_clean.con",                "isoph_225_clean.con", "isoph_230_clean.con", "isoph_235_clean.con", "isoph_240_clean.con",                "isoph_245_clean.con", "isoph_250_clean.con", "isoph_255_clean.con", "isoph_260_clean.con", "isoph_264_clean.con"]

# magnitudes
magnitudes = ["20.2", "21.0", "21.5", "22.0", "22.5", "23.0", "23.5", "24.0", "24.5", "25.0", "25.5", "26.0", "26.4"]
#fitcolors  = ['darksalmon', 'salmon', 'orangered', 'firebrick', 'crimson', 'mediumvioletred', 'darkviolet', 'slateblue', 'darkblue', 'dodgerblue', 'ightseagreen', 'darkgreen', 'goldenrod']

# plotting isophotes, fittings
#files = [0,1,2,3,4,5,6,7,8,9,10,11,12]
files = [0,1,2,3,5,8,10,12]

for item in files:
    raw   = raw_files[item]
    clean = clean_files[item]
    mag   =  magnitudes[item]
    ran_col = (round(np.random.random(),2), round(np.random.random(),2), round(np.random.random(),2))  # color selection
    
    # full isophotes points
    rx, ry = np.loadtxt(raw, unpack=True, usecols=(0,1)) # full isophotes
    plt.scatter(rx, ry, color=ran_col, s=2, marker='o', alpha=0.9, zorder=4)#, label=str(mag) + ' mag/"2  full isophote')
    #plt.scatter(rx, ry, color=fitcolor, s=3, marker='o', alpha=0.4, zorder=3, label=str(mag) + ' mag/"2  full isophote')
    
    # selection of points for fitting
    ix, iy   = np.loadtxt(clean, unpack=True, usecols=(0,1)) # clean isopohotes
    plt.scatter(ix, iy, color=ran_col, s=50, marker='o', alpha=0.6, zorder=5, facecolor='k', edgecolor=ran_col)#, label="Fitting points of " + str(mag) + ' mag/"2 isophote')
    #plt.scatter(ix, iy, color=fitcolor, s=15, marker='d', alpha=0.9, zorder=1, facecolor='None', edgecolor=ran_col, label="Fitting points of " + str(mag) + ' mag/"2 isophote')
    
    # elliptical fitting
    cenx, ceny, ang, a, b, rel_ax = isofit(clean)
    arc = np.pi
    R = np.arange(0,arc*np.pi, 0.01)
    cx = cenx + a*np.cos(R)*np.cos(ang) - b*np.sin(R)*np.sin(ang)
    cy = ceny + a*np.cos(R)*np.sin(ang) + b*np.sin(R)*np.cos(ang)
    
    # plot of elliptical fitting
    plt.plot(cx, cy, color=ran_col, linestyle='-', linewidth=3, zorder=5, label=str(mag) + ' mag/"2 fitting')
    plt.plot(cx, cy, color='k', linestyle='-', linewidth=5, zorder=2)
    #plt.plot(cx, cy, color=fitcolor, linestyle='-', linewidth=2, zorder=2, label=str(mag) + ' mag/"2 elliptical fit')
    plt.show()

# legend
plt.legend(loc=0, fontsize=14, bbox_to_anchor=(1,1))
# legend outside the plot, usage sample: ax.legend(bbox_to_anchor=(1.1, 1.05))

# save image
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_isophotes_fittings_'+str(round(np.random.random(),3))+'.png', transparent=True, bbox_inches='tight')

ping()


# In[75]:

plt.close()


                # Select direcytory
go_to("opt")

# figure structure
fig = plt.figure(figsize=(12,12))

# CIG96 coordinates in deg
RA  = 33.865231
DEC = 6.0020581
rad = 0.08  # zoom radius so, the larger rad, the larger FoV in the image

isoph = aplpy.FITSFigure('cig96_def_SB.fits', figure=fig, north=True, zorder=1)         # SB  version
isoph.set_tick_labels_font(size='small')
isoph.set_axis_labels_font(size='small')
isoph.set_title('Isophotes fitting', size=12)
isoph.show_colorscale(stretch='log', exponent=2, vmin=20, vmax=27.5, cmap="bone")         # SB version

# zoom
isoph.recenter(RA, DEC, radius=rad)

# Data selected
# image:            cig96_def_SB.fits
# isophotes: isoph_xxx.reg where xxx = 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265
#           meaning the 21.0, 21.5, 22.0 etc. magnitudes

mags = [21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0, 25.5, 26.0]
isoph_files = ["isoph_210.con", "isoph_215.con", "isoph_220.con", "isoph_225.con", \
               "isoph_230.con", "isoph_235.con", "isoph_240.con", "isoph_245.con", \
               "isoph_250.con", "isoph_255.con", "isoph_260.con"]


for elem, mag in zip(isoph_files, mags):
    cenX, cenY, phi, majax, minax, axrel = isofit(elem)
    print cenX, cenY, phi, majax, minax, axrel
    ran_col = round(np.random.random(),1)
    #isoph.show_regions(elem[:-4] + ".reg", zorder=2) # plots isophote points
    isoph.show_ellipses(cenX, cenY, majax, minax, angle=phi, linewidth=2, \
                        edgecolor=(ran_col,ran_col,ran_col), facecolor=(ran_col,ran_col,ran_col), zorder=2)
    
# save image
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_caha22m_', transparent=True, bbox_inches='tight')
                
                isoph.close()
plt.close()
                
# Background for the previous image.

# In[76]:

# Select direcytory
go_to("opt")

# figure structure
fig = plt.figure(figsize=(12,12))

# CIG96 coordinates in deg
RA  = 33.869231
DEC = 6.0020581
rad = 0.085  # zoom radius so, the larger rad, the larger FoV in the image

# make faint colormap: https://stackoverflow.com/questions/37327308/add-alpha-to-an-existing-matplotlib-colormap
# Choose colormap
cmap = plt.cm.Blues
# Get the colormap colors
my_cmap = cmap(np.arange(cmap.N))
# Set alpha
my_cmap[:,-1] = np.linspace(0, 0.6, cmap.N)
# Create new colormap
my_cmap = mpl.colors.ListedColormap(my_cmap)

isoph = aplpy.FITSFigure('cig96_def_SB.fits', figure=fig, north=True, zorder=1)         # SB  version
isoph.set_tick_labels_font(size='small')
isoph.set_axis_labels_font(size='small')
# inner contours
isoph.show_contour(levels=[19.4, 19.6, 19.8, 20.0, 20.2, 20.6, 21.0], cmap='Pastel2', zorder=1)#, 21.0, 21.4, 21.8, 22.2, 22.6, 23.0], cmap='Pastel2', zorder=1)
#isoph.set_title('Isophotes fitting', size=12)
vmin = 20
vmax = 27.5
isoph.show_colorscale(stretch='log', exponent=5, vmin=vmin, vmax=vmax, cmap=my_cmap)         # SB version
#isoph.set_nan_color('rosybrown') # for reddish maps
isoph.set_nan_color('lightslategrey') # for blueish maps

# labels
isoph.set_tick_labels_font(size=16)
isoph.axis_labels.set_font(size=16)
isoph.set_tick_size(0)
# ticks
isoph.tick_labels.set_xformat('hh:mm:ss')
isoph.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='out', bottom='on', top='off', left='on', right='off', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='out', bottom='on', top='off', left='on', right='off', length=6, width=1.5, color='black')

# zoom
isoph.recenter(RA, DEC, radius=rad)

# save image
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_caha22m_'+str(vmin)+'_'+str(vmax)+'.png', transparent=True, bbox_inches='tight')


# In[93]:

isoph.close()
plt.close()


# # 3. Ellipticity of the pseudoring (deprojected image)
# 
# Definition: $Ellipticity = (major - minor) / major$
# 
# **WARNING**: when using a de-projected image, you need to be careful with what coordinates you feed to the ellipse function. The distortion caused by the de-projection is relevant when taking with IRAF either the physical or image coordinates. Suggestion: take manual regions and save their local coordinates (XY, image).

# In[15]:

# Select direcytory
go_to("opt")

# pseudoring X,Y cooridnates file
pseudocoords = 'cig96_dep_xy_border.reg'

cenx, cey, ang, majax, minax, relax = isofit(pseudocoords)
print cenx, cey, ang, majax, minax, relax

e = (majax - minax)/majax

print 'Pseudoring ellipticity = ', np.round(e,3), " = ", 100*np.round(e,3), "%"

