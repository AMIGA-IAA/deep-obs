
# coding: utf-8

# # CIG96 - CIG96 - Asymmetry, optical plots, HI-pseudoring and cirrus
# 
# ##GENERAL INFORMATION
# 
# Notebook with the analysis of CIG96.
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

# In[2]:

import aplpy
import astropy.io.fits as fits
from   astropy import units as u
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


# ###b) Data and working directory

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
    'surf_bri' is the surface brightness of the isophote'''
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


# #1. HI system LSRK velocity and asymmetry recalculation
# 
# ##1.1 EVLA+VLA datacube LSRK velocity
# 
# Medium velocity at 20% of peak flux. As described in <a href="http://www.aanda.org/articles/aa/pdf/2011/08/aa16117-10.pdf">Espada et al. 2011, Fig.2</a>.

# In[4]:

# Select working directory
go_to('HI')

# read text file with:
# x column: radio velocity [km/s]
# y column: [Jy] Flux Density
data = open("cig96_circ_intflux_LSRK_fluxdens_Jy.txt")

# take columns to individual lists
vel = []
flux = []

for column in (raw.strip().split() for raw in data):
    vel.append(column[0])
    flux.append(column[1])

# make floats
flux = [float(elem) for elem in flux]
    
# cut at 20% that defines horizontal cut
cut = float(max(flux))*0.20

# figure structure
fig, ax = plt.subplots()

# plot and axes titles
#plt.title("CIG96 (VLA+EVLA) - System LSRK velocity")
plt.xlabel("LSRK radio velocity (km/s)")
plt.ylabel("Flux density (Jy)")
#ax.plot([min(vel),max(vel)], [0,0], color='blue', linestyle='-', linewidth=1)
plt.axhline(0.0, linestyle="--", color='k', alpha=0.3)
plt.grid(alpha=0.2)

# plot flux vs velocity
ax.plot(vel, flux, color='b', linestyle='-', linewidth=2, label="EVLA+VLA HI cube (this work)")

# line at 20% + medium velocity
#ax.plot([float(min(vel))+30, float(max(vel))-30], [cut, cut], color='r', linestyle=':', linewidth=2, zorder=1)
plt.axhline(cut, linestyle=":", color='r', alpha=0.5)

# visual selection of V(low) and V(high) - cuts of the horizontal F20 line
vl = 1425.6
vh = 1662.7

vmean = (vl + vh) / 2

print "20% of the flux highest peak\nF20 =", round(cut,4), "Jy"

# Markers for the system velocity plot
#ell = mpl.patches.Ellipse((vmean,cut), 2, 0.02, color='g', fill=True)
plt.plot(vmean, cut, 'go', markersize=4, label="LSRK central velocity")
plt.plot(vmean+5, cut, 'g|')
plt.plot(vmean-5, cut, 'g|')
#ax.annotate("V = " + str(round(vmean,0)) + " km/s",(vmean-50,cut/2))
ax.annotate("V$_{LSRK}$ = 1544 km/s",(vmean-50,cut/2))
#ax.add_artist(ell)

plt.legend()


# In[5]:

plt.close()


# ##1.2 Comparison between single-dish and EVLA+VLA velocities

# In[4]:

# Select working directory
go_to('HI')

# read text file with:
# x column: radio velocity [km/s]
# y column: [Jy] Flux Density
data = open("cig96-circ-radiovelLSRK-intflux.txt")
GB = open("CIG0096_HH98_sb.txt")
comp = open("cig96comp_integrated_flux.txt")

# declare individual lists
vel = []
flux = []
velGB = []
fluxGB = []
velc = []
fluxc = []

# extract column values
for column in (raw.strip().split() for raw in data):
    vel.append(column[0])
    flux.append(column[1])
    
for column in (raw.strip().split() for raw in GB):
    velGB.append(column[0])
    fluxGB.append(column[1])
    
for column in (raw.strip().split() for raw in comp):
    velc.append(column[0])
    fluxc.append(column[1])
    
# make floats
flux = [float(elem) for elem in flux]
fluxGB = [float(elem) for elem in fluxGB]
fluxc = [float(elem) for elem in fluxc]
maxflux = max(flux)
maxfluxc = max(fluxc)


# In[26]:

# figure structure
fig = plt.subplots(figsize=(9,8), sharex=True)
gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1.5], hspace=0.0) 
#fig.subplots_adjust(hspace=0)   # no space between plots

# axes with different sizes
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])

# plot and axes titles
#plt.title("CIG96 (VLA+EVLA) - System LSRK velocity")
ax0.set_xlabel("", fontsize=14)
ax0.set_ylabel("Flux density (Jy)", fontsize=18)
#ax.plot([min(vel),max(vel)], [0,0], color='blue', linestyle='-', linewidth=1)
ax0.axhline(0.0, linestyle="--", color='k', alpha=0.3)
ax0.grid(alpha=0.2)

# plot EVLA+VLA data
ax0.plot(vel, flux, color='b', linestyle='-', linewidth=2, alpha=0.7, label="CIG 96 EVLA+VLA (this work)")

# plot GB43m data
ax0.plot([float(velocity)-17 for velocity in velGB], fluxGB, color='m', linestyle='--', linewidth=2, alpha=0.5, label="Green Bank 43m")

# plot Companion data
ax0.plot([float(velocity)-17 for velocity in velc], fluxc, color='green', linewidth=2, alpha=0.7, label="CIG 96 companion")
ax1.plot([float(velocity)-17 for velocity in velc], fluxc, color='green', linewidth=2, alpha=0.7, label="CIG 96 companion")
ax1.axhline(0.0, linestyle="--", color='k', alpha=0.3)
ax1.set_xlabel("LSRK radio velocity (km/s)", fontsize=18)
ax1.set_ylabel("Flux density (Jy)", fontsize=18)
ax1.grid(alpha=0.2)

# tick sizes
ax0.set_yticklabels(np.arange(-0.1, 0.8, 0.1), fontsize=14)
ax1.set_xticklabels(np.arange(1300, 1850, 100), fontsize=14)
ax1.set_yticklabels(np.arange(-0.002, 0.009, 0.002), fontsize=14)

######## central velocities

# CIG96 and companion "cutc20c" cuts at 20% that defines horizontal cut
cut20 = float(max(flux))/5
cut30 = float(max(flux))/3.3333
cut40 = float(max(flux))/2.5
cut50 = float(max(flux))/2
cut20c = float(max(fluxc))/5 
cut50c = float(max(fluxc))/2 
cutlist = [cut20, cut30, cut40, cut50, cut20c]

# line at 20% + medium velocity
ax0.axhline(cut20, linestyle="--", linewidth=2, color='k', alpha=0.6)
ax1.axhline(cut20c, linestyle="--", linewidth=2, color='k', alpha=0.6)
#ax1.axhline(cut50c, linestyle="--", linewidth=2, color='k', alpha=0.6)

percent = [20, 30, 40, 50]
vl = [1425.6, 1429.7, 1432.2, 1434.6]
vh = [1662.7, 1658.8, 1656.2, 1653.6]
color = ['g', 'r', 'y', 'k']

percent = [20]
vl = [1425.6]
vh = [1662.7]
vlc = 1521.1
vhc = 1634.7
color = ['b']

# CIG 96 central velocities calculation
print "CIG 96 central velocity"
for cut, per, vmin, vmax, col in zip(cutlist, percent, vl, vh, color):
    vmean = float((vmin + vmax) / 2)
    verr = 0.23   # according to Fouque et al. 1990, the error is 3x sigma, where sigma is the channel width
    print per, "% of the flux highest peak\nF", per, "=", round(cut,4), "Jy"
    print "System velocity = ", vmean, "+/-", round(verr,2), "km/s", "\n"
    # Markers for the system velocity plot
    ax0.errorbar(vmean, cut, color='b', marker='o', linewidth='0', markersize=8, markeredgecolor='k', capsize=8,
                 label="V$_{LSRK\ (W" + str(per) + ")}$ = " + str(vmean) + "$\pm$" + str(verr) +" km/s")  # xerr=verr+0.5,
    ax0.plot([vmean,vmean], [0.0, 0.3], color=col, linestyle=":", alpha=0.7)
    
# Companion central velocities calculation
print "Companion central velocity"
vmean = float((vlc + vhc) / 2)
verr = 0.23   # according to Fouque et al. 1990, the error is 3x sigma, where sigma is the channel width
print per, "% of the flux highest peak\nF", per, "=", round(cut20c,4), "Jy"
print "System velocity = ", vmean, "+/-", round(verr,2), "km/s", "\n"
# Markers for the system velocity plot
ax1.errorbar(vmean, cut20c, color='g', marker='o', linewidth='0', markersize=8, markeredgecolor='k', capsize=8,
               label="V$_{LSRK\ (W" + str(per) + ")}$ = " + str(vmean) + "$\pm$" + str(verr) +" km/s")  # xerr=verr+0.5,
ax1.plot([vmean,vmean], [0.0, 0.0072], color='g', linestyle=":", alpha=0.7)


# plot limits
ax0.set_xlim(1300, 1810)
ax0.set_xticklabels([])
ax1.set_xlim(1300, 1810)
#ax[0].ylim(-0.05, 0.75)
#ax[1].xlim(1300, 1800)
#ax[1].ylim(-0.05, 0.02)

#plt.legend(fontsize=14)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96-lsrkvelocity.png', bbox_inches='tight')


# In[27]:

plt.close()


# ###SNR between noise level and highest peak (SNR_peak) for CIG 96 and Companion

# In[78]:

# noise level data
data = open("cig96_circ_intflux_noise_LSRK_fluxdens_Jy.txt")  # CIG 96

# take columns to individual lists
vel = []
flux = []

for column in (raw.strip().split() for raw in data):
    vel.append(column[0])
    flux.append(column[1])

# make floats
flux = [float(elem) for elem in flux]

print "CIG 96"
print "stdev =", np.std(flux)
print "rms =", np.sqrt(np.mean(np.square(flux)))

snr_peak = maxflux/np.std(flux)    # Companion maxflux is measured in the complete profile
print "snr_peak =", snr_peak


# In[79]:

# noise level data
data = open("cig96comp_integrated_flux_noise.txt")   # Compainon

# take columns to individual lists
vel = []
flux = []

for column in (raw.strip().split() for raw in data):
    vel.append(column[0])
    flux.append(column[1])

# make floats
flux = [float(elem) for elem in flux]

print "Companion"
print "stdev =", np.std(flux)
print "rms =", np.sqrt(np.mean(np.square(flux)))

#snr_peak = maxfluxc/np.std(flux)    # CIG 96 maxflux is measured in the complete profile
snr_peak = maxfluxc/np.std(flux)    # Companion maxflux is measured in the complete profile
print "snr_peak =", snr_peak


# ###NW feature integrated spectrum

# In[9]:

# Select working directory
go_to('HI')

# read text file with:
# x column: radio velocity [km/s]
# y column: [Jy] Flux Density
datanw = open("featNW.txt")
datase = open("featSE.txt")

# take columns to individual lists
velnw = []
fluxnw = []
velse = []
fluxse = []

for column in (raw.strip().split() for raw in datanw):
    velnw.append(column[0])
    fluxnw.append(column[1])
    
for column in (raw.strip().split() for raw in datase):
    velse.append(column[0])
    fluxse.append(column[1])

# make floats
fluxnw = [float(elem) for elem in fluxnw]
fluxse = [float(elem) for elem in fluxse]

# plot
fig = plt.figure(figsize=(14, 6))
# plot window design
vw = 0.45 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

nw = plt.subplot(121)
plt.xticks(np.arange(1250, 1850, 50), fontsize=12)
plt.yticks(np.arange(-0.004, 0.008, 0.001), fontsize=12)
#nw.axis('equal')
se = plt.subplot(122)
plt.xticks(np.arange(1250, 1850, 50), fontsize=12)
plt.yticks(np.arange(-0.004, 0.008, 0.001), fontsize=12)
#se.axis('equal')

# plot flux vs velocity
# NW and SE
nw.plot(velnw, fluxnw, color='g', linestyle='--', alpha=0.3, linewidth=2)
nw.plot(velse, fluxse, color='m', linestyle='--', alpha=0.3, linewidth=2)

# SE and NW
se.plot(velse, fluxse, color='m', linestyle='--', alpha=0.3,linewidth=2)
se.plot(velnw, fluxnw, color='g', linestyle='--', alpha=0.3, linewidth=2)

# NW feature only
vin  = 1480
vout = 1550
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velnw_feat  = velnw[begin:stop:1]
fluxnw_feat = fluxnw[begin:stop:1]
nw.plot(velnw_feat, fluxnw_feat, color='g', linestyle='-', alpha=0.9, linewidth=2)
# SE feature only
vin  = 1600
vout = 1640
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velse_feat  = velse[begin:stop:1]
fluxse_feat = fluxse[begin:stop:1]
se.plot(velse_feat, fluxse_feat, color='m', linestyle='-', alpha=0.9, linewidth=2)

# plot and axes titles
nw.set_xlabel("LSRK radio velocity (km/s)", fontsize=15)
nw.set_ylabel("Flux density (Jy)", fontsize=14)
nw.set_ylim(-0.0035, 0.006)
se.set_xlabel("LSRK radio velocity (km/s)", fontsize=15)
se.set_ylabel("", fontsize=0)
se.set_ylim(-0.0035, 0.006)

# baseline
nw.axhline(0.0, linestyle="--", color='k', alpha=0.3)
nw.grid(alpha=0.2)
se.axhline(0.0, linestyle="--", color='k', alpha=0.3)
se.grid(alpha=0.2)

# Markers for the system velocity plot
#ell = mpl.patches.Ellipse((vmean,cut), 2, 0.02, color='g', fill=True)
nw.plot(1544, 0.0, 'bo', markersize=8, label="CIG 96 LSRK central velocity")
se.plot(1544, 0.0, 'bo', markersize=8, label="CIG 96 LSRK central velocity")

# features location
nw.annotate("NW feature", (1380, 0.005), fontsize=14)
nw.annotate('', xy=(1495, 0.0044), xytext=(1460, 0.0048), arrowprops=dict(facecolor='g'))

se.annotate("SE feature", (1630, 0.0041), fontsize=14)
se.annotate('', xy=(1634, 0.0036), xytext=(1670, 0.0039), arrowprops=dict(facecolor='m'))

#nw.legend(loc=0, fontsize=11)
#se.legend(loc=0, fontsize=11)
plt.savefig("/home/prm/Desktop/paper_cig96/images/cig96_featsNWSE_intspectrum2.png", bbox_inches='tight')


# In[10]:

plt.close()


# In[ ]:

# Select working directory
go_to('HI')

# read text file with:
# x column: radio velocity [km/s]
# y column: [Jy] Flux Density
datanw = open("featNW.txt")
datase = open("featSE.txt")

# take columns to individual lists
velnw = []
fluxnw = []
velse = []
fluxse = []

for column in (raw.strip().split() for raw in datanw):
    velnw.append(column[0])
    fluxnw.append(column[1])
    
for column in (raw.strip().split() for raw in datase):
    velse.append(column[0])
    fluxse.append(column[1])

# make floats
fluxnw = [float(elem) for elem in fluxnw]
fluxse = [float(elem) for elem in fluxse]

# plot
fig = plt.figure(figsize=(20, 6))
# plot window design
vw = 0.45 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

# left

nw = plt.subplot(121)
plt.xticks(np.arange(1250, 1850, 50), fontsize=12)
plt.yticks(np.arange(-0.004, 0.008, 0.001), fontsize=12)
#nw.axis('equal')

# right

se = plt.subplot(122)
plt.xticks(np.arange(1250, 1850, 50), fontsize=12)
plt.yticks(np.arange(-0.004, 0.008, 0.001), fontsize=12)
#se.axis('equal')

# plot flux vs velocity
# NW and SE
nw.plot(velnw, fluxnw, color='g', linestyle='--', alpha=0.3, linewidth=2)
nw.plot(velse, fluxse, color='m', linestyle='--', alpha=0.3, linewidth=2)

# SE and NW
se.plot(velse, fluxse, color='m', linestyle='--', alpha=0.3,linewidth=2)
se.plot(velnw, fluxnw, color='g', linestyle='--', alpha=0.3, linewidth=2)

# NW feature only
vin  = 1480
vout = 1550
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velnw_feat  = velnw[begin:stop:1]
fluxnw_feat = fluxnw[begin:stop:1]
nw.plot(velnw_feat, fluxnw_feat, color='g', linestyle='-', alpha=0.9, linewidth=2)
# SE feature only
vin  = 1600
vout = 1640
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velse_feat  = velse[begin:stop:1]
fluxse_feat = fluxse[begin:stop:1]
se.plot(velse_feat, fluxse_feat, color='m', linestyle='-', alpha=0.9, linewidth=2)

# plot and axes titles
nw.set_xlabel("LSRK radio velocity (km/s)", fontsize=15)
nw.set_ylabel("Flux density (Jy)", fontsize=14)
nw.set_ylim(-0.0035, 0.006)
se.set_xlabel("LSRK radio velocity (km/s)", fontsize=15)
se.set_ylabel("", fontsize=0)
se.set_ylim(-0.0035, 0.006)

# baseline
nw.axhline(0.0, linestyle="--", color='k', alpha=0.3)
nw.grid(alpha=0.2)
se.axhline(0.0, linestyle="--", color='k', alpha=0.3)
se.grid(alpha=0.2)

# Markers for the system velocity plot
#ell = mpl.patches.Ellipse((vmean,cut), 2, 0.02, color='g', fill=True)
nw.plot(1544, 0.0, 'bo', markersize=8, label="CIG 96 LSRK central velocity")
se.plot(1544, 0.0, 'bo', markersize=8, label="CIG 96 LSRK central velocity")

# features location
nw.annotate("NW feature", (1380, 0.005), fontsize=14)
nw.annotate('', xy=(1495, 0.0044), xytext=(1460, 0.0048), arrowprops=dict(facecolor='g'))

se.annotate("SE feature", (1630, 0.0041), fontsize=14)
se.annotate('', xy=(1634, 0.0036), xytext=(1670, 0.0039), arrowprops=dict(facecolor='m'))

#nw.legend(loc=0, fontsize=11)
#se.legend(loc=0, fontsize=11)
plt.savefig("/home/prm/Desktop/paper_cig96/images/cig96_featsNWSE_intspectrum4.png", bbox_inches='tight')


# In[ ]:

plt.close()


# In[11]:

# Select working directory
go_to('HI')

# read text file with:
# x column: radio velocity [km/s]
# y column: [Jy] Flux Density
datanw = open("featNW.txt")
datase = open("featSE.txt")

# take columns to individual lists
velnw = []
fluxnw = []
velse = []
fluxse = []

for column in (raw.strip().split() for raw in datanw):
    velnw.append(column[0])
    fluxnw.append(column[1])
    
for column in (raw.strip().split() for raw in datase):
    velse.append(column[0])
    fluxse.append(column[1])

# make floats
fluxnw = [float(elem) for elem in fluxnw]
fluxse = [float(elem) for elem in fluxse]

# plot
fig = plt.figure(figsize=(8, 6))
# plot window design
vw = 0.45 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

nw = plt.subplot(111)

# plot flux vs velocity
# NW and SE
nw.plot(velnw, fluxnw, color='g', linestyle='--', alpha=0.4, linewidth=2)
nw.plot(velse, fluxse, color='m', linestyle='--', alpha=0.4, linewidth=2)
# NW feature only
vin  = 1480
vout = 1550
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velnw_feat  = velnw[begin:stop:1]
fluxnw_feat = fluxnw[begin:stop:1]
nw.plot(velnw_feat, fluxnw_feat, color='g', linestyle='-', alpha=0.9, linewidth=2)
# SE feature only
vin  = 1600
vout = 1640
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velse_feat  = velse[begin:stop:1]
fluxse_feat = fluxse[begin:stop:1]
nw.plot(velse_feat, fluxse_feat, color='m', linestyle='-', alpha=0.9, linewidth=2)

# plot and axes titles
plt.xticks(np.arange(1250, 1850, 50), fontsize=12)
plt.yticks(np.arange(-0.004, 0.006, 0.001), fontsize=12)
nw.set_xlabel("LSRK radio velocity (km/s)", fontsize=14)
nw.set_ylabel("Flux density (Jy)", fontsize=14)

nw.axhline(0.0, linestyle="--", color='k', alpha=0.3)
nw.grid(alpha=0.2)

# Markers for the system velocity plot
nw.plot(1544, 0.0, 'bo', markersize=8, label="CIG 96 LSRK central velocity")

# features location
nw.annotate("NW feature", (1405, 0.005), fontsize=12)
nw.annotate('', xy=(1495, 0.0044), xytext=(1460, 0.0048), arrowprops=dict(facecolor='g'))

nw.annotate("SE feature", (1630, 0.0041), fontsize=12)
nw.annotate('', xy=(1634, 0.0036), xytext=(1670, 0.0039), arrowprops=dict(facecolor='m'))

# legend
nw.legend(loc=4, fontsize=10)

# save figure
plt.savefig("/home/prm/Desktop/paper_cig96/images/cig96_featsNWSE_intspectrum1.png", bbox_inches='tight')


# In[12]:

plt.close()


# In[13]:

# Select working directory
go_to('HI')

# read text file with:
# x column: radio velocity [km/s]
# y column: [Jy] Flux Density
datanw = open("featNW.txt")
datase = open("featSE.txt")

# take columns to individual lists
velnw = []
fluxnw = []
velse = []
fluxse = []

for column in (raw.strip().split() for raw in datanw):
    velnw.append(column[0])
    fluxnw.append(column[1])
    
for column in (raw.strip().split() for raw in datase):
    velse.append(column[0])
    fluxse.append(column[1])

# make floats
fluxnw = [float(elem) for elem in fluxnw]
fluxse = [float(elem) for elem in fluxse]

# plot
fig = plt.figure(figsize=(8, 6))
# plot window design
vw = 0.45 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

nw = plt.subplot(111)

# plot flux vs velocity
# NW full
nw.plot(velnw, fluxnw, color='g', linestyle='--', alpha=0.4, linewidth=2)
plt.xticks(np.arange(1250, 1850, 50), fontsize=12)
plt.yticks(np.arange(-0.004, 0.008, 0.001), fontsize=12)

# NW feature only
vin  = 1480
vout = 1550
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velnw_feat  = velnw[begin:stop:1]
fluxnw_feat = fluxnw[begin:stop:1]
nw.plot(velnw_feat, fluxnw_feat, color='g', linestyle='-', alpha=0.8, linewidth=2)

# plot and axes titles
nw.set_xlabel("LSRK radio velocity (km/s)", fontsize=14)
nw.set_ylabel("Flux density (Jy)", fontsize=14)

nw.axhline(0.0, linestyle="--", color='k', alpha=0.3)
nw.grid(alpha=0.2)

# Markers for the system velocity plot
nw.plot(1544, 0.0, 'bo', markersize=8, label="CIG 96 LSRK central velocity")

# features location
nw.annotate("NW feature", (1405, 0.005), fontsize=12)
nw.annotate('', xy=(1495, 0.0044), xytext=(1460, 0.0048), arrowprops=dict(facecolor='g'))

# legend
nw.legend(loc=4, fontsize=10)

# save figure
plt.savefig("/home/prm/Desktop/paper_cig96/images/cig96_featNW_intspectrum.png", bbox_inches='tight')


# In[14]:

plt.close()


# In[15]:

# Select working directory
go_to('HI')

# read text file with:
# x column: radio velocity [km/s]
# y column: [Jy] Flux Density
datanw = open("featNW.txt")
datase = open("featSE.txt")

# take columns to individual lists
velnw = []
fluxnw = []
velse = []
fluxse = []

for column in (raw.strip().split() for raw in datanw):
    velnw.append(column[0])
    fluxnw.append(column[1])
    
for column in (raw.strip().split() for raw in datase):
    velse.append(column[0])
    fluxse.append(column[1])

# make floats
fluxnw = [float(elem) for elem in fluxnw]
fluxse = [float(elem) for elem in fluxse]

# plot
fig = plt.figure(figsize=(8, 6))
# plot window design
vw = 0.45 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

se = plt.subplot(111)

# plot flux vs velocity
# SE
se.plot(velse, fluxse, color='m', linestyle='--', alpha=0.4, linewidth=2)

# SE feature only
vin  = 1600
vout = 1640
begin = (vin - 1330)/10
stop  = (vout - 1330)/10 + 1
velse_feat  = velse[begin:stop:1]
fluxse_feat = fluxse[begin:stop:1]
se.plot(velse_feat, fluxse_feat, color='m', linestyle='-', alpha=0.8, linewidth=2)

# plot and axes titles
plt.xticks(np.arange(1250, 1850, 50), fontsize=12)
plt.yticks(np.arange(-0.004, 0.006, 0.001), fontsize=12)
se.set_xlabel("LSRK radio velocity (km/s)", fontsize=14)
se.set_ylabel("Flux density (Jy)", fontsize=14)

se.axhline(0.0, linestyle="--", color='k', alpha=0.3)
se.grid(alpha=0.2)

# Markers for the system velocity plot
se.plot(1544, 0.0, 'bo', markersize=8, label="CIG 96 LSRK central velocity")

# features location
se.annotate("SE feature", (1630, 0.0041), fontsize=12)
se.annotate('', xy=(1634, 0.0036), xytext=(1670, 0.0039), arrowprops=dict(facecolor='m'))

# legend
se.legend(loc=4, fontsize=10)

# save figure
plt.savefig("/home/prm/Desktop/paper_cig96/images/cig96_featSE_intspectrum.png", bbox_inches='tight')


# In[16]:

plt.close()


# #2. Asymmetry recalculation

# In[94]:

# Select working directory
go_to('HI')

data = open("cig96-circ-radiovelLSRK-intflux.txt")

# take columns to individual lists
vel = []
flux = []

for column in (raw.strip().split() for raw in data):
    vel.append(column[0])
    flux.append(column[1])

# make floats
vel = [float(elem) for elem in vel]
flux = [float(elem) for elem in flux]

# interpolation for better accuracy
vinterp = np.arange(min(vel),max(vel),0.5)
interpHI = np.interp(vinterp, vel, flux)

plt.clf()
plt.plot(vinterp, interpHI, 'o')
plt.show()

vmean = 1543.92

# conversion of scientific velocity and flux lists to float arrays
V = np.array(vinterp, dtype='float')
fl = np.array(interpHI, dtype='float')

# separate the two halves of the spectrum
vleft = V[V <= vmean]
vright = V[V >= vmean]

# trapezoid rule to calculate area below curve
L = np.trapz(fl[V <= vmean], x=vleft)
R = np.trapz(fl[V >= vmean], x=vright)

# asymmetry calculation
print "New asymmetry level =", L/R, "=", int(L/R*100-100), "% asymmetry level"
print "mean = ", (L - R) * 2. / (L + R)
print "L", L
print "R", R


# In[65]:

plt.close()


# #3. Optical SB contours
# 
# Dependent on the main reduction and calibration notebook: **Appendix A**.
# 
# - <a href="http://matplotlib.org/examples/pylab_examples/contour_demo.html">MatPlotLib contours</a>.
# - <a href="https://python4astronomers.github.io/plotting/aplpy.html">APLpy examples</a>.
# 
# ### Important: APLpy and A_0_0 keywords
# When working with APLpy, we need to have ONLY ONE set of dimension keywords ('A_0_0' and alike) in the fits file header. Should there be more than one set, APLpy will not plot the image.
# 
# ***How to do it ***
# Wipe the WCS and reset them, as explained in the main reduction and calibration notebook: <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96_CAHA2.2_jun16.ipynb#">**Appendix D**</a>.

                # Select working directory
go_to('opt')

image = "cig96_def_crop.fits"

im = aplpy.FITSFigure(image, figsize=(12,12), auto_refresh=True, north=True)

im.show_grayscale(vmin=49, vmax=53, stretch='linear')
im.add_colorbar()
im.colorbar.set_location('right')
#im.set_theme('publication')
#im.list_layers()

# contours
im.show_contour(image, levels=[50.1, 50.3, 50.5, 50.7, 50.912, 50.995])
                
                im.close()
                
# ###**Workaround**:
# 
# In case APLpy plot does not work or the image does not overlay correctly, we plot it again via matplotlib.

                # Select working directory
go_to('opt')

# FITS image selection
image = "cig96_def_crop.fits"
im = fits.open(image)
data = im[0].data
levels = [50.0,50.1,50.2,50.3,50.4,50.5,50.6,50.7,50.8,50.9,51]

plt.imshow(data.T, cmap='rainbow', vmin=49.5, vmax=58, origin='lower')
plt.gca().invert_xaxis()
plt.colorbar()
plt.contour(data.T,levels)
plt.show()
                
# ##3.1 Converting CIG96 from flux to SB

# In[429]:

# Select working directory
go_to('peris')

# read the image
image = fits.open('cig96_def_crop_wcs.fits', mode='update')
#image.info()

# select layer (in this case, it is unique)
im = image[0]
im.shape

# transform data to SB in magnitudes/arcsec2
im.data = 24.2823 + 0.5525 - 2.5*np.log10(im.data - 50.5017)

im.writeto('cig96_def_SB.fits', clobber=True)


# In[438]:

# SB conversion to counts (see reduction notebook, Appendix A)
SB = [20, 20.5, 21,22,23,24,25,26,27,28]
counts = [round(50.5017 + 10**((elem - 24.2823 - 0.5525)/(-2.5)),4) for elem in SB]
SB = [round(elem,2) for elem in SB]

print SB
print counts


# # 4. Sharpened image: large optical structures (smoothed image subtracted)
# 
# Image sharpened after subtracting the n-kernel-gaussian-smoothed image.
# 
# Dependent on the main reduction and calibration notebook: <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96_CAHA2.2_jun16.ipynb#">**Appendix E**</a>.

                # Select working directory
go_to('opt')

# Select image
image = "cig96_def_crop_smooth_gauss20.fits"

im = aplpy.FITSFigure(image, north=True, figsize=(12,12), auto_refresh=True)

#im.show_contour(image, levels=[-1.0, 0.0])
im.show_grayscale(vmin=-0.2, vmax=0.5, invert=False, stretch='linear')
#im.show_regions('ring-and-north-feat.reg')
im.show_regions('north-feat-wcs.reg')
#im.add_colorbar()
#im.colorbar.set_location('top')
                
                im.close()
                
# # 5. Deprojection and structures tracing
# 
# ##5.1 Deprojection visualization

                # Select working directory
go_to('opt')

# Select image
image = "cig96_dep.fits"

im = aplpy.FITSFigure(image, auto_refresh=True)

#im.show_contour(image, levels=[-1.0, 0.0])
im.show_grayscale(vmin=35, vmax=39, invert=True, stretch='log')
#im.show_regions('ring-and-north-feat.reg')
#im.add_colorbar()
#im.colorbar.set_location('top')
                
                im.close()
                
# ##5.2. Tracing the structures
# 
# ###5.2.1 Arms and internal features in the deprojected image
# 
# We trace the South and North structures by manually saving some positions via **DS9 regions** file as follows:
# 
# - format: ***XY***
# - coord: ***image***
# 
# Then we flip and mirror the Southern ones and overlay with the Northern ones to compare the spatial evolution.

                # Select working directory
go_to('opt')

# read positions text file:
dataN = open("arms_North_ONLY_DEPROJ_XYimage.reg")
dataS = open("arms_South_ONLY_DEPROJ_XYimage.reg")
dataC = open("arms_Center_ONLY_DEPROJ_XYimage.reg")
dataSt = open("arms_Star_ONLY_DEPROJ_XYimage.reg")

# center and star coordinates (in pixels)
cx = []
cy = []

for column in (raw.strip().split() for raw in dataC):
    cx.append(float(column[0])) # in pixels
    cy.append(float(column[1])) # in pixels
    
# center redefinition
xc = float(cx[0])
yc = float(cy[0])

# bright Eastern star
stx = []
sty = []

for column in (raw.strip().split() for raw in dataSt):
    stx.append((float(column[0]) - xc)*137.86/193.3) # in pixels
    sty.append((float(column[1]) - yc)*243.95/217.4) # in pixels
    
# star redefinition
xst = float(stx[0])
yst = float(sty[0])

# North coordinates (in pixels) to list, using center as reference
nx = []
ny = []

for column in (raw.strip().split() for raw in dataN):
    nx.append((float(column[0]) - xc)*137.86/193.3) # converted to arcsec
    ny.append((float(column[1]) - yc)*243.95/217.4) # converted to arcsec

# South coordinates (in pixels) to list, using center as reference
sx = []
sy = []

for column in (raw.strip().split() for raw in dataS):
    sx.append(float(column[0]) + xc) # in pixels
    sy.append(float(column[1]) + yc) # in pixels
    
# invert South region points and refer to center and correct to arcsec
sxinv = []
syinv = []

for coord in sx:
    sx_inv = (coord - 2*(coord - xc))*137.86/193.3 # converted to arcsec
    sxinv.append(sx_inv)

for coord in sy:
    sy_inv = (coord + 2*(yc - coord))*243.95/217.4 # converted to arcsec
    syinv.append(sy_inv)
    
# NOTE: on pix to arcsec conversion:
# x-axis: 137.86/193.3 = 137.86 arcsec spread in the original image through the 193.3 pix of the deprojected image
# y-axis: 243.95/217.4 = 243.95 arcsec spread in the original image through the 217.4 pix of the deprojected image
                
                # figure structure
fig, ax = plt.subplots()

# plot and axes titles
plt.title("Internal structures tracing - Deprojected image of CIG96", fontsize=14)
plt.axis("equal")
plt.xlabel("Distance from center (arcsec)", fontsize=13)
plt.xlim(-100,120)
plt.xticks(np.arange(-100, 120, 10), size=10, rotation=-40)
plt.ylabel("Distance from center (arcsec)", fontsize=13)
plt.yticks(np.arange(-150, 150, 10), size=8)
plt.ylim(-40,140)

center = ax.plot(xc-xc, yc-yc, 'yH', markersize=12, label='CIG96 center')
north = ax.plot(nx, ny, 'ro', label='Northern arms')
south = ax.plot(sxinv, syinv, 'bs', label='Southern arms')
star = ax.plot(xst, yst, 'r*', markersize=14, label='Star (as seen by Northern arms)')
star2 = ax.plot(-xst, -yst, 'b*', markersize=14, label='Star (as seen by South arms)')

fig.subplots_adjust(right=0.92)
plt.legend(numpoints=1, loc=3, fontsize=12, fancybox=True, bbox_to_anchor=(0.55, 0.72))
plt.grid(alpha=0.4)
plt.show()
# plot flux vs velocity
#ax.plot(vel, flux, color='b', linestyle='-', linewidth=2)
                
                plt.close()
                
# Inverting South and North regions files:

                # Select working directory
go_to('opt')

# read positions text file:
dataN = open("arms_North_ONLY_DEPROJ_XYimage.reg")
dataS = open("arms_South_ONLY_DEPROJ_XYimage.reg")
dataC = open("arms_Center_ONLY_DEPROJ_XYimage.reg")

# center and star coordinates (in pixels)
cx = []
cy = []

for column in (raw.strip().split() for raw in dataC):
    cx.append(float(column[0])) # in pixels
    cy.append(float(column[1])) # in pixels
    
# center redefinition
xc = float(cx[0])
yc = float(cy[0])

# flip North coordinates to plot inverted N-S markers on fits image
nx_flip = []
ny_flip = []

for col in (raw.strip().split() for raw in dataN):
    nx_flip.append(float(col[0])*(-1) + 2*xc) # in pixels
    ny_flip.append(float(col[1])*(-1) + 2*yc) # in pixels

# open new regions file
file = open('arms_North_INVERTED_XYimage.reg', 'w+')

nx_array = np.array(nx_flip)
ny_array = np.array(ny_flip)
data = np.array([nx_array, ny_array])

# here you transpose your data, so to have it in two columns
data = data.T

# save flipped North and South coordinates to regions file and close it
np.savetxt('arms_North_INVERTED_XYimage.reg', data, fmt=['%f','%f'])
file.close()


# flip South coordinates to plot inverted N-S markers on fits image
sx_flip = []
sy_flip = []

for col in (raw.strip().split() for raw in dataS):
    sx_flip.append(float(col[0])*(-1) + 2*xc) # in pixels
    sy_flip.append(float(col[1])*(-1) + 2*yc) # in pixels

# open new regions file
fileS = open('arms_South_INVERTED_XYimage.reg', 'w+')

sx_array = np.array(sx_flip)
sy_array = np.array(sy_flip)
datas = np.array([sx_array, sy_array])

# here you transpose your data, so to have it in two columns
datas = datas.T

# save flipped North and South coordinates to regions file and close it
np.savetxt('arms_South_INVERTED_XYimage.reg', datas, fmt=['%f','%f'])
fileS.close()
                
# New regions files (inverted North and inverted South) can be seen in ds9 over the deprojected image and fine tuned if necessary:
# - arms_North_INVERTED_XYimage.reg
# - arms_South_INVERTED_XYimage.reg
# 
# ###5.2.2 Pseudoring
# 
# We trace the pseudoring structure by manually saving some positions via **DS9 regions** file as follows:
# 
# - format: ***XY***
# - coord: ***image***

                # Select working directory
go_to('opt')

# read positions text file:
dataR = open("pseudoring_ONLY_DEPROJ_XYimage.reg")

rx = []
ry = []

for column in (raw.strip().split() for raw in dataR):
    rx.append(float(column[0]) - xc) # in pixels
    ry.append(float(column[1]) - yc) # in pixels
    
    
# center already defined in 5.2.1 as
#print "Center XY-image coordinates = ", xc, yc

# figure structure
fig, ax = plt.subplots()

ax.plot(rx, ry,'gd')
ax.axis('equal')

from numpy.linalg import eig, inv

rx = np.array(rx, dtype='float')
ry = np.array(ry, dtype='float')

fitEllipse(rx, ry)
                
                plt.close()
                
                # Select working directory
go_to('opt')

image = "cig96_dep.fits"

im = aplpy.FITSFigure(image, auto_refresh=True)#, north=True)

im.show_grayscale(vmin=35, vmax=45, stretch='linear', invert=True)
im.add_colorbar()
im.colorbar.set_location('right')
im.show_regions('pseudoring_ONLY_DEPROJ_wcs.reg')
                
                im.close()
                
# ###5.2.3 Arms and internal features in the non-deprojected image
# 
# In this case we will save the DS9 regions file as follows:
# 
# - format: ***ds9***
# - coord: ***fk5***

                # Select working directory
go_to('opt')

image = "cig96_def_crop.fits"

im = aplpy.FITSFigure(image, north=True, auto_refresh=True)

im.show_grayscale(vmin=49, vmax=100, stretch='log')#, invert=True)
im.add_colorbar()
im.colorbar.set_location('right')
im.show_regions('arms_ALL_def_wcs.reg')
#im.set_theme('publication')
#arms = np.loadtxt('arms_ALL_def_wcs.reg', dtype='float')
#xarms, yarms = arms[:, 0], arms[:, 1]
#im.show_markers(xarms, yarms, edgecolor='green', facecolor='none', marker='o', s=10, alpha=0.5)
                
                im.close()
                
# #6. R,G,B full calibration, PSF matching and combination with R deep
# 
# See <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96%20-%20RGB%20calibration%2C%20PSF%20matching%20and%20combination.ipynb">CIG96 - RGB calibration, PSF matching and combination</a> notebook.
# 
# R,G,B images alignment done via <a href="http://obswww.unige.ch/~tewes/alipy/">Alipy</a> module in Pique's notebook: <a href="http://amiga.iaa.es:9999/ipython/notebooks/pique/AlignFITS/AlignFITS-PabloApril.ipynb">Align FITS</a>.
# 
# ##6.1 R deep + R, G, B combined and final image

# In[20]:

# Select working directory
go_to('peris')

# select and open datacube image
image = "cig96_RdRGB_SB.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
cig96 = aplpy.FITSFigure(im, auto_refresh=True, north=True)
cig96.show_colorscale(stretch='linear', cmap='bone_r', vmin=24, vmax=28) # YlGnBu_r
cig96.set_nan_color('white')
cig96.recenter(RA, DEC, radius=rad)

# contour levels
#cig96.show_contour(levels=[21, 22, 23, 24, 25, 26, 26.5], cmap='RdYlBu_r', smooth=3)  # NOT valid since the combination is made with different bands

#save
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_allCAHA.png', bbox_inches='tight')


# In[21]:

cig96.close()


# In[33]:

phys = 10             # kpc
linea = phys/0.098417 # arcsec
linmin = linea/60     # arcmin

print linea, " arcsec"
print linmin, " arcmin"


# In[19]:

# Select working directory
go_to('peris')

# select and open datacube image
image = "cig96_RdRGB_SB_oriented.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
cig96 = aplpy.FITSFigure(im, auto_refresh=True, north=True)
cig96.show_colorscale(stretch='linear', cmap='bone_r', vmin=24, vmax=27.8) # YlGnBu_r
cig96.set_nan_color('white')
cig96.recenter(RA, DEC, radius=rad)

# features
feats = "cig96_CAHA_feats.reg"   # ds9 regions must be in FK5 + degrees
cig96.show_regions(feats)
plt.annotate('', xy=(280, 530), xytext=(220, 500), arrowprops=dict(facecolor='y'))

# scale bar
cig96.add_scalebar(0.0254, str(phys)+' kpc', color='orange', linewidth=2, fontsize=16)

# contour levels
#cig96.show_contour(levels=[21, 22, 23, 24, 25, 26, 26.5], cmap='RdYlBu_r', smooth=3)  # NOT valid since the combination is made with different bands

#save
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_allCAHASDSS_marks.png', bbox_inches='tight')


# In[92]:

cig96.close()


# In[36]:

# Select working directory
go_to('peris')

# select and open datacube image
image = "cig96_RdRGB_SB_oriented.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
fig = plt.figure(figsize=(9,9))
cig96 = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True)
cig96.show_colorscale(stretch='linear', cmap='bone', vmin=24, vmax=27.8) # YlGnBu_r
cig96.set_nan_color('white')
cig96.recenter(RA, DEC, radius=rad)
# labels
cig96.set_tick_labels_font(size=16)
cig96.axis_labels.set_font(size=16)
cig96.axis_labels.set_xposition('top')
# ticks
cig96.tick_labels.set_xposition('top')
cig96.tick_labels.set_xformat('hh:mm:ss')
cig96.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')
#plt.yticks(rotation=90)

# features
feats = "cig96_CAHA_feats.reg"   # ds9 regions must be in FK5 + degrees
cig96.show_regions(feats)
plt.annotate('', xy=(280, 530), xytext=(220, 500), arrowprops=dict(facecolor='y'))

# scale bar
cig96.add_scalebar(0.0254, str(phys)+' kpc', color='blue', linewidth=2, fontsize=16)

# contour levels
#cig96.show_contour(levels=[21, 22, 23, 24, 25, 26, 26.5], cmap='RdYlBu_r', smooth=3)  # NOT valid since the combination is made with different bands

#save
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_allCAHASDSS_marks_invert.png', bbox_inches='tight')


# In[37]:

cig96.close()


# ###6.1.1 VST image to SB

# In[6]:

# Select working directory
go_to('opt')

# flux to SB
image = fits.open("cig96_VST_companion_14arcmin_cropped.fits", mode='update')
# transform data to SB in magnitudes/arcsec2
im = image[0]
im.data = -2.5*np.log10(np.abs(im.data)*0.2) - 3.4948
# write image in SB
#im.writeto("cig96_VST_companion_14arcmin_cropped_SB.fits", clobber=True)
image.close()


                # Select working directory
go_to('opt')

# select and open datacube image
image = "cig96_VST_companion_SB.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
comp = aplpy.FITSFigure(im, auto_refresh=True, north=True)
comp.show_colorscale(stretch='linear', cmap='YlGnBu_r', vmin=24, vmax=28.5)
comp.set_nan_color('white')
#cig96.recenter(RA, DEC, radius=rad)
plt.annotate('', xy=(435, 2230), xytext=(435, 2580), arrowprops=dict(facecolor='cyan', shrink=0.05))

# contour levels
comp.show_contour(levels=[27.0, 27.5, 28.0, 28.4], cmap='RdYlBu_r', smooth=3)

# save figure
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_SB.png', bbox_inches='tight')
                
# In[ ]:

plt.close()


# ###6.1.2 VST image of companion

                # Select working directory
go_to('opt')

# select and open datacube image
image  = "cig96_VST_companion_20arcmin_cropped_SB.fits"
im     = fits.open(image)
imageF = "cig96_VST_companion_20arcmin_cropped.fits"
mom0   = "cig96.w.circ.HI.Filt.image.3.5sigma.mom0.fits"

#image = "cig96_VST_coadd_2arcmin_cropped_SB.fits"   # test image
#im = fits.open(image)                               # test image
#imageF = "cig96_VST_coadd_2arcmin_cropped.fits"     # test image

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.25 # degrees I think

# plot
fig = plt.figure(figsize=(9,9))
comp = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True)
comp.show_colorscale(stretch='linear', cmap='bone', vmin=24, vmax=29)
comp.set_nan_color('white')
comp.axis_labels.set_font(size=16)
comp.set_tick_labels_font(size=16)
# ticks
comp.tick_labels.set_xformat('hh:mm:ss')
comp.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')
#plt.yticks(rotation=90)

comp.recenter(RA, DEC, radius=rad)
plt.annotate('', xy=(435, 2530), xytext=(435, 2980), arrowprops=dict(facecolor='fuchsia', shrink=0.05))

# optical contours (converts SB into flux)
smooth = 11
# contours (converts SB into flux)
#lev = [24, 25, 26, 27.0, 27.5, 28.0]
#levs = []
#for elem in lev:
#    flux_lev = 10**((-0.4)*(elem+3.4948))
#    levs.append(flux_lev)
#    levs.sort()
    
#comp.show_contour(imageF, levels=levs, smooth=smooth, cmap="RdYlBu", linewidths=2, returnlevels=True)
comp.show_contour(levels=[24, 25, 26, 27.0, 27.5, 28.0, 28.4], cmap='RdYlBu', smooth=smooth) # over the SB image

# HI faintest contours
levels_border = list(np.arange(0.0040, 0.0150, 0.0001, dtype='float'))
comp.show_contour(data=mom0, levels=levels_border, colors='purple', linewidths=0.25)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_comp_'+str(smooth)+'pix_sm_large.png', bbox_inches='tight')
                
                comp.close()
                
# In[9]:

# Select working directory
go_to('opt')

# select and open datacube image
image  = "cig96_VST_companion_14arcmin_cropped_SB.fits"
im     = fits.open(image)
imageF = "cig96_VST_companion_14arcmin_cropped.fits"
mom0   = "cig96.w.circ.HI.Filt.image.3.5sigma.mom0.fits"

#image = "cig96_VST_coadd_2arcmin_cropped_SB.fits"   # test image
#im = fits.open(image)                               # test image
#imageF = "cig96_VST_coadd_2arcmin_cropped.fits"     # test image

# CIG96 coordinates in deg
RA = 33.98216
DEC = 5.9565033
rad = 0.22

# plot
fig = plt.figure(figsize=(9,9))
comp2 = aplpy.FITSFigure(im, figure=fig, auto_refresh=True, north=True)
comp2.show_colorscale(stretch='linear', cmap='bone', vmin=24, vmax=29)
comp2.set_nan_color('white')
comp2.axis_labels.set_font(size=16)
comp2.set_tick_labels_font(size=16)
# ticks
comp2.tick_labels.set_xformat('hh:mm:ss')
comp2.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')
#plt.yticks(rotation=90)

#comp2.recenter(RA, DEC, radius=rad)
#plt.annotate('', xy=(1070, 2300), xytext=(1070, 2700), arrowprops=dict(facecolor='fuchsia', shrink=0.05))

# optical contours (converts SB into flux)
smooth = 11
# contours (converts SB into flux)
#lev = [24, 25, 26, 27.0, 27.5, 28.0]
#levs = []
#for elem in lev:
#    flux_lev = 10**((-0.4)*(elem+3.4948))
#    levs.append(flux_lev)
#    levs.sort()
    
#comp2.show_contour(imageF, levels=levs, smooth=smooth, cmap="RdYlBu", linewidths=2, returnlevels=True)
comp2.show_contour(levels=[24, 25, 26, 27.0, 27.5, 28.0, 28.4], cmap='RdYlBu', smooth=smooth) # over the SB image

# HI faintest contours
levels_border = list(np.arange(0.0040, 0.0150, 0.0001, dtype='float'))
comp2.show_contour(data=mom0, levels=levels_border, colors='purple', linewidths=0.25)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_comp28x28_'+str(smooth)+'pix_sm.png', bbox_inches='tight')


# In[11]:

comp2.close()


# ###6.1.3 VST image: HI Western feature

# In[6]:

# Select working directory
go_to('opt')

# flux to SB
image = fits.open("cig96_VST_HIfeats.fits", mode='update')
# transform data to SB in magnitudes/arcsec2
im = image[0]
im.data = -2.5*np.log10(np.abs(im.data)*0.2) - 3.4948
# write image in SB
im.writeto("cig96_VST_HIfeats_SB.fits", clobber=True)
image.close


# In[15]:

# select and open datacube image
image = "cig96_VST_HIfeats_SB.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
feat = aplpy.FITSFigure(im, auto_refresh=True, north=True)
feat.show_colorscale(stretch='log', cmap='YlGnBu_r', vmin=24, vmax=28)
feat.set_nan_color('white')
#cig96.recenter(RA, DEC, radius=rad)
plt.annotate('', xy=(2100, 1350), xytext=(2100, 1200), arrowprops=dict(facecolor='blue'))
plt.annotate('', xy=(1000, 2130), xytext=(800, 2130), arrowprops=dict(facecolor='blue'))

# contour levels
feat.show_contour(levels=[26, 27.5], cmap='RdYlBu_r', smooth=3)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_HIfeat.png', bbox_inches='tight')


# In[ ]:

plt.close()


# ###6.1.4 VST image

# In[19]:

# Select working directory
go_to('opt')

# select and open datacube image
image = "cig96_VST_coadd_6arcmin_cropped_SB.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
cig96 = aplpy.FITSFigure(im, auto_refresh=True, north=True)
cig96.axis_labels.set_font(size=18)
cig96.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')
stretch='linear'
vmin=26
vmax=29
cig96.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone')

# zoom
#cig96.recenter(RA, DEC, rad)

# colorbar and contours
cig96.add_colorbar()
cig96.colorbar.set_font(size=16)
cig96.colorbar.set_axis_label_text("mag arcsec$^{-2}$")
cig96.colorbar.set_axis_label_font(size=14)
lev = [19, 20, 21, 22, 23, 24, 25.0, 25.5, 26.0, 26.5, 27.0, 27.5, 28.0, 28.5]
#cig96.show_contour(levels=lev, smooth=5, cmap="spectral_r", linewidths=2, returnlevels=True)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_SDSS.png', bbox_inches='tight')


# In[20]:

cig96.close()


# In[6]:

# Select working directory
go_to('opt')

# select and open datacube image
image = "cig96_VST_coadd_6arcmin_cropped_SB.fits"
im = fits.open(image)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.08

# plot
cig96 = aplpy.FITSFigure(im, auto_refresh=True, north=True)
cig96.axis_labels.set_font(size=18)
cig96.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')
stretch='linear'
vmin=24
vmax=29
cig96.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone')

# zoom
cig96.recenter(RA, DEC, rad)

# colorbar and contours
cig96.add_colorbar()
cig96.colorbar.set_font(size=16)
cig96.colorbar.set_axis_label_text("mag arcsec$^{-2}$")
cig96.colorbar.set_axis_label_font(size=14)
lev = [20.2, 21, 22, 23, 24, 25.0, 25.5, 26.0, 26.5, 27.0, 27.5, 28.0, 28.4]
cig96.show_contour(levels=lev, smooth=7, cmap="spectral_r", linewidths=2, returnlevels=True)

# save figure
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_'+stretch+'_'+str(vmin)+'_'+str(vmax)+'.png', bbox_inches='tight')


# In[5]:

cig96.close()


# In[39]:

# Select working directory
go_to('opt')

# select and open datacube image
imageSB = "cig96_VST_coadd_10arcmin_cropped_SB.fits"
imageF = "cig96_VST_coadd_10arcmin_cropped.fits"
imSB = fits.open(imageSB)
#imF = fits.open(imageF)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
fig = plt.figure(figsize=(9,9))
cig96 = aplpy.FITSFigure(imSB, figure=fig, auto_refresh=True, north=True)
cig96.axis_labels.set_font(size=16)
cig96.set_tick_labels_font(size=16)
#cig96.set_nan_color('white')
stretch='linear'
vmin=26
vmax=29
#vmin=-1e-12
#vmax=8e-12
#vmid=-3e-12
cig96.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone')

# zoom
cig96.recenter(RA, DEC, rad)

# ticks
cig96.ticks.show()
cig96.axis_labels.set_xposition('top')
cig96.tick_labels.set_xposition('top')
cig96.tick_labels.set_xformat('hh:mm:ss')
cig96.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', bottom='on', top='on', left='on', right='on', length=15, width=1.5, color='black')
plt.tick_params(axis=u'both', which='minor', direction='in', bottom='on', top='on', left='on', right='on', length=6, width=1.5, color='black')
#plt.yticks(rotation=90)

# colorbar
cig96.add_colorbar()
cig96.colorbar.set_font(size=17)
cig96.colorbar.set_axis_label_text("$\mu_{r}\ $(mag arcsec$^{-2}$)")
cig96.colorbar.set_axis_label_font(size=17)
cig96.colorbar.set_location(location="bottom")

# scale bar
cig96.add_scalebar(0.024, '10 kpc', color='blue', linewidth=2, fontsize=18)

# contours (converts SB into flux)
lev_num = 26.5
lev = [26.5]
levs = []

for elem in lev:
    flux_lev = 10**((-0.4)*(elem+3.5))
    levs.append(flux_lev)
    levs.sort()
smooth = 3    
cig96.show_contour(imageF, levels=levs, smooth=smooth, cmap="Reds_r", linewidths=2, returnlevels=True)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_'+stretch+str(vmin)+'_'+str(vmax)+'_265SBcont'+str(smooth)+'sm.png', bbox_inches='tight')


# In[40]:

plt.close()


                # Select working directory
go_to('opt')

# select and open datacube image
imageSB = "cig96_VST_coadd_2arcmin_cropped_SB.fits" # "cig96_VST_coadd_10arcmin_cropped_SB.fits"
imageF = "cig96_VST_coadd_2arcmin_cropped.fits" # "cig96_VST_coadd_10arcmin_cropped.fits"
imSB = fits.open(imageSB)
#imF = fits.open(imageF)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# plot
cig96 = aplpy.FITSFigure(imSB, auto_refresh=True, north=True)
cig96.axis_labels.set_font(size=18)
cig96.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')
stretch='linear'
vmin=26
vmax=29
#vmin=-1e-12
#vmax=8e-12
#vmid=-3e-12
cig96.show_colorscale(stretch=stretch, vmin=vmin, vmax=vmax, cmap='bone')

# zoom
#cig96.recenter(RA, DEC, rad)

# colorbar
cig96.add_colorbar()
cig96.colorbar.set_font(size=16)
cig96.colorbar.set_axis_label_text("mag arcsec$^{-2}$")
cig96.colorbar.set_axis_label_font(size=14)

# scale bar
cig96.add_scalebar(0.0254, '10 kpc', color='blue', linewidth=2, fontsize=18)

# contours: indicate SB and they'll be converted to flux
#lev = [20, 21, 22, 23, 24, 25.0, 25.5, 26.0, 26.5, 27.0, 27.5, 28.0, 28.4]
lev = [28.5]
levs = []
for elem in lev:
    flux_lev = 10**((-0.4)*(elem+3.5))
    levs.append(flux_lev)
    levs.sort()
    
#lev = [1.5848932e-13]    # 28.5 mag arcsec2
#lev = [1e-12]            # 26.5 mag arcsec2
cig96.show_contour(imageF, levels=levs, smooth=3, cmap="Reds_r", linewidths=2, returnlevels=True)

# save figure
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_'+stretch+str(vmin)+'_'+str(vmax)+'.png', bbox_inches='tight')
                
                cig96.close()
                
# ###6.1.3 Gaussian kernel over VST CIG96 image

# In[6]:

# Select working directory
go_to('opt')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.10

# select and open datacube image
image = "cig96_VST_coadd_20_cropped.fits"
im = fits.open(image, mode='update')
data = im[0].data

# smoothing
for i in [3,5]:
    im[0].data = scipy.ndimage.filters.gaussian_filter(data, i, mode='nearest')
    im.writeto(image[:-5]+"_"+str(i)+"sm.fits", clobber=True)


                # Select working directory
go_to('opt')

# convert to SB
image = "cig96_VST_coadd_20_cropped_3sm.fits"
im = fits.open(image, mode='update')
im[0].data = -2.5*np.log10(im[0].data/5) - 3.4948
im.writeto(image[:-5]+"_SB.fits", clobber=True)

# convert to SB
image = "cig96_VST_coadd_20_cropped.fits"
im = fits.open(image, mode='update')
im[0].data = -2.5*np.log10(im[0].data/5) - 3.4948
im.writeto(image[:-5]+"_SB.fits", clobber=True)

# convert to SB
image = "cig96_VST_coadd_10arcmin_cropped.fits"
im = fits.open(image, mode='update')
im[0].data = -2.5*np.log10(im[0].data/5) - 3.4948
im.writeto(image[:-5]+"_SB.fits", clobber=True)

# convert to SB
image = "cig96_VST_coadd_7arcmin_cropped.fits"
im = fits.open(image, mode='update')
im[0].data = -2.5*np.log10(im[0].data/5) - 3.4948
im.writeto(image[:-5]+"_SB.fits", clobber=True)
                
# In[82]:

phys = 25             # kpc
linea = phys/0.098417 # arcsec
linmin = linea/60     # arcmin

print linea, " arcsec"
print linmin, " arcmin"


# In[84]:

# select and open datacube image
image = "cig96_VST_coadd_20_cropped_3sm_SB.fits"
#image = "cig96_RdRGB_SB.fits"
im = fits.open(image)

# plot
cmapcol = 'bone'
cirrus = aplpy.FITSFigure(im, auto_refresh=True, north=True)
cirrus.show_colorscale(stretch='log', cmap=cmapcol, vmin=28, vmax=32)
cirrus.set_axis_labels_size(14)
if cmapcol[:-2] == "_r":
    cirrus.set_nan_color('black')
else:
    cirrus.set_nan_color('white')
cirrus.add_colorbar()
cirrus.colorbar.set_axis_label_text("$\mu_{r}\ (mag\ arcsec^{-2})$")
cirrus.colorbar.set_axis_label_font(size=14)

# scale bar
cirrus.add_scalebar(0.06, '25 kpc', color='orange', linewidth=2, fontsize=16)

# contour levels
#cirrus.show_contour(levels=[22000, 310000], stretch='log', cmap='RdYlBu_r', smooth=9)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_cirrus_VST_3sm.png', bbox_inches='tight')


# In[85]:

plt.close()


# ## 6.2 RGB images
# 
# <a href="http://aplpy.readthedocs.io/en/stable/howto_rgb.html">APLpy can create RGB files</a> (FITS and png) easily.
# 
# ### R,G,B Peris
# 
# In case a certain band dominates, We may fine tune the image composition by individually modifying the percentiles or values for each band.

# In[40]:

# Select working directory
go_to('peris')

# RGB creation
images = ["CIG96_R_SB.fits", "CIG96_G_SB.fits", "CIG96_B_SB.fits"]
fitsRGB = "CIG96_RGB_cube.fits"
imRGB   = "cig96_RGB_SB.png"

# creates cube
aplpy.make_rgb_cube(images, fitsRGB)

# creates RGB png images
vmin = 27   # Percentile values used to determine the minimum pixel value to use for R, G and B
vmax = 24   # Percentile values used to determine the MAXIMUM pixel value to use for R, G and B
scale = "log"    # stretch type: 'linear', 'log', 'sqrt', 'arcsinh'

imRGB = "cig96_RGB_SB_" + scale + "_" + str(vmax) + "_" + str(vmin) + ".png"   # outfile name

aplpy.make_rgb_image(fitsRGB, imRGB, stretch_r=scale, stretch_g=scale, stretch_b=scale,
                     vmax_r=vmax, vmax_g=vmax, vmax_b=vmax,
                     vmin_r=vmin, vmin_g=vmin, vmin_b=vmin)

plt.imshow(matplotlib.image.imread(imRGB))


# In[41]:

plt.close()


# ## 6.3 Flux ratios: R over G, G over B, B over R and all inverted
# 
# After dividing the fluxes, we plot the resulting images.
# 
# ***NOTE***: *CIG96_RG.fits* indicates "**R over G**". So on for the rest.

                # Select working directory
go_to('peris')

## divisions made in IRAF via imarith
## sources: CIG96_R_koords-skysub-norm.fits, CIG96_G_koords-skysub-norm.fits, CIG96_B_koords-skysub-norm.fits

# figure structure
fig = plt.figure(figsize=(12, 10))

#*subplot*: [ list of four floats ]
#    If specified, a subplot will be added at this position. The list should contain [xmin, ymin, dx, dy] 
#    where xmin and ymin are the position of the bottom left corner of the subplot, and dx and dy are the 
#    width and height of the subplot respectively. These should all be given in units of the figure width and height.
#    For example, [0.1, 0.1, 0.8, 0.8] will almost fill the entire figure, leaving a 10 percent margin on all sides.

# "recenter" is the subroutine to make a zoom on image: recenter (RA in deg, DEC in deg, radius in deg)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.07

# max and min pixel values to plot
valmax = 1.7
valmin = 0.5

# plot window design
vw = 0.395 # vertical width of the subplots in %
hw = 0.28 # horizontal width of the subplots in % 

# R over G
fRG = aplpy.FITSFigure('CIG96_RG.fits', figure=fig, subplot=[0.1, 0.52, hw, vw], north=True)
fRG.set_tick_labels_font(size='x-small')
fRG.set_axis_labels_font(size='small')
fRG.set_title("R over G", size=10)
fRG.show_grayscale(vmin=valmin, vmax=valmax)
fRG.recenter(RA, DEC, radius=rad)
# to hide
fRG.hide_xtick_labels()
fRG.hide_xaxis_label()

# G over R
fGR = aplpy.FITSFigure('CIG96_GR.fits', figure=fig, subplot=[0.1, 0.1, hw, vw], north=True)
fGR.set_tick_labels_font(size='x-small')
fGR.set_axis_labels_font(size='small')
fGR.set_title("G over R", size=10)
fGR.show_grayscale(vmin=valmin, vmax=valmax)
fGR.recenter(RA, DEC, radius=rad)

# G over B
fGB = aplpy.FITSFigure('CIG96_GB.fits', figure=fig, subplot=[0.4, 0.52, hw, vw], north=True)
fGB.set_tick_labels_font(size='x-small')
fGB.set_axis_labels_font(size='small')
fGB.set_title("G over B", size=10)
fGB.show_grayscale(vmin=valmin, vmax=valmax)
fGB.recenter(RA, DEC, radius=rad)
# to hide
fGB.hide_xtick_labels()
fGB.hide_xaxis_label()

# B over G
fBG = aplpy.FITSFigure('CIG96_BG.fits', figure=fig, subplot=[0.4, 0.1, hw, vw], north=True)
fBG.set_tick_labels_font(size='x-small')
fBG.set_axis_labels_font(size='small')
fBG.set_title("B over G", size=10)
fBG.show_grayscale(vmin=valmin, vmax=valmax)
fBG.recenter(RA, DEC, radius=rad)
# to hide
fGB.hide_yaxis_label()
fGB.hide_ytick_labels()
fBG.hide_yaxis_label()
fBG.hide_ytick_labels()

# B over R
fBR = aplpy.FITSFigure('CIG96_BR.fits', figure=fig, subplot=[0.7, 0.52, hw, vw], north=True)
fBR.set_tick_labels_font(size='x-small')
fBR.set_axis_labels_font(size='small')
fBR.set_title("B over R", size=10)
fBR.show_grayscale(vmin=valmin, vmax=valmax)
fBR.recenter(RA, DEC, radius=rad)
# to hide
fBR.hide_xtick_labels()
fBR.hide_xaxis_label()

# R over B
fRB = aplpy.FITSFigure('CIG96_RB.fits', figure=fig, subplot=[0.7, 0.1, hw, vw], north=True)
fRB.set_tick_labels_font(size='x-small')
fRB.set_axis_labels_font(size='small')
fRB.set_title("R over B", size=10)
fRB.show_grayscale(vmin=valmin, vmax=valmax)
fRB.recenter(RA, DEC, radius=rad)
# to hide
fBR.hide_yaxis_label()
fBR.hide_ytick_labels()
fRB.hide_yaxis_label()
fRB.hide_ytick_labels()

fig.canvas.draw()
                
# In[ ]:

plt.close()


# # 7. Pseudoring elliptical fit
# 
# Source: <a href="http://stackoverflow.com/questions/13635528/fit-a-ellipse-in-python-given-a-set-of-points-xi-xi-yi">Fit an ellipse in python given a set of x,y points</a>
# 
# ##7.1 Points (x,y) extraction via IRAF

                In IRAF:

imexam - "x" to extract x,y coordinates data into text file
    - folder: "deep_peris"

x,y resulting file: pseudoring-xy-----.dat

Works in ds9 as an xy-image region.
                
# ##7.2 Elliptical fittings
# 
# ###7.2.1 Fitting over the normal image
# 
# *fitEllipse* functions definitions:

# In[59]:

# Data selected
# image:            cig96_RGB.fits

# data files
x0, y0 = np.loadtxt('pseudoring-xy.dat', unpack=True, usecols=(0,1))
x1, y1 = np.loadtxt('pseudoring-xy-RdeepRGB.dat', unpack=True, usecols=(0,1))

# ellipses fitting
f0 = fitEllipse(x0, y0)
center0 = ellipse_center(f0)
phi0 = ellipse_angle_of_rotation(f0)
axes0 = ellipse_axis_length(f0)
    
print "center = ",  center0
print "angle of rotation = ",  phi0
print "axes = ", axes0

f1 = fitEllipse(x1, y1)
center1 = ellipse_center(f1)
phi1 = ellipse_angle_of_rotation(f1)
axes1 = ellipse_axis_length(f1)
    
print "center = ", center1
print "angle of rotation = ",  phi1
print "axes = ", axes1

# pi
arc = np.pi
# radius
R = np.arange(0, arc*np.pi, 0.01)

# ellipses points
a0, b0 = axes0
elx0 = center0[0] + a0*np.cos(R)*np.cos(phi0) - b0*np.sin(R)*np.sin(phi0)
ely0 = center0[1] + a0*np.cos(R)*np.sin(phi0) + b0*np.sin(R)*np.cos(phi0)

a1, b1 = axes1
elx1 = center1[0] + a1*np.cos(R)*np.cos(phi1) - b1*np.sin(R)*np.sin(phi1)
ely1 = center1[1] + a1*np.cos(R)*np.sin(phi1) + b1*np.sin(R)*np.cos(phi1)

# plot
fig = plt.figure(figsize=(12, 6))
# plot window design
vw = 0.45 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

ax0 = plt.subplot(121)
ax0.axis('equal')
ax0.plot(x0, y0, 'mo')
ax0.plot(elx0, ely0, color='red')
    
ax1 = plt.subplot(122)
ax1.axis('equal')
ax1.plot(x1, y1, 'bo')
ax1.plot(elx1, ely1, color='green')


# In[60]:

plt.close()


# ###7.2.2 Fitting over the de-projected image

# In[ ]:

# Data selected
# image:   CIG96_R_SB_dep.fits

# data files
x0, y0 = np.loadtxt('pseudoring-xy-dep.dat', unpack=True, usecols=(0,1))

# ellipses fitting
f0 = fitEllipse(x0, y0)
center0 = ellipse_center(f0)
phi0 = ellipse_angle_of_rotation(f0)
axes0 = ellipse_axis_length(f0)
    
print "center = ",  center0
print "angle of rotation = ",  phi0
print "axes = ", axes0

f1 = fitEllipse(x1, y1)
center1 = ellipse_center(f1)
phi1 = ellipse_angle_of_rotation(f1)
axes1 = ellipse_axis_length(f1)
    
print "center = ", center1
print "angle of rotation = ",  phi1
print "axes = ", axes1

# pi
arc = np.pi
# radius
R = np.arange(0, arc*np.pi, 0.01)

# ellipses points
a0, b0 = axes0
elx0 = center0[0] + a0*np.cos(R)*np.cos(phi0) - b0*np.sin(R)*np.sin(phi0)
ely0 = center0[1] + a0*np.cos(R)*np.sin(phi0) + b0*np.sin(R)*np.cos(phi0)

a1, b1 = axes1
elx1 = center1[0] + a1*np.cos(R)*np.cos(phi1) - b1*np.sin(R)*np.sin(phi1)
ely1 = center1[1] + a1*np.cos(R)*np.sin(phi1) + b1*np.sin(R)*np.cos(phi1)

# plot
fig = plt.figure(figsize=(12, 6))
# plot window design
vw = 0.45 # vertical width of the subplots in %
hw = 0.50 # horizontal width of the subplots in % 

ax0 = plt.subplot(121)
ax0.axis('equal')
ax0.plot(x0, y0, 'mo')
ax0.plot(elx0, ely0, color='red')
    
ax1 = plt.subplot(122)
ax1.axis('equal')
ax1.plot(x1, y1, 'bo')
ax1.plot(elx1, ely1, color='green')


# In[ ]:

plt.close()


# ##7.3 Elliptical fitting over pseudoring

# In[18]:

# Select working directory
go_to('peris')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.08   # zoom radius so, the larger rad, the larger FoV in the image

# max and min pixel values to plot
valmin = 0
valmax = 40

# plot window design
vw = 0.80 # vertical width of the subplots in %
hw = 0.38 # horizontal width of the subplots in % 


# RGB Peris combined + ellipse
# Ellipse fit data
center = [33.867654, 5.9998111]  # center position in degrees
dim = [0.032752181, 0.040832538] # width, height in deg
rot = 21.05                   # angle of rotation: 0.367401544695 rad = 21,050562277 deg


# In[21]:

# figure structure
fig = plt.figure(figsize=(12, 6))

fRGB = aplpy.FITSFigure('cig96_RGB_notcalibrated.fits', figure=fig, subplot=[0.1, 0.1, hw, vw], north=True)
fRGB.set_tick_labels_font(size='x-small')
fRGB.set_axis_labels_font(size='small')
fRGB.set_title("CIG96 - RGB Peris + Pseudoring fit", size=12)
fRGB.show_grayscale(vmin=valmin, vmax=valmax, invert=True)
# elliptical fit
fRGB.show_ellipses(center[0], center[1], 2.55*dim[0], 2.55*dim[1], angle=rot, linestyle=":", linewidth=4, edgecolor="magenta", facecolor="None")
fRGB.show_ellipses(center[0], center[1], 0.002, 0.002, angle=0.0, linewidth=1, edgecolor="magenta", facecolor="magenta")
# zoom
fRGB.recenter(RA, DEC, radius=rad)
# CIG96 center
fRGB.show_markers(RA, DEC, marker="x", facecolor="yellow")


# R deep + RGB Peris combined + ellipse
# Ellipse fit data
center = [33.8669, 6.0001361]
dim = [0.039267303, 0.048368206] # width, height in deg
rot = 20.95                      #  angle of rotation: 0.365746529512 rad = 20,955736872 deg

fRdRGB = aplpy.FITSFigure('cig96_RdRGB_SB.fits', figure=fig, subplot=[0.5, 0.1, hw, vw], north=True)
fRdRGB.set_tick_labels_font(size='x-small')
fRdRGB.set_axis_labels_font(size='small')
fRdRGB.set_title("CIG96 - R deep + RGB Peris + Pseudoring fit", size=12)
fRdRGB.show_grayscale(vmin=valmin, vmax=valmax, invert=True)

# elliptical fit
fRdRGB.show_ellipses(center[0], center[1], 2.1*dim[0], 2.1*dim[1], angle=rot, linestyle=":", linewidth=4, edgecolor="green", facecolor="None")
fRdRGB.show_ellipses(center[0], center[1], 0.002, 0.002, angle=0.0, linewidth=1, edgecolor="green", facecolor="green")
# zoom
fRdRGB.recenter(RA, DEC, radius=rad)
# show CIG96 center
fRdRGB.show_markers(RA, DEC, marker="x", facecolor="yellow")
# to hide
fRdRGB.hide_yaxis_label()
fRdRGB.hide_ytick_labels()

plt.legend()


# In[22]:

fRdRGB.close()
fRGB.close()
plt.close()


# #8. Disc isophotes fit and comparison with ring elliptical fit in R-Deep image
# 
# **See <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96%20-%20Isophotes.ipynb">isophotes notebook</a>.**
# 
# #9. HI moments and channels map
# 
# **See <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96%20-%20HI%20maps.ipynb">HI moments and channels map notebook</a>.**
# 
# #10. Deep optical image and HI contours

# In[9]:

# Select working directory
go_to('peris')

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.09

# open optical image
image = "cig96_RdRGB.fits"
im = fits.open(image)

# open HI image - used in contours
moment0 = "cig96-circ-HI-Filt-3.5sigma-mom0.fits"
m0 = fits.open(moment0)

# figure structure
fig = plt.figure(figsize=(14, 7))

# plot window design
vw = 0.80 # vertical width of the subplots in %
hw = 0.415 # horizontal width of the subplots in % 

# LEFT PLOT: only optical
cig96 = aplpy.FITSFigure(im, auto_refresh=True, north=True, figure=fig, subplot=[0.08, 0.1, hw, vw])
cig96.show_colorscale(stretch='log', cmap='gray_r', vmin=1, vmax=75)
cig96.set_nan_color('black')
cig96.recenter(RA-0.01, DEC-0.01, radius=rad)
cig96.set_title("CIG96 deep optical image", size=18)
cig96.add_grid()
cig96.grid.set_alpha(0.5)

# optical contours
#lev = np.linspace(12, 25, num=2)
#cig96.show_contour(im, colors="orange", linewidth=1.5, alpha=0.9, levels=lev)

# pseudoring indications
plt.annotate('Pseudoring', fontsize=14, xy=(501, 711), xytext=(681, 750), arrowprops=dict(facecolor='orange', shrink=0.05))
plt.annotate('', xy=(632, 542), xytext=(681, 750), arrowprops=dict(facecolor='orange', shrink=0.05))
plt.annotate('Pseudoring', fontsize=14, xy=(609, 395), xytext=(650, 330), arrowprops=dict(facecolor='orange', shrink=0.05))
plt.annotate('', xy=(520, 342), xytext=(650, 330), arrowprops=dict(facecolor='orange', shrink=0.05))

# legend
gray_patch = mpl.patches.Patch(color='gray', label='Faintest SB(R) = ~27.5 mag/arcsec$^2$')
orange_patch = mpl.patches.Patch(color='orange', label='Pseudoring SB(R) = 25.5 - 25.7 mag/arcsec$^2$')
plt.legend(handles=[gray_patch, orange_patch], fontsize=14)


# RIGHT PLOT: optical + HI
cig96HI = aplpy.FITSFigure(im, auto_refresh=True, north=True, figure=fig, subplot=[0.5, 0.1, hw, vw])
cig96HI.show_colorscale(stretch='log', cmap='gray_r', vmin=1, vmax=75)
cig96HI.set_nan_color('black')
cig96HI.recenter(RA-0.01, DEC-0.01, radius=rad)
cig96HI.set_title("CIG96 deep optical image + HI contours", size=18)

# HI contours
#lev = np.linspace(0.006, 0.7, num=15)
cig96HI.show_contour(m0, cmap="YlOrRd_r", linewidths=1.7, alpha=1.0, smooth=1, 
                     levels=[0.007,0.02,0.04,0.06,0.1,0.2,0.3,0.4,0.5,0.6,0.7])

# beam
cig96HI.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor="red")

#legend
red_line = mpl.patches.Patch(color='brown', label="N$_{HI}$ limit ~ 8.9 x 10$^{18}$ cm$^{-2}$")
plt.legend(handles=[red_line], fontsize=14)

# to hide
cig96HI.hide_yaxis_label()
cig96HI.hide_ytick_labels()

cig96HI.add_grid()
cig96HI.grid.set_alpha(0.5)


# In[10]:

cig96.close()
cig96HI.close()


# #11. Pseudoring over HI

# In[26]:

# Select working directory
go_to('peris')

# select and open datacube image
#HI = "cig96_HI_mom0_3.5rms.fits"
HI = "cig96_HI_mom0_3rms_coldens.fits"
im = fits.open(HI)


# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.085
#rad = 0.15

# figure
fig = plt.figure(figsize=(12, 12))

# plot
optHI = aplpy.FITSFigure(im, auto_refresh=True, figure=fig)
optHI.show_colorscale(stretch='linear', cmap='RdYlBu_r', vmin=0.0, vmax=12)
optHI.add_colorbar()#cmap='Paired', vmin=100, vmax=5000)
optHI.set_nan_color('white')
optHI.recenter(RA, DEC, radius=rad)
optHI.add_grid()
optHI.grid.set_alpha(0.2)
optHI.set_axis_labels_size(22)
optHI.set_tick_labels_font(size=18)
#optHI.set_title("CIG96 - Pseudoring footprint over HI (0$^{th}$ moment map and contours)", size=18)

# colorbar
optHI.colorbar.set_axis_label_text("$N_{HI}$ (10$^{20}$ at cm$^{-2})$")
optHI.colorbar.set_axis_label_font(size=22)
optHI.colorbar.set_font(size=20)

# beam
optHI.add_beam(28*u.arcsec, 28*u.arcsec, 0.0, facecolor='y', edgecolor='k')

# contours
#lowlev = np.linspace(1, 15, num=15)
#highlev = np.linspace(500, 4000, num=10)
#lev = np.ndarray.tolist(lowlev) + np.ndarray.tolist(highlev)
lev = np.arange(1,15, 1)
print lev
optHI.show_contour(levels=lev, cmap="gist_earth", label="HI contours", alpha=0.7)

# pseudoring region file
display_color = "black"  # green, red, black, white, orange, yellow
#reg = "pseudoring-WCS-RdeepRGB-" + display_color + ".reg"
reg = "cig96_pseudo_reg_dep_BandR_wcs_crosses.reg"
optHI.show_regions(reg)

# southern feature
feat = "cig96_south_feat.reg"
optHI.show_regions(feat)

# legend
#mpl.markers()
HIem   = mpl.patches.Patch(color=(0.2, 0.7, 0.4), label="HI emission")
HIcont = mpl.lines.Line2D([], [], color='purple', linewidth=4, linestyle='-', label='HI contours')
pseudo = mpl.lines.Line2D([], [], color='black', marker='+', markeredgecolor='black', markeredgewidth=2,
                          markerfacecolor='none', linestyle='none', markersize=12, label='Pseudoring')
#plt.legend(handles=[HIem, HIcont, pseudo], loc=0, numpoints=2)

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_HI_pseudo_marks.png', bbox_inches='tight')


# In[27]:

plt.close()


# #12. Pseudoring B - R color gradient
# 
# See notebook: <a href="http://localhost:8888/notebooks/DeSIGN/CAHA/cig96/CIG96%20-%20Pseudoring%20colors.ipynb">Pseudoring colors</a>.

# #13. Cirrus: Planck, WISE and optical images

# WISE3 image from SkyView

# In[84]:

# Select working directory
go_to('opt')

# select and open datacube image
image = "cig96_40arcmin_WISE3.fits"
im = fits.open(image)

print np.amax(im[0].data)
print np.amin(im[0].data)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.335

# plot
wise = aplpy.FITSFigure(im, auto_refresh=True, north=True)
wise.show_colorscale(stretch='log', cmap='bone', vmin=790, vmax=794)
wise.set_nan_color('white')
wise.recenter(RA, DEC, radius=rad)

# contour levels
levs = [789.9, 790.3, 790.7,  791.1, 791.2, 791.3, 792]
#wise.show_contour(levels=levs, cmap='Spectral_r')


# In[ ]:

wise.close()


# In[36]:

# Select working directory
go_to('peris')

opt = "cig96_RdRGB_SB_oriented_cropped.fits"
wise3 = "CIG96.W3.clean.fits"
planck857 = "cig96_planck857.fits"
planck857L = "cig96_planck857_1deg.fits"

# plot
fig = plt.figure(figsize=(16,6))

planckL = aplpy.FITSFigure(planck857L, figure=fig, subplot=[0.07, 0.15, 0.29, 0.77], zorder=1)
planckL.show_colorscale(cmap="bone")
planckL.show_rectangles(33.859974, 6.0017605, 0.25, 0.25, color='white', linestyle='--', linewidth=2, zorder=2)
planckL.axis_labels.set_font(size="12")
plt.tick_params(axis=u'both', which='both', direction='in', top='off', right='off')
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

planck = aplpy.FITSFigure(planck857, figure=fig, subplot=[0.41, 0.15, 0.29, 0.77])
planck.show_colorscale(cmap="bone")
planck.show_contour(opt, smooth=1, levels=[21, 22, 23, 24, 25, 26], cmap="winter", linewidths=0.8, alpha=0.8)
planck.axis_labels.set_font(size="12")
planck.axis_labels.hide_y()
plt.xticks(fontsize=8, rotation=-20)
plt.yticks(fontsize=8)
plt.tick_params(axis=u'both', which='both', direction='in', top='off', right='off')

wise = aplpy.FITSFigure(wise3, figure=fig, subplot=[0.70,  0.15, 0.29, 0.77])
wise.show_colorscale(cmap="hot", stretch='log', vmin=417, vmid=400, vmax=460)
wise.show_contour(opt, smooth=1, levels=[21, 22, 23, 24, 25, 26], cmap="winter", linewidths=0.8, alpha=0.8)
wise.axis_labels.set_font(size="12")
wise.hide_yaxis_label()
wise.tick_labels.hide_y()
plt.xticks(fontsize=8, rotation=-20)
plt.tick_params(axis=u'both', which='both', direction='in', top='off', right='off')

fig.canvas.draw()

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_planck_wise.png', bbox_inches='tight')


# In[40]:

plt.close()


# In[183]:

# Select working directory
go_to('peris')

opt = "cig96_RdRGB_SB_oriented_cropped.fits"
wise3 = "CIG96.W3.clean.fits"
planck857 = "cig96_planck857.fits"
planck857L = "cig96_planck857_1deg.fits"

# plot
fig = plt.figure(figsize=(16,6))

planckL = aplpy.FITSFigure(planck857L, figure=fig, subplot=[0.08, 0.1, 0.31, 0.84], zorder=1)
planckL.show_colorscale(cmap="bone")
planckL.show_colorbar()
planckL.colorbar.set_location(location="top")
planckL.show_rectangles(33.859974, 6.0017605, 0.25, 0.25, color='white', linestyle='--', linewidth=2, zorder=2)
planckL.axis_labels.set_font(size="12")

planck = aplpy.FITSFigure(planck857, figure=fig, subplot=[0.38, 0.1, 0.31, 0.84])
planck.show_colorscale(cmap="bone")
planck.show_colorbar()
planck.colorbar.set_location(location="top")
planck.show_contour(opt, smooth=1, levels=[20, 22, 24, 26], cmap="winter", linewidths=0.8, alpha=0.8)
planck.hide_yaxis_label()
planck.tick_labels.hide_y()
planck.axis_labels.set_font(size="12")

wise = aplpy.FITSFigure(wise3, figure=fig, subplot=[0.68, 0.1, 0.31, 0.84])
wise.show_colorscale(cmap="hot", stretch='log', vmin=417, vmid=400, vmax=460)
wise.show_colorbar()
wise.colorbar.set_location(location="top")
wise.show_contour(opt, smooth=1, levels=[20, 22, 24, 26], cmap="winter", linewidths=0.8, alpha=0.8)
wise.axis_labels.set_font(size="12")
wise.hide_yaxis_label()
wise.tick_labels.hide_y()

fig.canvas.draw()


# In[184]:

plt.close()


# ##WARNING
# The next plot takes 30 minutes to run!

                # Select working directory
go_to('opt')

opticalL = "cig96_VST_coadd_20_cropped_3sm_SB.fits"
planck857L = "cig96_planck857_20arcmin.fits"

opt = "cig96_VST_coadd_10arcmin_cropped_SB.fits"
planck857 = "cig96_planck857_12arcmin.fits"
wise3 = "CIG96.W3.clean.fits"

# plot
fig = plt.figure(figsize=(12,12))

# % widths
hw = 0.38
vw = 0.38

#### top row

optL = aplpy.FITSFigure(opticalL, figure=fig, subplot=[0.1, 0.50, hw, vw], zorder=1)
optL.show_colorscale(stretch='log', cmap="bone", vmin=28, vmax=32)
#optL.show_rectangles(33.871973, 5.99991, 0.27, 0.27, color='blue', linestyle='--', linewidth=2, zorder=2)
optL.axis_labels.set_font(size="12")
optL.hide_xaxis_label()
plt.tick_params(axis=u'both', which='both', direction='in', top='on', right='on')
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
optL.show_colorbar()
optL.colorbar.set_location(location="top")

planckL = aplpy.FITSFigure(planck857L, figure=fig, subplot=[0.55, 0.50, hw, vw], zorder=1)
planckL.show_colorscale(stretch='log', cmap="bone", vmin=1.8, vmax=2.7)
planckL.show_rectangles(33.871973, 5.99991, 0.25, 0.25, color='white', linestyle='--', linewidth=2, zorder=2)
planckL.axis_labels.set_font(size="12")
planckL.hide_xaxis_label()
planckL.hide_yaxis_label()
plt.tick_params(axis=u'both', which='both', direction='in', top='on', right='on')
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)


#### bottom row

planck = aplpy.FITSFigure(planck857, figure=fig, subplot=[0.1, 0.1, hw, vw], zorder=1)
planck.show_colorscale(cmap="bone")
planck.show_contour(opt, smooth=7, levels=[25, 26, 27, 27.5], cmap="winter", linewidths=0.8, alpha=0.8)
planck.axis_labels.set_font(size="12")
plt.yticks(fontsize=8)
plt.tick_params(axis=u'both', which='both', direction='in', top='on', right='on')
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

wise = aplpy.FITSFigure(wise3, figure=fig, subplot=[0.55, 0.1, hw, vw], zorder=1)
wise.show_colorscale(cmap="hot", stretch='log', vmin=417, vmid=400, vmax=460)
wise.show_contour(opt, smooth=7, levels=[25, 26, 27, 27.5], cmap="winter", linewidths=0.8, alpha=0.8)
wise.axis_labels.set_font(size="12")
wise.hide_yaxis_label()
plt.tick_params(axis=u'both', which='both', direction='in', top='on', right='on')
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

fig.canvas.draw()

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_planck_wise.png', bbox_inches='tight')
                
                plt.close()
                
# In[26]:

# Select working directory
go_to('cig')

opticalL = "cig96_VST_coadd_20_cropped_3sm_SB.fits"
planck857L = "cig96_planck857_20arcmin.fits"
opt = "cig96_VST_coadd_10arcmin_cropped_SB.fits"
planck857 = "cig96_planck857_12arcmin.fits"
conts = "cig96_mosaicW3_40_21sm.fits"

#opticalL = "cig96_VST_coadd_2arcmin_cropped_SB.fits"  # test image
#planck857L = "cig96_planck857_20arcmin.fits"          # test image
#opt = "cig96_VST_coadd_2arcmin_cropped_SB.fits"       # test image
#planck857 = "cig96_planck857_12arcmin.fits"           # test image
#conts = "cig96_HI_mom0_5rms_Msol.fits"                # test image

# plot
fig = plt.figure(figsize=(12,12))

# % widths
hw = 0.38
vw = 0.38

# fontsize for axes ticks
ftsize = 11

#### top row

optL = aplpy.FITSFigure(opticalL, figure=fig, subplot=[0.09, 0.50, hw+0.02, vw+0.02], zorder=1)
optL.show_colorscale(stretch='log', cmap="bone", vmin=28, vmax=32)
#optL.show_rectangles(33.871973, 5.99991, 0.27, 0.27, color='blue', linestyle='--', linewidth=2, zorder=2)
optL.axis_labels.set_font(size="14")
optL.hide_xaxis_label()
optL.tick_labels.set_xformat('hh:mm:ss')
optL.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', top='on', right='on', length=12, width=1.5, color='white')
plt.tick_params(axis=u'both', which='minor', direction='in', top='on', right='on', length=6, width=1.5, color='white')
plt.xticks(fontsize=ftsize)
plt.yticks(fontsize=ftsize)
optL.show_colorbar()
optL.colorbar.set_location(location="top")
optL.colorbar.set_font(size=14)

vstwise = aplpy.FITSFigure(opticalL, figure=fig, subplot=[0.50, 0.50, hw, vw], zorder=1)
vstwise.show_colorscale(stretch='log', cmap="bone", vmin=28, vmax=32)
#vstwise.show_contour(opt, smooth=7, levels=[25, 26, 27, 27.5], cmap="winter", linewidths=0.8, alpha=0.8)
vstwise.axis_labels.set_font(size="14")
vstwise.hide_xaxis_label()
vstwise.hide_yaxis_label()
#vstwise.axis_labels.set_yposition('right')
vstwise.tick_labels.set_yposition('right')
vstwise.tick_labels.set_xformat('hh:mm:ss')
vstwise.tick_labels.set_yformat('dd:mm')
vstwise.show_contour(conts, levels=[419.275, 419.5, 419.8, 420], filled=True, cmap="spring", linewidths=0.8, alpha=0.5)
plt.tick_params(axis=u'both', which='major', direction='in', top='on', right='on', length=12, width=1.5, color='white')
plt.tick_params(axis=u'both', which='minor', direction='in', top='on', right='on', length=6, width=1.5, color='white')
plt.xticks(fontsize=ftsize)
plt.yticks(fontsize=ftsize)

#### bottom row

planckL = aplpy.FITSFigure(planck857L, figure=fig, subplot=[0.10, 0.10, hw, vw], zorder=1)
planckL.show_colorscale(stretch='log', cmap="bone", vmin=1.8, vmax=2.7)
planckL.show_rectangles(33.871973, 5.99991, 0.25, 0.25, color='white', linestyle='--', linewidth=2, zorder=2)
planckL.axis_labels.set_font(size="14")
planckL.tick_labels.set_xformat('hh:mm:ss')
planckL.tick_labels.set_yformat('dd:mm')
plt.tick_params(axis=u'both', which='major', direction='in', top='on', right='on', length=12, width=1.5, color='white')
plt.tick_params(axis=u'both', which='minor', direction='in', top='on', right='on', length=6, width=1.5, color='white')
plt.xticks(fontsize=ftsize)
plt.yticks(fontsize=ftsize)

planck = aplpy.FITSFigure(planck857, figure=fig, subplot=[0.50, 0.10, hw, vw], zorder=1)
planck.show_colorscale(cmap="bone")
planck.show_contour(opt, smooth=7, levels=[25, 26, 27, 27.5], cmap="winter", linewidths=0.8, alpha=0.7)
planck.axis_labels.set_font(size="14")
planck.hide_yaxis_label()
planck.tick_labels.set_xformat('hh:mm:ss')
planck.tick_labels.set_yformat('dd:mm')
#planck.axis_labels.set_yposition('right')
planck.tick_labels.set_yposition('right')
plt.tick_params(axis=u'both', which='major', direction='in', top='on', right='on', length=12, width=1.5, color='white')
plt.tick_params(axis=u'both', which='minor', direction='in', top='on', right='on', length=6, width=1.5, color='white')
plt.xticks(fontsize=ftsize)
plt.yticks(fontsize=ftsize)

fig.canvas.draw()

# save figure
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_planck_wise_new.png', bbox_inches='tight')
#plt.savefig('/home/prm/Desktop/paper_cig96/images/test_test.png', bbox_inches='tight')   # test save


# In[27]:

plt.close()


# ###VST image with WISE contours

# In[ ]:

# Select working directory
go_to('opt')

# select and open images for background and contours
back  = "cig96_mosaicW3_40_13sm.fits"
im = fits.open(back)
conts = "cig96_VST_coadd_20_cropped_SB.fits"
#conts = "cig96_VST_coadd_2arcmin_cropped_SB.fits" # dummy image for testing (faster loading)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.35

# plot
cig96 = aplpy.FITSFigure(im, auto_refresh=True, north=True)
cig96.axis_labels.set_font(size=18)
cig96.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')
cig96.show_colorscale(stretch='linear', vmin=418, vmax=421, cmap='bone')

# zoom
cig96.recenter(RA, DEC, rad)

# colorbar
cig96.add_colorbar()
cig96.colorbar.set_font(size=16)
cig96.colorbar.set_axis_label_text("mag arcsec$^{-2}$")
cig96.colorbar.set_axis_label_font(size=14)

# scale bar
cig96.add_scalebar(0.031, '10 kpc', color='orange', linewidth=2, fontsize=18)

# contours
cig96.show_contour(conts, smooth=7, levels=[25, 26, 27, 27.5, 28.5, 29], cmap="winter", linewidths=0.8, alpha=0.8)

# save figure
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_'+stretch+str(vmin)+'_'+str(vmax)+'.png', bbox_inches='tight')
plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_WISE.png', bbox_inches='tight')


# In[7]:

cig96.close()


# ###VST image with WISE contours

# In[ ]:

# Select working directory
go_to('opt')

# select and open images for background and contours
back  = "cig96_VST_coadd_20_cropped_3sm_SB.fits"
#back  = "cig96_VST_coadd_2arcmin_cropped_SB.fits" # dummy image for testing (faster loading)
im = fits.open(back)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.33

# gaussian kernel pixel radius
sm = [21, 19, 17, 15, 13, 11]    # 3, 5, 7, 9, 11, 13, 15, 17, 19, 21 pixels for smoothing kernel size and available images

# plot images
for elem in sm:
    
    # plot
    cig96wise = aplpy.FITSFigure(im, auto_refresh=True, north=True)

    cig96wise.axis_labels.set_font(size=18)
    cig96wise.set_tick_labels_font(size=12)
    #cig96.set_nan_color('white')
    cig96wise.show_colorscale(stretch='log', cmap="bone", vmin=28, vmax=32)
    # zoom
    #cig96wise.recenter(RA, DEC, rad)
    # colorbar
    cig96wise.add_colorbar()
    cig96wise.colorbar.set_font(size=16)
    cig96wise.colorbar.set_axis_label_text("$\mu$ (mag arcsec$^{-2}$)")
    cig96wise.colorbar.set_axis_label_font(size=14)

    # scale bar
    cig96wise.add_scalebar(0.06, '20 kpc', color='blue', linewidth=2, fontsize=18)

    # contours
    conts = "cig96_mosaicW3_40_"+str(elem)+"sm.fits"
    #cig96wise.show_contour(conts, levels=[419.275, 450, 470], filled=True, cmap="spring", linewidths=0.8, alpha=0.5)
    #cig96wise.show_contour(conts, levels=[419.275, 419.5, 419.8, 420], filled=True, cmap="spring", linewidths=0.8, alpha=0.5)
    cig96wise.show_contour(conts, levels=[419.275, 419.5, 419.8, 420], filled=True, cmap="spring", linewidths=0.8, alpha=0.5)
    
    # save figure
    plt.savefig("/home/prm/Desktop/paper_cig96/images/cig96_VST_28_32_WISE_"+str(elem)+"sm.png", bbox_inches='tight')


# In[6]:

cig96wise.close()


# ###WISE contrast image

# In[112]:

# gaussian smoothing
def smooth(fnamein, kern_px):
    """
    Smooth the image in the first HDU of the input file. 
    kern_px: FWHM of kernel in pxiels
    """
    f1=pyfits.open(fnamein)
    f1[0].data=scipy.ndimage.gaussian_filter(f1[0].data, kern_px/(2*math.sqrt(2*math.log(2))))
    fnameout = fnamein[:-5]+"_"+str(kern_px)+"sm.fits"
    f1.writeto(fnameout, clobber=True)


# In[114]:

# Select working directory
go_to('opt')

smooth("cig96_mosaicW3_40.fits", 3)
smooth("cig96_mosaicW3_40.fits", 5)
smooth("cig96_mosaicW3_40.fits", 7)
smooth("cig96_mosaicW3_40.fits", 9)
smooth("cig96_mosaicW3_40.fits", 11)
smooth("cig96_mosaicW3_40.fits", 13)
smooth("cig96_mosaicW3_40.fits", 15)
smooth("cig96_mosaicW3_40.fits", 17)
smooth("cig96_mosaicW3_40.fits", 19)
smooth("cig96_mosaicW3_40.fits", 21)


                In IRAF via imarith:

standard image over smoothed image with gaussian kernels of 3, 5, 7, 9, 11, 13, 15, 17, 19 and 21 pixels radius

imarith:

operand1= cig96_mosaicW3_40.fits              Operand image or numerical constant
op      =                    /                Operator
operand2= cig96_mosaic_W3_40_3sm.fits        Operand image or numerical constant
result  = cig96_mosaicW3_40_st_over_3sm.fits  Resultant image

In a cl file:
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_3sm.fits result=cig96_mosaicW3_40_st_over_3sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_5sm.fits result=cig96_mosaicW3_40_st_over_5sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_7sm.fits result=cig96_mosaicW3_40_st_over_7sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_9sm.fits result=cig96_mosaicW3_40_st_over_9sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_9sm.fits result=cig96_mosaicW3_40_st_over_11sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_9sm.fits result=cig96_mosaicW3_40_st_over_13sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_9sm.fits result=cig96_mosaicW3_40_st_over_15sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_9sm.fits result=cig96_mosaicW3_40_st_over_17sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_9sm.fits result=cig96_mosaicW3_40_st_over_19sm.fits
imarith operand1=cig96_mosaicW3_40.fits op=/ operand2=cig96_mosaicW3_40_9sm.fits result=cig96_mosaicW3_40_st_over_21sm.fits
                
# In[34]:

# Select working directory
go_to('opt')

# select and open images for background and contours
back  = "cig96_mosaicW3_40_19sm.fits"
im = fits.open(back)

# CIG96 coordinates in deg
RA = 33.865231
DEC = 6.0020581
rad = 0.35

# plot
cig96 = aplpy.FITSFigure(im, auto_refresh=True, north=True)
cig96.axis_labels.set_font(size=18)
cig96.set_tick_labels_font(size=12)
#cig96.set_nan_color('white')
cig96.show_colorscale(stretch='linear', vmin=418, vmax=421, cmap='bone_r')

# zoom
cig96.recenter(RA, DEC, rad)

# colorbar
#cig96.add_colorbar()
#cig96.colorbar.set_font(size=16)
#cig96.colorbar.set_axis_label_text("mag arcsec$^{-2}$")
#cig96.colorbar.set_axis_label_font(size=14)

# scale bar
cig96.add_scalebar(0.06, '20 kpc', color='blue', linewidth=2, fontsize=18)

# contours
cig96.show_contour(levels=[419.25, 419.5, 419.8, 420], filled=True, cmap="spring", linewidths=0.8, alpha=0.5)

# save figure
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_VST_'+stretch+str(vmin)+'_'+str(vmax)+'.png', bbox_inches='tight')
#plt.savefig('/home/prm/Desktop/paper_cig96/images/cig96_cirrus_conts_WISE.png', bbox_inches='tight')


# In[32]:

cig96.close()


# In[ ]:



