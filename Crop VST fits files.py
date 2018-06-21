
# coding: utf-8

# #Crop VST fits files

# In[4]:

import pyfits as fits
import numpy as np
import os
import matplotlib.pyplot as plt

# renders interactive figures in the Notebook
get_ipython().magic(u'matplotlib nbagg')

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

#############################################################################

def crop_to_field(fitsfile, ra, dec, c_r=20.):
        '''
        Crops the noisy outskirts
        of the input image and returns
        a fits file which is cropped to
        the fds field.
        
        Example: (crops the 20*20 arcmin area in the center of field 11)
        >>> crop_to_field('Sci-Coadd.fits',fdsfield_coords[10][0],fdsfield_coords[10][1],c_r=10.)
        '''
        hdulist2 = fits.open(fitsfile)
        #parameters
        amin_to_pix    = 60.*abs((1./ (3600.*hdulist2[0].header['CD1_1'])))
        pix_to_amin    = 1./amin_to_pix
  
        NAXIS1 = hdulist2[0].header['NAXIS1']
        NAXIS2 = hdulist2[0].header['NAXIS2']
        print 'Image dimensions (x,y):',NAXIS1,NAXIS2

        #BORDERS:

        x_center = int(-np.cos(dec*2.*np.pi/360.)*(ra-hdulist2[0].header['CRVAL1'])*60.*amin_to_pix) + hdulist2[0].header['CRPIX1']
        y_center = int((dec-hdulist2[0].header['CRVAL2'])*60.*amin_to_pix + hdulist2[0].header['CRPIX2'])
        upper_border = int(y_center+c_r*amin_to_pix) 
        lower_border = int(y_center-c_r*amin_to_pix) 
        right_border = int(x_center+c_r*amin_to_pix)  # /np.cos(dec*2.*np.pi/360.)  #THIS WAS LEFT OUT TO GET SQUARES
        left_border  = int(x_center-c_r*amin_to_pix)  # /np.cos(dec*2.*np.pi/360.)
        matrix = hdulist2[0].data[lower_border:upper_border,left_border:right_border]
        print 'Field center in:',x_center,y_center,'(xpix,ypix)'
        print 'Cropped area:',left_border,right_border,'(xmin,xmax)',lower_border,upper_border,'(ymin,ymax)'
        hdulist2[0].data = matrix
        hdulist2[0].header['NAXIS2'] = upper_border-lower_border 
        hdulist2[0].header['NAXIS1'] = right_border-left_border
        hdulist2[0].header['CRPIX2'] -= lower_border
        hdulist2[0].header['CRPIX1'] -= left_border
        hdulist2.writeto(fitsfile[:-5]+'_'+str(c_r)+'arcmin_cropped.fits', clobber=True)
        
#######################################

def ping():
    os.system('spd-say "Ping"')


# In[5]:

ping()


# # VST
# ##CIG11

# In[3]:

os.chdir("/home/prm/Desktop/VST/cig11")

image = "cig11_vdef.fits"
RA = 3.633
DEC = -0.736

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# ##CIG33

# In[4]:

os.chdir("/home/prm/Desktop/VST/cig33")

image = "cig33_vdef.fits"
RA = 10.866497
DEC = -0.12441277

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# ##CIG59

# In[5]:

os.chdir("/home/prm/Desktop/VST/cig59")

image = "cig59_vdef.fits"
RA = 24.587031
DEC = 7.5346537

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# ##CIG96

# In[2]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd.fits"
RA = 33.865231
DEC = 6.0020581

crop_to_field(image, RA, DEC, 10)
crop_to_field(image, RA, DEC, 10)


# In[2]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd.fits"
RA = 33.865231
DEC = 6.0020581

crop_to_field(image, RA, DEC, 6)


# In[2]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd.fits"
RA = 33.865231
DEC = 6.0020581

crop_to_field(image, RA, DEC, 12.5)


# In[2]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd.fits"
RA = 33.89182
DEC = 6.043247

crop_to_field(image, RA, DEC, 2)


# ### - CIG96 companion

# In[3]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd_20_cropped.fits"
RA = 33.967818
DEC = 5.9838688

crop_to_field(image, RA, DEC, 10)


# In[3]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd.fits"
RA = 33.931723
DEC = 6.0175497

crop_to_field(image, RA, DEC, 14)


# In[ ]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd.fits"
RA = 33.967818
DEC = 5.9838688

crop_to_field(image, RA, DEC, 20)


# ### - CIG96 HI features

# In[4]:

os.chdir("/home/prm/Desktop/VST/cig96")
image = "cig96_VST_coadd_20_cropped.fits"
RA = 33.793461
DEC = 6.0468212

crop_to_field(image, RA, DEC, 4.5)


# ##CIG152

# In[6]:

os.chdir("/home/prm/Desktop/VST/cig152")

image = "cig152_vdef.fits"
RA = 67.748514
DEC = -2.0033558

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# ##CIG154

# In[7]:

os.chdir("/home/prm/Desktop/VST/cig154")

image = "cig154_vdef.fits"
RA = 71.9385
DEC = 1.8191

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# ##CIG279

# In[8]:

os.chdir("/home/prm/Desktop/VST/cig279")

image = "cig279_vdef.fits"
RA = 130.38419
DEC = 4.981089

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# In[7]:

os.chdir("/home/prm/Desktop/VST/cig279")

image = "cig279_vdef.fits"
RA  = 130.42283
DEC = 4.9569477

crop_to_field(image, RA, DEC, c_r=4.5)


# ##CIG568

# In[9]:

os.chdir("/home/prm/Desktop/VST/cig568")

image = "cig568_vdef.fits"
RA = 196.06271
DEC = 9.223596

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# In[8]:

os.chdir("/home/prm/Desktop/VST/cig568")

image = "cig568_vdef.fits"
RA = 196.06411
DEC = 9.2320021

crop_to_field(image, RA, DEC, c_r=3)


# ##CIG 1002

# In[10]:

os.chdir("/home/prm/Desktop/VST/cig1002")

image = "cig1002_vdef.fits"
RA = 345.17068
DEC = 8.4676962

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# ##CIG 1027

# In[11]:

os.chdir("/home/prm/Desktop/VST/cig1027")

image = "cig1027_vdef.fits"
RA = 353.859
DEC = 7.322

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# ##CIG 1047

# In[12]:

os.chdir("/home/prm/Desktop/VST/cig1047")

image = "cig1047_vdef.fits"
RA = 359.198
DEC = 1.3550 

crop_to_field(image, RA, DEC, c_r=2)
crop_to_field(image, RA, DEC, c_r=3)
#crop_to_field(image, RA, DEC, c_r=4)
#crop_to_field(image, RA, DEC, c_r=7)
#crop_to_field(image, RA, DEC, c_r=8)
#crop_to_field(image, RA, DEC, c_r=10)
#crop_to_field(image, RA, DEC, c_r=20)


# In[ ]:



