#!/usr/bin/env python

'''
CASA script that combines 2 datasets of CIG96: one from EVLA data and another from VLA data.
-
It can be adapted to any other pair of datasets.
Tasks performed: mstranform (same velocity frame), splitting of spws of interest, concatenation of ms files and clean.
'''
################################################
# HI LINE: splitting, concatenation and clean. #
################################################

# Setting same reference frame for both datasets before concatenation.
# mstransform(vis='../cig96-vla/vla.av276-282/cig96_vla.5kHz.corr.ms.contsub', outputvis= '../cig96-vla/vla.av276-282/cig96_vla.corr.ms.LSRK.contsub', regridms=True, datacolumn='data', outframe='LSRK')
# mstransform(vis='../evla/CDconf/cig96_CDconf.1kHz.corr.ms.line.contsub', outputvis= '../evla/CDconf/cig96_evla.corr.ms.line.LSRK.contsub', regridms=True, datacolumn='data', outframe='LSRK')

##### NOTE: no need to repeat this last 2 transformations. If this script must be repeated, start right BELOW this line.

# Same split command for VLA data. Constrained EVLA data data split around the mentioned channels
split(vis='../cig96-vla/vla.av276-282/cig96_vla.corr.ms.LSRK.contsub', outputvis='cig96.vla.spw_HI.ms', datacolumn='data', field='0', spw='2')
split(vis='../evla/CDconf/cig96_evla.corr.ms.line.LSRK.contsub', outputvis='cig96.evla.spw_HI.ms', datacolumn='data', field='0', spw='0~4:820~1220') # covers a bandwidth of ~1300 km/s (~650 km.s on each side of the HI line)
# 400 channels for the EVLA data.

# Concatenation. Weight calculation: w(i) = 1/rms(i)^2, where 'i' is each cube to be concatenated. The rms has been measured after making a quick HI image of each file: rms(VLA)~0.31 mJy/beam, rms(EVLA)~0.84 mJy/beam.
concat(vis=['cig96.vla.spw_HI.ms','cig96.evla.spw_HI.ms'], concatvis='cig96.all.weight.HI.ms', visweightscale=[10.40,1.42])

# Clean HI.
clean(vis='cig96.all.weight.HI.ms', imagename='cig96.all.weight.HI', mode='velocity', start='1330km/s', width='10km/s', nchan=48, interactive=False, phasecenter='J2000 02h15m27.6 +06d00m09', imsize=[600,600], niter=3000, cell='4.0arcsec', restfreq='1.42040575177GHz', weighting='natural', usescratch=T)

# Clean HI with the most regular resolution: circular beam from the largest beam axis of the weighted image.
clean(vis='cig96.all.weight.HI.ms', imagename='cig96.all.weight.circbeam.HI', mode='velocity', start='1330km/s', width='10km/s', nchan=48, phasecenter='J2000 02h15m27.6 +06d00m09', interactive=False, imsize=[600,600], niter=3000, cell='4.0arcsec', restfreq='1.42040575177GHz', weighting='natural', restoringbeam =['28arcsec'], usescratch=T)

# Clean HI with the most regular resolution and a velocity smoothing of 20 km/s.
clean(vis='cig96.all.weight.HI.ms', imagename='cig96.all.weight.circbeam.20kms.HI', mode='velocity', start='1330km/s', width='20km/s', nchan=24, phasecenter='J2000 02h15m27.6 +06d00m09', interactive=False, imsize=[600,600], niter=3000, cell='4.0arcsec', restfreq='1.42040575177GHz', weighting='natural', restoringbeam =['28arcsec'], usescratch=T)

# 0th, 1st and 2nd moments of the cubes
immoments(imagename='cig96.all.weight.HI.image', outfile='cig96.all.weight.HI.image.5sig.mom', moments=[0,1,2], includepix = [0.001245,1000])
immoments(imagename='cig96.all.weight.circbeam.HI.image', outfile='cig96.all.weight.circbeam.HI.image.5sig.mom', moments=[0,1,2], includepix = [0.00126,1000])
immoments(imagename='cig96.all.weight.circbeam.20kms.HI.image', outfile='cig96.all.weight.circbeam.20kms.HI.image.5sig.mom', moments=[0,1,2], includepix = [0.00124,1000])

# WT. First we transform to fits files the cubes.
exportfits(imagename='cig96.all.weight.HI.image', fitsimage='cig96.all.weight.HI.fits', velocity=True)
exportfits(imagename='cig96.all.weight.circbeam.HI.image', fitsimage='cig96.all.weight.circbeam.HI.fits', velocity=True)
exportfits(imagename='cig96.all.weight.circbeam.20kms.HI.image', fitsimage='cig96.all.weight.circbeam.20kms.HI.image.fits', velocity=True)

# Then we execute the WT script with SNR = 5.0 and image noise of 0.250 mJy/beam for VLA and 0.252 mJy/beam for EVLA.
execfile("wt.line.py")
# Filtered fits:
# cig96.all.weight.HI.fits.Filt.fits
# cig96.all.weight.circbeam.HI.fits.Filt.fits
# cig96.all.weight.circbeam.20kms.HI.image.fits.Filt.fits

# From fits to CASA image again via importfits.foptical
importfits(fitsimage='cig96.all.weight.HI.fits.Filt.fits', imagename='cig96.w.HI.filt.image')
importfits(fitsimage='cig96.all.weight.circbeam.HI.fits.Filt.fits', imagename='cig96.w.circ.HI.filt.image')
importfits(fitsimage='cig96.all.weight.circbeam.20kms.HI.image.fits.Filt.fits', imagename='cig96.w.circ.20kms.HI.filt.image')

# Axis switch: exportfits task changes stokes and velocity axii so after WT and importfits we have to reverse these two again to keep using blanking.py
imtrans(imagename='cig96.w.HI.filt.image', order='0132', outfile='cig96.w.HI.Filt.image')
imtrans(imagename='cig96.w.circ.HI.filt.image', order='0132', outfile='cig96.w.circ.HI.Filt.image')
imtrans(imagename='cig96.w.circ.20kms.HI.filt.image', order='0132', outfile='cig96.w.circ.20kms.HI.Filt.image')

# Remove wrongly-ordered images to avoid mistakes
!rm -Rf cig96.w.HI.filt.image cig96.w.circ.HI.filt.image cig96.w.circ.20kms.HI.filt.image

# DONE (Continuum calibration and imaging follows)
