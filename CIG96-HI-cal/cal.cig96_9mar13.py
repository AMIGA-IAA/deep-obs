# -*- coding: iso-8859-1 -*-
### EVLA DATA REDUCTION
### Project: 13A-341
### Dataset date: 9mar13
### Original dataset name: 13A-341.sb19189214.eb19331286.56361.695559375
### Renamed as: cig96_01.ms
###
### Configuration: D (1/3)

# ===============================================================================

### Import of EVLA data from SDM format:
########################################

# importevla(asdm='cig96_01', vis='cig96_01.ms')

# ===============================================================================

### Listobs inspection:
#######################

listobs(vis='cig96_01.ms')

# Listobs output summary: 27 antennae, 9 spectral windows, 2 polarizations (RR, LL), 1 dummy scan, 3 sources (2 calibrators + 1 target), 2048 channels of 7.8 kHz each for spw=0

# Spectral Window 0 (line):
# SpwID  Name          #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) BBC Num  Corrs  
# 0	 A0C0#0        2048    TOPO    1404.995         7.812     16000.0      12  RR  LL

#  Sources: 3
# ID   Code Name                RA               Decl           Epoch   SrcId      nRows
# 0    NONE 0137+331=3C48       01:37:41.299431 +33.09.35.13299 J2000   0          28431
# 1    K    0137+331=3C48       01:37:41.299431 +33.09.35.13299 J2000   1         338013
# 2    D    J0238+1636          02:38:38.930108 +16.36.59.27470 J2000   2         202176
# 3    NONE CIG 96              02:15:27.600000 +06.00.09.00001 J2000   3        1658475       

### BANDPASS and FLUX calibrator: 	1 0137+331 = 3C48
### PHASE calibrator: 			2 J0238+1636
### TARGET: 				3 CIG 96

# ===============================================================================

### Data inspection via plotms:
###############################

# Prior to anything else, there is dummy scan, field=0, so we remove it:

flagdata(vis='cig96_01.ms', mode='manual', field='0')

# We run plotms to inspect the general aspect of the data:

# plotms(vis='cig96_01.ms', spw='0', xaxis='time', yaxis='amp', coloraxis='field', avgchannel='1e9')

# ===============================================================================

# Antenna position correction:
##############################

gencal(vis='cig96_01.ms', caltable='cig96_01.ms.antpos', caltype='antpos')

# ===============================================================================

### Plot with the spatial distribution of the antennae:
#######################################################

plotants(vis='cig96_01.ms',figfile='plotants_cig96_01.png')

# ===============================================================================

### Flagging:
#############

# Log file does not provide with any insight of the antennae status, corruption, etc.

### Shadowing correction over the whole .ms:

flagdata(vis='cig96_01.ms', mode='shadow',  flagbackup=True)

# Percentage of data flagged in table selection: 16.6475%

### Zero clipping correction over the whole .ms:

flagdata(vis='cig96_01.ms', mode='clip', clipzeros=True, flagbackup=True)

# Percentage of data flagged in table selection: 0.00842345%

### Field 1 (bandpass and flux calibrator): good behaviour in spw=0, except for antenna ea19(17) and ea25(23) which 
# are showing low amplitudes in both RR and LL and for all baselines, so we remove them:

flagdata(vis='cig96_01.ms', mode='manual', field='1', spw='0', antenna='ea19', flagbackup=True)
flagdata(vis='cig96_01.ms', mode='manual', field='1', spw='0', antenna='ea25', flagbackup=True)

# apart from this antennas, the rest of the field looks good in both polarizations.

### Field 2 (phase calib.): there is an intense RFI located around channel 385 (freq.=1.408 GHz):
# NRAO known RFIs: https://science.nrao.edu/facilities/vla/observing/RFI/L-Band
# Known RFI in 1.408 GHz: http://www.vla.nrao.edu/astro/RFI/rfi_plots/SUB4.gif
# "The narrow feature at 1408 MHz is internally generated."
# Galactic plane emission in 1.420 GHz: confirm with Min.

# we flag it in field 2:

flagdata(vis='cig96_01.ms', mode='manual', field='2', spw='0:360~410', flagbackup=True)
flagdata(vis='cig96_01.ms', mode='manual', field='2', spw='0:1950~1975', flagbackup=True)

# Again, antenna ea25(23) is showing low amplitudes for all baselines so we remove it completely:

flagdata(vis='cig96_01.ms', mode='manual', field='2', spw='0', antenna='ea25', flagbackup=True)

### Field 3 (CIG96): shows a few RFIs: one in channels ~470, one in 164 and also the same RFIs as in field 2, one in 
# channel ~385 (freq.~1.408 GHz) and the same intense emission in channel ~1980 (freq.~1.4203 GHz). We remove 
# them all:

flagdata(vis='cig96_01.ms', mode='manual', field='3', spw='0:465~480', flagbackup=True)
flagdata(vis='cig96_01.ms', mode='manual', field='3', spw='0:360~410', flagbackup=True)
flagdata(vis='cig96_01.ms', mode='manual', field='3', spw='0:1950~1975', flagbackup=True)
flagdata(vis='cig96_01.ms', mode='manual', field='3', spw='0:1604', flagbackup=True)

### Quacking: specially fields 2 and 3 need quacking at the beginning but we do it over the 3 fields 
# to be safe:

flagdata(vis='cig96_01.ms', mode='quack', field='', spw='0', quackinterval=15.0, quackmode='beg', flagbackup=True)

### And we add a autoflagging for each of the calibrators (do NOT mix them in one command, use two to 
### avoid flagging the calibrator with highest amplitudes):

# flagdata(vis='cig96_01.ms',  mode='rflag', field='1', spw='0', ntime='scan', action='apply')
# flagdata(vis='cig96_01.ms',  mode='rflag', field='2', spw='0', ntime='scan', action='apply')

# Enough flagging for now.

# ===============================================================================

### Calibration:
### ============

### Reference antenna selection:
################################

# After checking plotants, we select one of the central-most and has shown no problems:

refant='ea22'

# ===============================================================================

### Opacity and gaincurve corrections:
######################################

# Should we do opacities correction?
# myTau = plotweather(vis='cig96_01.ms', doPlot=T)
# display cig96_01.ms.plotweather.png
# gencal(vis='cig96_01.ms', caltype='opac', caltable='cig96_01.ms.opcac', parameter=myTau)
# ANSWER: no, not needed for HI observations, the atmosphere barely has any influence with HI. Only 
# ionosphere in very flat elevations and at dawn/sunset.

# Should we do gaincurve correction?
# Takes into account the elevations of each antenna, see elevations plot: CIG96 and the phase calibrator 
# have a large elevation range.
#
# gencal(vis='cig96_01.ms', caltype='gceff', caltable='cig96_01.ms.gaincurve')
#
# ANSWER: no, it might even harm the observations due to the possible crosstalk at low elevations 
# in C and D configs. and the model it is based on.

# ===============================================================================

### Delays correction:
######################

# Prior to the bandpass correction, we need to correct the phase variations with time (they are large). 
# We use the bandpass calibrator in field 1 with the antenna position correction table:

gaincal(vis='cig96_01.ms', caltable='cig96_01.ms.delays', field='1', refant='ea22', gaintype='K', gaintable=['cig96_01.ms.antpos'])

# ===============================================================================

### Flux density:
#################

# Our flux calibrator is 3C48, in field = 0.
# First, we check the list of available models:

setjy(vis='cig96_01.ms', listmodels=T)

# Candidate modimages (*) at pepino (dae66) in path:
# /Applications/CASA.app/Contents/data/nrao/VLA/CalModels/
#
# Candidate modimages (*) at NRAO local workstations in path:
# /home/casa/packages/RHEL5/release/casapy-42.1.29047-001-1-64b/data/nrao/VLA/CalModels/
# 
# The model chosen has to be in accordance with the calibrator and band selected: 3C48 in L band:

setjy(vis='cig96_01.ms', field='1', modimage='/mnt/scops/data/data/paramimo/casapy-42.2.30986-1-64b/data/nrao/VLA/CalModels/3C48_L.im')

# ===============================================================================

### Bandpass calibration:
#########################

# For the bandpass and flux (field=1) and phase (field=2) calibrators we use solution interval 
# time of solint='5s':

gaincal(vis='cig96_01.ms', caltable='cig96_01.ms.bpphase5s', field='1', refant='ea22', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_01.ms.antpos', 'cig96_01.ms.delays'])

# The solution interval of 5 seconds has been calculated following with the VLA exp.time 
# calculator using the parameters:
#
# Freq. = 1.42 GHz
# Medium elevation (summer time)
# Bandwith freq. = 3,076 KHz (see listobs)
# RMS noise = 20 mJy (more than 3.0 in SNR)
# 
# The calculator estimates less than one second (0.2s) is enough to get such SNR or even higher 
# so we set solint=5s since 5s is the shortest integration time of our data. This should mean that 
# there should be a solution for all the intervals.

# We see the phase VS time figures to see that the phase is now waaayyy flatter:

# plotcal(caltable='cig96_01.ms.bpphase5s', xaxis='time', yaxis='phase')

# Apply phase solutions on the fly:

bandpass(vis='cig96_01.ms', caltable='cig96_01.ms.bandpass5s', field='1', refant='ea22', solint='inf', solnorm=T, minsnr=3.0, minblperant=3, gaintable=['cig96_01.ms.bpphase5s', 'cig96_01.ms.antpos', 'cig96_01.ms.delays'], interp=['nearest'])

# We check again the solutions:

# plotcal(caltable='cig96_01.ms.bpphase5s', field='1', xaxis='time', yaxis='phase')  

# Difference in phase is less, phase shows a flatter behaviour now.

# ===============================================================================


### Phase and amplitude calibration:
####################################

# Phase calibration for the calibrators, fields 1 and 2 (spw 0, where the line is):

gaincal(vis='cig96_01.ms',caltable='cig96_01.ms.intphase', field='1,2', refant='ea22', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_01.ms.antpos', 'cig96_01.ms.delays', 'cig96_01.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, 
# to use FOR THE TARGET later on:

# time.sleep(5)

gaincal(vis='cig96_01.ms', caltable='cig96_01.ms.scanphase', field='1,2', refant='ea22', calmode='p', solint='inf', minsnr=5.0, gaintable=['cig96_01.ms.antpos', 'cig96_01.ms.delays', 'cig96_01.ms.bandpass5s'])

# Derive amplitude solutions:

# time.sleep(5)

gaincal(vis='cig96_01.ms', caltable='cig96_01.ms.amp', field='1,2',  refant='ea22', calmode='ap', solint='inf', minsnr=5.0, gaintable=['cig96_01.ms.antpos', 'cig96_01.ms.delays', 'cig96_01.ms.bandpass5s', 'cig96_01.ms.intphase'])

# I check the tables with plotcal and some things do not look as expected: the first data show very 
# low intensity compares to the rest and some graphs show only one point:

# plotcal(caltable='cig96_01.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)
# plotcal(caltable='cig96_01.ms.amp', xaxis='time', yaxis='amp', iteration='baseline', subplot=331)

# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal 
# in terms of future application of the calibration to the data, (i.e., it's not an incremental table) 
# UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_01.ms', caltable='cig96_01.ms.amp', fluxtable='cig96_01.ms.flux', reference='1', transfer='2', incremental=False)

# Result:
# Flux density for J0238+1636 in SpW=0 (freq=1.40499e+09 Hz) is: 0.827323 +/- 0.0127157 (SNR = 65.0632, N = 50)

# ===============================================================================

### Application of the calibration to the .ms file:
###################################################

# Note: In all applycal steps we set calwt=F. It is very important to turn off this parameter which 
# determines if the weights are calibrated along with the data. Data from antennae with better receiver 
# performance and/or longer integration times should have higher weights, and it can be advantageous to 
# factor this information into the calibration. During the VLA era, meaningful weights were available for 
# each visibility. However, at the time of this observation, the VLA was not yet recording the information 
# necessary to calculate meaningful weights. Since these data weights are used at the imaging stage you 
# can get strange results from having calwt=T when the input weights are themselves not meaningful, 
# especially for self-calibration on resolved sources (your flux calibrator and target, for example).

applycal(vis='cig96_01.ms', field='1', gaintable=['cig96_01.ms.antpos', 'cig96_01.ms.delays', 'cig96_01.ms.bandpass5s', 'cig96_01.ms.intphase', 'cig96_01.ms.flux'], gainfield=['','','','1',''], calwt=F)

# time.sleep(5)

applycal(vis='cig96_01.ms', field='2', gaintable=['cig96_01.ms.antpos', 'cig96_01.ms.delays', 'cig96_01.ms.bandpass5s', 'cig96_01.ms.intphase', 'cig96_01.ms.flux'], gainfield=['','','','2',''], calwt=F)

# time.sleep(5)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_01.ms', field='3', gaintable=['cig96_01.ms.antpos', 'cig96_01.ms.delays', 'cig96_01.ms.bandpass5s', 'cig96_01.ms.scanphase', 'cig96_01.ms.flux'], gainfield=['','','','2',''], calwt=F)

# ===============================================================================

### Regridding of the .ms to a new frame:
#########################################

cvel(vis='cig96_01.ms', outputvis='cig96_01.ms.cvel', mode='velocity', field='', spw='0', restfreq='1.42040575177GHz', outframe='LSRK', veltype='radio')


# ===============================================================================

### Splitting of CIG96:
#######################

# Splitting of field = 3 (CIG96) calibrated data, in the spw=0:

split(vis='cig96_01.ms.cvel', datacolumn='corrected', outputvis='cig96_01.cvel.corr.ms', field='3', spw='0:400~1700')

# Now, for the new corr.ms file:
# field = 3 is stored as: field = 0  
# spw = 0 keeps being spw = 0

# Split of the whole spw:

# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.spw0.ms', field='3', spw='0')

# Emission in the higher frequency:
# 
# The emission line is seen in 1.4203 GHz. We split the fields in spw=0 from channel 780 to 1600:

# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.emission.ms', field='3', spw='0:780~1600')

# Splitting CIG96 fields in spw=0 with time binning of 20s and channel width of 10 channels:

# split(vis='cig96_01.ms.cvel', datacolumn='corrected', outputvis='cig96_01.corr.20s.10chan.ms', width=10, timebin='20s', field='3', spw='0:451~1700')

# ===============================================================================

### UV continuum subtraction:
#############################

# Since we have chopped out the first 390 channels where the RFI region was (see spw definition in the 
# split command) and the last 350 channels, the total number of channels has changed: from 2048 we 
# now have ~1350.
# 
# Now, we remove the continuum from source by subtracting the channels indicated:

uvcontsub(vis='cig96_01.cvel.corr.ms', field='0', fitspw='0:50~350;900~1050', fitorder=1)
# uvcontsub(vis='cig96_01.corr.20s.10chan.ms', field='0', fitspw='0:15~30;90~105', fitorder=1)

# Also, for the second emission in 1.4203 GHz, the emission is in a very noisy region so the continuum 
# subtraction won't be perfect, likely more noisy on higher frequencies since the continuum is
# selected in channels 780~1000:

# uvcontsub(vis='cig96_01.corr.emission.ms', field='0', fitspw='0:0~510;710~1300')

# For the time-channel binned file, we subtract as follows:

# uvcontsub(vis='cig96_01.corr.20s.10chan.ms', field='0', fitspw='0:45~93;110~181')

# For the whole spw 0 we use an as extended as possible continuum:

# uvcontsub(vis='cig96_01.corr.spw0.ms', field='0', fitspw='0:300~365;396~460;490~920;1110~1520')

# ===============================================================================

### Clean of the continuum subtracted data:
########################################### 

# Weighting natural with the factor of 6.25064 for the channel smoothing that will be necessary to combine with VLA data:

clean(vis='cig96_01.cvel.corr.ms.contsub', imagename='cig96_01.cvel.corr.contsub.natural.line', field='0', spw='0:450~800', mode='frequency', start=450, nchan=351, niter=10000, width='48.830kHz', threshold='1.5mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='8.0arcsec', restfreq='1.42040575177GHz', weighting='natural', usescratch=T)

# Beam size:  arcsec2

# Weighting uniform:

# clean(vis='cig96_01.corr.ms.contsub', imagename='cig96_01.corr.v3.contsub.uniform.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='uniform', usescratch=T)

# Beam size: 

# Weighting Briggs rob=0.0:

# clean(vis='cig96_01.corr.ms.contsub', imagename='cig96_01.corr.v3.contsub.rob0.0.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='briggs', robust='0.0', usescratch=T)

# Beam size: 

# Collapse of cube (moments):

immoments(imagename='cig96_01.cvel.corr.contsub.natural.line.image', axis='spectral', moments=[0,1], outfile='cig96_01.cvel.corr.contsub.natural.line.image.mom')


################ END
#
#
# IN CASE YOU WISH TO RESET THE WHOLE CALIBRATION, TAKE THE FOLLOWING STEPS:
#
# clearcal(vis='xxx')
#
# flagdata(vis='xxx',mode=unflag)
#
# Ready to redo calibration!





######################### OLD #######################




### Splitting of CIG96:
#######################

# Splitting of field = 3 (CIG96) calibrated data, in the spw=0:

# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.v3.ms', field='3', spw='0:399~2047')

# Now, for the new corr.ms file:
# field = 3 is stored as: field = 0  
# spw = 0 keeps being spw = 0

# Split of the whole spw:

# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.spw0.ms', field='3', spw='0')

# Emission in the higher frequency:
# 
# The emission line is seen in 1.4203 GHz. We split the fields in spw=0 from channel 780 to 1600:

# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.emission.ms', field='3', spw='0:780~1600')

# Splitting CIG96 fields in spw=0 with time binning of 20s and channel width of 10 channels:

# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.20s.10chan.ms', width=10, timebin='20s', field='3', spw='0')

# ===============================================================================

### UV continuum subtraction:
#############################

# Since we have chopped out the first 400 channels where the RFI region was (see spw definition in the 
# split command) the total number of channels have changed now: from 2048 we now have 1648.
# Now, we remove the continuum from source by subtracting the channels from 500to910 and from 730to1300.

# uvcontsub(vis='cig96_01.corr.v3.ms', field='0', fitspw='0:0~510;730~1300')

# Also, for the second emission in 1.4203 GHz, the emission is in a very noisy region so the continuum 
# subtraction won't be perfect, likely more noisy on higher frequencies since the continuum is
# selected in channels 780~1000:

# uvcontsub(vis='cig96_01.corr.emission.ms', field='0', fitspw='0:0~510;730~1300')

# For the time-channel binned file, we subtract as follows:

# uvcontsub(vis='cig96_01.corr.20s.10chan.ms', field='0', fitspw='0:45~93;110~181')

# For the whole spw 0 we use an as extended as possible continuum:

# uvcontsub(vis='cig96_01.corr.spw0.ms', field='0', fitspw='0:300~365;396~460;490~920;1110~1520')

# ===============================================================================

### Clean of the continuum subtracted data:
########################################### 

# Weighting natural:

# clean(vis='cig96_01.corr.ms.contsub', imagename='cig96_01.corr.v3.contsub.natural.line', field='0', spw='0:510~730', mode='channel', start=510, nchan=220, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)

# Beam size: 

# WARN MFCleanImageSkyModel   Clean not converging   
# Successfully deconvolved image
# Beam used in restoration: 99.5838 by 47.044 (arcsec) at pa -56.6383 (deg) 

# Weighting uniform:

# clean(vis='cig96_01.corr.ms.contsub', imagename='cig96_01.corr.v3.contsub.uniform.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='uniform', usescratch=T)

# Beam size: 

# Weighting Briggs rob=0.0:

# clean(vis='cig96_01.corr.ms.contsub', imagename='cig96_01.corr.v3.contsub.rob0.0.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='briggs', robust='0.0', usescratch=T)

# Beam size: 



### Noise level via viewer task:
################################

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# (cig96_01.corr.contsub.line.image)
#      Frequency       Velocity         Stokes BrightnessUnit       BeamArea 
#  1.41305e+09Hz    1553.49km/s              I        Jy/beam        4.44663 
#           Npts            Sum    FluxDensity           Mean            Rms 
#           1155   2.966822e-02   6.672063e-03   2.568677e-05   4.262555e-04 
#        Std dev        Minimum        Maximum   region count 
#   4.256651e-04  -1.101821e-03   1.550763e-03              1 
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

# For the time-channel binned file, we clean as follows:

# Weighting natural:

# clean(vis='cig96_01.corr.20s.10chan.ms.contsub', imagename='cig96_01.corr.20s.10chan.ms.contsub.natural.line.image', field='0', spw='0:93~110', mode='channel', start=0, niter=1000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)

# For the spw 0 split dataset:

# clean(vis='cig96_01.corr.spw0.ms.contsub', imagename='cig96_01.corr.spw0.ms.contsub.natural.line', field='0', spw='0:920~1110', mode='channel', start=920, nchan=190, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='10.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)

# ===============================================================================

### Corrected data inspection via msview:
#########################################

# We inspect the baselines that present solar interferences via msview

# Field 1 shows interference patterns as well as some noisy baselines:
# 
# 2-8, 5-26, 5-20, 8-13, 8-12, 2-12, 2-25, 2-20, 22-24, 0-26, 12-13, 0-1, 12-25, 6-8, 3-15, 1-7, 20-25, 7-14, 
# 3-19, 2-24, 2-6, 2-5, 11-19, 0-5, 12-22, 9-25, 4-21, 11-16, 24-25, 5-25, 1-26, 5-8, 2-26, 2-9, 5-24, 5-6, 
# 25-26, 1-5, 20-22, 0-2, 10-24, 5-9, 7-26, 0-14, 5-7, 1-2, 5-15, 1-25, 8-17, 17-22, 12-17, 13-17, 17-25, 6-17, 
# 17-26, 18-24, 3-17, 
#  
# Antennas 0=ea01, 1=ea02, 4=ea05, 6=ea07 and 18=ea20 always shows very large amp differences, we flag it:

# flagdata(vis='cig96_01.ms', mode='manual', field='1', spw='0', antenna='ea07;ea01;ea02;ea05;ea20')

# We flag the baselines that have shown interferences:

# flagdata(vis='cig96_01.ms', mode='manual', field='1', spw='0', antenna='2&8;5&26;5&20;8&13;8&12;2&12;2&25;2&20;22&24;0&26;12&13;0&1;12&25;6&8;3&15;1&7;20&25;7&14;3&19;2&24;2&6;2&5;11&19;0&5;12&22;9&25;4&21;11&16;24&25;5&25;1&26;5&8;2&26;2&9;5&24;5&6;25&26;1&5;20&22;0&2;10&24;5&9;7&26;0&14;5&7;1&2;5&15;1&25;8&17;17&22;12&17;13&17;17&25;6&17;17&26;18&24;3&17', flagbackup=True)

# We now delete the calibration to perform it again with the new flagged data:
# 
# To do so, we restart from the beginning, right after the flagging.
# 
# New flux value:
# 
# Flux density for J0238+1636 in SpW=0 (freq=1.40499e+09 Hz) is: 0.818094 +/- 0.0143805 (SNR = 56.8892, N = 44)
# 
# After applying the whole calibration again, we split field = 3 (CIG96) re-calibrated data, in the spw=0:

# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.v2.ms', field='3', spw='0:399~2047')
# split(vis='cig96_01.ms', datacolumn='corrected', outputvis='cig96_01.corr.allfields.ms', field='1~3', spw='0')

# and redo uvcontsub:

# Since we have chopped out the first 400 channels where the RFI region was (see spw definition in the 
# split command) the total number of channels have changed now: from 2048 we now have 1648.
# Now, we remove the continuum from source by subtracting the channels from 500to910 and from 730to1300.

# uvcontsub(vis='cig96_01.corr.v2.ms', field='0', fitspw='0:0~510;730~1300')

# Cleaning: weighting natural:

# clean(vis='cig96_01.corr.v2.ms.contsub', imagename='cig96_01.corr.v2.contsub.natural.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[300,300], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)




################ END
#
#
# IN CASE YOU WISH TO RESET THE WHOLE CALIBRATION, TAKE THE FOLLOWING STEPS:
#
# clearcal(vis='xxx')
#
# flagdata(vis='xxx',mode=unflag)
#
# Ready to redo calibration!

