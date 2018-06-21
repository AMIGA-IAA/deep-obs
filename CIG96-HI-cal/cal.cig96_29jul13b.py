# -*- coding: iso-8859-1 -*-
### EVLA DATA REDUCTION
### Project: 13A-341
### Dataset date: 29jul13b
### Original dataset name: 13A-341.sb24028522.eb24167588.56502.374954375.ms
### Renamed as: cig96_13.ms
###
### Configuration: C (8/10)

# ===============================================================================

###   CALIBRATION   ###

# ===============================================================================

### Import of EVLA data from SDM format:
########################################

# importevla(asdm='cig96_13', vis='cig96_13.ms')

# Original data file gets now CASA format and a new name: "cig96_13.ms"

# ===============================================================================

### Listobs inspection:
#######################

listobs(vis='cig96_13.ms')

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

# flagdata(vis='cig96_13.ms', mode='manual', field='0')

# We run plotms to inspect the general aspect of the data:

#plotms(vis='cig96_13.ms',xaxis='channel',yaxis='amp',field='1',spw='0',coloraxis='field',avgtime='1e9')
#plotms(vis='cig96_13.ms',xaxis='time',yaxis='amp',field='1',spw='0',coloraxis='field',avgchannel='1e9')

# ===============================================================================

# Antenna position correction:
##############################

gencal(vis='cig96_13.ms', caltable='cig96_13.ms.antpos', caltype='antpos')

# ===============================================================================

### Plot with the spatial distribution of the antennae:
#######################################################

plotants(vis='cig96_13.ms',figfile='cig96_13.ms.plotants.png')

# ===============================================================================

### Flagging:
#############

# Log file provides with some insight of the antennae status, corruption, etc.:

# Antenna(s) 19 (Data: Corrupted):
# L-band visibility amplitudes below normal due to receiver LNA problem.

# so we flag it.

flagdata(vis='cig96_13.ms', mode='manual', field='', antenna='ea19', flagbackup=True)

### Shadowing correction over the whole .ms:

flagdata(vis='cig96_13.ms', mode='shadow', flagbackup=True)

# Percentage of data flagged in table selection: 


### Zero clipping correction over the whole .ms: 

flagdata(vis='cig96_13.ms', mode='clip', clipzeros=True, flagbackup=True)

# Percentage of data flagged in table selection: 4.31437e-07%



### Field 1 (bandpass and flux calibrator): quacking of the first 30 seconds (all antennas showing low amplitudes):

flagdata(vis='cig96_13.ms', mode='quack', field='1', quackinterval=30.0, quackmode='beg')

# bad antennas:

flagdata(vis='cig96_13.ms', mode='manual', scan='2', antenna='3&5', correlation='LL', flagbackup=True)

# there is a baseline, 5&15 in LL correllation, that shows a higher Amp in Amp vs channel plot, we remove it:

flagdata(vis='cig96_13.ms', mode='manual', scan='2', antenna='5&15', correlation='LL', flagbackup=True)


### Field 2 (phase calib.): needs quacking, and anteannae 3 and 5 in corr LL need removal, they show some high values:

# We could remove the upper amplitude values:
# flagdata(vis='cig96_13.ms', mode='manual', field='2', antenna='3&5', correlation='LL', flagbackup=True)
# Decission: we let them be.

# Quacking needed:

flagdata(vis='cig96_13.ms', mode='quack', field='2', quackinterval=5.0, quackmode='beg')

# bad antennas: baseline 3&5 always shows the highest amplitudes in LL pol., quite above the average values of 0.015 in Amp

flagdata(vis='cig96_13.ms', mode='manual', field='2', antenna='3&5', correlation='LL', flagbackup=True)

# it has the RFI in channel ~350:

flagdata(vis='cig96_13.ms', mode='manual', field='2', spw = '0:330~380')

# and a small one in ~1250:

flagdata(vis='cig96_13.ms', mode='manual', field='2', antenna='4&8', spw = '0:1253')

# and the MW HI emission in ~1960:

flagdata(vis='cig96_13.ms', mode='manual', field='2', spw = '0:1940~1975')


### Field 3 (CIG96) flagging:

# quacking of the first 10 seconds of the rest of the scans:

flagdata(vis='cig96_13.ms', mode='quack', field='3', quackinterval=10.0, quackmode='beg')

# same RFI in ~350 and MW emission as in field 2; no RFI in ~1250:

flagdata(vis='cig96_13.ms', mode='manual', field='3', spw = '0:330~380')
flagdata(vis='cig96_13.ms', mode='manual', field='3', spw = '0:1940~1975')

# Enough flagging for now.

# ===============================================================================

### Calibration:
### ============

### Reference antenna selection:
################################

# After checking plotants, we select this one because it is quite central:

refant='ea10'

# ===============================================================================

### Opacity and gaincurve corrections:
######################################

# Should we do opacities correction?
# myTau = plotweather(vis='cig96_13.ms', doPlot=T)
# display cig96_13.ms.plotweather.png
# gencal(vis='cig96_13.ms', caltype='opac', caltable='cig96_13.ms.opcac', parameter=myTau)
# ANSWER: no, not needed for HI observations, the atmosphere barely has any influence with HI. Only ionosphere in very flat elevations and at dawn/sunset.

# Should we do gaincurve correction?
# Takes into account the elevations of each antenn, see elevations plot: CIG96 and the phase calibrator have a large elevation range.
# gencal(vis='cig96_13.ms', caltype='gceff', caltable='cig96_13.ms.gaincurve')
# ANSWER: no, it might even harm the observations due to the possible crosstalk at low elevations and the model it is based on.

# ===============================================================================

### Delays correction:
######################

# Prior to the bandpass correction, we need to correct the phase variations with time (they are large). 
# We use the bandpass calibrator in field 1 with the antenna position correction table:

gaincal(vis='cig96_13.ms', caltable='cig96_13.ms.delays', field='1', refant='ea10', gaintype='K', gaintable=['cig96_13.ms.antpos'])

# ===============================================================================

### Flux density:
#################

# Our flux calibrator is 3C48, in field = 1.
# 
# First, we check the list of available models:

setjy(vis='cig96_13.ms', listmodels=T)

# Candidate modimages (*) at pepino (dae66) in path:
# /Applications/CASA.app/Contents/data/nrao/VLA/CalModels/
#
# Candidate modimages (*) at NRAO local workstations in path:
# /home/casa/packages/RHEL5/release/casapy-42.1.29047-001-1-64b/data/nrao/VLA/CalModels/
# 
# The model chosen has to be in accordance with the calibrator and band selected: 3C48 in L band:

setjy(vis='cig96_13.ms', field='1', modimage='/mnt/scops/data/data/paramimo/casapy-42.2.30986-1-64b/data/nrao/VLA/CalModels/3C48_L.im')

# ===============================================================================

### Bandpass calibration:
#########################

# For the bandpass and flux (field=1) and phase (field=2) calibrators we use solution interval 
# time of solint='5s':

gaincal(vis='cig96_13.ms', caltable='cig96_13.ms.bpphase5s', field='1', refant='ea10', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays'])

# The solution interval of 5 seconds has been calculated following with the VLA exp.time 
# calculator using the parameters:
#
# Freq. = 1.42 GHz
# Medium elevation (summer time)
# Bandwith freq. = 16 MHz (see listobs)
# RMS noise = 1.5 mJy
# 
# The calculator estimates less than one second (0.2s) is enough to get such SNR or even higher 
# so we set solint=5s since 5s is the shortest integration time of our data. This should mean that 
# there should be a solution for all the intervals.

# We see the phase VS time figures to see that the phase is now waaayyy flatter:

#plotcal(caltable='cig96_13.ms.bpphase5s', xaxis='time', yaxis='phase')

# Apply phase solutions on the fly:

bandpass(vis='cig96_13.ms', caltable='cig96_13.ms.bandpass5s', field='1', refant='ea10', solint='inf', solnorm=T, minsnr=3.0, minblperant=3, gaintable=['cig96_13.ms.bpphase5s', 'cig96_13.ms.antpos', 'cig96_13.ms.delays'], interp=['nearest'])

# We check again the solutions:

#plotcal(caltable='cig96_13.ms.bpphase5s', field='1', xaxis='time', yaxis='phase')  

# Difference in phase is less, phase shows a flatter behaviour now.

# ===============================================================================


### Phase and amplitude calibration:
####################################

# Phase calibration for the calibrators, fields 1 and 2 (spw 0, where the line is):

gaincal(vis='cig96_13.ms',caltable='cig96_13.ms.intphase', field='1,2', refant='ea10', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays', 'cig96_13.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, 
# to use FOR THE TARGET later on:

# time.sleep(5)

gaincal(vis='cig96_13.ms', caltable='cig96_13.ms.scanphase', field='1,2', refant='ea10', calmode='p', solint='inf', minsnr=5.0, gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays', 'cig96_13.ms.bandpass5s'])

# Derive amplitude solutions:

# time.sleep(5)

gaincal(vis='cig96_13.ms', caltable='cig96_13.ms.amp', field='1',  refant='ea10', calmode='ap', solint='inf', minsnr=5.0, gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays', 'cig96_13.ms.bandpass5s', 'cig96_13.ms.intphase'], gainfield=['','','','1'],append=True)

gaincal(vis='cig96_13.ms', caltable='cig96_13.ms.amp', field='2',  refant='ea10', calmode='ap', solint='inf', minsnr=5.0, gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays', 'cig96_13.ms.bandpass5s', 'cig96_13.ms.intphase'], gainfield=['','','','2'],append=True)

# I check the tables with plotcal and some things do not look as expected: the first data show very 
# low intensity compares to the rest and some graphs show only one point:

# #######plotcal(caltable='cig96_13.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)
# #######plotcal(caltable='cig96_13.ms.amp', xaxis='time', yaxis='amp', iteration='baseline', subplot=331)

# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal 
# in terms of future application of the calibration to the data, (i.e., it's not an incremental table) 
# UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_13.ms', caltable='cig96_13.ms.amp', fluxtable='cig96_13.ms.flux', reference='1', transfer='2', incremental=False)

# Result:
# Flux density for J0238+1636 in SpW=0 (freq=1.40523e+09 Hz) is: 0.738149 +/- 0.00258917 (SNR = 285.09, N = 52)

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

applycal(vis='cig96_13.ms', field='1', gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays', 'cig96_13.ms.bandpass5s', 'cig96_13.ms.intphase', 'cig96_13.ms.flux'], gainfield=['','','','1',''], calwt=F)

# time.sleep(5)

applycal(vis='cig96_13.ms', field='2', gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays', 'cig96_13.ms.bandpass5s', 'cig96_13.ms.intphase', 'cig96_13.ms.flux'], gainfield=['','','','2',''], calwt=F)

# time.sleep(5)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_13.ms', field='3', gaintable=['cig96_13.ms.antpos', 'cig96_13.ms.delays', 'cig96_13.ms.bandpass5s', 'cig96_13.ms.scanphase', 'cig96_13.ms.flux'], gainfield=['','','','2',''], calwt=F)

# ===============================================================================

### Regridding of the .ms to a new frame:
#########################################

# cvel(vis='cig96_13.ms', outputvis='cig96_13.ms.cvel', mode='velocity', field='', spw='0', restfreq='1.42040575177GHz', outframe='LSRK', veltype='radio')
# 
# 
# # ===============================================================================
# 
# ### Splitting of CIG96:
# #######################
# 
# # Splitting of field = 3 (CIG96) calibrated data, in the spw=0:
# 
# split(vis='cig96_13.ms.cvel', datacolumn='corrected', outputvis='cig96_13.cvel.corr.ms', field='3', spw='0:400~1700')
# 
# # Now, for the new corr.ms file:
# # field = 3 is stored as: field = 0  
# # spw = 0 keeps being spw = 0
# 
# # Split of the whole spw:
# 
# # split(vis='cig96_13.ms', datacolumn='corrected', outputvis='cig96_13.corr.spw0.ms', field='3', spw='0')
# 
# # Emission in the higher frequency:
# # 
# # The emission line is seen in 1.4203 GHz. We split the fields in spw=0 from channel 780 to 1600:
# 
# # split(vis='cig96_13.ms', datacolumn='corrected', outputvis='cig96_13.corr.emission.ms', field='3', spw='0:780~1600')
# 
# # Splitting CIG96 fields in spw=0 with time binning of 20s and channel width of 10 channels:
# 
# split(vis='cig96_13.ms.cvel', datacolumn='corrected', outputvis='cig96_13.corr.20s.10chan.ms', width=10, timebin='20s', field='3', spw='0:451~1700')
# 
# # ===============================================================================
# 
# ### UV continuum subtraction:
# #############################
# 
# # Since we have chopped out the first 390 channels where the RFI region was (see spw definition in the 
# # split command) and the last 350 channels, the total number of channels has changed: from 2048 we 
# # now have ~1350.
# # 
# # Now, we remove the continuum from source by subtracting the channels indicated:
# 
# uvcontsub(vis='cig96_13.cvel.corr.ms', field='0', fitspw='0:50~350;900~1050', fitorder=1)
# uvcontsub(vis='cig96_13.corr.20s.10chan.ms', field='0', fitspw='0:15~30;90~105', fitorder=1)
# 
# # Also, for the second emission in 1.4203 GHz, the emission is in a very noisy region so the continuum 
# # subtraction won't be perfect, likely more noisy on higher frequencies since the continuum is
# # selected in channels 780~1000:
# 
# # uvcontsub(vis='cig96_13.corr.emission.ms', field='0', fitspw='0:0~510;710~1300')
# 
# # For the time-channel binned file, we subtract as follows:
# 
# # uvcontsub(vis='cig96_13.corr.20s.10chan.ms', field='0', fitspw='0:45~93;110~181')
# 
# # For the whole spw 0 we use an as extended as possible continuum:
# 
# # uvcontsub(vis='cig96_13.corr.spw0.ms', field='0', fitspw='0:300~365;396~460;490~920;1110~1520')
# 
# # ===============================================================================
# 
# ### Clean of the continuum subtracted data:
# ########################################### 
# 
# # Weighting natural with the factor of 6.25064 for the channel smoothing that will be necessary to combine with VLA data:
# 
# clean(vis='cig96_13.cvel.corr.ms.contsub', imagename='cig96_13.cvel.corr.contsub.natural.line', field='0', spw='0:450~800', mode='frequency', start=450, nchan=351, niter=10000, width='48.830kHz', threshold='1.5mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='8.0arcsec', restfreq='1.42040575177GHz', weighting='natural', usescratch=T)
# 
# # Beam size: 26.28 x 18.45 arcsec2
# 
# # Weighting uniform:
# 
# clean(vis='cig96_13.corr.ms.contsub', imagename='cig96_13.corr.v3.contsub.uniform.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='uniform', usescratch=T)
# 
# # Beam size: 
# 
# # Weighting Briggs rob=0.0:
# 
# clean(vis='cig96_13.corr.ms.contsub', imagename='cig96_13.corr.v3.contsub.rob0.0.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='briggs', robust='0.0', usescratch=T)
# 
# # Beam size: 
# 
# # Collapse of cube (moments):
# 
# immoments(imagename='cig96_13.cvel.corr.contsub.natural.line.image', axis='spectral', moments=[0,1], outfile='cig96_13.cvel.corr.contsub.natural.line.image.mom')
# 
# 
# ### Noise level via viewer task:
# ################################
# 
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# (cig96_13.corr.contsub.line.image)
#      Frequency       Velocity         Stokes BrightnessUnit       BeamArea 
#  1.41305e+09Hz    1553.49km/s              I        Jy/beam        4.44663 
#           Npts            Sum    FluxDensity           Mean            Rms 
#           1155   2.966822e-02   6.672063e-03   2.568677e-05   4.262555e-04 
#        Std dev        Minimum        Maximum   region count 
#   4.256651e-04  -1.101821e-03   1.550763e-03              1 
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# 
# # For the time-channel binned file, we clean as follows:
# 
# # Weighting natural:
# 
# clean(vis='cig96_13.corr.20s.10chan.ms.contsub', imagename='cig96_13.corr.20s.10chan.ms.contsub.natural.line.image', field='0', spw='0:93~110', mode='channel', start=0, niter=1000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)
# 
# # For the spw 0 split dataset:
# 
# clean(vis='cig96_13.corr.spw0.ms.contsub', imagename='cig96_13.corr.spw0.ms.contsub.natural.line', field='0', spw='0:920~1110', mode='channel', start=920, nchan=190, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[200,200], phasecenter='', cell='10.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)
# 
# # ===============================================================================
# 
# ### Corrected data inspection via msview:
# #########################################
# 
# # We inspect the baselines that present solar interferences via msview
# 
# # Field 1 shows interference patterns as well as some noisy baselines:
# # 
# # 2-8, 5-26, 5-20, 8-13, 8-12, 2-12, 2-25, 2-20, 22-24, 0-26, 12-13, 0-1, 12-25, 6-8, 3-15, 1-7, 20-25, 7-14, 
# # 3-19, 2-24, 2-6, 2-5, 11-19, 0-5, 12-22, 9-25, 4-21, 11-16, 24-25, 5-25, 1-26, 5-8, 2-26, 2-9, 5-24, 5-6, 
# # 25-26, 1-5, 20-22, 0-2, 10-24, 5-9, 7-26, 0-14, 5-7, 1-2, 5-15, 1-25, 8-17, 17-22, 12-17, 13-17, 17-25, 6-17, 
# # 17-26, 18-24, 3-17, 
# #  
# # Antennas 0=ea01, 1=ea02, 4=ea05, 6=ea07 and 18=ea20 always shows very large amp differences, we flag it:
# 
# flagdata(vis='cig96_13.ms', mode='manual', field='1', spw='0', antenna='ea07;ea01;ea02;ea05;ea20')
# 
# # We flag the baselines that have shown interferences:
# 
# flagdata(vis='cig96_13.ms', mode='manual', field='1', spw='0', antenna='2&8;5&26;5&20;8&13;8&12;2&12;2&25;2&20;22&24;0&26;12&13;0&1;12&25;6&8;3&15;1&7;20&25;7&14;3&19;2&24;2&6;2&5;11&19;0&5;12&22;9&25;4&21;11&16;24&25;5&25;1&26;5&8;2&26;2&9;5&24;5&6;25&26;1&5;20&22;0&2;10&24;5&9;7&26;0&14;5&7;1&2;5&15;1&25;8&17;17&22;12&17;13&17;17&25;6&17;17&26;18&24;3&17', flagbackup=True)
# 
# # We now delete the calibration to perform it again with the new flagged data:
# # 
# # To do so, we restart from the beginning, right after the flagging.
# # 
# # New flux value:
# # 
# # Flux density for J0238+1636 in SpW=0 (freq=1.40499e+09 Hz) is: 0.818094 +/- 0.0143805 (SNR = 56.8892, N = 44)
# # 
# # After applying the whole calibration again, we split field = 3 (CIG96) re-calibrated data, in the spw=0:
# 
# split(vis='cig96_13.ms', datacolumn='corrected', outputvis='cig96_13.corr.v2.ms', field='3', spw='0:399~2047')
# split(vis='cig96_13.ms', datacolumn='corrected', outputvis='cig96_13.corr.allfields.ms', field='1~3', spw='0')
# 
# # and redo uvcontsub:
# 
# # Since we have chopped out the first 400 channels where the RFI region was (see spw definition in the 
# # split command) the total number of channels have changed now: from 2048 we now have 1648.
# # Now, we remove the continuum from source by subtracting the channels from 500to910 and from 730to1300.
# 
# uvcontsub(vis='cig96_13.corr.v2.ms', field='0', fitspw='0:0~510;730~1300')
# 
# # Cleaning: weighting natural:
# 
# clean(vis='cig96_13.corr.v2.ms.contsub', imagename='cig96_13.corr.v2.contsub.natural.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[300,300], phasecenter='', cell='15.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)
# 
# 
# 
# 
# ################ END
# #
# #
# # IN CASE YOU WISH TO RESET THE WHOLE CALIBRATION, TAKE THE FOLLOWING STEPS:
# #
# # clearcal(vis='xxx')
# #
# # flagdata(vis='xxx',mode=unflag)
# #
# # Ready to redo calibration!
