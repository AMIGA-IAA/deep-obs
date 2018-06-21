# -*- coding: iso-8859-1 -*-
### EVLA DATA REDUCTION
### Project: 13A-341
### Dataset date: 10mar13 (first observation)
### Original dataset name: 13A-341.sb19189214.eb19331286.56361.695559375
### Renamed as: cig96_02.ms
###
### Configuration: D (2/3)

# ===============================================================================

### Import of EVLA data from SDM format:
########################################

# importevla(asdm='cig96_02', vis='cig96_02.ms')

# Original data file "AV282_B050723.xp1" gets now CASA format and a new name: "cig96_02.ms"

# ===============================================================================

### Listobs inspection:
#######################

listobs(vis='cig96_02.ms')

# Listobs output summary: 27 antennae, 9 spectral windows, 2 polarizations (RR, LL), 1 dummy scan, 
# 3 fields (2 calibrators + 1 target), 2048 channels of 7.8 kHz each for spw=0:

# Spectral Window 0 (line):
# SpwID  Name          #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) BBC Num  Corrs  
# 0	 A0C0#0        2048    TOPO    1404.995         7.812     16000.0      12  RR  LL

#  Fields: 4
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    NONE 0137+331=3C48       01:37:41.299431 +33.09.35.13299 J2000   0          31590
#  1    K    0137+331=3C48       01:37:41.299431 +33.09.35.13299 J2000   1         334854
#  2    D    J0238+1636          02:38:38.930108 +16.36.59.27470 J2000   2         199017
#  3    NONE CIG 96              02:15:27.600000 +06.00.09.00001 J2000   3        1655316

### BANDPASS and FLUX calibrator: 	1 0137+331 = 3C48
### PHASE calibrator: 			2 J0238+1636
### TARGET: 				3 CIG 96

# ===============================================================================

### Data inspection via plotms:
###############################

# Prior to anything else, there is dummy scan, field=0, so we remove it:

flagdata(vis='cig96_02.ms', mode='manual', field='0')

# We run plotms to inspect the general aspect of the data:

# plotms(vis='cig96_02.ms', xaxis='time', yaxis='amp', field='1', spw='0', coloraxis='field', avgchannel='1e9')

# ===============================================================================

# Antenna position correction:
##############################

gencal(vis='cig96_02.ms', caltable='cig96_02.ms.antpos', caltype='antpos')

# ===============================================================================

### Plot with the spatial distribution of the antennae:
#######################################################

plotants(vis='cig96_02.ms',figfile='plotants_cig96_02.png')

# ===============================================================================

### Flagging:
#############

# Log file provides some insight of the antennae status, corruption, etc.:

# 09Jul 10:30:26   09Jul 11:30:18   FRONT END         C130665       1.00   59.9
# Antenna(s) 19 (Data: Corrupted):
# All IFs weak at L-band due to a long-term receiver LNA issue.

# 09Jul 10:30:26   09Jul 11:30:18   FRONT END         C130665       0.50   29.9
# Antenna(s) 25 (Data: Corrupted):
# IFs C and D very weak at L-band due to long-term receiver LNA issue.

# Just in field 1 we could already observe weird behaviour of these antennas so we flag ea19 
# and ea25 in all fields:

flagdata(vis='cig96_02.ms', mode='manual', field='', spw='', antenna='ea25;ea19', flagbackup=True)

### Shadowing correction over the whole .ms:

flagdata(vis='cig96_02.ms', mode='shadow',  flagbackup=True)

# Percentage of data flagged in table selection: 25.0291%

### Zero clipping correction over the whole .ms:

flagdata(vis='cig96_02.ms', mode='clip', clipzeros=True, flagbackup=True)

# Percentage of data flagged in table selection: 0.574026%

### Field 1 (bandpass and flux calibrator): good behaviour in spw=0, except for antenna ea19(17) which is 
# showing low amplitudes in both RR and LL and for all baselines, but it has already been removed.

# Also, the baseline 6&18 shows bad behaviour, we remove it:

flagdata(vis='cig96_02.ms', mode='manual', field='1', spw='', antenna='ea07&ea20', flagbackup=True)

# Apart from this antenna,there is an intense RFI located around channel 385 (freq.=1.408 GHz):
# NRAO known RFIs: https://science.nrao.edu/facilities/vla/observing/RFI/L-Band
# Known RFI in 1.408 GHz: http://www.vla.nrao.edu/astro/RFI/rfi_plots/SUB4.gif
# "The narrow feature at 1408 MHz is internally generated."
# Galactic plane emission in 1.420 GHz: confirm with Min.

flagdata(vis='cig96_02.ms', mode='manual', field='1', spw='0:370~410', flagbackup=True)

# and it needs quacking of the first minutes:

flagdata(vis='cig96_02.ms', mode='quack', field='1', spw='0', quackinterval=180.0, quackmode='beg', flagbackup=True)

### Field 2 (phase calib.): there is also the intense RFI located around channel 385 (freq.=1.408 GHz) so
###  we flag it:

flagdata(vis='cig96_02.ms', mode='manual', field='2', spw='0:370~410', flagbackup=True)

# also, there is some Galactic emission around channel 1955 (freq.=1.4203 GHz), we remove it:

flagdata(vis='cig96_02.ms', mode='manual', field='2', spw='0:1945~1965', flagbackup=True)

# Quacking is needed in each of the 3 scans 3, 7 and 10:

flagdata(vis='cig96_02.ms', mode='quack', field='2', spw='0', scan='3', quackinterval=60.0, quackmode='beg', flagbackup=True)
flagdata(vis='cig96_02.ms', mode='quack', field='2', spw='0', scan='7', quackinterval=15.0, quackmode='beg', flagbackup=True)
flagdata(vis='cig96_02.ms', mode='quack', field='2', spw='0', scan='10', quackinterval=15.0, quackmode='beg', flagbackup=True)

# there is algo a baseline between ea07 and ea27 that is not working well in RR in scan 7

flagdata(vis='cig96_02.ms', mode='manual', field='2', spw='0', scan='7', antenna='ea07&ea27', correlation='RR', flagbackup=True)

### Field 3 (CIG96): shows the same RFI in channel ~385 (freq.~1.408 GHz) as in field 2; also the same 
# intense emission in channel ~1980 (freq.~1.4203 GHz).

# We flag the RFI and the galactic emission:

flagdata(vis='cig96_02.ms', mode='manual', field='3', spw='0:370~410', flagbackup=True)
flagdata(vis='cig96_02.ms', mode='manual', field='3', spw='0:1945~1965', flagbackup=True)

### Quacking: fields 2 and 3 need quacking at the beginning:

flagdata(vis='cig96_02.ms', mode='quack', field='2,3', spw='0', quackinterval=15.0, quackmode='beg', flagbackup=True)

### And we add a autoflagging for each of the calibrators (do NOT mix them in one command, use two to 
### avoid flagging the calibrator with highest amplitudes):

# flagdata(vis='cig96_02.ms',  mode='rflag', field='1', spw='0', ntime='scan', action='apply')
# flagdata(vis='cig96_02.ms',  mode='rflag', field='2', spw='0', ntime='scan', action='apply')

# Enough flagging for now.

# ===============================================================================

### Calibration:
### ============

### Reference antenna selection:
################################

# After checking plotants, we select one of the central-most and has shown no problems:

refant='ea26'

# ===============================================================================

### Opacity and gaincurve corrections:
######################################

# Should we do opacities correction?
# myTau = plotweather(vis='cig96_02.ms', doPlot=T)
# display cig96_02.ms.plotweather.png
# gencal(vis='cig96_02.ms', caltype='opac', caltable='cig96_02.ms.opcac', parameter=myTau)
# ANSWER: no, not needed for HI observations, the atmosphere barely has any influence with HI. Only 
# ionosphere in very flat elevations and at dawn/sunset.

# Should we do gaincurve correction?
# Takes into account the elevations of each antenna, see elevations plot: CIG96 and the phase calibrator 
# have a large elevation range.
#
# gencal(vis='cig96_02.ms', caltype='gceff', caltable='cig96_02.ms.gaincurve')
#
# ANSWER: no, it might even harm the observations due to the possible crosstalk at low elevations 
# in C and D configs. and the model it is based on.

# ===============================================================================

### Delays correction:
######################

# Prior to the bandpass correction, we need to correct the phase variations with time (they are large). 
# We use the bandpass calibrator in field 1 with the antenna position correction table:

gaincal(vis='cig96_02.ms', caltable='cig96_02.ms.delays', field='1', refant='ea26', gaintype='K', gaintable=['cig96_02.ms.antpos'])

# ===============================================================================

### Flux density:
#################

# Our flux calibrator is 3C48, in field = 0.
# First, we check the list of available models:

setjy(vis='cig96_02.ms', listmodels=T)

# Candidate modimages (*) at pepino (dae66) in path:
# /Applications/CASA.app/Contents/data/nrao/VLA/CalModels/
#
# Candidate modimages (*) at NRAO local workstations in path:
# /home/casa/packages/RHEL5/release/casapy-42.1.29047-001-1-64b/data/nrao/VLA/CalModels/
# 
# The model chosen has to be in accordance with the calibrator and band selected: 3C48 in L band:

setjy(vis='cig96_02.ms', field='1', modimage='/mnt/scops/data/data/paramimo/casapy-42.2.30986-1-64b/data/nrao/VLA/CalModels/3C48_L.im')

# ===============================================================================

### Bandpass calibration:
#########################

# For the bandpass and flux (field=1) and phase (field=2) calibrators we use solution interval 
# time of solint='5s':

gaincal(vis='cig96_02.ms', caltable='cig96_02.ms.bpphase5s', field='1', refant='ea26', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_02.ms.antpos','cig96_02.ms.delays' ])

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

# plotcal(caltable='cig96_02.ms.bpphase5s', xaxis='time', yaxis='phase')

# Apply phase solutions on the fly:

bandpass(vis='cig96_02.ms', caltable='cig96_02.ms.bandpass5s', field='1', refant='ea26', solint='inf', solnorm=T, minsnr=10.0, minblperant=3, gaintable=['cig96_02.ms.bpphase5s', 'cig96_02.ms.antpos'], interp=['nearest'])

# We check again the solutions:

# plotcal(caltable='cig96_02.ms.bpphase5s', field='1', xaxis='time', yaxis='phase')  

# Difference in phase is less, phase shows a flatter behaviour now.

# ===============================================================================


### Phase and amplitude calibration:
####################################

# Phase calibration for the calibrators, fields 1 and 2 (spw 0, where the line is):

gaincal(vis='cig96_02.ms',caltable='cig96_02.ms.intphase', field='1,2', refant='ea26', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_02.ms.antpos', 'cig96_02.ms.delays', 'cig96_02.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, 
# to use FOR THE TARGET later on:

gaincal(vis='cig96_02.ms', caltable='cig96_02.ms.scanphase', field='1,2', refant='ea26', calmode='p', solint='inf', minsnr=5.0, gaintable=['cig96_02.ms.antpos', 'cig96_02.ms.delays', 'cig96_02.ms.bandpass5s'])

# Derive amplitude solutions:

gaincal(vis='cig96_02.ms', caltable='cig96_02.ms.amp', field='1,2',  refant='ea26', calmode='ap', solint='inf', minsnr=5.0, gaintable=['cig96_02.ms.antpos', 'cig96_02.ms.delays', 'cig96_02.ms.bandpass5s', 'cig96_02.ms.intphase'])

# I check the tables with plotcal and some things do not look as expected: the first data show very 
# low intensity compares to the rest and some graphs show only one point:

# plotcal(caltable='cig96_02.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)
# plotcal(caltable='cig96_02.ms.amp', xaxis='time', yaxis='amp', iteration='baseline', subplot=331)

# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal 
# in terms of future application of the calibration to the data, (i.e., it's not an incremental table) 
# UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_02.ms', caltable='cig96_02.ms.amp', fluxtable='cig96_02.ms.flux', reference='1', transfer='2', incremental=False)

# Result:
# 
# Flux density for J0238+1636 in SpW=0 (freq=1.405e+09 Hz) is: 0.811138 +/- 0.00458405 (SNR = 176.948, N = 46)

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

applycal(vis='cig96_02.ms', field='1', gaintable=['cig96_02.ms.antpos', 'cig96_02.ms.delays', 'cig96_02.ms.bandpass5s', 'cig96_02.ms.intphase', 'cig96_02.ms.flux'], gainfield=['','','','1',''], calwt=F)

# time.sleep(5)

applycal(vis='cig96_02.ms', field='2', gaintable=['cig96_02.ms.antpos', 'cig96_02.ms.delays', 'cig96_02.ms.bandpass5s', 'cig96_02.ms.intphase', 'cig96_02.ms.flux'], gainfield=['','','','2',''], calwt=F)

# time.sleep(5)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_02.ms', field='3', gaintable=['cig96_02.ms.antpos', 'cig96_02.ms.delays', 'cig96_02.ms.bandpass5s', 'cig96_02.ms.scanphase', 'cig96_02.ms.flux'], gainfield=['','','','2',''], calwt=F)

# ===============================================================================

### Splitting of CIG96:
#######################

# Splitting of field = 3 (CIG96) calibrated data, in the spw=0:

split(vis='cig96_02.ms', datacolumn='corrected', outputvis='cig96_02.corr.ms', field='3', spw='0:399~2047')

# Now, for the new corr.ms file:
# field = 3 is stored as: field = 0  
# spw = 0 keeps being spw = 0

# Emission in the higher frequency:
# 
# The emission line is seen in 1.4203 GHz. We split the fields in spw=0 from channel 780 to 1600:

split(vis='cig96_02.ms', datacolumn='corrected', outputvis='cig96_02.corr.emission.ms', field='3', spw='0:780~1600')

# Splitting CIG96 fields in spw=0 with time binning of 20s and channel width of 10 channels:

split(vis='cig96_02.ms', datacolumn='corrected', outputvis='cig96_02.corr.20s.10chan.ms', width=10, timebin='20s', field='3', spw='0')

# ===============================================================================

### UV continuum subtraction:
#############################

# Since we have chopped out the first 400 channels where the RFI region was (see spw definition in the 
# split command) the total number of channels have changed now: from 2048 we now have 1648.
# Now, we remove the continuum from source by subtracting the channels from 500to910 and from 730to1300.

uvcontsub(vis='cig96_02.corr.ms', field='0', fitspw='0:0~510;730~1300')

# Also, for the second emission in 1.4203 GHz, the emission is in a very noisy region so the continuum 
# subtraction won't be perfect, likely more noisy on higher frequencies since the continuum is
# selected in channels 780~1000:

uvcontsub(vis='cig96_02.corr.emission.ms', field='0', fitspw='0:0~510;730~1300')

# For the time-channel binned file, we subtract as follows:

uvcontsub(vis='cig96_02.corr.20s.10chan.ms', field='0', fitspw='0:45~93;110~181')

# ===============================================================================

### Clean of the continuum subtracted data:
########################################### 

# Weighting natural:

clean(vis='cig96_02.corr.ms.contsub', imagename='cig96_02.corr.contsub.natural.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[300,300], phasecenter=0, cell='15.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)

# Beam size: 

# WARN MFCleanImageSkyModel   Clean not converging   
# Successfully deconvolved image
# Beam used in restoration: 99.5838 by 47.044 (arcsec) at pa -56.6383 (deg) 

# Weighting uniform:

clean(vis='cig96_02.corr.ms.contsub', imagename='cig96_02.corr.contsub.uniform.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[300,300], phasecenter=0, cell='15.0arcsec', restfreq='1.4204GHz', weighting='uniform', usescratch=T)

# Beam size: 

# Weighting Briggs rob=0.0:

clean(vis='cig96_02.corr.ms.contsub', imagename='cig96_02.corr.contsub.rob0.0.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[300,300], phasecenter=0, cell='15.0arcsec', restfreq='1.4204GHz', weighting='briggs', robust='0.0', usescratch=T)

# Beam size: 

# Weighting Briggs rob=2.0:

clean(vis='cig96_02.corr.ms.contsub', imagename='cig96_02.corr.contsub.rob2.0.line', field='0', spw='0:510~730', mode='channel', start=0, niter=10000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[300,300], phasecenter=0, cell='15.0arcsec', restfreq='1.4204GHz', weighting='briggs', robust='2.0', usescratch=T)

# Beam size: 

### Noise level via viewer task:
################################

---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(cig96_02.corr.contsub.line.image)
     Frequency       Velocity         Stokes BrightnessUnit       BeamArea 
 1.41305e+09Hz    1553.49km/s              I        Jy/beam        4.44663 
          Npts            Sum    FluxDensity           Mean            Rms 
          1155   2.966822e-02   6.672063e-03   2.568677e-05   4.262555e-04 
       Std dev        Minimum        Maximum   region count 
  4.256651e-04  -1.101821e-03   1.550763e-03              1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

# For the time-channel binned file, we clean as follows:

# Weighting natural:

clean(vis='cig96_02.corr.20s.10chan.ms.contsub', imagename='cig96_02.corr.20s.10chan.ms.contsub.natural.line.image', field='0', spw='0:93~110', mode='channel', start=0, niter=7000, threshold='1.2mJy', interactive=T, npercycle=100, imsize=[300,300], phasecenter=0, cell='15.0arcsec', restfreq='1.4204GHz', weighting='natural', usescratch=T)



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
