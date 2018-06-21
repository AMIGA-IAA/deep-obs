### VLA DATA REDUCTION
### Project: AV276
### Dataset date: 19jul04
### Configuration: D

# ===============================================================================

### Import of VLA data from xp format:
######################################

importvla(archivefiles='AV282_B050723.xp1', vis='cig96_vla02.ms',bandname='L')

# Original data file "AV276_A040719.xp2" gets now CASA format and a new name: "cig96_vla02.ms"

# ===============================================================================

### Listobs inspection:
#######################

listobs(vis='cig96_vla02.ms')

# Listobs output summary:

# SpwID  Name                                   #Chans   Frame   Ch1(MHz)  ChanWid(kHz)  TotBW(kHz)  Corrs  
# 0      63*48.8 kHz channels @ 1.41 GHz (BARY)     63   BARY    1408.078        48.830      3076.3  RR  LL
# 1      63*48.8 kHz channels @ 1.4 GHz (BARY)      63   BARY    1397.036        48.828      3076.2  RR  LL
# 2      63*48.8 kHz channels @ 1.41 GHz (BARY)     63   BARY    1411.483        48.824      3075.9  RR  LL
# 3      63*48.8 kHz channels @ 1.39 GHz (BARY)     63   BARY    1391.391        48.824      3075.9  RR  LL
# 4      63*48.8 kHz channels @ 1.39 GHz (BARY)     63   BARY    1387.969        48.824      3075.9  RR  LL      

# Fields = 10
# 0 1400+621
# 1 0834+555 
# 2 CIG240  
# 3 CIG551  
# 4 0137+331
# 5 0318+164
# 6 CIG0096
# 7 CIG0123 
# 8 0841+708 
# 9 CIG0421  

### Fields of interest:
# 4 0137+331
# 5 0318+164
# 6 CIG0096

### BANDPASS and FLUX calibrator: 	4 0137+331 = 3C48
### PHASE calibrator: 				5 0318+164

# ===============================================================================

### Data inspection via plotms:
###############################

plotms(vis='cig96_vla02.ms',xaxis='time',yaxis='amp',field='4~6',coloraxis='field',acgchannel='1000',timerange='09:10:00~18:50:00')

# Field 4 (bandpass and flux caibrator) shows normal amplitudes in the 64 channels except for baseline: 13&19
# Fields 5 (phase calib.) and 6 (CIG96) seem to have a lot of noise.

# Antenna position correction: gencal-caltype='antpos' does not work with VLA data. Instead, we need to use option 'antposvla'.
# 
# Example:
# 
# gencal(vis='test.ms',caltable='test.G',caltype='antposvla',
#              antenna='ea09,ea10',
#             parameter=[0.01,0.02,0.03, -0.03,-0.01,-0.02])
# 
#          --> Antenna position corrections (in the traditional VLA-centric
#               frame) will be introduced in meters for
#               antenna ea09 (dBx=0.01, dBy=0.02, dBz=0.03) and for
#               antenna ea10 (dBx=-0.03, dBy=-0.01, dBz=-0.02)
#               These offsets will be rotated to the ITRF frame before
#               storing them in the caltable.
#               See the above example for caltype='ph' for details
#               of the sign convention adopted when applying antpos 
#               corrections.
#
# In the VLA baseline corrections website:
# http://www.vla.nrao.edu/astro/archive/baselines/
# TO BE CHECKED
       
# gencal(vis='cig96_vla02.ms',
# 	caltable='cig96_vla02.ms.antposvla',
#	caltype='antposvla',
#	antenna='')

# ===============================================================================

### Plot with the spatial distribution of the antennae:
#######################################################

plotants(vis='cig96_vla02.ms',figfile='cig96_vla02.ms.plotants.png')

# ===============================================================================

### Flagging:
### =========

# Log file does not provide with any insight of the antennae status, corruption, etc.

# There is no dummy scan, no scan to remove.

### Shadowing correction:

flagdata(vis='cig96_vla02.ms', mode='shadow',  flagbackup=True)
	
# INFO Shadow	=> Percentage of data flagged in table selection: 1.40264%

### Zero clipping correction:

flagdata(vis='cig96_vla02.ms', mode='clip', clipzeros=True, flagbackup=True)

# INFO Clip	=> Percentage of data flagged in table selection: 2.01485e-06%

### Field 4 (bandpass and flux calibrator) flagging of the initial data of the scans 14, 15 and 16:

flagdata(vis='cig96_vla02.ms', mode='manual', field='4', spw='2', scan='14', timerange='09:21:40~09:21:58')

flagdata(vis='cig96_vla02.ms', mode='manual', field='4', spw='3', scan='15', timerange='09:32:00~09:32:20')

flagdata(vis='cig96_vla02.ms', mode='manual', field='4', spw='4', scan='16', timerange='09:37:20~09:37:40')

# flagging of weird behaviour baselines:

flagdata(vis='cig96_vla02.ms', mode='manual', field='4', spw='2', scan='14', antenna='VA19&VA28', correlation='RR')

flagdata(vis='cig96_vla02.ms', mode='manual', field='4', spw='2', scan='14', antenna='VA07&VA21', correlation='LL')

flagdata(vis='cig96_vla02.ms', mode='manual', field='4', spw='3,4', scan='15,16', antenna='VA21&VA26', correlation='LL')

# and we get rid of data above 4 Jy in LL polarization:

flagdata(vis='cig96_vla02.ms', mode='clip', field='4', clipminmax=[-100,4], clipoutside=T, correlation='LL')


### Field 5 (phase calib.) flagging:

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='2', scan='43', timerange='14:56:44~14:56:46')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='2', scan='43', timerange='14:56:54~14:56:56')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='3', scan='44', timerange='14:57:34~14:57:36')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='3', scan='44', timerange='14:57:44~14:57:46')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='3', scan='44', timerange='14:57:54~14:57:56')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='3', scan='44', timerange='14:58:04~14:58:06')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='3', scan='44', timerange='14:58:14~14:58:16')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='3', scan='44', timerange='14:58:24~14:58:26')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='3', scan='44', timerange='14:58:14~34:58:36')

# flagdata(vis='cig96_vla02.ms', mode='manual', field='5', spw='2', scan='37', timerange='14:58:14~34:58:36')


### all of the above seem to be a consequence of antennae VA16 (3) and VA26 (23) causing problems in scans 37, 43 and 44. And more specifically, only in the LL polarization. 


### QUESTION: can the flagging be done just by the complete removal of VA16 and VA26 in LL from the respective scans? Or as a clipping of the data above 4 Jy? The lower intensity data seem to be fine so, is this an option?

flagdata(vis='cig96_vla02.ms', mode='manual', field='5', correlation='LL', antenna='VA16,VA26')

# or:
# flagdata(mode='clip', field='5', clipminmax=[-100,4], clipoutside=T, correlation='LL')

# Quacking of the first seconds of all the scans of field 5:

flagdata(vis='cig96_vla02.ms', mode='quack', field='5', quackinterval=10.0, quackmode='beg')

flagdata(vis='cig96_vla02.ms', mode='manual', field='5', timerange='10:51:20~10:51:30')

#############
# NOTE: both calibrators show scans with different average amplitudes that depend on the spectral window (frequency) at which they have been observed. Our source has been observed in spectral window 2, so we stick to spw = 2 for the calibrators too (where the freq. for spw=2 is 1.411 GHz). The other 2 spw's have lower frequencies.
#############


### Field 6 (CIG96) flagging: antennae VA16(13) and VA26(23) give problems in polarization LL, not RR.

# Same as for field 5, is it optimal to just remove the antennae in LL polarization?

flagdata(vis='cig96_vla02.ms', mode='manual', field='6', correlation='LL', antenna='VA16,VA26')

# QUESTION: possible data clipping above 0.14 Jy?

# flagdata(mode='clip', field='6', clipminmax=[-100,0.14], clipoutside=T, correlation='LL')

# Quacking of the first seconds of all the scans of field 6:

flagdata(vis='cig96_vla02.ms', mode='quack', field='6', quackinterval=5.0, quackmode='beg')

# enough flagging for now.

# ===============================================================================

### Calibration:
### ============

### Reference antenna selection:
################################

# After checking plotants, we select VA25 because it is quite central (located on the N2 position):

refant='VA25'

# ===============================================================================

### Opacity and gaincurve corrections:
######################################

# Should we do opacities correction?
# myTau = plotweather(vis='cig96_vla02.ms', doPlot=T)
# display cig96_vla02.ms.plotweather.png
# gencal(vis='cig96_vla02.ms', caltype='opac', caltable='cig96_vla02.ms.opcac', parameter=myTau)
# ANSWER: no, not needed for HI observations, the atmosphere barely has any influence with HI. Only ionosphere in very flat elevations and at dawn/sunset.

# Should we do gaincurve correction?
# Takes into account the elevations of each antenn, see elevations plot: CIG96 and the phase calibrator have a large elevation range.
# gencal(vis='cig96_vla02.ms', caltype='gceff', caltable='cig96_vla02.ms.gaincurve')
# ANSWER: no, it might even harm the observations due to the possible crosstalk at low elevations and the model it is based on.

# ===============================================================================

### Delays correction:
######################

# Prior to the bandpass correction, we need to correct the phase variations with time (they are large). We use the bandpass calibrator in field 4.

gaincal(vis='cig96_vla02.ms', caltable='cig96_vla02.ms.delays', field='4', refant='VA25', gaintype='G')

# we have added no gaintables to the gaincal command because we have NOT calculated any yet: "antpos" is not defined yet, "gaincurve" is not needed, and "opacities" is not needed either.

# The answer of the previous command is:

# 2014-05-05 16:02:59 INFO gaincal	##########################################
# 2014-05-05 16:02:59 INFO gaincal	##### Begin Task: gaincal            #####
# 2014-05-05 16:02:59 INFO gaincal	gaincal(vis="cig96_vla02.ms",caltable="cig96_vla02.ms.delays",field="4",spw="",
# 2014-05-05 16:02:59 INFO gaincal	        intent="",selectdata=True,timerange="",uvrange="",antenna="",
# 2014-05-05 16:02:59 INFO gaincal	        scan="",observation="",msselect="",solint="inf",combine="",
# 2014-05-05 16:02:59 INFO gaincal	        preavg=-1.0,refant="VA25",minblperant=4,minsnr=3.0,solnorm=False,
# 2014-05-05 16:02:59 INFO gaincal	        gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,
# 2014-05-05 16:02:59 INFO gaincal	        npointaver=3,phasewrap=180.0,gaintable=[''],gainfield=[''],interp=[''],
# 2014-05-05 16:02:59 INFO gaincal	        spwmap=[],gaincurve=False,opacity=[],parang=False)
# 2014-05-05 16:02:59 INFO gaincal	Opening MS: cig96_vla02.ms for calibration.
# 2014-05-05 16:03:00 INFO gaincal	Initializing nominal selection to the whole MS.
# 2014-05-05 16:03:00 INFO gaincal	Beginning selectvis--(MSSelection version)-------
# 2014-05-05 16:03:00 INFO gaincal	Reseting solve/apply state
# 2014-05-05 16:03:00 INFO gaincal	Performing selection on MeasurementSet
# 2014-05-05 16:03:00 INFO gaincal	 Selecting on field: '4'
# 2014-05-05 16:03:00 INFO gaincal	By selection 393900 rows are reduced to 40300
# 2014-05-05 16:03:01 INFO gaincal	Frequency selection: Selecting all channels in all spws.
# 2014-05-05 16:03:01 INFO gaincal	Beginning setsolve--(MSSelection version)-------
# 2014-05-05 16:03:01 INFO gaincal	Arranging to SOLVE:
# 2014-05-05 16:03:01 INFO gaincal	.   G Jones: table=cig96_vla02.ms.delays append=false solint=inf refant='VA25' minsnr=3 apmode=AP solnorm=false
# 2014-05-05 16:03:01 INFO gaincal	Beginning solve-----------------------------
# 2014-05-05 16:03:01 INFO gaincal	The following calibration terms are arranged for apply:
# 2014-05-05 16:03:01 INFO gaincal	.   (None)
# 2014-05-05 16:03:01 INFO gaincal	The following calibration term is arranged for solve:
# 2014-05-05 16:03:01 INFO gaincal	.   G Jones: table=cig96_vla02.ms.delays append=false solint=inf refant='VA25' minsnr=3 apmode=AP solnorm=false
# 2014-05-05 16:03:01 INFO gaincal	Solving for G Jones
# 2014-05-05 16:03:01 INFO gaincal	For solint = inf, found 3 solution intervals.
# 2014-05-05 16:03:04 INFO gaincal	  Found good G Jones solutions in 3 slots.
# 2014-05-05 16:03:04 INFO gaincal	Applying refant: VA25
# 2014-05-05 16:03:04 INFO gaincal	Writing solutions to table: cig96_vla02.ms.delays
# 2014-05-05 16:03:04 INFO gaincal	Finished solving.
# 2014-05-05 16:03:04 INFO gaincal	Calibration solve statistics per spw:  (expected/attempted/succeeded):
# 2014-05-05 16:03:04 INFO gaincal	  Spw 0: 0/0/0
# 2014-05-05 16:03:04 INFO gaincal	  Spw 1: 0/0/0
# 2014-05-05 16:03:04 INFO gaincal	  Spw 2: 1/1/1
# 2014-05-05 16:03:04 INFO gaincal	  Spw 3: 1/1/1
# 2014-05-05 16:03:04 INFO gaincal	  Spw 4: 1/1/1
# 2014-05-05 16:03:04 INFO gaincal	##### End Task: gaincal              #####
# 2014-05-05 16:03:04 INFO gaincal	##########################################


# ===============================================================================

### Flux density:
#################

# Our flux calibrator is 3C48, in field = 4.
# First, we check the list of available models:

setjy(vis='cig96_vla02.ms', listmodels=T)

# Candidate modimages (*) at Pablo's computer in path:
#
# /Applications/CASA.app/Contents/data/nrao/VLA/CalModels/
#
# The model chosen has to be in accordance to the band selected, in this case: L band:

setjy(vis='cig96_vla02.ms', field='4', modimage='/Applications/CASA.app/Contents/data/nrao/VLA/CalModels/3C48_L.im')

# ===============================================================================

### Bandpass calibration:
#########################

# Real bandpass calibration. For the BANDPASS and FLUX(field=4) and PHASE(field=5) calibrators we use solution interval time of solint='5s':

gaincal(vis='cig96_vla02.ms',caltable='cig96_vla02.ms.bpphase5s', field='4', spw='2', refant='VA25', calmode='p', solint='5s', minsnr=5.0)
	
# The solution interval of 5 seconds has been calculated following with the VLA exp.time calculator using the parameters:
# Freq. = 1.42 GHz
# Medium elevation (summer time)
# Bandwith freq. = 3,076 KHz (see listobs)
# RMS noise = 20 mJy (more than 3.0 in SNR)
# 
# The calculator estimates less than one second (0.2s) is enough to get such SNR or even higher so we set solint=5s since 5s is the shortest integration time of our data. This should mean that there should be a solution for all the intervals.

# We see the phase VS time figures to see that the phase is now waaayyy flatter:

plotcal(caltable='cig96_vla02.ms.bpphase5s', xaxis='time', yaxis='phase')  

# Apply phase solutions on the fly:

bandpass(vis='cig96_vla02.ms', caltable='cig96_vla02.ms.bandpass5s', field='4', spw='2', refant='VA25', solint='inf', solnorm=T, minsnr=10.0, minblperant=3, gaintable=['cig96_vla02.ms.bpphase5s'], interp=['nearest'])

# We check again the solutions:

plotcal(caltable='cig96_vla02.ms.bpphase5s', field='4', spw='2', xaxis='time', yaxis='phase')  

# Difference in phase is less, phase shows a much flatter behaviour now.

# ===============================================================================

### Phase and amplitude calibration:
####################################

# Phase calibration for the calibrators, fields 4 and 5 (spw 2, where the line is):

gaincal(vis='cig96_vla02.ms',caltable='cig96_vla02.ms.intphase', field='4,5', spw='2', refant='VA25', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_vla02.ms.delays', 'cig96_vla02.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, to use FOR THE TARGET later on:

# time.sleep(5)

gaincal(vis='cig96_vla02.ms', caltable='cig96_vla02.ms.scanphase', field='4,5', spw='2', refant='VA25', calmode='p', solint='inf', minsnr=5.0, gaintable=['cig96_vla02.ms.delays', 'cig96_vla02.ms.bandpass5s'])

# Derive amplitude solutions:

# time.sleep(5)

gaincal(vis='cig96_vla02.ms', caltable='cig96_vla02.ms.amp', field='4,5', spw='2', refant='VA25', calmode='ap', solint='inf', minsnr=3.0, gaintable=['cig96_vla02.ms.delays', 'cig96_vla02.ms.bandpass5s', 'cig96_vla02.ms.intphase'])

# I check the tables with plotcal and some things do not look as expected: the first data show very low intensity compares to the rest and some graphs show only one point:

# plotcal(caltable='cig96_vla02.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)
# plotcal(caltable='cig96_vla02.ms.amp', xaxis='time', yaxis='amp', iteration='baseline', subplot=331)


# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal in terms of future application of the calibration to the data, (i.e., it's not an incremental table) UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_vla02.ms', caltable='cig96_vla02.ms.amp', fluxtable='cig96_vla02.ms.flux', reference='4', incremental=False)

# Result:
# Flux density for 0318+164 in SpW=2 is: 7.79475 +/- 0.00296097 (SNR = 2632.5, N = 48)

# ===============================================================================

### Application of the calibration to the .ms file:
###################################################

# Note: In all applycal steps we set calwt=F. It is very important to turn off this parameter which determines if the weights are calibrated along with the data. Data from antennae with better receiver performance and/or longer integration times should have higher weights, and it can be advantageous to factor this information into the calibration. During the VLA era, meaningful weights were available for each visibility. However, at the time of this observation, the VLA was not yet recording the information necessary to calculate meaningful weights. Since these data weights are used at the imaging stage you can get strange results from having calwt=T when the input weights are themselves not meaningful, especially for self-calibration on resolved sources (your flux calibrator and target, for example).

applycal(vis='cig96_vla02.ms', field='4', gaintable=['cig96_vla02.ms.delays', 'cig96_vla02.ms.bandpass5s', 'cig96_vla02.ms.intphase', 'cig96_vla02.ms.flux'], calwt=F)

# time.sleep(5)

applycal(vis='cig96_vla02.ms', field='4', gaintable=['cig96_vla02.ms.delays', 'cig96_vla02.ms.bandpass5s', 'cig96_vla02.ms.intphase', 'cig96_vla02.ms.flux'], calwt=F)

# time.sleep(5)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_vla02.ms', field='6', gaintable=['cig96_vla02.ms.delays', 'cig96_vla02.ms.bandpass5s', 'cig96_vla02.ms.scanphase', 'cig96_vla02.ms.flux'], gainfield=['','','5','5'], calwt=F)


# ===============================================================================

### Imaging of CIG96:
#####################

# Splitting of field = 6 (CIG96) calibrated data, in the spw=2:

split(vis='cig96_vla02.ms', datacolumn='corrected', outputvis='cig96_vla02.corr.ms', field='6', spw='2')

# Now, for the new corr.ms file:
# Field = 6 is stored as: field = 0  
# spw = 2 is stored as: spw = 0.

# Splitting CIG96 fields in spw=2 with channel binning of 8 and time binning of 20s:

split(vis='cig96_vla02.ms', datacolumn='corrected', outputvis='cig96_vla02.corr.8chann.20s.ms', timebin='20s', width=8, field='6', spw='2')

# Splitting CIG96 fields in spw=2 with time binning of 20s:

split(vis='cig96_vla02.ms', datacolumn='corrected', outputvis='cig96_vla02.corr.20s.ms', timebin='20s', field='6', spw='2')


### UV continuum subtraction:
#############################

# With plotms we can see
# Remove continuum from source by subtracting the channels from 0to900 and from 1150to1700 (the fact they are dividided by 10 is because they have been binned in channel in the split command above):

uvcontsub(vis='cig96_vla02.corr.ms', field='0', fitspw='0:4~18;46~56')


### Cleaning:
#############

# We do different cleanings changing the weighting, number of iterations, cell size

# Configuration 1: weighting='natural', niter=1000, cell='15.0arcsec'

clean(vis='cig96_vla02.corr.ms.contsub', imagename='cig96_vla02.corr.contsub.conf1.line', field='0', spw='0:18~46', mode='channel', niter=1000, start=0, interactive=T, imsize=[150,150], phasecenter=0, restfreq='1.42GHz', cell='15.0arcsec')

# Restorign beam: 65,31" 52,04"

# Image moments 0 and 1:
########################

immoments(imagename='cig96_vla02.corr.contsub.conf1.line.image', moments=[0, 1], axis='spectral', chans='18~46', stokes='I')

# Configuration 1 noise level via viewer task:

---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(cig96_vla02.corr.contsub.conf1.line.image.integrated)
        Stokes       Velocity          Frame        Doppler      Frequency 
             I    1612.58km/s           BARY          RADIO    1.41236e+09 
BrightnessUnit       BeamArea           Npts            Sum    FluxDensity 
  Jy/beam.km/s        17.1167           1325   1,112999e+01   6,502434e-01 
          Mean            Rms        Std dev        Minimum        Maximum 
  8,399993e-03   5,917827e-02   5,860119e-02  -1,459343e-01   1,921714e-01 
  region count 
             1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



# Configuration 2: weighting='natural', niter=10000, cell='13.0arcsec'

clean(vis='cig96_vla02.corr.ms.contsub', imagename='cig96_vla02.corr.contsub.conf2.line', field='0', spw='0:18~46', mode='channel', niter=5000, start=0, interactive=T, imsize=[150,150], phasecenter=0, restfreq='1.42GHz', cell='13.0arcsec')
		
# Restorign beam: 66.0728" 51.7902"


# Image moments 0 and 1:
########################

immoments(imagename='cig96_vla02.corr.contsub.conf2.line.image', moments=[0, 1], axis='spectral', chans='18~46', stokes='I')
	
	
# Configuration 2 noise level via viewer task:

---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(cig96_vla02.corr.contsub.conf2.line.image.integrated)
        Stokes       Velocity          Frame        Doppler      Frequency 
             I    1612.58km/s           BARY          RADIO    1.41236e+09 
BrightnessUnit       BeamArea           Npts            Sum    FluxDensity 
  Jy/beam.km/s        22.9429           2040  -2,646981e+01  -1,153726e+00 
          Mean            Rms        Std dev        Minimum        Maximum 
 -1,297540e-02   5,827867e-02   5,682980e-02  -1,680215e-01   1,572564e-01 
  region count 
             1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


# Configuration 3: weighting='briggs' (aka robust), niter=10000, cell='13.0arcsec'

clean(vis='cig96_vla02.corr.ms.contsub',
	imagename='cig96_vla02.corr.contsub.conf3.line',
	field='0',
	spw='0:18~46',
	mode='channel',
	niter=10000,
	start=0,
	interactive=T,
	imsize=[150,150],
	phasecenter=0,
	weighting='briggs',
	restfreq='1.42GHz',
	cell='13.0arcsec')

# Restorign beam: 44.96" x 41.54"

# Image moments 0 and 1:
########################

immoments(imagename='cig96_vla02.corr.contsub.conf3.line.image', moments=[0, 1], axis='spectral', chans='18~46', stokes='I')

# Configuration 3 noise level via viewer task:

---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(cig96_vla02.corr.contsub.conf3.line.image.integrated)
        Stokes       Velocity          Frame        Doppler      Frequency 
             I    1612.58km/s           BARY          RADIO    1.41236e+09 
BrightnessUnit       BeamArea           Npts            Sum    FluxDensity 
  Jy/beam.km/s        12.5254           1135   6,712193e+00   5,358875e-01 
          Mean            Rms        Std dev        Minimum        Maximum 
  5,913826e-03   6,557953e-02   6,534113e-02  -1,930735e-01   2,045056e-01 
  region count 
             1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


# Configuration 4: weighting='briggs' (aka robust), niter=10000, cell='9.0arcsec'

clean(vis='cig96_vla02.corr.ms.contsub',
	imagename='cig96_vla02.corr.contsub.conf4.line',
	field='0',
	spw='0:18~46',
	mode='channel',
	niter=10000,
	start=0,
	interactive=T,
	imsize=[150,150],
	phasecenter=0,
	weighting='briggs',
	restfreq='1.42GHz',
	cell='9.0arcsec')

# Restorign beam: 43.95" x 39.16" 

# Image moments 0 and 1:
########################

immoments(imagename='cig96_vla02.corr.contsub.conf4.line.image', moments=[0, 1], axis='spectral', chans='18~46', stokes='I')

# Configuration 4 noise level via viewer task:

---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(cig96_vla02.corr.contsub.conf4.line.image.integrated)
        Stokes       Velocity          Frame        Doppler      Frequency 
             I    1612.58km/s           BARY          RADIO    1.41236e+09 
BrightnessUnit       BeamArea           Npts            Sum    FluxDensity 
  Jy/beam.km/s        24.0748            987  -1,471879e+00  -6,113764e-02 
          Mean            Rms        Std dev        Minimum        Maximum 
 -1,491266e-03   5,903086e-02   5,904193e-02  -1,842583e-01   1,673624e-01 
  region count 
             1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

################ END
#
#
#
