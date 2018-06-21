### VLA DATA REDUCTION
### Project: AV282
### Dataset date: 23jul05
### Configuration: C

# ===============================================================================

### Import of VLA data from xp format:
######################################

importvla(archivefiles='AV282_B050723.xp1', vis='cig96_vla282.ms',bandname='L')

# Original data file "AV282_B050723.xp1" gets now CASA format and a new name: "cig96_vla282.ms"

# ===============================================================================

### Listobs inspection:
#######################

listobs(vis='cig96_vla282.ms')

# Listobs output summary: 24 antennae, one spectral window (spw='0'), 3 sources (2 calibrators + 1 target), 63 channels of 48.8 kHz each

# Spectral Windows:  (1 unique spectral windows and 1 unique polarization setups)
# SpwID  Name                                   #Chans   Frame   Ch1(MHz)  ChanWid(kHz)  TotBW(kHz)  Corrs  
#  0     63*48.8 kHz channels @ 1.41 GHz (BARY)     63   BARY    1411.483        48.824      3075.9  RR LL

#  Sources: 3
# 	  ID   Name                SpwId RestFreq(MHz)  SysVel(km/s) 
# 	  0    0137+331            any   1420.40575     1572         
# 	  1    0318+164            any   1420.40575     1572         
# 	  2    CIG0096             any   1420.40575     1572         

### Fields of interest:
# 0 0137+331
# 1 0318+164
# 2 CIG0096

### BANDPASS and FLUX calibrator: 	0 0137+331 = 3C48
### PHASE calibrator: 				1 0318+164
### TARGET: 						2 CIG0096

# ===============================================================================

### Data inspection via plotms:
###############################

plotms(vis='cig96_vla282.ms',xaxis='time',yaxis='amp',field='4~6',coloraxis='field',acgchannel='1000',timerange='09:10:00~18:50:00')

# Field 0 (bandpass and flux calibrator) shows normal amplitudes in both of its scans: 1 and 13.
# Field 1 (phase calibrator) looks good.
# Field 2 (CIG96) has noise.

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
       
# gencal(vis='cig96_vla282.ms', caltable='cig96_vla282.ms.antposvla', caltype='antposvla', antenna='')

# ===============================================================================

### Plot with the spatial distribution of the antennae:
#######################################################

plotants(vis='cig96_vla282.ms',figfile='cig96_vla282.ms.plotants.png')

# ===============================================================================

### Flagging:
#############

# Log file does not provide with any insight of the antennae status, corruption, etc.

# There is no dummy scan, no scan to remove.

### Shadowing correction over the whole .ms:

flagdata(vis='cig96_vla282.ms',
	mode='shadow', 
	flagbackup=True)
	
# INFO Shadow	=> Percentage of data flagged in table selection: 0.513585%

### Zero clipping correction over the whole .ms:

flagdata(vis='cig96_vla282.ms', 
	mode='clip', 
	clipzeros=True,
	flagbackup=True)

# No flagging done, no zeros

### Field 0 (bandpass and flux calibrator) looks fine in both polarizations; it only needs quacking of the first 10 seconds:

flagdata(vis='cig96_vla282.ms', mode='quack', field='0', quackinterval=10.0, quackmode='beg')


### Field 1 (phase calib.): same as with field = 0, it looks good, its scans only need quacking:

flagdata(vis='cig96_vla282.ms', mode='quack', field='1', quackinterval=10.0, quackmode='beg')


### Field 2 (CIG96) flagging: polarization LL and RR show different problems:

# in RR we only have two problems: for all the channels in scan = 11, antenna VA08(7) seems to be giving problems in all baselines so we remove it:

flagdata(vis='cig96_vla282.ms', mode='manual', field='2', scan='11', correlation='RR', antenna='VA08')

# also, baseline VA03&VA22 (2&17) does show high intensity values sowe remove it:

flagdata(vis='cig96_vla282.ms', mode='manual', field='2', scan='11', correlation='RR', antenna='VA03&VA22')

# in LL correlation, in scan 3, antenna VA24(19) is not working properly 

flagdata(vis='cig96_vla282.ms', mode='manual', field='2', scan='3', correlation='LL', antenna='VA24')

# as well as in scan 14, in which antenna VA24(19) shows a lot of intensity for all the baselines:

flagdata(vis='cig96_vla282.ms', mode='manual', field='2', scan='14', correlation='LL', antenna='VA24')

# we remove the first seconds of scan 16, no more quacking needed apart from this scan:

flagdata(vis='cig96_vla282.ms', mode='quack', field='2', scan='16', quackinterval=5.0, quackmode='beg')

# enough flagging for now.

# ===============================================================================

### Calibration:
### ============

### Reference antenna selection:
################################

# After checking plotants, we select VA15(12) because it is quite central (located on the E4 position):

refant='VA15'

# ===============================================================================

### Opacity and gaincurve corrections:
######################################

# Should we do opacities correction?
# myTau = plotweather(vis='cig96_vla282.ms', doPlot=T)
# display cig96_vla282.ms.plotweather.png
# gencal(vis='cig96_vla282.ms', caltype='opac', caltable='cig96_vla282.ms.opcac', parameter=myTau)
# ANSWER: no, not needed for HI observations, the atmosphere barely has any influence with HI. Only ionosphere in very flat elevations and at dawn/sunset.

# Should we do gaincurve correction?
# Takes into account the elevations of each antenn, see elevations plot: CIG96 and the phase calibrator have a large elevation range.
# gencal(vis='cig96_vla282.ms', caltype='gceff', caltable='cig96_vla282.ms.gaincurve')
# ANSWER: no, it might even harm the observations due to the possible crosstalk at low elevations and the model it is based on.

# ===============================================================================

### Delays correction:
######################

# Prior to the bandpass correction, we need to correct the phase variations with time (they are large). We use the bandpass calibrator in field 4.

gaincal(vis='cig96_vla282.ms', caltable='cig96_vla282.ms.delays', field='0', refant='VA15', gaintype='G')

# we have added no gaintables to the gaincal command because we have NOT calculated any yet: "antpos" is not defined yet, "gaincurve" is not needed, and "opacities" is not needed either.

# The answer of the previous command is:

# 2014-05-05 16:02:59 INFO gaincal	##########################################
# 2014-05-05 16:02:59 INFO gaincal	##### Begin Task: gaincal            #####
# 2014-05-05 16:02:59 INFO gaincal	gaincal(vis="cig96_vla282.ms",caltable="cig96_vla282.ms.delays",field="4",spw="",
# 2014-05-05 16:02:59 INFO gaincal	        intent="",selectdata=True,timerange="",uvrange="",antenna="",
# 2014-05-05 16:02:59 INFO gaincal	        scan="",observation="",msselect="",solint="inf",combine="",
# 2014-05-05 16:02:59 INFO gaincal	        preavg=-1.0,refant="VA25",minblperant=4,minsnr=3.0,solnorm=False,
# 2014-05-05 16:02:59 INFO gaincal	        gaintype="G",smodel=[],calmode="ap",append=False,splinetime=3600.0,
# 2014-05-05 16:02:59 INFO gaincal	        npointaver=3,phasewrap=180.0,gaintable=[''],gainfield=[''],interp=[''],
# 2014-05-05 16:02:59 INFO gaincal	        spwmap=[],gaincurve=False,opacity=[],parang=False)
# 2014-05-05 16:02:59 INFO gaincal	Opening MS: cig96_vla282.ms for calibration.
# 2014-05-05 16:03:00 INFO gaincal	Initializing nominal selection to the whole MS.
# 2014-05-05 16:03:00 INFO gaincal	Beginning selectvis--(MSSelection version)-------
# 2014-05-05 16:03:00 INFO gaincal	Reseting solve/apply state
# 2014-05-05 16:03:00 INFO gaincal	Performing selection on MeasurementSet
# 2014-05-05 16:03:00 INFO gaincal	 Selecting on field: '4'
# 2014-05-05 16:03:00 INFO gaincal	By selection 393900 rows are reduced to 40300
# 2014-05-05 16:03:01 INFO gaincal	Frequency selection: Selecting all channels in all spws.
# 2014-05-05 16:03:01 INFO gaincal	Beginning setsolve--(MSSelection version)-------
# 2014-05-05 16:03:01 INFO gaincal	Arranging to SOLVE:
# 2014-05-05 16:03:01 INFO gaincal	.   G Jones: table=cig96_vla282.ms.delays append=false solint=inf refant='VA15' minsnr=3 apmode=AP solnorm=false
# 2014-05-05 16:03:01 INFO gaincal	Beginning solve-----------------------------
# 2014-05-05 16:03:01 INFO gaincal	The following calibration terms are arranged for apply:
# 2014-05-05 16:03:01 INFO gaincal	.   (None)
# 2014-05-05 16:03:01 INFO gaincal	The following calibration term is arranged for solve:
# 2014-05-05 16:03:01 INFO gaincal	.   G Jones: table=cig96_vla282.ms.delays append=false solint=inf refant='VA15' minsnr=3 apmode=AP solnorm=false
# 2014-05-05 16:03:01 INFO gaincal	Solving for G Jones
# 2014-05-05 16:03:01 INFO gaincal	For solint = inf, found 3 solution intervals.
# 2014-05-05 16:03:04 INFO gaincal	  Found good G Jones solutions in 3 slots.
# 2014-05-05 16:03:04 INFO gaincal	Applying refant: VA25
# 2014-05-05 16:03:04 INFO gaincal	Writing solutions to table: cig96_vla282.ms.delays
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

# Our flux calibrator is 3C48, in field = 0.
# First, we check the list of available models:

setjy(vis='cig96_vla282.ms', listmodels=T)

# Candidate modimages (*) at Pablo's computer in path:
#
# /Applications/CASA.app/Contents/data/nrao/VLA/CalModels/
#
# The model chosen has to be in accordance to the band selected, in this case: L band:

setjy(vis='cig96_vla282.ms', field='0', modimage='/Applications/CASA.app/Contents/data/nrao/VLA/CalModels/3C48_L.im')

# ===============================================================================

### Bandpass calibration:
#########################

# Real bandpass calibration. For the BANDPASS and FLUX(field=4) and PHASE(field=5) calibrators we use solution interval time of solint='5s':

gaincal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.bpphase5s',
	field='0',
	refant='VA15',
	calmode='p',
	solint='5s',
	minsnr=5.0)
	
# The solution interval of 5 seconds has been calculated following with the VLA exp.time calculator using the parameters:
# Freq. = 1.42 GHz
# Medium elevation (summer time)
# Bandwith freq. = 3,076 KHz (see listobs)
# RMS noise = 20 mJy (more than 3.0 in SNR)
# 
# The calculator estimates less than one second (0.2s) is enough to get such SNR or even higher so we set solint=5s since 5s is the shortest integration time of our data. This should mean that there should be a solution for all the intervals.

# We see the phase VS time figures to see that the phase is now waaayyy flatter:

plotcal(caltable='cig96_vla282.ms.bpphase5s',
	xaxis='time',
	yaxis='phase')  

# Apply phase solutions on the fly:

bandpass(vis='cig96_vla282.ms', caltable='cig96_vla282.ms.bandpass5s', field='0', refant='VA15', solint='inf', solnorm=T, minsnr=10.0, minblperant=3, gaintable=['cig96_vla282.ms.bpphase5s'], interp=['nearest'])

# We check again the solutions:

plotcal(caltable='cig96_vla282.ms.bpphase5s', field='0', xaxis='time', yaxis='phase')  

# Difference in phase is less, phase shows a much flatter behaviour now.

# ===============================================================================

### Phase and amplitude calibration:
####################################

# Phase calibration for the calibrators, fields 4 and 5 (spw 2, where the line is):

gaincal(vis='cig96_vla282.ms',caltable='cig96_vla282.ms.intphase', field='0,1', refant='VA15', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, to use FOR THE TARGET later on:

time.sleep(5)

gaincal(vis='cig96_vla282.ms', caltable='cig96_vla282.ms.scanphase', field='0,1', refant='VA15', calmode='p', solint='inf', minsnr=5.0, gaintable=['cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s'])

# Derive amplitude solutions:

time.sleep(5)

gaincal(vis='cig96_vla282.ms', caltable='cig96_vla282.ms.amp', field='0,1',  refant='VA15', calmode='ap', solint='inf', minsnr=5.0, gaintable=['cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.intphase'])

# I check the tables with plotcal and some things do not look as expected: the first data show very low intensity compares to the rest and some graphs show only one point:

plotcal(caltable='cig96_vla282.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)
plotcal(caltable='cig96_vla282.ms.amp', xaxis='time', yaxis='amp', iteration='baseline', subplot=331)


# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal in terms of future application of the calibration to the data, (i.e., it's not an incremental table) UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_vla282.ms', caltable='cig96_vla282.ms.amp', fluxtable='cig96_vla282.ms.flux', reference='0', incremental=False)

# Result:
# Flux density for 0318+164 in SpW=0 is: 7.82742 +/- 0.00806651 (SNR = 970.36, N = 48)

# ===============================================================================

### Application of the calibration to the .ms file:
###################################################

# Note: In all applycal steps we set calwt=F. It is very important to turn off this parameter which determines if the weights are calibrated along with the data. Data from antennae with better receiver performance and/or longer integration times should have higher weights, and it can be advantageous to factor this information into the calibration. During the VLA era, meaningful weights were available for each visibility. However, at the time of this observation, the VLA was not yet recording the information necessary to calculate meaningful weights. Since these data weights are used at the imaging stage you can get strange results from having calwt=T when the input weights are themselves not meaningful, especially for self-calibration on resolved sources (your flux calibrator and target, for example).

applycal(vis='cig96_vla282.ms', field='0', gaintable=['cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.intphase', 'cig96_vla282.ms.flux'], calwt=F)

time.sleep(5)

applycal(vis='cig96_vla282.ms', field='0', gaintable=['cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.intphase', 'cig96_vla282.ms.flux'], calwt=F)

time.sleep(5)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_vla282.ms', field='2', gaintable=['cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.scanphase', 'cig96_vla282.ms.flux'], gainfield=['','','1','1'], calwt=F)


# ===============================================================================

### Imaging of CIG96:
#####################

# Splitting of field = 0 (CIG96) calibrated data, in the spw=0:

split(vis='cig96_vla282.ms', datacolumn='corrected', outputvis='cig96_vla282.corr.ms', field='2')

# Now, for the new corr.ms file:
# Field = 2 is stored as: field = 0  
# spw = 0 is not modified since it was the only spw.

# Splitting CIG96 fields in spw=0 with channel binning of 8 and time binning of 20s:

split(vis='cig96_vla282.ms', datacolumn='corrected', outputvis='cig96_vla282.corr.8chann.20s.ms', timebin='20s', width=8, field='2')

# Splitting CIG96 fields in spw=0 with time binning of 20s:

split(vis='cig96_vla282.ms', datacolumn='corrected', outputvis='cig96_vla282.corr.20s.ms', timebin='20s', field='2')



### UV continuum subtraction:
#############################

# With plotms we can see
# Remove continuum from source by subtracting the channels from 0to900 and from 1150to1700 (the fact they are dividided by 10 is because they have been binned in channel in the split command above):

uvcontsub(vis='cig96_vla282.corr.ms', field='0', fitspw='0:3~18;46~57')

# and now we make the cleaning over the continuum subtracted:

clean(vis='cig96_vla282.corr.ms.contsub',
	imagename='cig96_vla282.corr.contsub.line',
	field='0',
	spw='0:18~46',
	mode='channel',
	niter=1000,
	interactive=T,
	imsize=[150,150],
	phasecenter=0,
	cell='15.0arcsec')
	
# Beam size: bmaj: 40.3626", bmin: 21.8761", bpa: 47.2371 deg

### Noise level via viewer task:
################################

---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
(cig96_vla282.corr.contsub.line.image)
     Frequency       Velocity         Stokes BrightnessUnit       BeamArea 
 1.41305e+09Hz    1553.49km/s              I        Jy/beam        4.44663 
          Npts            Sum    FluxDensity           Mean            Rms 
          1155   2.966822e-02   6.672063e-03   2.568677e-05   4.262555e-04 
       Std dev        Minimum        Maximum   region count 
  4.256651e-04  -1.101821e-03   1.550763e-03              1 
---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----





################ END
#
#
#
#
#
################  EVLA SCRIPT USED AS REFERENCE  ################

# First of all, in case we don't have the CASA measurement set file (.ms) yet but Science Data Model observation, we have to convert it:

# importevla(asdm='cig96_16', vis='cig96_vla282.ms')

### WORK WITH .MS FILE ###

# We import 'time' to introduce time delays between commands (to allow complete writing):

import time

# Example, 5s delay:   time.sleep(5)

# Summary of the observation and identification of sources and calibrators:

listobs(vis='cig96_vla282.ms')

# From listobs we extract the following information:
#
# We have 9 spectral windows (spw): from 0 to 8, both included.
#
# ID/Field - Object - Type
# ------------------------
# 1	-	0137+331=3C48	-	flux calibrator
# 2	-	J0238+1636	-	phase calibrator
# 3	-	CIG96		-	target

# Creation of antenna position table for corrections (automatic task): caltype='antpos'.
# From obs.log we can see which antennae have wrong baseline positions but if we set: antenna=’’,
# this triggers an automated lookup of antenna positions for EVLA.
# If no corrections are found, it will throw an exception with the message, 'no offsets found. no caltable
# created'. The task may terminate with a SEVERE error message and may sound alarming, but it simply means
# that it cannot produce the caltable.

gencal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.antpos',
	caltype='antpos',
	antenna='')

# Some positions of certain antennas corrected!

# We see the band A0/C0 dedicated to the line (spw = 0) before doing anything else:

plotms

# E.g.:
# plotms(vis='cig96_vla282.ms', xaxis='time', yaxis='amp', spw='0', field='', avgchannel='200', coloraxis='field')
# plotms(vis='cig96_vla282.ms', xaxis='freq', yaxis='amp', spw='0', field='3', avgtime='1e9', coloraxis='field')

# There is dummy scan, so we need to remove it:

flagdata(vis='cig96_vla282.ms', 
	mode='manual',
	scan='1')

# Flagging of the shadowed antennas:

flagdata(vis='cig96_vla282.ms',
	mode='shadow', 
	flagbackup=True)

# We write out the zero values with the clip mode:

flagdata(vis='cig96_vla282.ms', 
	mode='clip', 
	clipzeros=True,
	flagbackup=True)

# Flagging: quack mode, manual mode.

# We remove corrupted antenna ea19 according to the log so I remove it completely:

flagdata(vis='cig96_vla282.ms', 
	mode='manual', 
	antenna='ea19')

# We remove the antenna ea04 and baseline ea06-ea16 in field = 1:

flagdata(vis='cig96_vla282.ms', mode='manual', spw='0', field='1', antenna='ea04')
flagdata(vis='cig96_vla282.ms', mode='manual', spw='0', field='1', antenna='ea06&ea16')

# quack = remove/keep specific time range at scan beginning/end
# Quacking needed:
# I remove the initial seconds of some scans manually, flagging the time range in which I see the low amplitude values:

flagdata(vis='cig96_vla282.ms', mode='manual', spw='0', scan='3', timerange='08:53:25~08:53:30')
flagdata(vis='cig96_vla282.ms', mode='manual', spw='0', scan='4', timerange='08:55:12~08:55:13')
flagdata(vis='cig96_vla282.ms', mode='manual', spw='0', scan='8', timerange='09:23:16~09:23:19')

# We remove the RFI around channels 350-360 visible in spw=0 and fields=2,3:

flagdata(vis='cig96_vla282.ms', mode='manual', field='2~3', spw = '0:345~360')

# E.g.:

# flagdata(vis='cig96_vla282.ms', mode='quack', quackmode='beg', scan = '3' , spw = '0' , quackinterval=51.0)
# flagdata(vis='cig96_vla282.ms', mode='quack', quackmode='beg', scan = '7,10' , spw = '0' , quackinterval=16.0)
# flagdata(vis='cig96_vla282.ms', mode='manual', spw = '0' , antenna = '23', correlation = 'LL')
# flagdata(vis='cig96_vla282.ms', mode='manual', spw = '0:383~387')

# Create and check the weather report. Plots are written in a png file and myTau variable saves the opacities.

myTau = plotweather(vis='cig96_vla282.ms', doPlot=T)

# Plot can be seen typing on the console the following command:

### !display cig96_vla282.ms.plotweather.png

# Create calibration table for the opacities:

gencal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.opac',
	caltype='opac',
	spw='0~8',
	parameter=myTau)

# Create calibration table for the gaincurve (behavior of each antenna as a function of elevation):

gencal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.gaincurve',
	caltype='gceff')

# We apply these 2 calibrations (a-priori calibration) to avoid carrying on with data that in the end will result unnecessary

applycal(vis='cig96_vla282.ms',
	field='1',
	gaintable=['cig96_vla282.ms.antpos','cig96_vla282.ms.gaincurve', 'cig96_vla282.ms.opac'],
	calwt=F)

# Antennae position plot:

plotants(vis='cig96_vla282.ms',
	figfile='cig96_vla282.ms.ant.png')

# I select "ea13" since it is centered and has not presented any problems so far.

# After flagdata command flagging, you have to force a complete reload of the cache to look at the same plot again with the new flags applied. To do this, either check the "force reload" box in the lower left, or do Shift+Plot.

# If the plotms tool hangs during a plot, try clicking the cancel button on the  load progress GUI, and/or if you see a "table locked" message try typing "clearstat" on the CASA command line.

# For now, enough flagging.

# Setting the flux density. Our flux calibrator is 3C48, ID=1. First we see the list of available models:

setjy(vis='cig96_vla282.ms', listmodels=T)

# Candidate modimages (*) at dae66/pepino in: /usr/local/bin/casapy-41.0.prerelease-3-64b/data/nrao/VLA/CalModels
# The model chosen is in accordance to the band selected, in this case: L band.

setjy(vis='cig96_vla282.ms',
	field='1',
	modimage='/usr/local/bin/casapy-41.0.prerelease-3-64b/data/nrao/VLA/CalModels/3C48_L.im')

# Checking the bandpass (phase/amp VS freq). First, we have to select the reference antenna checking the plotants figure.

# Before solving for the bandpass, we have to correct the phase variations with time, because they are 
# not small. To do so, first we have to correct the delay (the slope of phase across frequency) setting the 
# reference antenna:

gaincal(vis='cig96_vla282.ms', caltable='cig96_vla282.ms.delays', field='1', refant='VA15', gaintype='G', gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.gaincurve', 'cig96_vla282.ms.opac'])


#      gaintype -- Type of gain solution (G, T, or GSPLINE)  
#                   default: ’G’;
#	            example: gaintype=’GSPLINE’  
#              ’G’ means determine gains for each polarization and sp_wid  
#              ’T’ obtains one solution for both polarizations;  Hence. their phase offset must be first removed using a prior G.
#              ’GSPLINE’ makes a spline fit to the calibrator data.  It is  
#                   useful for noisy data and fits a smooth curve through the  
#                   calibrated amplitude and phase.  However,  
#                   at present GSPLINE is somewhat experimental.  Use with  
#                   caution and check solutions.  
#              ’K’ solves for simple antenna-based single-band delays  
#                   via FFTs of the spectra on baselines to the  
#                   reference antenna.  (This is not global fringe-fitting.)  
#              ’KCROSS’ solves for a global cross-hand  
#                   delay.  Use parang=T and apply prior gain and  
#                   bandpass solutions.  

# Real bandpass calibration. For the FLUX(ID=1) and PHASE(ID=2) calibrators we use solution interval time of: solint='5s':

gaincal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.bpphase5s',
	field='1',
	spw='0~8',
	refant='VA15',
	calmode='p',
	solint='5s',
	minsnr=3.0)

# We see the phase VS time figures to see that the phase is now waaayyy flatter:

plotcal(caltable='cig96_vla282.ms.bpphase5s',
	xaxis='time',
	yaxis='phase')  

# Check solutions:
# plotcal(caltable='cig96_vla282.ms.bpphase5s', xaxis='time', yaxis='phase', iteration='antenna', subplot=331, plotrange=[0,0,-180, 180])

# Apply phase solutions on the fly:

bandpass(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.bandpass5s',
	field='1',
	spw='0~8',
	refant='VA15',
	solint='inf',
	solnorm=T,
	minsnr=10.0,
	minblperant=3,
	gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.delays', 'cig96_vla282.ms.bpphase5s'],
	interp=['nearest'])

# Check solutions:
# plotcal(caltable='cig96_vla282.ms.bandpass5s', xaxis='chan', yaxis='amp', iteration='antenna', subplot=331)
# plotcal(caltable='cig96_vla282.ms.bandpass5s', xaxis='chan', yaxis='phase', iteration='antenna', subplot=331, plotrange=[0,0,-10,10])

# plotcal(caltable='cig96_vla282.ms.bandpass5s',xaxis='freq',yaxis='amp',iteration='antenna',subplot=331)

# plotcal(caltable='cig96_vla282.ms.bandpass5s',xaxis='freq',yaxis='phase',iteration='antenna',subplot=331)


# Apply calibration to the flux cal:

# Not necessary because it is applied below:
# applycal(vis='cig96_vla282.ms', field='1', gaintable=['cig96_vla282.ms.gaincurve', 'cig96_vla282.ms.opac','cig96_vla282.ms.delays','cig96_vla282.ms.bandpass5s'], gainfield=['','','1','1'], calwt=F)

# As the log says:
# 2013-06-21 16:26:02 INFO applycal	Adding CORRECTED_DATA column(s).
# 2013-06-21 16:26:03 INFO applycal	Initializing CORRECTED_DATA (to DATA).

# Plotms show now a good behavoir in the bandpass.

# Next step is gain calibration (phase/amp vs time).
# Start solving for the phases separately, since phases change on a much shorter timescale than amplitude.


######     PHASE AND AMPLITUDE CALIBRATION     ######


# First for the calibrators, using solint='int' or ='5s':

gaincal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.intphase',
	field='1,2',
	spw='0~8',
	refant='VA15',
	calmode='p',
	solint='5s',
	minsnr=3.0,
	gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, to use for the targets later on:

time.sleep(5)

gaincal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.scanphase',
	field='1,2',
	spw='0~8',
	refant='VA15',
	calmode='p',
	solint='inf',
	minsnr=3.0,
	gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s'])

# Derive amplitude solutions:

time.sleep(5)

gaincal(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.amp',
	field='1,2',
	spw='0~8',
	refant='VA15',
	calmode='ap',
	solint='inf',
	minsnr=3.0,
	gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.intphase'])

# I check the tables with plotcal and everything looks really good!:
# plotcal(caltable='cig96_vla282.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)

# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal
# in terms of future application of the calibration to the data, (i.e., it's not an incremental table)
# UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_vla282.ms',
	caltable='cig96_vla282.ms.amp',
	fluxtable='cig96_vla282.ms.flux',
	reference='1',
	incremental=False)

# We obtain:

# 2013-08-08 09:51:45 INFO fluxscale	 Fitted spectrum for J0238+1636 with fitorder=1: Flux density = 0.745584 +/- 0.0502467 (freq=1.3827 GHz) spidx=0.162341 +/- 0.445142


#####   APPLYING CALIBRATION TO THE .ms   #####

# I apply the solutions to the measurement set.

# Note: In all applycal steps we set calwt=F. It is very important to turn off this parameter which determines if the weights are calibrated along with the data. Data from antennas with better receiver performance and/or longer integration times should have higher weights, and it can be advantageous to factor this information into the calibration. During the VLA era, meaningful weights were available for each visibility. However, at the time of this observation, the VLA was not yet recording the information necessary to calculate meaningful weights. Since these data weights are used at the imaging stage you can get strange results from having calwt=T when the input weights are themselves not meaningful, especially for self-calibration on resolved sources (your flux calibrator and target, for example).

applycal(vis='cig96_vla282.ms',
	field='1',
	gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.intphase', 'cig96_vla282.ms.flux'],
	calwt=F)

time.sleep(5)

applycal(vis='cig96_vla282.ms',
	field='2',
	gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.intphase', 'cig96_vla282.ms.flux'],
	calwt=F)

time.sleep(5)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_vla282.ms',
	field='3',
	gaintable=['cig96_vla282.ms.antpos', 'cig96_vla282.ms.delays', 'cig96_vla282.ms.bandpass5s', 'cig96_vla282.ms.scanphase', 'cig96_vla282.ms.flux'],
	gainfield=['','','2','2'],
	calwt=F)
 
# Some bad data may reveal now. If so, we use flagdata again like this:
# flagdata(vis='cig96_vla282.ms', mode='manual', spw='xx~xx', antenna='eaXX, ...', correlation='')

# We make some plots 

# Repeat the calibration process until everything looks fine.


#######################
# Imaging of cig96_16 #  
#######################

# Splitting just CIG96 data in the spw=0:

split(vis='cig96_vla282.ms',
	datacolumn='corrected',
	outputvis='cig96_16.corr.ms',
	field='3',
	spw='0')

time.sleep(5)

# Splitting CIG96 fields in spw=0 with channel and time binning:

split(vis='cig96_vla282.ms',
	datacolumn='corrected',
	outputvis='cig96_16.corr.8chann.20s.ms',
	timebin='20s',
	width=8,
	field='3',
	spw='0')

# We are only getting one field (in this case it is number 3, CIG96), so in the new .ms it will be renamed as field 0. Since it is only one, we do not need any specification command like: field = '0' because it is unique.

# Now we should concatenate the different .ms from cig96 (after calibrating all the datasets:
# concat(vis=['cig96_7.corr.ms','cig96_160marA.ms','cig96_160marB.ms'], concatvis='cig96_comb.ms')
# but first we have to calibrate the 2 observations (called A and B) from 10th of march.

# Imaging the target:

# Continuum:

clean(vis='cig96_16.corr.10chann.20s.ms',
	imagename='cig96_16.corr.cont',
	field='0',
	spw='0:40~90;115~170',
	mode='mfs',
	interpolation='linear',
	niter=1000,
	interactive=True,
	imsize=[250,250],
	cell='10.0arcsec')

# Line:

clean(vis='cig96_16.corr.10chann.20s.ms',
	imagename='cig96_16.corr.line',
	field='0',
	spw='0:90~115',
	mode='channel',
	niter=1000,
	interactive=T,
	imsize=[250,250],
	cell='15.0arcsec')

# Remove continuum from source by subtracting the channels from 0to900 and from 1150to1700 (the fact they are dividided by 10 is because they have been binned in channel in the split command above):

uvcontsub(vis='cig96_16.corr.10chann.20s.ms',
	field='0',
	fitspw='0:40~90;115~170')

# and now we make the cleaning over the continuum subtracted:

clean(vis='cig96_16.corr.10chann.20s.ms.contsub',
	imagename='cig96_16.corr.contsub.line',
	field='0',
	spw='0:92~115',
	mode='channel',
	niter=1000,
	interactive=T,
	imsize=[150,150],
	phasecenter=0,
	cell='15.0arcsec')




