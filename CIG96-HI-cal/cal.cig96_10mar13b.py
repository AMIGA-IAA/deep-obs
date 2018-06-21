# -*- coding: iso-8859-1 -*-
### EVLA DATA REDUCTION
### Project: 13A-341
### Dataset date: 10mar13 (first observation)
### Original dataset name: 13A-341.sb19189214.eb19331287.56361.737114351854
### Renamed as: cig96_03.ms
###
### Configuration: D (3/3)

# ===============================================================================

### Import of EVLA data from SDM format:
########################################

# importevla(asdm='cig96_03', vis='cig96_03.ms')

# ===============================================================================

### Listobs inspection:
#######################

listobs(vis='cig96_03.ms')

# Listobs output summary: 27 antennae, 9 spectral windows, 2 polarizations (RR, LL), 1 dummy scan, 
# 3 fields (2 calibrators + 1 target), 2048 channels of 7.8 kHz each for spw=0:

# Spectral Window 0 (line):
# SpwID  Name          #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) BBC Num  Corrs  
# 0	 A0C0#0        2048    TOPO    1404.995         7.812     16000.0      12  RR  LL

#  Fields: 3
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
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

# flagdata(vis='cig96_03.ms', mode='manual', field='0')

# We run plotms to inspect the general aspect of the data:

# plotms(vis='cig96_03.ms', xaxis='time', yaxis='amp', field='1', spw='0', coloraxis='field', avgchannel='1e9')

# ===============================================================================

# Antenna position correction:
##############################

gencal(vis='cig96_03.ms', caltable='cig96_03.ms.antpos', caltype='antpos')

# ===============================================================================

### Plot with the spatial distribution of the antennae:
#######################################################

plotants(vis='cig96_03.ms',figfile='plotants_cig96_03.png')

# ===============================================================================

### Flagging:
#############

# Log file provides some insight of the antennae status, corruption, etc.:

# 10Mar 17:41:22 10Mar 18:41:15 FRONT END C132657 0.50 29.9 
# Antenna(s) 25 (Data: Lost):
# LCP visibility amplitudes near zero at L-band due to a receiver problem.

# We remove corrupted antenna ea25(23) according to the log so I remove it completely:

flagdata(vis='cig96_03.ms', mode='manual', field='', spw='', antenna='ea25', flagbackup=True)

### Shadowing correction over the whole .ms:

flagdata(vis='cig96_03.ms', mode='shadow',  flagbackup=True)

# Percentage of data flagged in table selection: 25.0291%

### Zero clipping correction over the whole .ms:

flagdata(vis='cig96_03.ms', mode='clip', clipzeros=True, flagbackup=True)

# Percentage of data flagged in table selection: 0.574026%

### Field 1 (flux and bandpass calibr.):

# Shows an RFI around channel ~385:

flagdata(vis='cig96_03.ms', mode='manual', field='1', spw = '0:370~400')

# Antenna ea19 shows very low amplitudes, we remove it:

flagdata(vis='cig96_03.ms', mode='manual', field='1', spw='', antenna='ea19', flagbackup=True)

# also, baseline 14&15 and 2&14 has some weird and high amplitude points, we remove them:

flagdata(vis='cig96_03.ms', mode='manual', field='1', spw='', antenna='14&15', flagbackup=True)
flagdata(vis='cig96_03.ms', mode='manual', field='1', spw='', antenna='2&14', flagbackup=True)

# no quacking needed.

### Field 2 (phase calib.): needs some quacking in the 3 scans

flagdata(vis='cig96_03.ms', mode='quack', field='2', quackinterval=15.0, quackmode='beg', flagbackup=True)

# Antenna ea19 shows very low amplitudes in LL polarization, we remove it:

flagdata(vis='cig96_03.ms', mode='manual', field='2', antenna='ea19', flagbackup=True)

# We remove the RFI around channel ~385 visible in spw=0 and fields=2,3:

flagdata(vis='cig96_03.ms', mode='manual', field='2', spw = '0:370~400')

# also the Galactic emission in ~1950:

flagdata(vis='cig96_03.ms', mode='manual', field='2', spw = '0:1945~1965')

### Field 3 (CIG96):

# Same RFIs and Galactic emission as in field 2:

flagdata(vis='cig96_03.ms', mode='manual', field='3', spw = '0:370~400')
flagdata(vis='cig96_03.ms', mode='manual', field='3', spw = '0:1945~1965')

# and some quacking at the beginning of the scans:

flagdata(vis='cig96_03.ms', mode='quack', field='3', quackinterval=15.0, quackmode='beg', flagbackup=True)

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
# myTau = plotweather(vis='cig96_03.ms', doPlot=T)
# display cig96_03.ms.plotweather.png
# gencal(vis='cig96_03.ms', caltype='opac', caltable='cig96_03.ms.opcac', parameter=myTau)
# ANSWER: no, not needed for HI observations, the atmosphere barely has any influence with HI. Only 
# ionosphere in very flat elevations and at dawn/sunset.

# Should we do gaincurve correction?
# Takes into account the elevations of each antenna, see elevations plot: CIG96 and the phase calibrator 
# have a large elevation range.
#
# gencal(vis='cig96_03.ms', caltype='gceff', caltable='cig96_03.ms.gaincurve')
#
# ANSWER: no, it might even harm the observations due to the possible crosstalk at low elevations 
# in C and D configs. and the model it is based on.

# ===============================================================================

### Delays correction:
######################

# Prior to the bandpass correction, we need to correct the phase variations with time (they are large). 
# We use the bandpass calibrator in field 1 with the antenna position correction table:

gaincal(vis='cig96_03.ms', caltable='cig96_03.ms.delays', field='1', refant='ea26', gaintype='K', gaintable=['cig96_03.ms.antpos'])

# ===============================================================================

### Flux density:
#################

# Our flux calibrator is 3C48, in field = 0.
# First, we check the list of available models:

setjy(vis='cig96_03.ms', listmodels=T)

# Candidate modimages (*) at pepino (dae66) in path:
# /Applications/CASA.app/Contents/data/nrao/VLA/CalModels/
#
# Candidate modimages (*) at NRAO local workstations in path:
# /home/casa/packages/RHEL5/release/casapy-42.1.29047-001-1-64b/data/nrao/VLA/CalModels/
# 
# The model chosen has to be in accordance with the calibrator and band selected: 3C48 in L band:

setjy(vis='cig96_03.ms', field='1', modimage='/mnt/scops/data/data/paramimo/casapy-42.2.30986-1-64b/data/nrao/VLA/CalModels/3C48_L.im')

# ===============================================================================

### Bandpass calibration:
#########################

# For the bandpass and flux (field=1) and phase (field=2) calibrators we use solution interval 
# time of solint='5s':

gaincal(vis='cig96_03.ms', caltable='cig96_03.ms.bpphase5s', field='1', refant='ea26', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_03.ms.antpos','cig96_03.ms.delays' ])

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

# plotcal(caltable='cig96_03.ms.bpphase5s', xaxis='time', yaxis='phase')

# Apply phase solutions on the fly:

bandpass(vis='cig96_03.ms', caltable='cig96_03.ms.bandpass5s', field='1', refant='ea26', solint='inf', solnorm=T, minsnr=10.0, minblperant=3, gaintable=['cig96_03.ms.bpphase5s', 'cig96_03.ms.antpos'], interp=['nearest'])

# We check again the solutions:

# plotcal(caltable='cig96_03.ms.bpphase5s', field='1', xaxis='time', yaxis='phase')  

# Difference in phase is less, phase shows a flatter behaviour now.

# ===============================================================================


### Phase and amplitude calibration:
####################################

# Phase calibration for the calibrators, fields 1 and 2 (spw 0, where the line is):

gaincal(vis='cig96_03.ms',caltable='cig96_03.ms.intphase', field='1,2', refant='ea26', calmode='p', solint='5s', minsnr=5.0, gaintable=['cig96_03.ms.antpos', 'cig96_03.ms.delays', 'cig96_03.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, 
# to use FOR THE TARGET later on:

gaincal(vis='cig96_03.ms', caltable='cig96_03.ms.scanphase', field='1,2', refant='ea26', calmode='p', solint='inf', minsnr=5.0, gaintable=['cig96_03.ms.antpos', 'cig96_03.ms.delays', 'cig96_03.ms.bandpass5s'])

# Derive amplitude solutions:

gaincal(vis='cig96_03.ms', caltable='cig96_03.ms.amp', field='1,2',  refant='ea26', calmode='ap', solint='inf', minsnr=5.0, gaintable=['cig96_03.ms.antpos', 'cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.intphase'])

# I check the tables with plotcal and some things do not look as expected: the first data show very 
# low intensity compares to the rest and some graphs show only one point:

# plotcal(caltable='cig96_03.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)
# plotcal(caltable='cig96_03.ms.amp', xaxis='time', yaxis='amp', iteration='baseline', subplot=331)

# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal 
# in terms of future application of the calibration to the data, (i.e., it's not an incremental table) 
# UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_03.ms', caltable='cig96_03.ms.amp', fluxtable='cig96_03.ms.flux', reference='1', transfer='2', incremental=False)

# Result:
# 
# Flux density for J0238+1636 in SpW=0 (freq=1.405e+09 Hz) is: 0.828242 +/- 0.00983464 (SNR = 84.2169, N = 50)

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

applycal(vis='cig96_03.ms', field='1', gaintable=['cig96_03.ms.antpos', 'cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.intphase', 'cig96_03.ms.flux'], gainfield=['','','','1',''], calwt=F)

# time.sleep(5)

applycal(vis='cig96_03.ms', field='2', gaintable=['cig96_03.ms.antpos', 'cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.intphase', 'cig96_03.ms.flux'], gainfield=['','','','2',''], calwt=F)

# time.sleep(5)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_03.ms', field='3', gaintable=['cig96_03.ms.antpos', 'cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.scanphase', 'cig96_03.ms.flux'], gainfield=['','','','2',''], calwt=F)

# ===============================================================================












# E.g.:

# flagdata(vis='cig96_03.ms', mode='quack', quackmode='beg', scan = '3' , spw = '0' , quackinterval=51.0)

# flagdata(vis='cig96_03.ms', mode='quack', quackmode='beg', scan = '7,10' , spw = '0' , quackinterval=16.0)

# flagdata(vis='cig96_03.ms', mode='manual', spw = '0' , antenna = '23', correlation = 'LL')

# flagdata(vis='cig96_03.ms', mode='manual', spw = '0:383~387')

# Create and check the weather report. Plots are written in a png file and myTau variable saves the opacities.

myTau = plotweather(vis='cig96_03.ms', doPlot=T)

# Plot can be seen typing on the console the following command:

### !display cig96_03.ms.plotweather.png

# Create calibration table for the opacities:

gencal(vis='cig96_03.ms',
	caltable='cig96_03.ms.opac',
	caltype='opac',
	spw='0~8',
	parameter=myTau)

# Create calibration table for the gaincurve (behavior of each antenna as a function of elevation):

gencal(vis='cig96_03.ms',
	caltable='cig96_03.ms.gaincurve',
	caltype='gceff')

# We apply these 2 calibrations (a-priori calibration) to avoid carrying on with data that in the end will result unnecessary

applycal(vis='cig96_03.ms',
	field='1',
	gaintable=['cig96_03.ms.gaincurve', 'cig96_03.ms.opac'],
	calwt=F)

# Antennae position plot:

plotants(vis='cig96_03.ms',
	figfile='cig96_03.ms.ant.png')

# I select "ea13" (ID=12) since it is centered and has not presented any problems so far.

# After flagdata command flagging, you have to force a complete reload of the cache to look at the same plot again with the new flags applied. To do this, either check the "force reload" box in the lower left, or do Shift+Plot.

# If the plotms tool hangs during a plot, try clicking the cancel button on the  load progress GUI, and/or if you see a "table locked" message try typing "clearstat" on the CASA command line.

# For now, enough flagging.

# Setting the flux density. Our flux calibrator is 3C48, ID=1. First we see the list of available models:

setjy(vis='cig96_03.ms', listmodels=T)

# Candidate modimages (*) at dae66/pepino in: /usr/local/bin/casapy-41.0.prerelease-3-64b/data/nrao/VLA/CalModels
# The model chosen is in accordance to the band selected, in this case: L band.

setjy(vis='cig96_03.ms',
	field='1',
	modimage='/usr/local/bin/casapy-41.0.prerelease-3-64b/data/nrao/VLA/CalModels/3C48_L.im')

# Checking the bandpass (phase/amp VS freq). First, we have to select the reference antenna checking the plotants figure.

# Before solving for the bandpass, we have to correct the phase variations with time, because they are 
# not small. To do so, first we have to correct the delay (the slope of phase across frequency) setting the 
# reference antenna:

gaincal(vis='cig96_03.ms', caltable='cig96_03.ms.delays', field='1', refant='ea13', gaintype='G', gaintable=['cig96_03.ms.gaincurve', 'cig96_03.ms.opac'])


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

# NOTE: normally, in this and all the following 'gaincal' tasks, the tables list of 'gaintable' should include a 
# 'cig96_03.ms.antpos' table, but we don't have it in this case for the reasons stated above: there are no
# antennas position correction to be made, hence, no such table exists in March 9th 2013 case.

# Real bandpass calibration. For the FLUX(ID=1) and PHASE(ID=2) calibrators we use solution interval time of: solint='5s':

gaincal(vis='cig96_03.ms',
	caltable='cig96_03.ms.bpphase5s',
	field='1',
	spw='0~8',
	refant='ea13',
	calmode='p',
	solint='5s',
	minsnr=3.0)

# We see the phase VS time figures to see that the phase is now waaayyy flatter:

plotcal(caltable='cig96_03.ms.bpphase5s',
	xaxis='time',
	yaxis='phase')  

# Check solutions:
# plotcal(caltable='cig96_03.ms.bpphase5s', xaxis='time', yaxis='phase', iteration='antenna', subplot=331, plotrange=[0,0,-180, 180])

# Apply phase solutions on the fly:

bandpass(vis='cig96_03.ms',
	caltable='cig96_03.ms.bandpass5s',
	field='1',
	spw='0~8',
	refant='ea13',
	solint='inf',
	solnorm=T,
	minsnr=10.0,
	minblperant=3,
	gaintable=['cig96_03.ms.delays', 'cig96_03.ms.bpphase5s'],
	interp=['nearest'])

# Check solutions:
# plotcal(caltable='cig96_03.ms.bandpass5s', xaxis='chan', yaxis='amp', iteration='antenna', subplot=331)
# plotcal(caltable='cig96_03.ms.bandpass5s', xaxis='chan', yaxis='phase', iteration='antenna', subplot=331, plotrange=[0,0,-10,10])

# plotcal(caltable='cig96_03.ms.bandpass5s',xaxis='freq',yaxis='amp',iteration='antenna',subplot=331)

# plotcal(caltable='cig96_03.ms.bandpass5s',xaxis='freq',yaxis='phase',iteration='antenna',subplot=331)


# Apply calibration to the flux cal:

# Not necessary because it is applied below:
# applycal(vis='cig96_03.ms', field='1', gaintable=['cig96_03.ms.gaincurve', 'cig96_03.ms.opac','cig96_03.ms.delays','cig96_03.ms.bandpass5s'], gainfield=['','','1','1'], calwt=F)

# As the log says:
# 2013-06-21 16:26:02 INFO applycal	Adding CORRECTED_DATA column(s).
# 2013-06-21 16:26:03 INFO applycal	Initializing CORRECTED_DATA (to DATA).

# Plotms show now a good behavoir in the bandpass.

# Next step is gain calibration (phase/amp vs time).
# Start solving for the phases separately, since phases change on a much shorter timescale than amplitude.


######     PHASE AND AMPLITUDE CALIBRATION     ######


# First for the calibrators, using solint='int' or ='5s':

gaincal(vis='cig96_03.ms',
	caltable='cig96_03.ms.intphase',
	field='1,2',
	spw='0~8',
	refant='ea13',
	calmode='p',
	solint='5s',
	minsnr=3.0,
	gaintable=['cig96_03.ms.delays', 'cig96_03.ms.bandpass5s'])

# We create another calibration table, using solint='inf', i.e., finding one solution over the whole scan, to use for the targets later on:

gaincal(vis='cig96_03.ms',
	caltable='cig96_03.ms.scanphase',
	field='1,2',
	spw='0~8',
	refant='ea13',
	calmode='p',
	solint='inf',
	minsnr=3.0,
	gaintable=['cig96_03.ms.delays', 'cig96_03.ms.bandpass5s'])

# Derive amplitude solutions:

gaincal(vis='cig96_03.ms',
	caltable='cig96_03.ms.amp',
	field='1,2',
	spw='0~8',
	refant='ea13',
	calmode='ap',
	solint='inf',
	minsnr=3.0,
	gaintable=['cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.intphase'])

# I check the tables with plotcal and everything looks really good!:
# plotcal(caltable='cig96_03.ms.amp', xaxis='time', yaxis='amp', iteration='antenna', subplot=331)

# Now I derive the flux for the rest of the sources. Note that the flux table REPLACES the amp.gcal
# in terms of future application of the calibration to the data, (i.e., it's not an incremental table)
# UNLESS we set incremental=T (as of CASA 4.0):

myflux = fluxscale(vis='cig96_03.ms',
	caltable='cig96_03.ms.amp',
	fluxtable='cig96_03.ms.flux',
	reference='1',
	incremental=False)

# We obtain:

# 2013-08-08 09:51:45 INFO fluxscale	 Fitted spectrum for J0238+1636 with fitorder=1: Flux density = 0.745584 +/- 0.0502467 (freq=1.3827 GHz) spidx=0.162341 +/- 0.445142


#####   APPLYING CALIBRATION TO THE .ms   #####

# I apply the solutions to the measurement set.

# Note: In all applycal steps we set calwt=F. It is very important to turn off this parameter which determines if the weights are calibrated along with the data. Data from antennas with better receiver performance and/or longer integration times should have higher weights, and it can be advantageous to factor this information into the calibration. During the VLA era, meaningful weights were available for each visibility. However, at the time of this observation, the VLA was not yet recording the information necessary to calculate meaningful weights. Since these data weights are used at the imaging stage you can get strange results from having calwt=T when the input weights are themselves not meaningful, especially for self-calibration on resolved sources (your flux calibrator and target, for example).

applycal(vis='cig96_03.ms',
	field='1',
	gaintable=['cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.intphase', 'cig96_03.ms.flux'],
	calwt=F)

applycal(vis='cig96_03.ms',
	field='2',
	gaintable=['cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.intphase', 'cig96_03.ms.flux'],
	calwt=F)

# For the target sources we use the scanphase.gcal table:

applycal(vis='cig96_03.ms',
	field='3',
	gaintable=['cig96_03.ms.delays', 'cig96_03.ms.bandpass5s', 'cig96_03.ms.scanphase', 'cig96_03.ms.flux'],
	gainfield=['','','2','2'],
	calwt=F)

# This last calibration shows this message:

# The following MS spws have no corresponding cal spws: 3 5 
# 2014-01-13 17:45:17	WARN	applycal::Calibrater::correct	Spectral window(s) 3, 5, 
# 2014-01-13 17:45:17	WARN	applycal::Calibrater::correct+	  could not be corrected due to missing (pre-)calibration
# 2014-01-13 17:45:17	WARN	applycal::Calibrater::correct+	    in one or more of the specified tables.
# 2014-01-13 17:45:17	WARN	applycal::Calibrater::correct+	    Please check your results carefully!

# plotms shows NO problem.
 
# Some bad data may reveal now. If so, we use flagdata again like this:
# flagdata(vis='cig96_03.ms', mode='manual', spw='xx~xx', antenna='eaXX, ...', correlation='')

# We make some plots 

# Repeat the calibration process until everything looks fine.


#######################
# Imaging of cig96_03 #  
#######################

# Splitting just CIG96 data in the spw=0:

split(vis='cig96_03.ms',
	datacolumn='corrected',
	outputvis='cig96_03.corr.ms',
	field='3',
	spw='0')

# Splitting CIG96 fields in spw=0 with channel and time binning:

split(vis='cig96_03.ms',
	datacolumn='corrected',
	outputvis='cig96_03.corr.10chann.20s.ms',
	timebin='20s',
	width=10,
	field='3',
	spw='0')

# We are only getting one field (in this case it is number 3, CIG96), so in the new .ms it will be renamed as field 0. Since it is only one, we do not need any specification command like: field = '0' because it is unique.

# Now we should concatenate the different .ms from cig96 (after calibrating all the datasets:
# concat(vis=['cig96_7.corr.ms','cig96_030marA.ms','cig96_030marB.ms'], concatvis='cig96_comb.ms')
# but first we have to calibrate the 2 observations (called A and B) from 10th of march.

# Imaging the target:

# Continuum:

clean(vis='cig96_03.corr.10chann.20s.ms',
	imagename='cig96_03.corr.cont',
	field='0',
	spw='0:40~90;115~170',
	mode='mfs',
	interpolation='linear',
	niter=1000,
	interactive=True,
	imsize=[250,250],
	cell='10.0arcsec')

# Line:

clean(vis='cig96_03.corr.10chann.20s.ms',
	imagename='cig96_03.corr.line',
	field='0',
	spw='0:90~115',
	mode='channel',
	niter=1000,
	interactive=T,
	imsize=[250,250],
	cell='15.0arcsec')

# Remove continuum from source by subtracting the channels from 0to900 and from 1150to1700 (the fact they are dividided by 10 is because they have been binned in channel in the split command above):

uvcontsub(vis='cig96_03.corr.10chann.20s.ms',
	field='0',
	fitspw='0:40~90;115~170')

# and now we make the cleaning over the continuum subtracted:

clean(vis='cig96_03.corr.10chann.20s.ms.contsub',
	imagename='cig96_03.corr.contsub.line',
	field='0',
	spw='0:92~115',
	mode='channel',
	niter=1000,
	interactive=T,
	imsize=[150,150],
	phasecenter=0,
	cell='15.0arcsec')




