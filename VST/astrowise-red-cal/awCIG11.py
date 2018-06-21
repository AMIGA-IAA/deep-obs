# NECESSARY MODULES
import matplotlib.pyplot as plt
import urllib

# SET PROJECT, PRIVILEGES AND IMPORT NECESSARY FEEDING FUNCTIONS

context.set_project('AMIGADEEP')
context.set_privileges(2)
exec(open('/data/users/ramirez/support_scripts/singlejobsender3.py').read())    
exec(open('/data/users/ramirez/support_scripts/select_frames3.py').read())

# CIG0011 DATA

target = 'CIG0011'
cigra, cigdec = 3.63, -0.73

# PHOTOMETRY FOR ALL NIGHTS
# 7 observations from NLT observations progress:
#  	1418782 	2/3-Nov-2016 	04:03:59 > 05:01:31 	CIG0011_1 	A	weather comments
#	1418786 	2/3-Dec-2016 	02:33:12 > 03:32:58 	CIG0011_2 	A 	weather comments 
#	1418792 	28/29-Jun-2017 	08:51:29 > 09:49:15 	CIG0011_4 	A 	weather comments
#	1418795 	1/2-Jul-2017 	08:08:01 > 09:07:10 	CIG0011_5 	A 	weather comments
#	1418798 	1/2-Jul-2017 	09:07:11 > 10:03:23 	CIG0011_6 	A 	weather comments
#	1418789 	19/20-Aug-2017 	08:48:51 > 09:47:29 	CIG0011_3 	A 	weather comments
#	1418801 	15/16-Oct-2017 	01:42:57 > 02:41:24 	CIG0011_7 	A 	weather comments 
# 6 nights in total: 2/3-Nov-2016, 2/3-Dec-2016, 28/29-Jun-2017, 1/2-Jul-2017, 19/20-Aug-2017, 15/16-Oct-2017

date_start, date_end = [datetime.datetime(2016, 11, 2),  datetime.datetime(2016, 11, 3)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2016, 12, 2),  datetime.datetime(2016, 12, 3)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2017, 6, 28),  datetime.datetime(2017, 6, 29)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2017, 7, 1),  datetime.datetime(2017, 7, 2)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2017, 8, 19),  datetime.datetime(2017, 8, 20)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2017, 10, 15),  datetime.datetime(2017, 10, 16)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

# if there is a night change (tasks last for more than the present day):
#date = datetime.datetime(2017, 6, X)      # Replace the X with the day you started producing these calibrations
#query = ((ReducedScienceFrame.creation_date >= date) & (ReducedScienceFrame.OBJECT=='STD,ZEROPOINT')).project_only()

########################################################

# BACKGROUND, REDUCTION, ASTROMETRY, REGRIDDING

query_raw = (RawScienceFrame.OBJECT==target)
names=[i.filename for i in query_raw]
len(names)

########################################################

# BACKGROUND MODELS
singlejobsender('Background', names, cigra, cigdec) 

# backgrounds check

raw = RawScienceFrame.OBJECT==target
backs, noback = [], []
for i in raw:
  q = (BackgroundFrame.raw == i)
  if len(q):
     backs.append(q[0])
  else:
     noback.append(i)

len(raw), len(backs), len(noback)

# if noback != 0 then:
backfiles = [i.filename for i in noback]
singlejobsender('Background', backfiles, cigra, cigdec)

########################################################

# REDUCTION

query_raw = (RawScienceFrame.OBJECT==target)
names=[i.filename for i in query_raw]
len(names)

singlejobsender('Reduce', names, cigra, cigdec)

# reduced frames check

raw = RawScienceFrame.OBJECT==target
red, nored = [], []
for i in raw:
  q = (ReducedScienceFrame.raw == i)
  if len(q):
    if q[0].back!=None:
      red.append(q[0])
    else:
      nored.append(i)
  else:
    nored.append(i)

len(raw), len(red), len(nored)

# if nored != 0 then:
redfiles = [i.filename for i in nored]
singlejobsender('Reduce', redfiles, cigra, cigdec)

########################################################

# ASTROMETRY

query_red = select_frames('ReducedScienceFrame', 'OCAM_r_SDSS', cigra-0.8, cigra+0.8, cigdec-0.8, cigdec+0.8)
names=[i.filename for i in query_red]

singlejobsender('Astrometry', names, cigra, cigdec)

# astrometry check

red = ReducedScienceFrame.OBJECT==target
astro, noastro=[ ], [ ]

for i in red:
   q = AstrometricParameters.reduced==i
   if len(q):
        astro.append(q[0])
   else:
        noastro.append(i)

len(red), len(astro), len(noastro)

# if noastro != 0 then:
astrofiles = [i.filename for i in noastro]
singlejobsender('Astrometry', astrofiles, cigra, cigdec)


########################################################

# REGRIDDING

singlejobsender('Regrid', names, cigra, cigdec) 

# regridded frames check

red = ReducedScienceFrame.OBJECT==target
regr,noregr=[ ],[ ]
for i in red:
   q = RegriddedFrame.reduced==i
   if len(q):
        regr.append(q[0])
   else:
        noregr.append(i)

len(red), len(regr), len(noregr)

# if noregr != 0 then
regfiles = [i.filename for i in noregr]
singlejobsender('Regrid', regfiles, cigra, cigdec)

########################################################

# FRAMES VISUALIZATION

#psf = [i.psf_radius for i in query_red]
#plt.hist(psf)

#medians = [i.imstat.median for i in query_red]
#plt.hist(medians)

########################################################

# PSF AND GOOD BACKGROUND MODEL SELECTION

query_all = select_frames('RegriddedFrame', 'OCAM_r_SDSS', cigra-0.8, cigra+0.8, cigdec-0.8, cigdec+0.8)
len(query_all)
goodseeing = [ ]
for i in query_all:
      if (i.psf_radius < 2.0) & (i.reduced.back!= None):goodseeing.append(i.filename)


len(goodseeing)

########################################################

# COADDITION

dpu.run('Coadd', instrument='OMEGACAM',reg_filenames=goodseeing, dpu_aweversion='develop', dpu_time = 8*60*60, C=1)

########################################################

# DPU STATUS

dpu.get_status()

########################################################

# FINAL IMAGE

query_final = (CoaddedRegriddedFrame.OBJECT==target).project_only().max('creation_date')
query_final.retrieve()

# DOWNLOAD LINK

url = 'http://%s:%s/%s' % (Env['data_server'], Env['data_port'], urllib.parse.quote(query_final.filename))
print(url)

########################################################

# REGRIDDING RE-CORRECTION (DONE ONLY ONCE)

q=(CoaddedRegriddedFrame.filename=='Sci-PRAMIREZMORETA-OMEGACAM-------OCAM_r_SDSS---Coadd-Astrom-Sci-58073.8510212-7311bc9b5a3c5c060b52a62bfbacd15ef8818131.fits')[0]

reg= q.regridded_frames

cigra,cigdec = reg[0].astrom.CRVAL1,reg[0].astrom.CRVAL2

names = [i.reduced.filename for i in reg]

singlejobsender('Regrid', names, cigra, cigdec)

query_all = select_frames('RegriddedFrame', 'OCAM_r_SDSS', cigra-0.8, cigra+0.8, cigdec-0.8, cigdec+0.8)

goodseeing = [ ]
for i in query_all:
   if (i.psf_radius < 2.0) & (i.reduced.back!= None):goodseeing.append(i.filename)


dpu.run('Coadd', instrument='OMEGACAM',reg_filenames=goodseeing, dpu_aweversion='develop', dpu_time = 12*60*60, C=1)

# FINAL IMAGE

query_final = (CoaddedRegriddedFrame.OBJECT==target).project_only().max('creation_date')
query_final.retrieve()



