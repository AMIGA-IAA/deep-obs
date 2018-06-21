# NECESSARY MODULES
import matplotlib.pyplot as plt
import urllib

# SET PROJECT, PRIVILEGES AND IMPORT NECESSARY FEEDING FUNCTIONS

context.set_project('AMIGADEEP')
context.set_privileges(2)
exec(open('/data/users/ramirez/support_scripts/singlejobsender3.py').read())    
exec(open('/data/users/ramirez/support_scripts/select_frames3.py').read())

# CIG33 DATA

target = 'CIG0033'
cigra, cigdec = 10.9, -0.1

# PHOTOMETRY FOR ALL NIGHTS
# 7 observations from NLT observations progress:
# 	1418809		2/3-Nov-2016 	05:02:26 > 06:00:49 	CIG0033_2 	A 	weather comments
#	1418812 	3/4-Nov-2016 	05:22:59 > 06:22:14 	CIG0033_3 	B 	weather comments
#	1418815 	4/5-Nov-2016 	01:09:35 > 02:07:32 	CIG0033_4 	B 	weather comments
#	1418818 	4/5-Nov-2016 	02:07:34 > 03:03:17 	CIG0033_5 	B 	weather comments
#	1418821 	1/2-Dec-2016 	00:44:17 > 01:47:05 	CIG0033_6 	A 	weather comments
#	1418824 	1/2-Dec-2016 	01:47:11 > 02:44:20 	CIG0033_7 	B 	weather comments
#	1418805 	1/2-Dec-2016 	03:01:40 > 03:59:46 	CIG0033_1 	A

# 4 nights in total: 2/3-Nov-2016, 3/4-Nov-2016, 4/5-Nov-2016, 1/2-Dec-2016.

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

date_start, date_end = [datetime.datetime(2016, 11, 3),  datetime.datetime(2016, 11, 4)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2016, 11, 4),  datetime.datetime(2016, 11, 5)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2016, 12, 1),  datetime.datetime(2016, 12, 2)]
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
#date = datetime.datetime(2017,6,X)      # Replace the X with the day you started producing these calibrations
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

len(raw), len(regr), len(noregr)

# if noregr != 0 then
regfiles = [i.filename for i in noregr]
singlejobsender('Regrid', regfiles, cigra, cigdec)

########################################################

# FRAMES VISUALIZATION

psf = [i.psf_radius for i in query_red]
plt.hist(psf)

medians = [i.imstat.median for i in query_red]
plt.hist(medians)

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

dpu.run('Coadd', instrument='OMEGACAM',reg_filenames=goodseeing, dpu_aweversion='current', dpu_time = 8*60*60, C=1)

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

#


# REGRIDDING CORRECTION (DONE ONLY ONCE)

q=(CoaddedRegriddedFrame.filename=='Sci-PRAMIREZMORETA-OMEGACAM-------OCAM_r_SDSS---Coadd-Astrom-Sci-57940.4344097-bc9b1cfb514e8217e3b716ada320ca330f51e4c3.fits')[0]

reg= q.regridded_frames

cigra,cigdec = reg[0].astrom.CRVAL1,reg[0].astrom.CRVAL2

names = [i.reduced.filename for i in reg]

singlejobsender('Regrid', names, cigra, cigdec)

query_all = select_frames('RegriddedFrame', 'OCAM_r_SDSS', cigra-0.8, cigra+0.8, cigdec-0.8, cigdec+0.8)

goodseeing = [ ]
for i in query_all:
   if (i.psf_radius < 2.0) & (i.reduced.back!= None):goodseeing.append(i.filename)


dpu.run('Coadd', instrument='OMEGACAM',reg_filenames=goodseeing, dpu_aweversion='develop', dpu_time = 12*60*60, C=1)


