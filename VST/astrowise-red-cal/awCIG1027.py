# NECESSARY MODULES
import matplotlib.pyplot as plt
import urllib

# SET PROJECT, PRIVILEGES AND IMPORT NECESSARY FEEDING FUNCTIONS

context.set_project('AMIGADEEP')
context.set_privileges(2)
exec(open('/data/users/ramirez/support_scripts/singlejobsender3.py').read())    
exec(open('/data/users/ramirez/support_scripts/select_frames3.py').read())

# CIG1027 DATA

target = 'CIG1027'
cigra, cigdec = 353.8, 7.3

# PHOTOMETRY FOR ALL NIGHTS
# 7 observations from NLT observations progress:
# 	1227585 	12/13-Sep-2015 	06:26:58 > 07:26:26 	CIG1027_3 	B 	weather comments 
#	1418723 	25/26-Oct-2016 	00:20:20 > 01:21:57 	CIG1027_2 	B 	weather comments
#	1418726 	25/26-Oct-2016 	02:26:03 > 03:27:54 	CIG1027_4 	A 	weather comments
#	1418729 	1/2-Nov-2016 	03:01:48 > 04:00:47 	CIG1027_5 	A 	weather comments
#	1418732 	1/2-Nov-2016 	04:00:57 > 04:57:29 	CIG1027_6 	B 	weather comments
#	1418735 	2/3-Nov-2016 	03:02:18 > 04:03:58 	CIG1027_7 	A 	weather comments
#	1418719 	3/4-Nov-2016 	01:08:56 > 02:10:04 	CIG1027_1 	A 	weather comments 
# 5 nights in total: 12/13-Sep-2015, 25/26-Oct-2016, 1/2-Nov-2016, 2/3-Nov-2016, 3/4-Nov-2016

date_start, date_end = [datetime.datetime(2015, 9, 12),  datetime.datetime(2015, 9, 13)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2016, 10, 25),  datetime.datetime(2016, 10, 26)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

date_start, date_end = [datetime.datetime(2016, 11, 1),  datetime.datetime(2016, 11, 2)]
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))
len(date_query)                   
names = [i.filename for i in date_query]
len(names)
singlejobsender('Reduce', names, cigra, cigdec)
query = ((ReducedScienceFrame.DATE_OBS>date_start)&(ReducedScienceFrame.DATE_OBS<date_end)&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

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

len(red), len(regr), len(noregr)

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

#



#cod = (CoaddedRegriddedFrame.OBJECT=='CIG0279')
#cod = cod[3]     # or whatever the number, the last image in any case
#red = cod.regridded_frames
#names = [i.reduced.filename for i in red]
#singlejobsender('Regrid', names, cigra, cigdec)
#query_all = select_frames('RegriddedFrame', 'OCAM_r_SDSS', cigra-0.8, cigra+0.8, cigdec-0.8, cigdec+0.8)
#goodseeing = [ ]
#for i in query_all:
#   if (i.psf_radius < 2.0) & (i.reduced.back!= None):goodseeing.append(i.filename)
#dpu.run('Coadd', instrument='OMEGACAM',reg_filenames=goodseeing, dpu_aweversion='current', dpu_time = 8*60*60, C=1)


