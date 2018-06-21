####################################
### REDUCTION AND CALIBRATION OF ###
###            CIG96             ###
####################################
#
# CIG96 OBSERVING NIGHTS
# date_start, date_end = [datetime.datetime(2016, 10, 6),  datetime.datetime(2016, 10, 7)]
# date_start, date_end = [datetime.datetime(2016, 10, 9),  datetime.datetime(2016, 10, 10)]
# date_start, date_end = [datetime.datetime(2016, 10, 20), datetime.datetime(2016, 10, 21)]
# date_start, date_end = [datetime.datetime(2016, 11, 1),  datetime.datetime(2016, 11, 2)]
# date_start, date_end = [datetime.datetime(2016, 11, 2),  datetime.datetime(2016, 11, 3)]
# date_start, date_end = [datetime.datetime(2016, 12, 2),  datetime.datetime(2016, 12, 3)]
# date_start, date_end = [datetime.datetime(2016, 12, 3),  datetime.datetime(2016, 12, 4)]
# date_start, date_end = [datetime.datetime(2016, 12, 20), datetime.datetime(2016, 12, 21)]
#
#-----------------------------------------------------------------------------

### PROJECT SETTING
context.set_project('AMIGADEEP')

### PRIVILEGES SETTING (1=PROPRIETARY, 2=SHARED)
context.set_privileges(2)

### IMPORT FUNCTIONS
exec(open('/data/users/ramirez/support_scripts/singlejobsender3.py').read())    
exec(open('/data/users/ramirez/support_scripts/select_frames3.py').read())

#-----------------------------------------------------------------------------

### TARGET AND COORDINATES DEFINITION
target = 'CIG0096'
cigra, cigdec = 33.8, 6.0        # must have only 1 decimal; define grid center (in degrees) for later use in 'singlejobsender'

#-----------------------------------------------------------------------------

### PHOTOMETRIC SOLUTION FROM ALL NIGHTS WITH DATA FROM SELECTED TARGET

### Night 1
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 10, 6),  datetime.datetime(2016, 10, 7)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#############

### Night 2
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 10, 9),  datetime.datetime(2016, 10, 10)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#############

### Night 3
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 10, 20), datetime.datetime(2016, 10, 21)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#############

### Night 4
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 11, 1),  datetime.datetime(2016, 11, 2)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#############

### Night 5
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 11, 2),  datetime.datetime(2016, 11, 3)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#############

### Night 6
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 12, 2),  datetime.datetime(2016, 12, 3)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#############

### Night 7
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 12, 3),  datetime.datetime(2016, 12, 4)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#############

### Night 8
# date-hour selection
date_start, date_end = [datetime.datetime(2016, 12, 20), datetime.datetime(2016, 12, 21)]

# query raw data in previous date-hour selection
date_query = ((RawScienceFrame.DATE_OBS>date_start)&(RawScienceFrame.DATE_OBS<date_end)&(RawScienceFrame.OBJECT=='STD,ZEROPOINT'))

# size of the date_query (number of files)
len(date_query)                   

# name of files
names = [i.filename for i in date_query]
len(names)

# reduce first query
singlejobsender('Reduce', names, cigra, cigdec)

# query reduced data
query = ((ReducedScienceFrame.creation_date>datetime.date.today())&(ReducedScienceFrame.OBJECT=='STD,ZEROPOINT'))
names=[i.filename for i in query]

# astrometry and photometry of reduced data
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Photom', names, cigra, cigdec)

#-----------------------------------------------------------------------------

### REDUCTION

query_raw = (RawScienceFrame.OBJECT==target)
names=[i.filename for i in query_raw]
singlejobsender('Background', names, cigra, cigdec) 
singlejobsender('Reduce', names, cigra, cigdec)


query_red = select_frames('ReducedScienceFrame', 'OCAM_r_SDSS', cigra-0.3, cigra+0.3, cigdec-0.3, cigdec+0.3)
names=[i.filename for i in query_red]
singlejobsender('Astrometry', names, cigra, cigdec)
singlejobsender('Regrid', names, cigra, cigdec) 

#-----------------------------------------------------------------------------

### FRAMES VISUALIZATION
# Now you have frames to be stacked as a final 'coadd' but first, it is good to discard any frames with outlying intensity levels or psf fwhms:
psf = [i.psf_radius for i in query_red]
import matplotlib.pyplot as plt
plt.hist(psf)
medians = [i.imstat.median for i in query_red]
plt.hist(medians)

#-----------------------------------------------------------------------------

### GOOD FRAMES SELECTION
query_all = select_frames('RegriddedFrame', cigra-0.3, cigra+0.3, cigdec-0.3, cigdec+0.3)   # all files within 0.3 radius from the specified coordinates

goodseeing = [ ]

for i in query_all:
      if i.psf_radius < 1.2:goodseeing.append(i.filename)   # defines

#-----------------------------------------------------------------------------

### COADDITION
dpu.run('Coadd', instrument='OMEGACAM',reg_filenames=goodseeing, dpu_aweversion='current', dpu_time = 8*60*60, C=1)

#-----------------------------------------------------------------------------

### DPU STATUS

# You can check the status of the reduction with:
dpu.get_status()

### FINAL IMAGE
# When your reduction is ready you can retrieve it as follows:
query_final = (CoaddedRegriddedFrame.OBJECT==target).project_only().max('creation_date')
query_final.retrieve()

### DOWNLOAD LINK
# If you would rather get a download link you to be used e.g. with wget:
import urllib
url = 'http://%s:%s/%s' % (Env['data_server'], Env['data_port'], urllib.parse.quote(query_final.filename))
print(url)
