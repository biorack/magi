"""
This script prepares and sends a summary table of MAGI web jobs submitted
in the past week.

There is some serious timezone confusion because of the way I get file
timestamps, how python datetime works, etc. Make sure to run this script
before 2am Sundays to ensure there isn't a timezone mismatch within one
report when daylight savings time changes.

"""

import utils
import pandas as pd
import datetime
import os

magi_task_root = '/project/projectdirs/metatlas/projects/magi_tasks'
job_table = os.path.join(magi_task_root, 'submission_log.txt')

today = datetime.datetime.today()
today = today.replace(hour=0, minute=0, second=0, microsecond=0)
one_week = today - datetime.timedelta(days=7)
y = one_week.strftime('%Y')
m = one_week.strftime('%m')
d = one_week.strftime('%d')

jobs = utils.retrieve_jobs(sift=[
	('year_gte', y),
	('month_gte', m),
	('day_gte', d)
	])

# load submission table and only keep jobs submitted in the last week
df = pd.read_csv(job_table)
df['submit_time'] = pd.to_datetime(df['submit_time'])
df = df[df['submit_time'] >= one_week]

import pickle
with open('./tmp.pkl', 'w') as f:
	pickle.dump(jobs, f)

# make the report table
data = {
    'magi_web_pk': [],
    'email': [],
    'upload_time': [],
    'job_dir': [],
    'start_time': [],
    'end_time': []
}
for j in jobs:
    data['magi_web_pk'].append(j['pk'])
    data['email'].append(j['fields']['email'])
    uptime = j['fields']['uploaded_at']
    data['upload_time'].append(pd.to_datetime(uptime))
    y = uptime.split('-')[0]
    m = uptime.split('-')[1]
    jdir = os.path.join(magi_task_root, y, m, j['pk'])
    data['job_dir'].append(jdir)
    
    # job start time
    fname = os.path.join(jdir, 'admin', 'start_time.txt')
    if os.path.isfile(fname):
        with open(fname, 'r') as f:
            start_time = f.read()
        start_time = datetime.datetime.strptime(start_time, '%a %b %d %H:%M:%S UTC %Y ')
        data['start_time'].append(start_time)
    else:
        data['start_time'].append(pd.np.nan)
    
    # job end time
    fname = os.path.join(jdir, 'magi_results.csv')
    if os.path.isfile(fname):
        timestamp = os.path.getmtime(fname)
        end_time = pd.to_datetime(datetime.datetime.fromtimestamp(timestamp))
        data['end_time'].append(end_time)
    else:
        data['end_time'].append(pd.np.nan)

a = pd.DataFrame(data)

# convert time zones
# convert utc time to pacific
a['upload_time'] = a['upload_time'].apply(lambda x: x.tz_localize('UTC').tz_convert('US/Pacific'))
a['start_time'] = a['start_time'].apply(lambda x: x.tz_localize('UTC').tz_convert('US/Pacific'))
# convert start and end times to pacific
a['end_time'] = a['end_time'].apply(lambda x: x.tz_localize('US/Pacific'))

# merge tables and calculate job time
b = df.merge(a, on='magi_web_pk')
b['job_runtime'] = b['end_time'] - b['start_time']

# email formatted table
