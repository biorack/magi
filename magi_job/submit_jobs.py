"""
This script submits appropriate jobs when its on the appropriate machine
"""

import sys
# load local utils
import utils
import os
import subprocess
import socket
import re
import datetime

# get machine name to infer how to submit the job
host = socket.gethostname()

magi_task_root = '/project/projectdirs/metatlas/projects/magi_tasks'

all_jobs = utils.retrieve_jobs(sift=[('runflag', False)])
if all_jobs is None:
	sys.exit()

os.umask(002)

# submission log is list of lists in order:
# [magi web pk, nersc submission id, cori/genepool]
submission_log = []
admin_msg = ''
for job_data in all_jobs:
    submit = False
    uptime = job_data['fields']['uploaded_at']
    y = uptime.split('-')[0]
    m = uptime.split('-')[1]
    pk = job_data['pk']
    script_dir = os.path.join(magi_task_root, y, m, pk, 'admin')
    listdir = os.listdir(script_dir)
    try:
        job_script = [x for x in listdir if 'job_script' in x][0]
    except Exception as e:
        print script_dir
        print e.args
        continue
    script_path = os.path.join(script_dir, job_script)
    submit_protocol = job_script.split('.')[1]
    if submit_protocol not in ['sbatch', 'qsub']:
        print 'Could not determine which machine to submit the following script to: %s' % (script_path)
        continue

    if job_script.split('.')[1] == 'sbatch' and 'cori' in host:
        cmd = ['sbatch']
        cmd.append(script_path)
        submit = True
    elif job_script.split('.')[1] == 'qsub' and 'genepool' in host:
        cmd = ['qsub']
        cmd.append(script_path)
        submit = True

    if submit: 
        # submit the job
        submit_time = datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S.%f')
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        
        # change the runflag
        utils.change_params(pk, 'runflag', 'True')

        # successful submission
        if 'submitted' in out.lower():
            pattern = 'job ([\d]*)'
            match = re.search(pattern, out)
            nersc_id = match.groups()[0]
            submission_log.append([pk, nersc_id, host, submit_time])

            # email the user a success email
            subj = 'MAGI Job submitted!'
            msg = 'Hello, your MAGI job was just submitted to NERSC:\n'
            msg += 'https://magi-dev.nersc.gov/jobs/?id=%s\n' % (pk)
            msg += 'You should receive an email from NERSC soon when your job starts and ends.\n'
            msg += '\n'
            msg += "When your job starts at NERSC, you can monitor your job's progress by looking at the log files\n"
            msg += 'OUTPUT: https://magi-dev.nersc.gov/files//processed/%s/%s/%s/log_out.txt\n' % (y, m, pk)
            msg += 'ERROR: https://magi-dev.nersc.gov/files//processed/%s/%s/%s/log_err.txt\n\n' % (y, m, pk)
            msg += 'If you have any questions or if your job fails at NERSC (you will get an email), please contact us by replying to this email.\n'

        # unsuccessful submission
        else:
            admin_msg += 'ERROR submitting magi web job %s\n' % (pk)
            admin_msg += 'STDOUT: %s\n' % (out)
            admin_msg += 'STDERR: %s\n' % (err)
            admin_msg += 'submit command: %s\n\n' %(cmd)

            subj = 'Could not submit MAGI job'
            msg = 'Sorry, there was a problem submitting your MAGI job.\n'
            msg += 'The MAGI team will try to figure out what is wrong ASAP!\n'

        # send the email
        utils.email_user(job_data['fields']['email'], subj, msg)
# if needed, notify admin
if admin_msg != '':
    # for dev, just print to stdout
    print admin_msg
    # for prod, email magi_web
    # subj = 'ERROR SUBMITTING MAGI WEB JOB(S)'
    # utils.email_user('magi_web@lbl.gov', subj, msg)

# finally, write to the logfile
os.umask(007)
with open(os.path.join(magi_task_root, 'submission_log.txt'), 'a') as f:
    for j in submission_log:
        f.write(','.join(j) + '\n')
