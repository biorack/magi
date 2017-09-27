import sys
import os
import pandas as pd
import socket
import datetime
import subprocess
import re
# utils is in same dir
import utils

base_url = 'https://magi.nersc.gov/'
magi_task_root = '/project/projectdirs/metatlas/projects/magi_tasks'

# get machine name to infer how to submit the job
host = socket.gethostname()

# set umask for process
os.umask(002)

# get jobs from magiweb
all_jobs = utils.retrieve_jobs()

# mirror
utils.mirror_inputs(all_jobs)

# adjust paths and remove jobs where the input files could not be found
all_jobs = utils.adjust_file_paths(all_jobs)

# keep jobs that havent been submitted yet
keep = []
for job in all_jobs:
	if not job['fields']['runflag']:
		keep.append(job)
if len(keep) == 0:
    sys.exit()
else:
	unrun_jobs = keep

# keep only jobs that need a job script made
script_jobs, mass_search = utils.jobs_to_script(all_jobs)

# load up compound dataframe if necessary
if mass_search:
    reference_compounds = pd.read_pickle(utils.my_settings.compounds_df)

# process each job
for job in script_jobs:
    # determine fasta language and translate if needed
    if job['fields']['fasta_file'] != '':
        job = utils.determine_fasta_language(job)
    
    # conduct accurate mass search if needed
    if job['fields']['metabolite_file'] != '':
        if job['fields']['is_mass_search']:
            proceed = utils.accurate_mass_checkpoint(job)
            if proceed:
                try:
                    job = utils.accurate_mass_search_wrapper(job, reference_compounds)
                except RuntimeError as e:
                    job_link = os.path.join(base_url, 'jobs/?id=%s' % (job['pk']))
                    subj = 'Error processing your MAGI job'
                    if e.args[0] == 'too many compounds':
                        fname = 'too many compounds'
                        msg = 'Your compound search resulted in too many '
                        msg += 'compounds, please reduce the number of adducts or '
                        msg += 'lower the ppm by editing your job here: %s. ' % (job_link)
                        msg += 'You can reply to this email for more help. '
                        msg += 'Thanks for using MAGI!' 
                    else:
                        raise e
                    utils.email_user(job['fields']['email'], subj, msg)
                    utils.save_job_params(job, fname=fname)
                    print 'emailed \n %s \n %s' % (job['fields']['email'], msg)
                    continue
            else:
                continue
        n_compounds = pd.read_csv(job['fields']['metabolite_file']).shape[0]
    else:
        n_compounds = 0
    # create job script
    utils.job_script(job, n_cpd=n_compounds)

# submission log is list of lists in order:
# [magi web pk, nersc submission id, cori/genepool]
submission_log = []
admin_msg = ''
for job_data in unrun_jobs:
    submit = False
    job_path = utils.get_job_dir(job_data)
    script_dir = os.path.join(magi_task_root, job_path, 'admin')
    listdir = os.listdir(script_dir)
    pk = job_data['pk']
    job_script = [x for x in listdir if 'job_script' in x]
    if len(job_script) == 0:
        continue
    else:
        job_script = job_script[0]
    script_path = os.path.join(script_dir, job_script)
    submit_protocol = job_script.split('.')[1]
    if submit_protocol not in ['sbatch']:
        print 'Could not determine which machine to submit the following script to: %s' % (script_path)
        continue

    if job_script.split('.')[1] == 'sbatch' and 'cori' in host:
        cmd = ['sbatch']
        cmd.append(script_path)
        submit = True

    if submit: 
        # submit the job
        submit_time = datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S.%f')
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except OSError as e:
            print cmd
            raise e

        out, err = p.communicate()
        
        # change the runflag
        utils.change_params(pk, 'runflag', 'True')

        # successful submission
        if 'submitted' in out.lower():
            pattern = 'job ([\d]*)'
            match = re.search(pattern, out)
            nersc_id = match.groups()[0]
            submission_log.append([pk, nersc_id, host, submit_time])
            job_link = os.path.join(base_url, 'jobs/?id=%s' % (pk))
            # email the user a success email
            subj = 'MAGI Job submitted!'
            msg = 'Hello, your MAGI job was just submitted to NERSC:\n'
            msg += job_link + '\n'
            msg += 'Your job should start soon depending on the queue.\n'
            msg += 'If you have any questions, please contact us by replying to this email.\n'

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