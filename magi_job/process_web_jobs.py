"""
This script processes magi web jobs and writes job scripts for them
"""

from __future__ import print_function
import pandas as pd
import os
import sys
# load local utils
import utils

# set umask for python process
os.umask(002)

# get jobs from magiweb that havent been submitted yet
all_jobs = utils.retrieve_jobs(sift=[('runflag', False)])
if all_jobs is None:
    sys.exit()

# adjust paths and remove jobs where the input files could not be found
all_jobs = utils.adjust_file_paths(all_jobs)

# keep only jobs that need a job script made
all_jobs, mass_search = utils.jobs_to_script(all_jobs)
# load up compound dataframe if necessary
if mass_search:
    reference_compounds = pd.read_pickle(utils.my_settings.compounds_df)

# process each job
for job in all_jobs:
    # determine fasta language and translate if needed
    job = utils.determine_fasta_language(job)
    # conduct accurate mass search if needed
    if job['fields']['is_mass_search']:
        proceed = utils.accurate_mass_checkpoint(job)
        if proceed:
            try:
                job = utils.accurate_mass_search_wrapper(job, reference_compounds)
            except RuntimeError as e:
                job_link = 'https://magi-dev.nersc.gov/jobs/?id=%s' % (job['pk'])
                subj = 'Error processing your MAGI job'
                if e.args[0] == 'too many compounds':
                    fname = 'too many compounds'
                    msg = 'Your compound search resulted in too many '
                    msg += 'compounds, please reduce the number of adducts or '
                    msg += 'lower the ppm by editing your job here: %s. ' % (job_link)
                    msg += 'You can reply to this email for more help. '
                    msg += 'Thanks for using MAGI!' 
                elif e.args[0] == 'original_compound not floatable':
                    fname = 'original_compound not floatable'
                    msg = 'You elected to conduct an accurate mass search, '
                    msg += 'but at least one of the values in the '
                    msg += 'original_compound column in your metabolite table '
                    msg += 'does not look like a number. \nPlease double-check '
                    msg += 'your input and edit/resubmit your job by '
                    msg += 'clicking here: %s.\n\n' % (job_link)
                    msg += 'You can reply to this email for more help. '
                    msg += 'Thanks for using MAGI!'
                else:
                    raise e
                utils.email_user(job['fields']['email'], subj, msg)
                utils.save_job_params(job, fname=fname)
                print('emailed \n %s \n %s' % (job['fields']['email'], msg))
                continue
        else:
            continue
    
    n_compounds = pd.read_csv(job['fields']['metabolite_file']).shape[0]
    # create job script
    utils.job_script(job, n_cpd=n_compounds)
