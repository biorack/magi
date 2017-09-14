"""
This script processes magi web jobs and writes job scripts for them
"""

from __future__ import print_function
import pandas as pd
import os
import sys
sys.path.insert(
    0,
    '/project/projectdirs/metatlas/projects/metatlas_reactions')
# load utils
import utils
# load local settings
from local_settings import local_settings as settings_loc
my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)

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
    reference_compounds = pd.read_pickle(my_settings.compounds_df)

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
                if e.args[0] == 'too many compounds':
                    job_link = 'https://magi-dev.nersc.gov/jobs/?id=%s' % (job['pk'])
                    msg = 'Your compound search resulted in too many compounds, please reduce the number of adducts or lower the ppm by editing your job here: %s. You can reply to this email for more help. Thanks for using MAGI!' % (job_link)
                    subj = 'Error processing your MAGI job'
                    utils.email_user(job['fields']['email'], subj, msg)
                    utils.save_job_params(job)
                    # need to stop this job from going again until fixed somehow
                    # idea 1: make a small txt file that is the json of their job
                    # and titled "email_sent"
                    # if "email_sent" is present, see if the new job params are the same
                    # if they are the same, see if their compound file is identical
                    # if either of those 2 conditionals arent met, delete the email_sent file
                    # and rerun accurate mass searching
                    continue
                else:
                    raise e
        else:
            continue
    
    n_compounds = pd.read_csv(job['fields']['metabolite_file']).shape[0]
    # create job script
    utils.job_script(job, n_cpd=n_compounds)
