"""
This script processes magi web jobs and writes job scripts for them
"""

from __future__ import print_function
import pandas as pd
import os
import sys
sys.path.insert(
    0,
    '/project/projectdirs/metatlas/projects/metatlas_reactions/')
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
        email = job['fields']['email']
        job = utils.accurate_mass_search_wrapper(job, reference_compounds)
        if job == 'too many compounds':
            msg = 'Your compound search resulted in too many compounds, please reduce the number of adducts or lower the ppm'
            utils.email_user(email, msg)
            continue
    # create job script
    utils.job_script(job)
