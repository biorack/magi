"""
This script mirrors magi web inputs to NERSC
"""

import sys
import os
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

# set umask for process
os.umask(002)

# get jobs from magiweb
all_jobs = utils.retrieve_jobs()

utils.mirror_inputs(all_jobs)

