"""
This script mirrors magi web inputs to NERSC
"""

import sys
import os
# utils is in same dir
import utils

# set umask for process
os.umask(002)

# get jobs from magiweb
all_jobs = utils.retrieve_jobs()

utils.mirror_inputs(all_jobs)

