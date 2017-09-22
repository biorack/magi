This directory holds scripts and utilities for interfacing with the MAGI website for running jobs at NERSC

A few things need to be changed when the URL for magi_web is changed/different/going from dev to prod:

* change `magiweburl` variable in `user_settings.py`
* change `base_url` variable in `process_web_jobs.py`; just under imports
* change `base_url` variable in `submit_jobs.py`; just under imports 

