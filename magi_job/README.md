# This directory holds scripts and utilities for interfacing with the MAGI website for running jobs at NERSC

* Be mindful of which localsettings file you are pointed to, in particular what the `magiweburl` variable is (magi or magi-dev?)
* `magi_cron_tab.txt` is a copy of the crontab currently on `pasteur@corigrid.nersc.gov`
* `magi_web_mirror_process_submit.sh` runs a python script that gets the jobs, processes jobs that need to be processed, makes job scripts if necessary, submits jobs if necessary, and sends emails as necessary about this process
* `start_end_emailer.sh` runs a python script that checks if a job ended with an error, whether it started, and whether it successfully ended. Emails user and admin (if necessary).
* `utils.py` holds a bunch of helper functions for the scripts above

## How to set up Pasteur MAGI successfully
1. Log on to pasteur on corigrid
2. Go into `~/repos/magi`; this should be on the `pasteur` branch
3. `checkout` `master` branch, `pull`
4. `checkout` `pasteur` branch, `git merge origin/master`
5. Deal with any merge conflicts. You only need to keep the `pasteur` version of things if the conflict is regarding some sort of path.
5. Run `python pasteur_setup.py`