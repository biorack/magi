1. log on as pasteur (on `cori17`)
2. navigate to repo: `~/repos/magi`. Should be checked out to `pasteur` branch
3. pause CRON jobs
4. check out `master` branch and pull from remote
5. deal with any merge conflicts (if any, should just be for `local_settings` path)
6. check out `pasteur` branch
7. merge in master: `git merge origin/master`
8. deal with any merge conflicts (if any, should just be for `local_settings` path)
9. run python script `pasteur_setup.py`
10. test
11. turn on CRON jobs
