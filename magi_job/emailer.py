import utils
import os
import subprocess

base_url = utils.my_settings.magiweburl
magi_task_root = utils.my_settings.magi_task_path
MAGI_EMAIL = utils.my_settings.admin_email

# get only jobs that have been submitted
run_jobs = utils.retrieve_jobs(sift=[('runflag', 'True')])

run_jobs = utils.adjust_file_paths(run_jobs)

for job in run_jobs:
    job_path = utils.get_job_dir(job)
    admin_path = os.path.join(magi_task_root, job_path, 'admin')
    job_link = os.path.join(base_url, 'jobs/?id=%s' % (job['pk']))
    # check if a job ended with error code
    if os.path.isfile(os.path.join(admin_path, 'incomplete')):
        if not os.path.isfile(os.path.join(admin_path, 'incomplete_email.txt')):
            subj = 'MAGI JOB ERROR!'
            msg = ''
            msg += '%s' % (job_link)

            utils.email_user(MAGI_EMAIL, subj, msg)

            subprocess.call(['touch', '%s/incomplete_email.txt' % (admin_path)])
            continue
    # check if a job started
    if os.path.isfile(os.path.join(admin_path, 'start_time.txt')):
        if not os.path.isfile(os.path.join(admin_path, 'start_email.txt')):
            subj = 'MAGI Job Started!'
            msg = ''
            msg += 'Your MAGI job has started!\n'
            msg += '%s\n\n' % (job_link)
            msg += "You can monitor your job's progress by looking at the log files:\n"
            msg += 'OUTPUT: %s\n' % (os.path.join(base_url, 'files//processed/%s/' % (job_path), 'log_out.txt'))
            msg += 'ERROR: %s\n\n' % (os.path.join(base_url, 'files//processed/%s/' % (job_path), 'log_err.txt'))
            msg += 'If you have any questions, please contact us by replying to this email.\n'

            utils.email_user(job['fields']['email'], subj, msg)

            subprocess.call(['touch', '%s/start_email.txt' % (admin_path)])
    # check if a job successfully ended
    if os.path.isfile(os.path.join(admin_path, 'end_time.txt')):
        if not os.path.isfile(os.path.join(admin_path, 'end_email.txt')):
            subj = 'MAGI Job is Done!'
            msg = ''
            msg += 'Your MAGI job is done!\n'
            msg += 'See your results here: %s\n\n' % (job_link)
            msg += 'If you have any questions, please contact us by replying to this email.\n'

            utils.email_user(job['fields']['email'], subj, msg)

            subprocess.call(['touch', '%s/end_email.txt' % (admin_path)])
