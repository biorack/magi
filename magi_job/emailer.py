import utils
import os
import subprocess
import glob

magi_task_root = '/project/projectdirs/metatlas/projects/magi_tasks'
base_url = base_url = 'https://magi.nersc.gov/'

# get only jobs that have been submitted
run_jobs = utils.retrieve_jobs(sift=[('runflag', 'True')])

run_jobs = utils.adjust_file_paths(run_jobs)

for job in run_jobs:
	job_path = utils.get_job_dir(job)
	admin_path = os.path.join(magi_task_root, job_path, 'admin')
	job_link = os.path.join(base_url, 'jobs/?id=%s' % (job['pk']))
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

			subprocess.call(['date', '>', 'start_email.txt'],
				cwd='%s' % (admin_path))
			
			print 'emailed %s start email for job %s' % (job['fields']['email'], job['pk'])
			continue
	result_files = glob.glob(os.path.join(magi_task_root, job_path, '*results.csv'))
	if len(result_files) > 0:
		if not os.path.isfile(os.path.join(admin_path, 'end_email.txt')):
			subj = 'MAGI Job is Done!'
			msg = ''
			msg += 'Your MAGI job is done!\n'
			msg += 'See your results here: %s\n\n' % (job_link)
			msg += 'If you have any questions, please contact us by replying to this email.\n'

			utils.email_user(utils.email_user(job['fields']['email'], subj, msg))

			subprocess.call(['date', '>', 'end_email.txt'],
				cwd='%s' % (admin_path))

			print 'emailed %s end email for job %s' % (job['fields']['email'], job['pk'])
			continue
