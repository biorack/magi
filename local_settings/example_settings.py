"""
Two files need to be created.  Neither of these will be in the repo.

A file called "local_settings.py" should be put in this directory with a single line in it.
This line will define a settings.py file that has all the variables declared to interact with the database

For example:
	SETTINGS_FILE = 'user_settings_bpb_laptop'

In this directory the 
For local setup, in the SETTINGS_FILE 'user_settings_bpb_laptop.py' you might find content like this:
	metatlas_rxn_path = '/Users/bpb/repos/metatlas_reactions/'
	db_password_file = None
	if db_password_file:
	     with open(db_passwd_file) as fid:
	         pw = fid.read().strip()
	db_username = None
	db_name = 'bpb_workspace.db'
	sql_program = 'sqlite'
	sql_program = 'mysql+pymysql'
	db_path = '%s:///%s'%(sql_program,db_name)

For production setup as database administrator, in the SETTINGS_FILE 'user_settings_admin_NERSC.py' you might find content like this:
	import os
	metatlas_rxn_path = '/Users/bpb/repos/metatlas_reactions/'
	db_password_file = os.path.join(metatlas_rxn_path,'db_admin_auth.txt')
	if db_password_file:
	     with open(db_passwd_file) as fid:
	         pw = fid.read().strip()
	db_username = 'admin'
	db_name = 'metatlas_rxn.db'
	sql_program = 'mysql+pymysql'
	#db_path = '%s://%s:%s@scidb1.nersc.gov/%s' % (sql_program,db_username,pw, db_name)

Finally, for production setup as database user, in the SETTINGS_FILE 'user_settings_user_NERSC.py' you might find content like this:
	import os
	metatlas_rxn_path = '/Users/bpb/repos/metatlas_reactions/'
	db_password_file = os.path.join(metatlas_rxn_path,'db_user_auth.txt')
	if db_password_file:
	     with open(db_passwd_file) as fid:
	         pw = fid.read().strip()
	db_username = 'metatlas_rxn_user'
	db_name = 'metatlas_rxn.db'
	sql_program = 'mysql+pymysql'
	#db_path = '%s://%s:%s@scidb1.nersc.gov/%s' % (sql_program,db_username,pw, db_name)

""" 