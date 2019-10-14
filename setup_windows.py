description="""This script will setup magi for windows. 
It will write local settings to a file in the local settings folder and it will extract all databases that are needed for magi.
"""

import os
import shutil
import zipfile
import argparse

def main():
    parser = argparse.ArgumentParser(description=description)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-f', help='full setup. This is the default', action='store_true')
    group.add_argument('-d', help='database only', action='store_true', default=False)
    group.add_argument('-p', help='write local settings only', action='store_true', default=False)
    args = parser.parse_args()

    repo_path = os.getcwd()
    
    if args.d:
        extract_database = True
        set_paths = False
    if args.p:
        extract_database = False
        set_paths = True
    if args.f:
        extract_database = True
        set_paths = True
    else:
        extract_database = True
        set_paths = True
    
    if set_paths:
        fname = raw_input('USER INPUT: Settings Name (leave blank for default): ')
    
    if extract_database:
        print("Starting to prepare database...")
    # step one: extract the db tarball
        db_dir = os.path.join(repo_path, 'workflow', 'database')
        # wipe out existing db directory
        if os.path.isdir(db_dir):
            print('Deleting old database files...')
            shutil.rmtree(db_dir)
        print 'Extracting database files...'
        os.makedirs(db_dir)
        for zip_file in ['blastdb.zip', 'pickle_files.zip']:
            zip_file = os.path.join(repo_path, 'workflow', zip_file)
            with zipfile.ZipFile(zip_file,"r") as zip_ref:
                zip_ref.extractall(os.path.join(repo_path, 'workflow','database'))
        print 'Done'

    if set_paths:
    # step two: make local_settings.py file
        if fname == '':
            fname = 'user_settings'

        lines = ["import os\n\n"]
        lines.append("repo_location = '%s'\n\n" % (repo_path))
        lines.append("blastbin =          os.path.join(repo_location, 'workflow','blastbin')\n")
        lines.append("\n")
        lines.append("refseq_path =       os.path.join(repo_location, 'workflow','database','mrs_reaction_filtered_refseq_db_newrxns_actinofixed.pkl')\n")
        lines.append("refseq_db =         os.path.join(repo_location, 'workflow','database','mrs_reaction_filtered_refseq_db_fasta_blastdb_actinofix')\n")
        lines.append("mrs_reaction_path = os.path.join(repo_location, 'workflow','database','mrs_reaction_newrxns_added_actinofix.pkl')\n")
        lines.append("compounds_df =      os.path.join(repo_location, 'workflow','database','unique_compounds_groups_magi.pkl')\n")
        lines.append("mst_path =          os.path.join(repo_location, 'workflow','database','graph.pkl')\n")
        lines.append("chemnet_pickle =    os.path.join(repo_location, 'workflow','database','compound_groups.pkl')\n")
        lines.append("c2r =               os.path.join(repo_location, 'workflow','database','c2r.pkl')\n")
        lines.append("\n")
        lines.append("magi_task_path = ''\n")
        lines.append("magiweburl = 'https://magi.nersc.gov'\n")
        lines.append("magiwebsuperuser = ''\n")
        lines.append("magiwebsuperuserpass = ''\n")
        lines.append("admin_email = ''\n")

        with open(os.path.join('local_settings',fname+'.py'), 'w') as f:
            for line in lines:
                f.write(line)
        with open(os.path.join('local_settings','local_settings.py'), 'a') as f:
            f.write("SETTINGS_FILE= '%s'\n" % (fname))
        print 'Successfully wrote local settings files'

    print 'Setup Done!'

if __name__ == '__main__':
    main()