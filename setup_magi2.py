import os
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='full setup', action='store_true')
    parser.add_argument('-d', help='database only', action='store_true', default=False)
    args = parser.parse_args()

    repo_path = os.getcwd()

    if args.d:
        args.f = False
    else:
        args.f = True

    if args.f:
        fname = input('USER INPUT: Settings Name (leave blank for default): ')

    if args.f or args.d:
        # step one: extract the db tarball
        print('Extracting database files...')
        db_dir = os.path.join(repo_path, 'workflow_2', 'database')
        # wipe out existing db directory
        if os.path.isdir(db_dir):
            subprocess.call(['rm', '-r', db_dir])
        os.makedirs(db_dir)
        for zipfile in ['reference_sequences.zip', 'MAGI_database.db.zip']:
            cmd = ['unzip', os.path.join("..",zipfile)]
            subprocess.call(
                cmd,
                cwd = os.path.join(repo_path, 'workflow_2','database')
                )
        print('Done')

    if args.f:
        # step two: make local_settings.py file

        if fname == '':
            fname = 'magi_2_user_settings'

        lines = ["import os\n\n"]
        lines.append("# Location where MAGI is stored locally \n")
        lines.append("repo_location = '%s'\n\n" % (repo_path))
        lines.append("# Location where NCBI BLAST tools are stored \n")
        lines.append("blastbin =          os.path.join(repo_location, 'workflow/blastbin')\n")
        lines.append("\n")
        lines.append("# Database with UniProt reference sequences of proteins that have a Rhea reation\n")
        lines.append("refseq_path =       os.path.join(repo_location, 'workflow_2','database','reference_sequences','reaction_to_gene_reference.csv')\n")
        lines.append("refseq_db =         os.path.join(repo_location, 'workflow_2','database','reference_sequences','rhea2uniprot.db')\n")

        with open('local_settings/%s.py' %(fname), 'w') as f:
            for line in lines:
                f.write(line)
        with open('local_settings/local_settings.py', 'a') as f:
            f.write("SETTINGS_FILE= '%s'\n" % (fname))
        print('Successfully wrote local settings files')

    print('Setup Done!')

if __name__ == '__main__':
    main()
