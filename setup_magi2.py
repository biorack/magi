import os
import shutil
from zipfile import ZipFile
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='full setup', action='store_true')
    parser.add_argument('-d', help='database only', action='store_true', default=False)
    parser.add_argument('-s', help='local settings only', action='store_true', default=False)
    args = parser.parse_args()
    # Get path to MAGI
    repo_path = os.getcwd()

    if args.d or args.s:
        args.f = False
    else:
        args.f = True

    if args.f or args.s:
        fname = input('USER INPUT: Settings Name (leave blank for default): ')

    if args.f or args.d:
        # step one: extract the database files
        print('Extracting database files...')
        db_dir = os.path.join(repo_path, 'workflow_2', 'database')
        # wipe out existing db directory
        if os.path.isdir(db_dir):
            shutil.rmtree(db_dir)
        os.mkdir(db_dir)
        for zipfile in ['rhea2uniprot.db.zip', 'reaction_to_gene_reference.zip', 'MAGI_database.zip']:
            with ZipFile(os.path.join(repo_path, "workflow_2", zipfile), "r") as zip_file:
                zip_file.extractall(db_dir)
        print('Done')

    if args.f or args.s:
        # step two: make local_settings.py file

        if fname == '':
            fname = 'magi_2_user_settings'

        lines = ["import os\n\n"]
        lines.append("# Location where MAGI is stored locally \n")
        lines.append("repo_location = r'%s'\n\n" % (repo_path))
        lines.append("# Location where NCBI BLAST tools are stored \n")
        lines.append("blastbin = os.path.join(repo_location, 'workflow','blastbin')\n")
        lines.append("\n")
        lines.append("# Location where MAGI database is stored \n")
        lines.append("magi_database = os.path.join(repo_location, 'workflow_2','database','MAGI_database.db')")
        lines.append("\n\n")
        lines.append("# Database with UniProt reference sequences of proteins that have a Rhea reation\n")
        lines.append("refseq_path = os.path.join(repo_location, 'workflow_2','database','reaction_to_gene_reference.csv')\n")
        lines.append("refseq_db = os.path.join(repo_location, 'workflow_2','database','rhea2uniprot.db')\n")

        with open('local_settings/%s.py' %(fname), 'w') as f:
            for line in lines:
                f.write(line)
        with open('local_settings/local_settings.py', 'a') as f:
            f.write("SETTINGS_FILE= '%s'\n" % (fname))
        print('Successfully wrote local settings files')

    print('Setup Done!')

if __name__ == '__main__':
    main()
