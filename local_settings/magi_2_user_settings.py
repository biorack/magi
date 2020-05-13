import os
# Location where MAGI is stored locally
repo_location = 'C:\\Users\\hanne\\Documents\\magi'

magi_database = "C:\\Users\\hanne\\Downloads\\database_tables\\MAGI_database.db"

# Location where NCBI BLAST tools are stored
blastbin =          os.path.join(repo_location, 'workflow/blastbin')

# Database with UniProt reference sequences of proteins that have a Rhea reation
refseq_path =       "/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/refseqs/reaction_to_gene_reference.csv"
refseq_db =         "/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/refseqs/rhea2uniprot.db"

#magi_results_storage = os.path.join(repo_location, 'outputs')
#magi_task_path = ''
#magiweburl = 'https://magi.nersc.gov'
#magiwebsuperuser = ''
#magiwebsuperuserpass = ''
#admin_email = ''
