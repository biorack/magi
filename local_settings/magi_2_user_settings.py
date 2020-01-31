import os
# Location where MAGI is stored locally
repo_location = '/Users/northenlab/Desktop/h_leegwater/magi'

magi_database = "/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/database/C2R_database_noHs_v20200128.db"

# Location where NCBI BLAST tools are stored
blastbin =          os.path.join(repo_location, 'workflow/blastbin')

# Database with UniProt reference sequences of proteins that have a Rhea reation
refseq_path =       "/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/refseqs/reaction_to_gene_reference.csv"
refseq_db =         "/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/refseqs/rhea2uniprot.db"
#mrs_reaction_path = os.path.join(repo_location, 'workflow/database/mrs_reaction_newrxns_added_actinofix.pkl')
#compounds_df =      os.path.join(repo_location, 'workflow/database/unique_compounds_groups_magi.pkl')
#mst_path =          os.path.join(repo_location, 'workflow/database/magi_cpd_similarity.graphml')
#chemnet_pickle =    os.path.join(repo_location, 'workflow/database/compound_groups.pkl')
#c2r =               os.path.join(repo_location, 'workflow/database/c2r.pkl')

#magi_results_storage = os.path.join(repo_location, 'outputs')
#magi_task_path = ''
#magiweburl = 'https://magi.nersc.gov'
#magiwebsuperuser = ''
#magiwebsuperuserpass = ''
#admin_email = ''
