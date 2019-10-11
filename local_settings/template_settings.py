import os

repo_location = ''
blastbin =          os.path.join(repo_location, 'workflow/blastbin')

refseq_path =       os.path.join(repo_location, 'workflow/database/mrs_reaction_filtered_refseq_db_newrxns_actinofixed.pkl')
refseq_db =         os.path.join(repo_location, 'workflow/database/mrs_reaction_filtered_refseq_db_fasta_blastdb_actinofix')
mrs_reaction_path = os.path.join(repo_location, 'workflow/database/mrs_reaction_newrxns_added_actinofix.pkl')
compounds_df =      os.path.join(repo_location, 'workflow/database/unique_compounds_groups_magi.pkl')
mst_path =          os.path.join(repo_location, 'workflow/database/graph.pkl')
chemnet_pickle =    os.path.join(repo_location, 'workflow/database/compound_groups.pkl')