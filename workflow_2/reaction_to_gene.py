# Reaction to gene

import sys
import os
import argparse
import pandas as pd
from multiprocessing import cpu_count as counting_cpus
import time
import blast_helpers as blast
import workflow_helpers as mg
import sqlite3

my_settings = mg.get_settings()

def prepare_compound_to_reaction_results(compound_to_reaction_path, output_dir, path_to_database=my_settings.magi_database):
    """ 
    This function opens the compound to reaction results file and extracts the retro rules reaction IDs.
    For these IDs, it will look up the Rhea reaction IDs, which are coupled to the MAGI reference genome database.
    It returns a data frame with rhea reaction IDs
    """
    # Get all retro rules reaction IDs from the c2r search
    reaction_ids = pd.read_csv(compound_to_reaction_path)
    reaction_ids = reaction_ids["reaction_ID"].to_frame()
    reaction_ids.drop_duplicates(inplace = True)
    # Add retro rules IDs
    with sqlite3.connect(path_to_database) as connection:
        query = "SELECT reaction_ID, retro_rules_ID FROM Retro_rules_reactions"
        reference_df = pd.read_sql_query(query, connection)
    reaction_ids = reaction_ids.merge(
                    reference_df,
                    how = "left", 
                    on = "reaction_ID"
                    )

    # Lookup matching rhea ID
    with sqlite3.connect(path_to_database) as connection:
        query = "SELECT * FROM Retro_rules_to_rhea_reactions"
        reference_df = pd.read_sql_query(query, connection)
    reaction_ids = reaction_ids.merge(
                    reference_df,
                    how = "left", 
                    on = "retro_rules_ID"
                    )
    no_match_ids = reaction_ids[pd.isna(reaction_ids["rhea_ID"])]["retro_rules_ID"]

    # Write missing IDs to file
    if len(no_match_ids) > 0:
        print("WARNING: No Rhea ID found for some retro rules reactions. See log_unmatched_reactions.txt for details.")
        with open(os.path.join(output_dir, "log_unmatched_reactions.txt"), "w") as logfile:
            for non_matching_id in no_match_ids:
                logfile.write("{}\n".format(non_matching_id))

    # Remove duplicates and return Rhea IDs.
    reaction_ids = reaction_ids[pd.notna(reaction_ids["rhea_ID"])]["rhea_ID"].drop_duplicates().apply(int).to_frame()
    return reaction_ids
    
def load_refseq_sequences_dataframe(refseq_path = my_settings.refseq_path):
    """
    Load table with reference sequences for which a reaction is known.
    It returns a pandas dataframe with two columns: a uniprot_ID and its corresponding sequence.
    """
    print( '!!! loading refseq and reaction tables...')
    # this table is only refseqs that are found in mrs-reaction
    print( '!!! Reference sequences in this file: {}'.format(refseq_path))
    refseq = mg.load_dataframe(refseq_path)
    refseq = refseq[refseq['sequence'].notnull()]
    print( '!!! {} reference sequences'.format(len(refseq)))
    refseq.set_index("uniprot_ID", inplace = True)
    return refseq
    
def get_uniprot_ids_for_rhea_reactions(rhea_reactions, path_to_database = my_settings.magi_database):
    """
    Get all protein identifiers that match reactions from compound to reactions
    """
    with sqlite3.connect(path_to_database) as connection:
        query = "SELECT rhea_ID, uniprot_ID FROM Proteins"
        rhea_to_uniprot_ids = pd.read_sql_query(query, connection)

    rhea_reactions = rhea_reactions.merge(
        rhea_to_uniprot_ids, 
        how = "left", 
        on = "rhea_ID")
    uniprot_ids = rhea_reactions["uniprot_ID"].drop_duplicates().dropna()
    return uniprot_ids

def workflow(compound_to_reaction_path, genome_db_path, blast_filter, output_dir, intermediate_files_dir, cpu_count):
    # reaction to gene search
    def keep_top_blast_helper(x, param=blast_filter):
        """
        # setup multiprocessing helpers based on input params
        x is the normal input
        param is the defined parameter
        """
        return blast.keep_top_blast(x, filt=param)
    
    # Load reference sequences
    refseq_sequences = load_refseq_sequences_dataframe()

    # Perform reaction to gene search
    print( '\n!@# Conducting reaction to gene search | TLOG %s' % (time.time()))
    sys.stdout.flush()
    start = time.time()

    ## Get all uniprot identifiers that match any of the proteins that we have
    rhea_reactions = prepare_compound_to_reaction_results(compound_to_reaction_path, output_dir)
    rhea_proteins = get_uniprot_ids_for_rhea_reactions(rhea_reactions)
    
    # rhea_proteins is the "query_list" for multi_blast()
    # refseq_sequences is the full refseq table #filter?
    # database_path is the path to the user's input genome's blast database
    print( '!!! {} reference sequences to search'.format(len(rhea_proteins)))
    sys.stdout.flush()

    # BLAST
    reaction_to_gene = blast.multi_blast(rhea_proteins, refseq_sequences, 
        genome_db_path, intermediate_files_dir, cpu=cpu_count, 
        raise_blast_error=False)

    # Add reaction info
    reaction_to_gene = blast.refseq_to_reactions(reaction_to_gene,
        'query acc.')

    reaction_to_gene.to_csv(os.path.join(intermediate_files_dir,
        'reaction_blast.csv'), index = False)
    print( '!!! r2g blast results saved to %s' \
            %(os.path.join(intermediate_files_dir, 'reaction_blast.csv')))

    reaction_groups = reaction_to_gene.groupby('query acc.')
    multidx = reaction_groups['e_score'].apply(keep_top_blast_helper).index
    del reaction_groups
    idx = multidx.levels[1]
    reaction_to_gene_top = reaction_to_gene.loc[idx]
    del reaction_to_gene
    reaction_to_gene_top_path = os.path.join(intermediate_files_dir, 
                                            'reaction_to_gene.csv')
    reaction_to_gene_top.to_csv(reaction_to_gene_top_path, index = False)
    print( '!@# reaction_to_gene table done in %s minutes'\
            %((time.time()-start)/60))
    print( '!!! reaction_to_gene table saved to %s'\
            % reaction_to_gene_top_path)
    
    return reaction_to_gene_top_path

def main():
    # Parse arguments and prepare for reaction to gene workflow
    magi_parameters = mg.general_magi_preparation()
    
    if magi_parameters["gene_to_reaction_only"] or magi_parameters["compound_to_reaction_only"]:
        print("Not performing MAGI reaction to gene workflow")
    else:
        # Set compound to reaction file path
        if magi_parameters["compound_to_reaction"] is not None:
            compound_to_reaction_path = magi_parameters["compound_to_reaction"]
        else:
            compound_to_reaction_path = mg.get_intermediate_file_path(magi_parameters["output_dir"], "compound_to_reaction_path")
        # Set genome database path
        if magi_parameters["genome_db"] is not None:
            genome_db_path = magi_parameters["genome_db"]
        else:
            genome_db_path = mg.get_intermediate_file_path(magi_parameters["output_dir"], "genome_db_path")

        # Start r2g workflow
        reaction_to_gene_path = workflow(compound_to_reaction_path = compound_to_reaction_path, 
            genome_db_path = genome_db_path, 
            blast_filter = magi_parameters["blast_filter"], 
            output_dir=magi_parameters["output_dir"],
            intermediate_files_dir = magi_parameters["intermediate_files_dir"], 
            cpu_count = magi_parameters["cpu_count"])
        mg.write_intermediate_file_path(magi_parameters["output_dir"], "reaction_to_gene_path", reaction_to_gene_path)
        print("Reaction to gene search finished")

if __name__ == "__main__":
    main()
