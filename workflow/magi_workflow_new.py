"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.0b workflow script

Lines beginning with "@@@" are input job parameters
Lines beginning with "!!!" are verbose log info
Lines beginning with "!@#" are checkpoints/announcements

Required inputs are a FASTA proteins file and a Compounds file.
The FASTA file should be in standard FASTA format, E.g.:
'''fasta
> UNIQUE_GENE_ID OTHER_INFORMATION
AMINO_ACID_SEQUENCE

> UNIQUE_GENE_ID OTHER_INFORMATION
AMINO_ACID_SEQUENCE
'''
*Please note the space in between UNIQUE_GENE_ID and OTHER_INFORMATION*

The FASTA file will be used to make a BLAST database, and to make a
gene table with these columns:
Gene_ID        | header                           | sequence
-----------------------------------------------------------------------
UNIQUE_GENE_ID | UNIQUE_GENE_ID OTHER_INFORMATION | AMINO_ACID_SEQUENCE

The Compounds file should be in some standard csv table format (CSV, tab-
delimited, etc.), and is required to have a column named 
"original_compound" (case-sensitive). In this column should be standard
InChI Keys representing a compound structure, or m/z values if the accurate 
mass search module will be used.

The Compounds table may also have a column named "compound_score" 
(case-sensitive), where the user can provide a score for each 
compound-row. If this column name does not exist, one will be created
and populated with 1.0.
"""
# Import modules
import sys
import os
import argparse
from multiprocessing import cpu_count as counting_cpus
import pandas as pd
import numpy as np
import time
import pickle
import datetime
import workflow_helpers_new as mg
import magi_workflow_gene_to_reaction


    
def perform_accurate_mass_search(compounds_file, adduct_file, polarity, accurate_mass_search_only, ppm_cutoff):
    import magi_workflow_accurate_mass_search
    mass_searched_compounds_filename = magi_workflow_accurate_mass_search.workflow(compounds_file, adduct_file, polarity, accurate_mass_search_only, ppm_cutoff)
    return mass_searched_compounds_filename

def compound_to_reaction_search(legacy, level, compound_to_reaction, compounds_file, main_start, cpu_count, fasta_file, output_dir, intermediate_files_dir, pactolus): 
    import magi_workflow_compound_to_reaction
    compounds = mg.load_compound_results(compounds_file, pactolus, output_dir)
    if compound_to_reaction is None:
        compound_to_reaction_path = magi_workflow_compound_to_reaction.workflow(compounds_to_search = compounds, tautomer_legacy = legacy, 
                                                                                neighbor_level = level, cpu_count = cpu_count, intermediate_files_dir = intermediate_files_dir)
        compound_to_reaction = pd.read_pickle(compound_to_reaction_path)
    else:
        compound_to_reaction = pd.read_pickle(compound_to_reaction)
        print( '\n!@# compound_to_reaction successfully loaded')
    return compounds, compound_to_reaction

def reaction_to_gene_search(compound_to_reaction, genome_db_path, blast_filter, reaction_to_gene, intermediate_files_dir, cpu_count):
    import magi_workflow_reaction_to_gene
    
    if reaction_to_gene is None:
        reaction_to_gene_top_path = magi_workflow_reaction_to_gene.workflow(compound_to_reaction, genome_db_path, blast_filter, intermediate_files_dir, cpu_count)
        reaction_to_gene_top = pd.read_pickle(reaction_to_gene_top_path)
    else: 
        reaction_to_gene_top = pd.read_pickle(reaction_to_gene)
        print( '\n!@# reaction_to_gene successfully loaded')
    return reaction_to_gene_top



def main():
    print_version_info()
    args = parse_arguments()
    args = print_parameters(args)
    
    #Set all MAGI parameters
    if args.accurate_mass_search is not None:
        accurate_mass_search = True
        polarity = args.accurate_mass_search
    else:
        accurate_mass_search = False
    accurate_mass_search_only = args.accurate_mass_search_only
    adduct_file = args.adduct_file
    annotations = args.annotations
    blast_filter = args.blast_filter
    chemnet_penalty = args.chemnet_penalty
    compound_to_reaction = args.compound_to_reaction
    compounds_file = args.compounds
    cpu_count = args.cpu_count
    fasta_file = args.fasta
    final_weights = args.final_weights
    gene_to_reaction = args.gene_to_reaction
    intermediate_files = args.intermediate_files
    legacy = args.legacy
    level = args.level
    merged_before_score = args.merged_before_score
    output_dir = args.output
    pactolus = args.pactolus
    ppm_cutoff = args.ppm_cutoff
    reaction_to_gene = args.reaction_to_gene
    reciprocal_closeness = args.reciprocal_closeness
    
    # start with magi workflow
    output_dir, intermediate_files_dir = mg.make_output_dirs(output_dir, fasta_file, compounds_file, intermediate_files) #TODO: rename this function
    main_start = time.time() # overall program timer
    
    if accurate_mass_search:
        compounds_file = perform_accurate_mass_search(compounds_file, adduct_file, polarity, accurate_mass_search_only, ppm_cutoff)
    # Gene to reaction search
    if fasta_file is not None:
        if gene_to_reaction is None:
            gene_to_reaction_path, genome_db_path = magi_workflow_gene_to_reaction.workflow(fasta_file, intermediate_files_dir, cpu_count, annotations, blast_filter)
            gene_to_reaction_top = pd.read_pickle(gene_to_reaction_path)
        else:
            gene_to_reaction_top = pd.read_pickle(gene_to_reaction)
            print( '\n!@# gene_to_reaction successfully loaded')
            
    if compounds_file is not None:
        compounds, compound_to_reaction = compound_to_reaction_search(legacy, level, compound_to_reaction, compounds_file, main_start, cpu_count, fasta_file, output_dir, intermediate_files_dir, pactolus)
        reaction_to_gene_top = reaction_to_gene_search(compound_to_reaction, genome_db_path, blast_filter, reaction_to_gene, intermediate_files_dir, cpu_count)
    merged_dataframe = merging_g2r_and_r2g_searches(compound_to_reaction, reaction_to_gene_top, gene_to_reaction_top, merged_before_score, intermediate_files_dir)
    final_dataframe, start = calculate_scores_and_format_tables(merged_dataframe, reciprocal_closeness, final_weights, chemnet_penalty)
    save_outputs(final_dataframe, compounds, main_start, start, output_dir)
    
if __name__ == "__main__":
    main()