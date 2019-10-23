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
import warnings
from multiprocessing import cpu_count as counting_cpus
import pandas as pd
import numpy as np
import time
import pickle
import datetime
import workflow_helpers_new as mg
import magi_workflow_gene_to_reaction

def print_version_info():
    """Print versions of modules that may be troublesome for magi."""
    print('!!! Python version:'+ sys.version)
    print('!!! numpy version: '+ np.__version__)
    print('!!! pandas version:'+ pd.__version__)
    #print('!!! pickle version:'+ pickle.__version__)
    print('#'*80)

def parse_arguments():
    def is_existing_file(filepath):
        """Checks if a file exists and return absolute path if it exists"""
        if not os.path.exists(filepath):
            msg = "{0} does not exist".format(filepath)
            raise argparse.ArgumentTypeError(msg)
        else:
            return os.path.abspath(filepath)
    def percentage_values_to_decimal(percentage):
        """Turns the blast filter and reciprocal closeness percentages 
        into decimal numbers"""
        try:
            percentage = int(percentage)
        except:
            msg = "Please enter an integer value"
            raise argparse.ArgumentTypeError(msg)        
        if percentage > 100:
            msg = "Max value is 100"
            raise argparse.ArgumentTypeError(msg)
        elif percentage < 0:
            msg = "Value cannot be negative"
            raise argparse.ArgumentTypeError(msg)
        else:
            decimal = percentage/100.
        return decimal
    
    def positive_number(number):
        """Checks if none of the number/numbers are negative"""
        try:
            number = float(number)
        except:
            msg = "Please enter a numeric value"
            raise argparse.ArgumentTypeError(msg)        
        if number < 0:
            msg = "Value cannot be negative"
            raise argparse.ArgumentTypeError(msg)
        else:
            return number
    
    def set_cpu_count(cpu_count):
        max_cpu = counting_cpus()  
        if cpu_count == 0:
            cpu_count = max_cpu    
        if cpu_count > max_cpu:
            msg = "ERROR: You have exceeded the cpus on this machine ({})".format(max_cpu)
            raise argparse.ArgumentTypeError(msg)
        return cpu_count
    try:
        """parse arguments"""
        parser = argparse.ArgumentParser()
        # required arguments
        required_args = parser.add_argument_group('Required arguments')
        required_args.add_argument('-f', '--fasta', type=is_existing_file,
            help='path to fasta file of genes in sample')
        required_args.add_argument('-c', '--compounds', type=is_existing_file,
            help='path to observed compounds file')
        
        # jump-start the script after certain computations
        start_halfway_args = parser.add_argument_group('Arguments to jump-start the script after certain computations')
        start_halfway_args.add_argument('--gene_to_reaction', type=is_existing_file,
            help='path to gene_to_reaction file, must be in pickle format')
        start_halfway_args.add_argument('--compound_to_reaction', type=is_existing_file,
            help='path to compound_to_reaction file, must be in pickle format')
        start_halfway_args.add_argument('--reaction_to_gene', type=is_existing_file,
            help='path to reaction_to_gene file, must be in pickle format')
        start_halfway_args.add_argument('--merged_before_score', type=is_existing_file,
            help='path to merged_before_score table, must be in hdf5 format,\
            with the key "merged_before_score"')
        
        # optional runtime variables
        optional_args = parser.add_argument_group("Optional runtime variables")
        optional_args.add_argument('-a', '--annotations', type=is_existing_file,
            help='path to annotation file for genes in sample', 
            default=None)
        optional_args.add_argument('-n', '--cpu_count', 
            help='number of cpus to use for multiprocessing. Default is to use max!', 
            type=int, default=0)
        optional_args.add_argument('-o', '--output', 
            help='path to a custom output', 
            type=str)
        optional_args.add_argument('-l', '--level', 
            help='how many levels deep to search the chemical network', 
            type=int, choices=[0,1,2,3], default=2)
        optional_args.add_argument('--legacy', dest='legacy', action='store_true',
            help='use legacy tautomer searching; default is no')
        optional_args.add_argument('--no-legacy', dest='legacy', action='store_false',
            help='use precomputed compound-to-reaction; default is yes')
        optional_args.set_defaults(legacy=False)
        optional_args.add_argument('--mute', 
            help='mutes pandas warnings', 
            action='store_true')
        optional_args.add_argument('--pactolus', 
            help='Flag to tell MAGI that the compounds input is a pactolus file', 
            action='store_true')
        optional_args.add_argument('--test', 
            help='TBD: run MAGI only on the first # of pactolus compounds', 
            type=int)
        optional_args.add_argument('--debug', 
            help='TBD: prints a lot of info', 
            action='store_true')
        optional_args.add_argument('--blast_filter', 
            help='How stringent to filter the top BLAST results, as percent;\
            default is 85 meaning that only BLAST results within 85%% of the top\
            result will be taken.', 
            type=percentage_values_to_decimal, default=0.85)
        optional_args.add_argument('--reciprocal_closeness', 
            help='Cutoff to call a reciprocal disagreement as "close", as percent;\
            default is 75 meaning that a reciprocal disagreement will be classified\
            as "close" if the lower blast score (e score) is within 75%% of the higher\
            score', 
            type=percentage_values_to_decimal, default=0.75)
        optional_args.add_argument('--final_weights', 
            help='Defined weights to weight the final scoring for the scores:\
            compound_score reciprocal_score homology_score reaction_connection', 
            type=positive_number, nargs=4, default=None)
        optional_args.add_argument('--chemnet_penalty', 
            help='Base factor in the chemical network search level penalty', 
            type=positive_number, default=4)
        optional_args.add_argument('--intermediate_files',
            help='What directory within --output to store intermediate files',
            type=str, default='intermediate_files')
        
        # Parameters for accurate mass search
        mass_search_args = parser.add_argument_group("Arguments for the optional accurate mass search")
        mass_search_args.add_argument('--accurate_mass_search',
                            type=str, choices=['pos','neg','neut'], default=None,
                            help = "Perform accurate mass search on m/z values in original_compound column in the input compounds file. \
                            Specify if the masses are measured in negative mode, positive mode or if they have been transformed to neutral masses."
                            )
        mass_search_args.add_argument('--accurate_mass_search_only',
                            action='store_true', default=False,
                            help = "If this parameter is used, MAGI will stop after performing the accurate mass search."
                            )
        mass_search_args.add_argument('--adduct_file',
                            type=str, default=None,
                            help="optionally specify which adducts to investigate. If not specified, it will search for M+,M+H,M+NH4,M+Na in positive mode or M-H,M+Cl,M+FA-H,M+Hac-H in negative mode."
                            )
        mass_search_args.add_argument('--ppm_cutoff', 
            help='The ppm cutoff for the accurate mass search. Default is 10 ppm.', 
            type=int, default=10)
        
        args = parser.parse_args()
        
        # Check parameters and set number of required CPUs
        if args.fasta is None and args.compounds is None:
            raise argparse.ArgumentTypeError('ERROR: either FASTA or metabolites file is required')
        args.cpu_count = set_cpu_count(args.cpu_count)
    except argparse.ArgumentTypeError as ex:
        print(ex.message)
        sys.exit(1)
    return args

def print_parameters(args):
    print('~~~~~PARAMETERS~~~~~~')

    # print your paths in stdout for logging
    print('@@@ FASTA file input: %s' %(args.fasta))
    print('@@@ Compound input: %s' %(args.compounds))
    if args.annotations is not None:
        print( '@@@ annotations input: %s' %(args.annotations))
    # print parameters
    if args.fasta is not None: 
        print( '@@@ BLAST filter: %s' % (args.blast_filter))
    if args.compounds is not None:
        print( '@@@ Using precomputed compound results: %s' % (not args.legacy))
        print( '@@@ Chemnet search to level %s' % (args.level))
        print( '@@@ Reciprocal closeness: %s' % (args.reciprocal_closeness))
        print( '@@@ Chemnet base penalty: %s' % (args.chemnet_penalty))
    print( '@@@ MAGI score weights: %s' % (args.final_weights))
    print( '@@@ Using %s CPUs' % (args.cpu_count))
    
    if args.gene_to_reaction is not None:
        print( '@@@ gene_to_reaction input: %s' %(args.gene_to_reaction))
    
    if args.compound_to_reaction is not None:
        print( '@@@ compound_to_reaction input: %s' %(args.compound_to_reaction))
    
    if args.reaction_to_gene is not None:
        print( '@@@ reaction_to_gene input: %s' %(args.reaction_to_gene))
    if args.mute:
        print( '!!! Warnings are muted')
        warnings.filterwarnings('ignore')
    return args
    
def perform_accurate_mass_search(compounds_file, adduct_file, polarity, accurate_mass_search_only, ppm_cutoff):
    import magi_workflow_accurate_mass_search
    mass_searched_compounds_filename = magi_workflow_accurate_mass_search.workflow(compounds_file, adduct_file, polarity, accurate_mass_search_only, ppm_cutoff)
    return mass_searched_compounds_filename

def compound_to_reaction_search(legacy, level, compound_to_reaction, compounds_file, main_start, cpu_count, fasta_file, output_dir, intermediate_files_dir, pactolus): 
    import magi_workflow_compound_to_reaction
    compounds = mg.load_compound_results(compounds_file, pactolus, output_dir)
    if compound_to_reaction is None:
        compound_to_reaction_path = magi_workflow_compound_to_reaction.workflow(legacy, level, compounds, main_start, cpu_count, output_dir, intermediate_files_dir, pactolus)
        compound_to_reaction = pd.read_pickle(compound_to_reaction_path)
    else:
        compound_to_reaction = pd.read_pickle(compound_to_reaction)
        print( '\n!@# compound_to_reaction successfully loaded')
    return compounds, compound_to_reaction

def reaction_to_gene_search(compound_to_reaction, genome_db_path, blast_filter, reaction_to_gene, intermediate_files_dir, cpu_count):
    import magi_workflow_reaction_to_gene_search
    reaction_to_gene_top = magi_workflow_reaction_to_gene_search.workflow(compound_to_reaction, genome_db_path, blast_filter, reaction_to_gene, intermediate_files_dir, cpu_count)
    return reaction_to_gene_top

def merging_g2r_and_r2g_searches(compound_to_reaction, reaction_to_gene_top, gene_to_reaction_top, merged_before_score, intermediate_files_dir):
    if merged_before_score is None:
    	print( '\n!@# Merging final table | TLOG %s' % (time.time()))
    	sys.stdout.flush()
    	start = time.time()
    
    	compound_to_gene = pd.merge(compound_to_reaction, reaction_to_gene_top, 
    								on='reaction_id', how='left')
    	del reaction_to_gene_top
    	del compound_to_reaction
    
    	compound_to_gene_small = compound_to_gene[['subject acc.', 'reaction_id',
    								'e_score', 'compound_score',
    								'original_compound', 'level', 'neighbor',
    								'note']]
    	del compound_to_gene
    
    	# okay to drop duplicates, because i only care about these columns 
    	# anyway; if these are duplicated then other information doesn't really 
    	# matter or can easily be re-expanded by joining 
    	compound_to_gene_small.drop_duplicates(inplace=True)
    
    	gene_to_reaction_small = gene_to_reaction_top[['query acc.', 'reaction_id',
    													'e_score']]
    	del gene_to_reaction_top
    	gene_to_reaction_small.drop_duplicates(inplace=True)
    
    	# Make an integrated dataframe, joining on the gene
    	df = pd.merge(compound_to_gene_small, gene_to_reaction_small, 
    		left_on='subject acc.', right_on='query acc.', 
    		suffixes=('_r2g', '_g2r'), how='outer')
    
    	df.reset_index(inplace=True, drop=True)
    	df.drop_duplicates(inplace=True)
    
    	# Clean up reaction_id_r2g column
    	idx = df[df['reaction_id_r2g'] == ''].index
    	df.loc[idx, 'reaction_id_r2g'] = np.nan
    	df['reaction_id_r2g'] = df['reaction_id_r2g'].astype(float)
    
    	# Clean up stupid NaNs in string columns
    	def check_str(x):
    	    if isinstance(x, str):
    	        return True
    	    else:
    	        return False
    
    	for c in df.columns:
    	    if len(df[c].apply(type).unique()) > 1:
    	        string_checked = df[c].apply(check_str)
    	        if string_checked.any():
    	            df[c].fillna('', inplace=True)
    
    	# Clean up neighbor column
    	df['neighbor'] = df['neighbor'].astype(str)
    
    	df.to_hdf(os.path.join(intermediate_files_dir, 'merged_before_score.h5'),
    		'merged_before_score', mode='w', format='table',
    		complib='blosc', complevel=9)
    
    	print( '!@# Final Merged table done in %s minutes'\
    		%((time.time() - start) / 60))
    	print( '!!! Final Merged table saved to %s'\
    			% (os.path.join(intermediate_files_dir, 'merged_before_score.h5')))
    else:
    	del reaction_to_gene_top
    	del compound_to_reaction
    	del gene_to_reaction_top
    	df = pd.read_hdf(merged_before_score, 'merged_before_score')
    	print( '\n!@# merged_before_score successfully loaded')
    return df

def calculate_scores_and_format_tables(df, reciprocal_closeness, final_weights, chemnet_penalty):
    print( '\n!@# Calculating final scores | TLOG %s' % (time.time()))
    start = time.time()
    sys.stdout.flush()
    
    # score reciprocal agreement
    df = mg.reciprocal_agreement(df, closeness_threshold=reciprocal_closeness)
    
    # calculate homology score
    score = mg.homology_score(df)
    # the nulls get a really low score
    score[pd.isnull(score)] = 1
    df['homology_score'] = score
    
    # reaction connection score says if the compound got connected to any
    # reaction in the database. Can't have zero because that messes up
    # geometric mean, so added a small number.
    df['reaction_connection'] = df[['reaction_id_r2g', 'reaction_id_g2r']]\
    								.apply(pd.notnull).sum(axis=1) + 0.01
    
    # calculate final MAGI integrated score
    scoring_data = ['compound_score', 'reciprocal_score', \
    				'homology_score', 'reaction_connection']
    scores = []
    to_score = df[scoring_data].values
    if final_weights is not None:
    	weights = np.asarray([final_weights] * to_score.shape[0])
    	data = mg.magi_score(to_score, weights)
    else:
    	data = mg.magi_score(to_score)
    scores.append(data)
    df['MAGI_score'] = scores[0] / (chemnet_penalty ** df['level'].values)
    
    # find gene ids that are floats, convert those to strings, without the decimal
    float_entries = df['subject acc.'].apply(lambda x: isinstance(x, float))
    df.loc[float_entries, 'subject acc.'] = df.loc[float_entries, \
    				'subject acc.'].apply(lambda x: "{:.0f}".format(x))
    
    print( '!@# Scoring done in %s minutes' %((time.time() - start) / 60))
    print( '\n!@# Formatting final table | TLOG %s' % (time.time()))
    start = time.time()
    sys.stdout.flush()
    
    # sort the final table and drop key duplicates
    df = df.sort_values(
    	['original_compound', 'MAGI_score'], 
    	ascending=[True, False]
    	).drop_duplicates(
    		['original_compound', 'level', 'neighbor', 'compound_score',
    		 'reciprocal_score', 'query acc.', 'reaction_id_r2g',
    		 'reaction_id_g2r']
    		 )
    mrs_reaction_path = mg.get_settings(); mrs_reaction_path = mrs_reaction_path.mrs_reaction_path
    mrs_reaction = mg.load_dataframe(mrs_reaction_path)
    df = df.merge(mrs_reaction[['database_id']],
    	left_on='reaction_id_r2g', right_index=True, how='left')
    df = df.merge(mrs_reaction[['database_id']],
    	left_on='reaction_id_g2r', right_index=True, how='left',
    	suffixes=('_r2g', '_g2r'))
    cols = df.columns.values
    idx = pd.np.argwhere(cols == 'query acc.')[0][0]
    cols[idx] = 'gene_id'
    df.columns = cols
    df = df[['MAGI_score','gene_id', 'original_compound', 'neighbor',
    	'note', 'compound_score','level','homology_score','reciprocal_score',
    	'reaction_connection', 'e_score_r2g','database_id_r2g', 'e_score_g2r',
    	'database_id_g2r']]
    return df, start

def save_outputs(df, compounds, main_start, start, output_dir):
    # save the full dataframe
    df.to_csv(os.path.join(output_dir, 'magi_results.csv'), index=False)
    print( 'full results saved to {}'.format(os.path.join(output_dir, 'magi_results.csv')))
    # save a compound-centric dataframe, where only the best row for each
    # original_compound was chosen (this is only for compound scoring, do
    # not use this for any kind of gene function analysis!)
    
    compound_centric = df[pd.notnull(df['original_compound'])]\
    					 .sort_values('MAGI_score', ascending=False)\
    					 .drop_duplicates(['original_compound', 'compound_score'])
    compound_centric = pd.merge(
    	compound_centric, compounds,
    	on=['original_compound', 'compound_score'],
    	how='right')
    compound_centric.to_csv(os.path.join(output_dir, 
    	'magi_compound_results.csv'), index=False)
    
    gene_centric = df.sort_values(['MAGI_score', 'e_score_g2r'], 
    	ascending=[False, False])\
    	.drop_duplicates(['gene_id', 'database_id_g2r'])
    gene_centric.to_csv(os.path.join(output_dir,
    	'magi_gene_results.csv'), index=False)
    
    print( '!@# MAGI Scoring done in %s minutes' %((time.time() - start) / 60))
    print( '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60))
    print( '!!! final results stored to %s' \
    		%(os.path.join(output_dir, 'magi_results.csv')))

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