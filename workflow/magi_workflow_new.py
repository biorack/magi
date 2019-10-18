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
import multiprocessing as mp
import pandas as pd
import numpy as np
import time
import pickle
import datetime
import workflow_helpers as mg

def print_version_info():
    """Print versions of modules that may be troublesome for magi."""
    print('!!! Python version:'+ sys.version)
    print('!!! numpy version: '+ np.__version__)
    print('!!! pandas version:'+ pd.__version__)
    #print('!!! pickle version:'+ pickle.__version__)
    print('#'*80)

def parse_arguments():
    """parse arguments"""
    parser = argparse.ArgumentParser()
    # required
    parser.add_argument('-f', '--fasta', 
        help='path to fasta file of genes in sample')
    parser.add_argument('-c', '--compounds', 
        help='path to observed compounds file')
    
    # jump-start the script after certain computations
    parser.add_argument('--gene_to_reaction', 
        help='path to gene_to_reaction file, must be in pickle format')
    parser.add_argument('--compound_to_reaction', 
        help='path to compound_to_reaction file, must be in pickle format')
    parser.add_argument('--reaction_to_gene', 
        help='path to reaction_to_gene file, must be in pickle format')
    parser.add_argument('--merged_before_score', 
        help='path to merged_before_score table, must be in hdf5 format,\
        with the key "merged_before_score"')
    
    # optional runtime variables
    parser.add_argument('-a', '--annotations', 
        help='path to annotation file for genes in sample', 
        default=None)
    parser.add_argument('-n', '--cpu_count', 
        help='number of cpus to use for multiprocessing. Default is to use max!', 
        type=int, default=0)
    parser.add_argument('-o', '--output', 
        help='path to a custom output', 
        type=str)
    parser.add_argument('-l', '--level', 
        help='how many levels deep to search the chemical network', 
        type=int, choices=[0,1,2,3], default=2)
    parser.add_argument('--legacy', dest='legacy', action='store_true',
        help='use legacy tautomer searching; default is no')
    parser.add_argument('--no-legacy', dest='legacy', action='store_false',
        help='use precomputed compound-to-reaction; default is yes')
    parser.set_defaults(legacy=False)
    parser.add_argument('--mute', 
        help='mutes pandas warnings', 
        action='store_true')
    parser.add_argument('--pactolus', 
        help='Flag to tell MAGI that the compounds input is a pactolus file', 
        action='store_true')
    parser.add_argument('--test', 
        help='TBD: run MAGI only on the first # of pactolus compounds', 
        type=int)
    parser.add_argument('--debug', 
        help='TBD: prints a lot of info', 
        action='store_true')
    parser.add_argument('--blast_filter', 
        help='How stringent to filter the top BLAST results, as percent;\
        default is 85 meaning that only BLAST results within 85%% of the top\
        result will be taken.', 
        type=int, choices=range(0, 101), default=85)
    parser.add_argument('--reciprocal_closeness', 
        help='Cutoff to call a reciprocal disagreement as "close", as percent;\
        default is 75 meaning that a reciprocal disagreement will be classified\
        as "close" if the lower blast score (e score) is within 75%% of the higher\
        score', 
        type=int, choices=range(0, 101), default=75)
    parser.add_argument('--final_weights', 
        help='Defined weights to weight the final scoring for the scores:\
        compound_score reciprocal_score homology_score reaction_connection', 
        type=float, nargs=4, default=None)
    parser.add_argument('--chemnet_penalty', 
        help='Base factor in the chemical network search level penalty', 
        type=float, default=4)
    parser.add_argument('--intermediate_files',
        help='What directory within --output to store intermediate files',
        type=str, default='intermediate_files')
    
    # Parameters for accurate mass search
    parser.add_argument('--accurate_mass_search',
                        type=str, choices=['pos','neg','neut'], default=None,
                        help = "Perform accurate mass search on m/z values in original_compound column in the input compounds file. \
                        Specify if the masses are measured in negative mode, positive mode or if they have been transformed to neutral masses."
                        )
    parser.add_argument('--accurate_mass_search_only',
                        action='store_true', default=False,
                        help = "If this parameter is used, MAGI will stop after performing the accurate mass search."
                        )
    parser.add_argument('--adduct_file',
                        type=str, default=None,
                        help="optionally specify which adducts to investigate. If not specified, it will search for M+,M+H,M+NH4,M+Na in positive mode or M-H,M+Cl,M+FA-H,M+Hac-H in negative mode."
                        )
    parser.add_argument('--ppm_cutoff', 
        help='The ppm cutoff for the accurate mass search. Default is 10 ppm.', 
        type=int, default=10)
    
    args = parser.parse_args()
    return args

def check_parameters(args):
    """Check if the input parameters are set correctly for MAGI"""
    if args.fasta is None and args.compounds is None:
        print('ERROR: either FASTA or metabolites file is required')
        sys.exit()
    if args.fasta is not None:
        if not os.path.isfile(args.fasta):
            raise IOError('%s does not exist!' %(args.fasta))
        else:
            args.fasta = os.path.abspath(args.fasta)
    if args.compounds is not None:
        if not os.path.isfile(args.compounds):
            raise IOError('%s does not exist!' %(args.compounds))
        else:
            args.compounds = os.path.abspath(args.compounds)    
    print('~~~~~PARAMETERS~~~~~~')

    # print your paths in stdout for logging
    print('@@@ FASTA file input: %s' %(args.fasta))
    print('@@@ Compound input: %s' %(args.compounds))
    if args.annotations is not None:
        args.annotations =  os.path.abspath(args.annotations)
        print( '@@@ annotations input: %s' %(args.annotations))
    # convert parameters
    args.blast_filter = args.blast_filter / 100.
    args.reciprocal_closeness = args.reciprocal_closeness / 100.
    
    # check parameters
    if args.blast_filter < 0:
        raise RuntimeError('argument blast_filter cannot be negative!')
    if args.reciprocal_closeness < 0:
        raise RuntimeError('argument reciprocal_closeness cannot be negative!')
    if args.final_weights is not None:
        for w in args.final_weights:
            if w < 0:
                raise RuntimeError('argument final_weights cannot be negative!')
    if args.chemnet_penalty < 0:
        raise RuntimeError('argument chemnet_penalty cannot be negative!')
    
    max_cpu = mp.cpu_count()
    
    if args.cpu_count > max_cpu:
        raise RuntimeError('You have exceeded the cpus on this machine (%s)!'
            %(max_cpu))
    if args.cpu_count == 0:
        args.cpu_count = max_cpu
    
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
        args.gene_to_reaction =  os.path.abspath(args.gene_to_reaction)
        if not os.path.isfile(args.gene_to_reaction):
            raise IOError('%s does not exist!' %(args.gene_to_reaction))
        else:
            print( '@@@ gene_to_reaction input: %s' %(args.gene_to_reaction))
    
    if args.compound_to_reaction is not None:
        args.compound_to_reaction =  os.path.abspath(args.compound_to_reaction)
        if not os.path.isfile(args.compound_to_reaction):
            raise IOError('%s does not exist!' %(args.compound_to_reaction))
        else:
            print( '@@@ compound_to_reaction input: %s' %(args.compound_to_reaction))
    
    if args.reaction_to_gene is not None:
        args.reaction_to_gene =  os.path.abspath(args.reaction_to_gene)
        if not os.path.isfile(args.reaction_to_gene):
            raise IOError('%s does not exist!' %(args.reaction_to_gene))
        else:
            print( '@@@ reaction_to_gene input: %s' %(args.reaction_to_gene))
    
    if args.mute:
        print( '!!! Warnings are muted')
        warnings.filterwarnings('ignore')
    return args
    
def make_output_dir_name(output_dir, fasta_file, compounds_file, intermediate_files):
    """set up where the results will be stored"""
    if output_dir is None:
        # autoname the directory based on fasta, or compound file
        # this will change eventually
        if fasta_file is not None:
            experiment_name = os.path.splitext(os.path.basename(fasta_file))[0]
            experiment_dir = os.path.abspath(os.path.dirname(fasta_file))
        else:
            experiment_name = os.path.splitext(os.path.basename(compounds_file))[0]
            experiment_dir = os.path.abspath(os.path.dirname(compounds_file))
        today = datetime.datetime.now()
        experiment_name += today.strftime('_%Y%m%d')
        experiment_path = os.path.join(experiment_dir, experiment_name)
    else:
        experiment_path = output_dir
    
    output_dir = os.path.abspath(experiment_path)
    intermediate_files_dir = os.path.join(output_dir, intermediate_files)
    
    print( '!!! Saving all results here:', output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    return output_dir, intermediate_files_dir
    
def load_fasta_genome(fasta_filename, intfile_path, annotation_file = None):
    """load genome"""
    print( '\n!!! LOADING GENOME')
    genome, genome_db_path = mg.load_genome(fasta_filename, intfile_path, 
                                        annotation_file)
    return genome, genome_db_path
    
def perform_accurate_mass_search(compounds_file, adduct_file, polarity, accurate_mass_search_only, ppm_cutoff):
        if compounds_file is None:
            raise RuntimeError("No compounds file specified. Exiting...")
        else:
            # Perform accurate mass search and set compounds file to mass-searched file.
            print("\n!!! Performing accurate mass search for {}".format(compounds_file))
            #Make list of adducts to search for
            if adduct_file is not None:
                adduct_file = os.path.abspath(adduct_file)
                print('@@@ Adduct file input: %s' %(adduct_file))
                with open(adduct_file) as adduct_file:
                    adducts = []
                    try:
                        for line in adduct_file:
                            adducts.append(line.rstrip())
                    except:
                        print("File cannot be converted to adducts list. Please specify one adduct per line.")
                        raise
            elif polarity == 'pos':
                adducts = ['M+', 'M+H', 'M+NH4', 'M+Na']
            elif polarity == 'neg':
                adducts = ['M-H', 'M+Cl', 'M+FA-H', 'M+Hac-H']
            elif polarity == 'neut':
                adducts = ['']
            else:
                raise RuntimeError('Could not understand polarity')
            mass_searched_compounds_filename = mg.accurate_mass_search(compounds_file, polarity, adducts, ppm_cutoff)
            print("\n!!! Accurate mass search done. Mass-searched file stored in {}".format(compounds_file))    
            if accurate_mass_search_only:
                sys.exit() # done with mass search. Exiting
        return mass_searched_compounds_filename

def load_compound_results(compounds_file, pactolus, output_dir):
    """ load compound results"""
    print( '\n!!! LOADING COMPOUNDS')
    compounds = mg.load_dataframe(compounds_file)
    # auto-rename pactolus columns
    if pactolus:
        compounds = mg.reformat_pactolus(compounds)
    # remove any missing compounds
    compounds = compounds[~pd.isnull(compounds['original_compound'])]
    compounds.fillna('', inplace=True)

    if 'original_compound' not in compounds.columns:
        raise RuntimeError('Could not find "original_compound" as a column, please\
            rename the column corresponding to inchi keys for the compounds')
    
    # remove compounds not in the database or network
    print( '!!! Scrubbing compounds')
    compounds['adj'] = compounds['original_compound'].apply(
        lambda x: '-'.join(x.split('-')[:2]))

    mg.compounds['adj'] = mg.compounds['inchi_key'].apply(
            lambda x: '-'.join(x.split('-')[:2]))

    filt = compounds.merge(mg.compounds, on='adj', how='left', suffixes=('', '_db'))
    # categorize their reason for not being searched
    not_in_db = filt[pd.isnull(filt['cpd_group'])]
    not_in_db['not_searched_reason'] = 'Not in metabolite database'
    not_in_net = filt[filt['cpd_group'] < 0]
    not_in_net['not_searched_reason'] = 'Not in similarity network yet'
    # combine into one table
    not_searched = pd.concat([not_in_db, not_in_net])
    # make the columns same as user input
    cols = compounds.columns[~compounds.columns.str.contains('adj')].tolist()
    cols.append('not_searched_reason')
    not_searched = not_searched[cols]
    # inform the user and save file
    if not_searched.shape[0] > 0:
        print( 'WARNING: some input compounds were not found in the metabolite database or chemical network; please report these compounds! (see log_unsearched_compounds.csv)')
        print( '!@#', not_searched['original_compound'].unique().shape[0],\
            'Compounds not being searched; see log_unsearched_compounds.csv')
        not_searched.to_csv(os.path.join(output_dir,
            'log_unsearched_compounds.csv'), index=False)

    to_search = filt[filt['cpd_group'] > 0]['original_compound'].unique()
    compounds = compounds[compounds['original_compound'].isin(to_search)]

    u_cpds = compounds['original_compound'].unique()
    print( '!@#', len(u_cpds), 'total input compounds to search\n')

    if 'compound_score' not in compounds.columns:
        print( 'WARNING: "compound_score" not found as a column; assuming that\
            there is no score for compounds, and setting the compound scores \
            to 1.0')
        compounds['compound_score'] = 1.0
    else:
        compounds['compound_score'] = compounds['compound_score'].apply(float)
    return compounds



def gene_to_reaction_search(blast_filter, gene_to_reaction, intermediate_files_dir, cpu_count, compounds_file, output_dir, legacy, level, genome, main_start):
    # Conduct gene to reaction search

    def keep_top_blast_helper(x, param=blast_filter):
        """
        # setup multiprocessing helpers based on input params
        x is the normal input
        param is the defined parameter
        """
        return mg.keep_top_blast(x, filt=param)
    #TODO: Maybe find a better way to set if something is already done?    
    if gene_to_reaction is None:
        print( '!@# Conducting gene to reaction search | TLOG %s' % (time.time()))
        start = time.time()
        gene_blast = mg.multi_blast(genome.index, genome, mg.refseq_dbpath, 
            intermediate_files_dir, raise_blast_error=False, cpu=cpu_count)

        print( '!@# Homology searching done in %s minutes' \
                %((time.time() - start) / 60))
        gene_blast.to_pickle(os.path.join(intermediate_files_dir, 'gene_blast.pkl'))
        print( '!!! g2r blast results saved to %s' \
                %(os.path.join(intermediate_files_dir, 'g2r_blast.pkl')))

        start = time.time()
        gene_to_reaction = mg.refseq_to_reactions(gene_blast, 'subject acc.')
        del gene_blast
        gene_groups = gene_to_reaction.groupby('query acc.')
        multidx = gene_groups['e_score'].apply(keep_top_blast_helper).index
        idx = multidx.levels[1]
        gene_to_reaction_top = gene_to_reaction.loc[idx]
        del gene_to_reaction
        print( '!@# gene_to_reaction table completed in %s minutes' \
                %((time.time() - start) / 60))
        # if not compounds file, then just quit
        if compounds_file is None:
            df = gene_to_reaction_top.merge(mg.mrs_reaction[['database_id']],
                left_on='reaction_id', right_index=True, how='left')
            df.to_csv(os.path.join(output_dir, 'magi_gene_results.csv'))
            print( '!!! gene to reaction results saved to %s' \
                    %(os.path.join(output_dir, 'magi_gene_results.csv')))
            print( '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60))
            sys.exit()
        else:
            gene_to_reaction_top.to_pickle(os.path.join(intermediate_files_dir, 
                                                    'gene_to_reaction.pkl'))
            print( '!!! gene to reaction results saved to %s' \
                    %(os.path.join(intermediate_files_dir, 'gene_to_reaction.pkl')))
    else:
        gene_to_reaction_top = pd.read_pickle(gene_to_reaction)
        print( '\n!@# gene_to_reaction successfully loaded')
    return gene_to_reaction_top

def compound_to_reaction_search(legacy, level, compound_to_reaction, compounds, main_start, cpu_count, fasta_file, output_dir, intermediate_files_dir): 
    """Perform compound to reaction search"""
    def connect_compound_to_reaction_mp_helper(inchikey, 
                                            tautomer=legacy, 
                                            neighbor_level=level):
        try:
            out = mg.connect_compound_to_reaction(inchikey, 
                                            tautomer=tautomer, 
                                            neighbor_level=neighbor_level)
        except Exception as e:
            print( inchikey )
            sys.stdout.flush()
            raise RuntimeError('offending inchikey: %s; error message: %s' \
                                %(inchikey, e.args))
        return out
    
    if compound_to_reaction is None:
        print( '\n!@# Conducting compound to reaction search | TLOG %s' % (time.time()))
        sys.stdout.flush()
        start = time.time()

        input_compounds = compounds['original_compound'].unique()
        if os.name == 'nt': #Check if the operating system is windows or linux/mac
            #To do: fix this
            print("!!! Operating system is Windows. No multiprocessing used for compound to reaction search")
            out = map(connect_compound_to_reaction_mp_helper, input_compounds)
        else:
            p = mp.Pool(cpu_count)
            out = p.map(connect_compound_to_reaction_mp_helper, input_compounds)
            p.close()
            p.terminate()
        
        compound_to_reaction = pd.concat(out)
        del out
        compound_to_reaction.reset_index(inplace=True, drop=True)
    
        # connect the compound score
        compound_to_reaction = pd.merge(compounds, compound_to_reaction, 
                                        on='original_compound', how='inner')
        print( '!@# compound_to_reaction table done in %s minutes'\
                    %((time.time()-start)/60))
        # if no fasta file then just save these results and quit
        if fasta_file is None:
            compound_to_reaction.to_csv(os.path.join(output_dir, 
                                                    'magi_compound_results.csv'))
            print( '!!! compound_reaction table saved to %s'\
                    % (os.path.join(output_dir, 'magi_compound_results.csv')))
            print( '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60))
            sys.exit()
        else:
            compound_to_reaction.to_pickle(os.path.join(intermediate_files_dir, 
                                                    'compound_to_reaction.pkl'))
    
            print( '!!! compound_reaction table saved to %s'\
                    % (os.path.join(intermediate_files_dir, 'compound_to_reaction.pkl')))
    else:
        compound_to_reaction = pd.read_pickle(compound_to_reaction)
        print( '\n!@# compound_to_reaction successfully loaded')
    return compound_to_reaction

def reaction_to_gene_search(compound_to_reaction, genome_db_path, blast_filter, reaction_to_gene, intermediate_files_dir, cpu_count):
    # reaction to gene search
    def keep_top_blast_helper(x, param=blast_filter):
        """
        # setup multiprocessing helpers based on input params
        x is the normal input
        param is the defined parameter
        """
        return mg.keep_top_blast(x, filt=param)
    
    if reaction_to_gene is None:
    	print( '\n!@# Conducting reaction to gene search | TLOG %s' % (time.time()))
    	sys.stdout.flush()
    	start = time.time()
    
    	# set up a list of reference sequences to blast against the genome
    	reactions = compound_to_reaction[compound_to_reaction\
    						['reaction_id'] != '']['reaction_id'].tolist()
    	reactions_refseqs = mg.mrs_reaction.loc[reactions, 'refseq_id']
    	del reactions
    	reactions_refseqs = reactions_refseqs[reactions_refseqs != '']
    	rseq_list = []
    	for reaction in reactions_refseqs:
    	    for rseq in reaction.split('|'):
    	        if rseq != '':
    	            rseq_list.append(rseq)
    	rseq_list = list(set(rseq_list))
    
    	# rseq_list is the "query_list" for multi_blast()
    	# query_full_table is the refseq table
    	# database_path is the path to the genome's blast database
    	print( '!!!', len(rseq_list), 'reference sequences to search')
    	sys.stdout.flush()
    
    	reaction_to_gene_blast = mg.multi_blast(rseq_list, mg.refseq, 
    		genome_db_path, intermediate_files_dir, cpu=cpu_count, 
    		raise_blast_error=False)
    
    	reaction_to_gene = mg.refseq_to_reactions(reaction_to_gene_blast,
    		'query acc.')
    	del reaction_to_gene_blast
    
    	reaction_to_gene.to_pickle(os.path.join(intermediate_files_dir,
    		'reaction_blast.pkl'))
    	print( '!!! r2g blast results saved to %s' \
    			%(os.path.join(intermediate_files_dir, 'r2g_blast.pkl')))
    
    	reaction_groups = reaction_to_gene.groupby('query acc.')
    	multidx = reaction_groups['e_score'].apply(keep_top_blast_helper).index
    	del reaction_groups
    	idx = multidx.levels[1]
    	reaction_to_gene_top = reaction_to_gene.loc[idx]
    	del reaction_to_gene
    	reaction_to_gene_top.to_pickle(os.path.join(intermediate_files_dir, 
    											'reaction_to_gene.pkl'))
    	print( '!@# reaction_to_gene table done in %s minutes'\
    			%((time.time()-start)/60))
    	print( '!!! reaction_to_gene table saved to %s'\
    			% (os.path.join(intermediate_files_dir, 'reaction_to_gene.pkl')))
    else:
    	reaction_to_gene_top = pd.read_pickle(reaction_to_gene)
    	print( '\n!@# reaction_to_gene successfully loaded')
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
    
    df = df.merge(mg.mrs_reaction[['database_id']],
    	left_on='reaction_id_r2g', right_index=True, how='left')
    df = df.merge(mg.mrs_reaction[['database_id']],
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
    print( 'full results saved to', os.path.join(output_dir, 'magi_results.csv'))
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
    args = check_parameters(args)
    
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
    output_dir, intermediate_files_dir = make_output_dir_name(output_dir, fasta_file, compounds_file, intermediate_files) #TODO: rename this function
    main_start = time.time() # overall program timer
    if fasta_file is not None:
        genome, genome_db_path = load_fasta_genome(fasta_file, intermediate_files_dir, annotations)
    if accurate_mass_search:
        args, compounds_file = perform_accurate_mass_search(compounds_file, adduct_file, polarity, accurate_mass_search_only, ppm_cutoff)
    if compounds_file is not None:
        compounds = load_compound_results(compounds_file, pactolus, output_dir) #TODO: rename this function
    if fasta_file is not None:
        gene_to_reaction_top = gene_to_reaction_search(blast_filter, gene_to_reaction, intermediate_files_dir, cpu_count, compounds_file, output_dir, legacy, level, genome, main_start)
        del genome
    if compounds_file is not None:
        compound_to_reaction = compound_to_reaction_search(legacy, level, compound_to_reaction, compounds, main_start, cpu_count, fasta_file, output_dir, intermediate_files_dir)
        reaction_to_gene_top = reaction_to_gene_search(compound_to_reaction, genome_db_path, blast_filter, reaction_to_gene, intermediate_files_dir, cpu_count)
    merged_dataframe = merging_g2r_and_r2g_searches(compound_to_reaction, reaction_to_gene_top, gene_to_reaction_top, merged_before_score, intermediate_files_dir)
    final_dataframe, start = calculate_scores_and_format_tables(merged_dataframe, reciprocal_closeness, final_weights, chemnet_penalty)
    save_outputs(final_dataframe, compounds, main_start, start, output_dir)
    
if __name__ == "__main__":
    main()