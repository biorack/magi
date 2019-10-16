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
        print '@@@ annotations input: %s' %(args.annotations)
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
        print '@@@ BLAST filter: %s' % (args.blast_filter)
    if args.compounds is not None:
        print '@@@ Using precomputed compound results: %s' % (not args.legacy)
        print '@@@ Chemnet search to level %s' % (args.level)
        print '@@@ Reciprocal closeness: %s' % (args.reciprocal_closeness)
        print '@@@ Chemnet base penalty: %s' % (args.chemnet_penalty)
    print '@@@ MAGI score weights: %s' % (args.final_weights)
    print '@@@ Using %s CPUs' % (args.cpu_count)
    
    if args.gene_to_reaction is not None:
        args.gene_to_reaction =  os.path.abspath(args.gene_to_reaction)
        if not os.path.isfile(args.gene_to_reaction):
            raise IOError('%s does not exist!' %(args.gene_to_reaction))
        else:
            print '@@@ gene_to_reaction input: %s' %(args.gene_to_reaction)
    
    if args.compound_to_reaction is not None:
        args.compound_to_reaction =  os.path.abspath(args.compound_to_reaction)
        if not os.path.isfile(args.compound_to_reaction):
            raise IOError('%s does not exist!' %(args.compound_to_reaction))
        else:
            print '@@@ compound_to_reaction input: %s' %(args.compound_to_reaction)
    
    if args.reaction_to_gene is not None:
        args.reaction_to_gene =  os.path.abspath(args.reaction_to_gene)
        if not os.path.isfile(args.reaction_to_gene):
            raise IOError('%s does not exist!' %(args.reaction_to_gene))
        else:
            print '@@@ reaction_to_gene input: %s' %(args.reaction_to_gene)
    
    if args.mute:
        print '!!! Warnings are muted'
        warnings.filterwarnings('ignore')
    return args
    
def make_output_dir_name(args):
    """set up where the results will be stored"""
    if args.output is None:
        # autoname the directory based on fasta, or compound file
        # this will change eventually
        if args.fasta is not None:
            experiment_name = os.path.splitext(os.path.basename(args.fasta))[0]
            experiment_dir = os.path.abspath(os.path.dirname(args.fasta))
        else:
            experiment_name = os.path.splitext(os.path.basename(args.compounds))[0]
            experiment_dir = os.path.abspath(os.path.dirname(args.compounds))
        today = datetime.datetime.now()
        experiment_name += today.strftime('_%Y%m%d')
        experiment_path = os.path.join(experiment_dir, experiment_name)
    else:
        experiment_path = args.output
    
    experiment_path = os.path.abspath(experiment_path)
    intfile_path = os.path.join(experiment_path, args.intermediate_files)
    
    print '!!! Saving all results here:', experiment_path
    if not os.path.isdir(experiment_path):
        os.makedirs(experiment_path)
    return args, intfile_path, experiment_path
    
def load_fasta_genome(fasta_filename, intfile_path, annotation_file = None):
    """load genome"""
    print '\n!!! LOADING GENOME'
    genome, genome_db_path = mg.load_genome(fasta_filename, intfile_path, 
                                        annotation_file)
    return genome, genome_db_path
    
def perform_accurate_mass_search(args):
        if args.compounds is None:
            raise RuntimeError("No compounds file specified. Exiting...")
        else:
            # Perform accurate mass search and set compounds file to mass-searched file.
            print("\n!!! Performing accurate mass search for {}".format(args.compounds))
            polarity = args.accurate_mass_search    
            #Make list of adducts to search for
            if args.adduct_file is not None:
                args.adduct_file = os.path.abspath(args.adduct_file)
                print('@@@ Adduct file input: %s' %(args.adduct_file))
                with open(args.adduct_file) as adduct_file:
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
            args.compounds = mg.accurate_mass_search(args.compounds, polarity, adducts, args.ppm_cutoff)
            print("\n!!! Accurate mass search done. Mass-searched file stored in {}".format(args.compounds))    
            if args.accurate_mass_search_only:
                sys.exit() # done with mass search. Exiting
        return args

def load_compound_results(args, experiment_path):
    """ load compound results"""
    print '\n!!! LOADING COMPOUNDS'
    compounds = mg.load_dataframe(args.compounds)
    # auto-rename pactolus columns
    if args.pactolus:
        compounds = mg.reformat_pactolus(compounds)
    # remove any missing compounds
    compounds = compounds[~pd.isnull(compounds['original_compound'])]
    compounds.fillna('', inplace=True)

    if 'original_compound' not in compounds.columns:
        raise RuntimeError('Could not find "original_compound" as a column, please\
            rename the column corresponding to inchi keys for the compounds')
    
    # remove compounds not in the database or network
    print '!!! Scrubbing compounds'
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
        print 'WARNING: some input compounds were not found in the metabolite database or chemical network; please report these compounds! (see log_unsearched_compounds.csv)'
        print '!@#', not_searched['original_compound'].unique().shape[0],\
            'Compounds not being searched; see log_unsearched_compounds.csv'
        not_searched.to_csv(os.path.join(experiment_path,
            'log_unsearched_compounds.csv'), index=False)

    to_search = filt[filt['cpd_group'] > 0]['original_compound'].unique()
    compounds = compounds[compounds['original_compound'].isin(to_search)]

    u_cpds = compounds['original_compound'].unique()
    print '!@#', len(u_cpds), 'total input compounds to search\n'

    if 'compound_score' not in compounds.columns:
        print 'WARNING: "compound_score" not found as a column; assuming that\
            there is no score for compounds, and setting the compound scores \
            to 1.0'
        compounds['compound_score'] = 1.0
    else:
        compounds['compound_score'] = compounds['compound_score'].apply(float)
    return compounds



def gene_to_reaction_search(args, genome, intfile_path, experiment_path, main_start):
    # Conduct gene to reaction search

    def keep_top_blast_helper(x, param=args.blast_filter):
        """
        # setup multiprocessing helpers based on input params
        x is the normal input
        param is the defined parameter
        """
        return mg.keep_top_blast(x, filt=param)
    #TODO: Maybe find a better way to set if something is already done?    
    if args.gene_to_reaction is None:
        print '!@# Conducting gene to reaction search | TLOG %s' % (time.time())
        start = time.time()
        gene_blast = mg.multi_blast(genome.index, genome, mg.refseq_dbpath, 
            intfile_path, raise_blast_error=False, cpu=args.cpu_count)

        print '!@# Homology searching done in %s minutes' \
                %((time.time() - start) / 60)
        gene_blast.to_pickle(os.path.join(intfile_path, 'gene_blast.pkl'))
        print '!!! g2r blast results saved to %s' \
                %(os.path.join(intfile_path, 'g2r_blast.pkl'))

        start = time.time()
        gene_to_reaction = mg.refseq_to_reactions(gene_blast, 'subject acc.')
        del gene_blast
        gene_groups = gene_to_reaction.groupby('query acc.')
        multidx = gene_groups['e_score'].apply(keep_top_blast_helper).index
        idx = multidx.levels[1]
        gene_to_reaction_top = gene_to_reaction.loc[idx]
        del gene_to_reaction
        print '!@# gene_to_reaction table completed in %s minutes' \
                %((time.time() - start) / 60)
        # if not compounds file, then just quit
        if args.compounds is None:
            df = gene_to_reaction_top.merge(mg.mrs_reaction[['database_id']],
                left_on='reaction_id', right_index=True, how='left')
            df.to_csv(os.path.join(experiment_path, 'magi_gene_results.csv'))
            print '!!! gene to reaction results saved to %s' \
                    %(os.path.join(experiment_path, 'magi_gene_results.csv'))
            print '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60)
            sys.exit()
        else:
            gene_to_reaction_top.to_pickle(os.path.join(intfile_path, 
                                                    'gene_to_reaction.pkl'))
            print '!!! gene to reaction results saved to %s' \
                    %(os.path.join(intfile_path, 'gene_to_reaction.pkl'))
    else:
        gene_to_reaction_top = pd.read_pickle(args.gene_to_reaction)
        print '\n!@# gene_to_reaction successfully loaded'
    return gene_to_reaction_top

def compound_to_reaction_search(args, compounds, experiment_path, intfile_path, main_start): 
    """Perform compound to reaction search"""
    def connect_compound_to_reaction_mp_helper(inchikey, 
                                            tautomer=args.legacy, 
                                            neighbor_level=args.level):
        try:
            out = mg.connect_compound_to_reaction(inchikey, 
                                            tautomer=tautomer, 
                                            neighbor_level=neighbor_level)
        except Exception as e:
            print inchikey
            sys.stdout.flush()
            raise RuntimeError('offending inchikey: %s; error message: %s' \
                                %(inchikey, e.args))
        return out
    
    if args.compound_to_reaction is None:
        print '\n!@# Conducting compound to reaction search | TLOG %s' % (time.time())
        sys.stdout.flush()
        start = time.time()

        input_compounds = compounds['original_compound'].unique()
        if os.name == 'nt': #Check if the operating system is windows or linux/mac
            #To do: fix this
            print("!!! Operating system is Windows. No multiprocessing used for compound to reaction search")
            out = map(connect_compound_to_reaction_mp_helper, input_compounds)
        else:
            p = mp.Pool(args.cpu_count)
            out = p.map(connect_compound_to_reaction_mp_helper, input_compounds)
            p.close()
            p.terminate()
        
        compound_to_reaction = pd.concat(out)
        del out
        compound_to_reaction.reset_index(inplace=True, drop=True)
    
        # connect the compound score
        compound_to_reaction = pd.merge(compounds, compound_to_reaction, 
                                        on='original_compound', how='inner')
        print '!@# compound_to_reaction table done in %s minutes'\
                    %((time.time()-start)/60)
        # if no fasta file then just save these results and quit
        if args.fasta is None:
            compound_to_reaction.to_csv(os.path.join(experiment_path, 
                                                    'magi_compound_results.csv'))
            print '!!! compound_reaction table saved to %s'\
                    % (os.path.join(experiment_path, 'magi_compound_results.csv'))
            print '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60)
            sys.exit()
        else:
            compound_to_reaction.to_pickle(os.path.join(intfile_path, 
                                                    'compound_to_reaction.pkl'))
    
            print '!!! compound_reaction table saved to %s'\
                    % (os.path.join(intfile_path, 'compound_to_reaction.pkl'))
    else:
        compound_to_reaction = pd.read_pickle(args.compound_to_reaction)
        print '\n!@# compound_to_reaction successfully loaded'
    return compound_to_reaction

def main(args):
    print_version_info()
    args = check_parameters(args)
    args, intfile_path, experiment_path = make_output_dir_name(args) #TODO: rename this function
    main_start = time.time() # overall program timer
    if args.fasta is not None:
        genome, genome_db_path = load_fasta_genome(args.fasta, intfile_path, args.annotations)
    if args.accurate_mass_search is not None:
        args = perform_accurate_mass_search(args)
    if args.compounds is not None:
        compounds = load_compound_results(args, experiment_path) #TODO: rename this function
    if args.fasta is not None:
        gene_to_reaction_top = gene_to_reaction_search(args, genome, intfile_path, experiment_path, main_start)
        del genome
    if args.compounds is not None:
        compound_to_reaction = compound_to_reaction_search(args, compounds, experiment_path, intfile_path, main_start)
    
if __name__ == "__main__":
    arguments = parse_arguments()
    main(arguments)