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
    

def main(args):
    print_version_info()
    args = check_parameters(args)
    args, intfile_path, experiment_path = make_output_dir_name(args) #TODO: rename this function
    
if __name__ == "__main__":
    arguments = parse_arguments()
    main(arguments)