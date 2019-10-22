#"Run gene to reaction search in separate script"
#
import sys
import os
#import argparse
#import warnings
#import multiprocessing as mp
import pandas as pd
#import numpy as np
import time
#import pickle
#import datetime
import blast_helpers as mg
import workflow_helpers_new as magi_settings
#
def workflow(compound_to_reaction, genome_db_path, blast_filter, reaction_to_gene, intermediate_files_dir, cpu_count):
    print( '!!! loading refseq and reaction tables...')
    # this table is only refseqs that are found in mrs-reaction
    my_settings = magi_settings.get_settings()
    refseq_path = my_settings.refseq_path
    print( '!!! Reference sequences in this file: {}'.format(refseq_path))
    refseq = magi_settings.load_dataframe(refseq_path)
    refseq.dropna(inplace=True)
    print( '!!! {} reference sequences'.format(len(refseq)))
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
    	print( '!!! {} reference sequences to search'.format(len(rseq_list)))
    	sys.stdout.flush()
    
    	reaction_to_gene_blast = mg.multi_blast(rseq_list, refseq, 
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

#def parse_arguments():
#    return "arguments are parsed"
#
#def format_output():
#    return "output is formatted"
#
##def main():
##    parse_stuff()
##    workflow()
##    format_output()
#
#if __name__ == "__main__":
#    #main()
#    print("Are you sure? This is not ready yet")