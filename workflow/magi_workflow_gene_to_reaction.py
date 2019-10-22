"Run gene to reaction search in separate script"

import sys
import os
import pandas as pd
import time
import blast_helpers as mg


def workflow(blast_filter, gene_to_reaction, intermediate_files_dir, cpu_count, compounds_file, output_dir, legacy, level, genome, main_start):
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

def parse_arguments():
    return "arguments are parsed"

def format_output():
    return "output is formatted"

#def main():
#    parse_stuff()
#    workflow()
#    format_output()

if __name__ == "__main__":
    #main()
    print("Are you sure? This is not ready yet")