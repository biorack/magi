"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.1 reaction to gene workflow script

Lines beginning with "@@@" are input job parameters
Lines beginning with "!!!" are verbose log info
Lines beginning with "!@#" are checkpoints/announcements

The reaction-to-gene search is conducted by querying a reaction's reference sequence 
against a BLAST database created by using the user's input protein FASTA file. 
The results of this search can be interpreted as "the input proteins that can probably catalyze this reaction", 
with the probability being represented by the E-score.

Required input is either a compound_to_reaction pickle file or
the directory in which intermediate files are stored. The input of this
script should have been created in the compound to reaction workflow.
"""

import sys
import os
import argparse
import pandas as pd
from multiprocessing import cpu_count as counting_cpus
import time
import blast_helpers as blast
import workflow_helpers_new as mg

def parse_arguments():
    def is_existing_file(filepath):
        """Checks if a file exists and return absolute path if it exists"""
        if not os.path.exists(filepath):
            msg = "{0} does not exist".format(filepath)
            raise argparse.ArgumentTypeError(msg)
        else:
            return os.path.abspath(filepath)
        
    def is_database(db_path):
        """
        Checks if the genome database exists.
        Three files should be present, a .phr file, a .pin file and a .psq file
        """
        for file_extension in [".phr", ".pin", ".psq"]:
            db_file = db_path + file_extension
            if not os.path.exists(db_file):
                msg = "{0} does not exist".format(db_file)
                raise argparse.ArgumentTypeError(msg)
        return os.path.abspath(db_path)
    
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
        parser.add_argument('--genome_db', help = "path to genome .db files", type=is_database, required = True)
        parser.add_argument('--compound_to_reaction', help = "compound to reaction pickle", type=is_existing_file, required = True)
        parser.add_argument('-o', '--output', 
            help='path to a custom output', required = True,
            type=str)
        # optional runtime variables
        parser.add_argument('-n', '--cpu_count', 
            help='number of cpus to use for multiprocessing. Default is to use max!', 
            type=int, default=0)

        parser.add_argument('--blast_filter', 
            help='How stringent to filter the top BLAST results, as percent;\
            default is 85 meaning that only BLAST results within 85%% of the top\
            result will be taken.', 
            type=percentage_values_to_decimal, default=0.85)
        parser.add_argument('--intermediate_files',
            help='What directory within --output to store intermediate files',
            type=str, default='intermediate_files')

        args = parser.parse_args()
        
        # Set number of required CPUs
        args.cpu_count = set_cpu_count(args.cpu_count)
        
    except argparse.ArgumentTypeError as ex:
        print(ex.message)
        sys.exit(1)
    return args
    
def workflow(compound_to_reaction, genome_db_path, blast_filter, intermediate_files_dir, cpu_count):
    # Read stuff
    print( '!!! loading refseq and reaction tables...')
    # this table is only refseqs that are found in mrs-reaction
    my_settings = mg.get_settings()
    refseq_path = my_settings.refseq_path
    print( '!!! Reference sequences in this file: {}'.format(refseq_path))
    refseq = mg.load_dataframe(refseq_path)
    refseq.dropna(inplace=True)
    print( '!!! {} reference sequences'.format(len(refseq)))
    # reaction to gene search
    def keep_top_blast_helper(x, param=blast_filter):
        """
        # setup multiprocessing helpers based on input params
        x is the normal input
        param is the defined parameter
        """
        return blast.keep_top_blast(x, filt=param)
    
    # Perform R2G
    print( '\n!@# Conducting reaction to gene search | TLOG %s' % (time.time()))
    sys.stdout.flush()
    start = time.time()

    # set up a list of reference sequences to blast against the genome
    reactions = compound_to_reaction[compound_to_reaction\
                        ['reaction_id'] != '']['reaction_id'].tolist()
    reactions_refseqs = blast.mrs_reaction.loc[reactions, 'refseq_id']
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

    reaction_to_gene_blast = blast.multi_blast(rseq_list, refseq, 
        genome_db_path, intermediate_files_dir, cpu=cpu_count, 
        raise_blast_error=False)

    reaction_to_gene = blast.refseq_to_reactions(reaction_to_gene_blast,
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
    reaction_to_gene_top_path = os.path.join(intermediate_files_dir, 
                                            'reaction_to_gene.pkl')
    reaction_to_gene_top.to_pickle(reaction_to_gene_top_path)
    print( '!@# reaction_to_gene table done in %s minutes'\
            %((time.time()-start)/60))
    print( '!!! reaction_to_gene table saved to %s'\
            % reaction_to_gene_top_path)
    
    return reaction_to_gene_top_path

def main():
    # Parse arguments and prepare for reaction to gene workflow
    magi_parameters = mg.general_magi_preparation()

    # Set compound to reaction file path
    if magi_parameters["compound_to_reaction"] is not None:
        compound_to_reaction_path = magi_parameters["compound_to_reaction"]
    else:
        compound_to_reaction_path = mg.get_intermediate_file_path(magi_parameters["intermediate_files_dir"], "compound_to_reaction_path")
    # Set genome database path
    if magi_parameters["genome_db"] is not None:
        genome_db_path = magi_parameters["genome_db"]
    else:
        genome_db_path = mg.get_intermediate_file_path(magi_parameters["intermediate_files_dir"], "genome_db_path")
    
    print("opening {}".format(magi_parameters["compound_to_reaction"]))
    compound_to_reaction = pd.read_pickle(compound_to_reaction_path)
    # Start r2g workflow
    reaction_to_gene_path = workflow(compound_to_reaction = compound_to_reaction, 
        genome_db_path = genome_db_path, 
        blast_filter = magi_parameters["blast_filter"], 
        intermediate_files_dir = magi_parameters["intermediate_files_dir"], 
        cpu_count = magi_parameters["cpu_count"])
    mg.write_intermediate_file_path(magi_parameters["intermediate_files_dir"], "reaction_to_gene_path", reaction_to_gene_path)
    print("Reaction to gene search successfully finished")

if __name__ == "__main__":
    main()