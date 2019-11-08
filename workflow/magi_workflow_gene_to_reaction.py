"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.1 gene to reaction workflow script

Lines beginning with "@@@" are input job parameters
Lines beginning with "!!!" are verbose log info
Lines beginning with "!@#" are checkpoints/announcements

The gene-to-reaction search queries one user input protein sequence against a database of all reactions' reference sequences. 
This can be interpreted as "the reactions that an input protein can probably catalyze".

Required input is a FASTA file.
The FASTA file should be in standard FASTA format, E.g.:
'''fasta
>UNIQUE_GENE_ID OTHER_INFORMATION
AMINO_ACID_SEQUENCE

>UNIQUE_GENE_ID OTHER_INFORMATION
AMINO_ACID_SEQUENCE
'''
*Please note the space in between UNIQUE_GENE_ID and OTHER_INFORMATION*

The FASTA file will be used to make a BLAST database, and to make a
gene table with these columns:
Gene_ID        | header                           | sequence
-----------------------------------------------------------------------
UNIQUE_GENE_ID | UNIQUE_GENE_ID OTHER_INFORMATION | AMINO_ACID_SEQUENCE

"""

import sys
import os
import argparse
import pandas as pd
from multiprocessing import cpu_count as counting_cpus
import time
import subprocess
import workflow_helpers_new as mg
import blast_helpers as blast

def parse_arguments():
    """
    This function parses the arguments that MAGI needs to perform a gene to reaction search.
    """
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
        parser.add_argument('-f', '--fasta', type=is_existing_file,required=True,
            help='path to fasta file of genes in sample')

        parser.add_argument('-a', '--annotations', type=is_existing_file,
            help='path to annotation file for genes in sample', 
            default=None)
        parser.add_argument('-n', '--cpu_count', 
            help='number of cpus to use for multiprocessing. Default is to use max!', 
            type=int, default=0)
        parser.add_argument('-o', '--output', 
            help='path to a custom output', 
            type=str)
        parser.add_argument('--blast_filter', 
            help='How stringent to filter the top BLAST results, as percent;\
            default is 85 meaning that only BLAST results within 85%% of the top\
            result will be taken.', 
            type=percentage_values_to_decimal, default=0.85)
        parser.add_argument('--intermediate_files',
            help='What directory within --output to store intermediate files',
            type=str, default='intermediate_files')
        parser.add_argument('--gene_to_reaction_only',
            help="Use this parameter if you are only interested in the gene to reaction search", 
            action='store_true', default=False)
        
        args = parser.parse_args()
        
        # Set number of required CPUs
        args.cpu_count = set_cpu_count(args.cpu_count)
    except argparse.ArgumentTypeError as ex:
        print(ex.message)
        sys.exit(1)
    return args

def load_genome(fasta, intfile_path, annotation_file=None):
    """
    This function will take some standard fasta a input and convert it
    to an appropriate genome dataframe.
    This function will also construct a blast database for the genome.

    Inputs
    ------
    fasta: path to standard fasta file of all the genes.
           FASTA sequence header format should be:
           >UNIQUE_GENE_IDENTIFIER OTHER_INFORMATION
           Note the space between the unique identifier and other info

    intfile_path: path to the temporary storage place of the database. If None, the
               BLAST database is not made, and the gene table is not
               stored.

    annotation_file: path to an optional .TAB file that contains
                     annotations for the genes in the fasta file. One
                     column must be named "Gene_ID" and correspond to
                     the UNIQUE_GENE_IDENTIFIER in the fasta headers.

    Outputs
    -------
    genome: gene sequence table that is merged with the annotation table
            if one is provided. This table is saved as a pickle file to
            intfile_path/gene_fastas/filename.pkl
    db_path: path to the genome's BLAST database
    """
    # TODO: find a way to handle windows text files (\n\r for new lines)
    # TODO: find a way to handle additional '>' characters in the header

    # process the fasta file
    with open(fasta, 'r') as f:
        genes = f.read()
    if genes.strip()[0] != '>':
        raise RuntimeError('%s does not appear to be a FASTA file' % (fasta))

    data = []
    for entry in genes.split('>')[1:]:
        header = '>' + entry.splitlines()[0]
        gene_id = header.split(' ')[0][1:]
        sequence = ''.join(entry.splitlines()[1:])
        data.append([gene_id, header, sequence])
    genome = pd.DataFrame(data, columns=['Gene_ID', 'header', 'sequence'])
    if genome['Gene_ID'].duplicated().any():
        first_dup = genome[genome['Gene_ID'].duplicated()].head(1)
        first_header = first_dup.iloc[0, 1]
        first_identifier = first_dup.iloc[0, 0]
        raise RuntimeError('There are duplicated Gene_ID fields! please check\
            your unique gene identifiers in the fasta headers. For the first\
            FASTA sequence of %s, I parsed the unique gene identifier to be \
            %s' % (repr(first_header), repr(first_identifier)))

    if annotation_file is not None:
        # Make gene info table
        annotation_table = mg.load_dataframe(annotation_file)
        # turn gene IDs into strings to keep things consistent
        annotation_table['Gene_ID'] = annotation_table['Gene_ID'].apply(str)

        # clean up the nans
        annotation_table.fillna(value='', inplace=True)

        # Really only if the annotation file was an IMG-derived table
        try:
            annotation_table['EC'] = annotation_table['Enzyme'].apply(mg.ec_parse)

            # move the new EC column to front for ease of visualization
            cols = annotation_table.columns.tolist()
            newcols = cols[-1:] + cols[:-1]
            annotation_table = annotation_table[newcols]
        except KeyError:
            print( 'Could not find a column corresponding to EC annotations, \
                    skipping EC parsing')

        genome = pd.merge(genome, annotation_table, on='Gene_ID', how='left')

    print( '!@# FASTA file loaded with {} genes'.format(len(genome)))

    genome.set_index('Gene_ID', inplace=True, drop=True)
    genome_name = os.path.splitext(os.path.basename(fasta))[0]
    if intfile_path is not None:
        gene_fasta_path = os.path.join(intfile_path, 'gene_fastas')
        if not os.path.isdir(gene_fasta_path):
            os.makedirs(gene_fasta_path)
        gene_seq_path = os.path.join(gene_fasta_path,
            '%s_sequences.pkl' % (genome_name))
        # gene_seq_path = '%s/gene_fastas/%s_sequences.pkl' \
        #                 % (intfile_path, genome_name)
        genome.to_pickle(gene_seq_path)
        print( '!!! saved gene_sequence table here: {}'.format(gene_seq_path))

        # Make the blast database of the fasta
        # this command makes the blast database for a given fasta file.
        blastbin = mg.get_settings(); blastbin = blastbin.blastbin
        makeblastdb_path = os.path.join(blastbin, 'makeblastdb')
        fasta_path = fasta
        db_path = os.path.join(intfile_path, 'BLAST_dbs',(os.path.splitext(os.path.basename(fasta_path))[0]+'.db'))
        print( '!!! blast database stored here: {}'.format(db_path))
        if os.name == 'nt': #Check if the operating system is windows or linux/mac
            make_db_command = '%s -in %s -out %s -dbtype prot' \
                % (makeblastdb_path+'.exe', fasta_path, db_path)
        else:
            make_db_command = '%s -in %s -out %s -dbtype prot' \
                % (makeblastdb_path, fasta_path, db_path)
        blastp = subprocess.Popen(
            make_db_command,
            shell=True, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if blastp.stdout is not None:
            print( blastp.stdout.read())
        if blastp.stderr is not None:
            print( blastp.stderr.read())
    else:
        db_path = None
    return genome, db_path

def workflow(fasta_file, intermediate_files_dir, cpu_count,
             annotations=None, blast_filter=0.85):
    """
    Conduct gene to reaction search. 
    The gene to reaction object will be stored as a picle file.
    
    Inputs
    ------
    fasta_file:             The fasta file to analyze
    annotations:        Optional annotations file. 
                            For more info, see the load_genome function.
    blast_filter:           How stringent to filter the top BLAST results. 
                            Should be a number between 0 and 1. Default is 0.85.
    intermediate_files_dir: Location to store the gene to reaction file
    cpu_count:              Amount of CPUs to use for BLAST searches

    Outputs
    -------
    gene_to_reaction_path:  path to location of gene to reaction pickle file
    genome_db_path:         path to the genome database
    """

    def keep_top_blast_helper(x, param=blast_filter):
        """
        # setup multiprocessing helpers based on input params
        x is the normal input
        param is the defined parameter
        """
        return blast.keep_top_blast(x, filt=param)
    
    print( '\n!!! LOADING GENOME')
    genome, genome_db_path = load_genome(fasta_file, intermediate_files_dir, 
                                        annotations)
    print("!!! Genome successfully loaded")
    print( '!@# Conducting gene to reaction search | TLOG %s' % (time.time()))
    start = time.time()
    gene_blast = blast.multi_blast(genome.index, genome, blast.refseq_dbpath, 
        intermediate_files_dir, raise_blast_error=False, cpu=cpu_count)

    print( '!@# Homology searching done in %s minutes' \
            %((time.time() - start) / 60))
    gene_blast.to_pickle(os.path.join(intermediate_files_dir, 'gene_blast.pkl'))
    print( '!!! g2r blast results saved to %s' \
            %(os.path.join(intermediate_files_dir, 'g2r_blast.pkl')))

    start = time.time()
    gene_to_reaction = blast.refseq_to_reactions(gene_blast, 'subject acc.')
    del gene_blast
    gene_groups = gene_to_reaction.groupby('query acc.')
    multidx = gene_groups['e_score'].apply(keep_top_blast_helper).index
    idx = multidx.levels[1]
    gene_to_reaction_top = gene_to_reaction.loc[idx]
    del gene_to_reaction
    print( '!@# gene_to_reaction table completed in %s minutes' \
            %((time.time() - start) / 60))

    #Write output file
    gene_to_reaction_top_path = os.path.join(intermediate_files_dir, 
                                            'gene_to_reaction.pkl')
    gene_to_reaction_top.to_pickle(gene_to_reaction_top_path)

    print( '!!! gene to reaction results saved to %s' \
            %gene_to_reaction_top_path)
    return gene_to_reaction_top_path, genome_db_path

def format_output(gene_to_reaction_top_path, output_dir, intermediate_files_dir):
    """
    Add information on the genes to the 
    Inputs
    ------
    gene_to_reaction_top_path: Path to the picle object that contains 
                               all reactions that match genes in the input file.
    output_dir:                Location to store the gene results dataframe
    main_start:                Time when MAGI started
    """
    gene_to_reaction_top = pd.read_pickle(gene_to_reaction_top_path)
    df = gene_to_reaction_top.merge(blast.mrs_reaction[['database_id']],
        left_on='reaction_id', right_index=True, how='left')
    df.to_csv(os.path.join(output_dir, 'magi_gene_results.csv'))
    print( '!!! gene to reaction results saved to %s' \
            %(os.path.join(output_dir, 'magi_gene_results.csv')))
    with open(os.path.join(intermediate_files_dir, "timer.txt"),"r") as timerfile:
        main_start = float(timerfile.read())
    print( '\n!@# MAGI analysis complete in %s minutes' %(
                                      (time.time() - main_start) / 60))

def main():
    # Parse arguments and prepare for gene to reaction workflow
    magi_parameters = mg.general_magi_preparation()
    
    #Run gene to reaction workflow
    gene_to_reaction_path, genome_db_path = workflow(fasta_file=magi_parameters["fasta"], 
             intermediate_files_dir=magi_parameters["intermediate_files_dir"], 
             cpu_count=magi_parameters["cpu_count"],
             annotations=magi_parameters["annotations"], 
             blast_filter=magi_parameters["blast_filter"])
    mg.write_intermediate_file_path(magi_parameters["intermediate_files_dir"], "gene_to_reaction_path", gene_to_reaction_path)
    mg.write_intermediate_file_path(magi_parameters["intermediate_files_dir"], "genome_db_path", genome_db_path)
    
    #Format output if this is the last step of the workflow
    if magi_parameters["gene_to_reaction_only"]:
        g2r_file = os.path.join(magi_parameters["intermediate_files_dir"], 
                                            'gene_to_reaction.pkl')
        format_output(g2r_file, magi_parameters["output_dir"], magi_parameters["intermediate_files_dir"])

if __name__ == "__main__":
    main()