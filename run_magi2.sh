#!/bin/bash

# Name of the run 
name=magi2_quicktest

# Full path to location where magi is installed
magi_path=.

# Input fasta file with protein sequences and unique identifiers in the header
# Example of a header:
# >gene_1 some description
fasta=./tests/full_workflow_test/s_coelicolor_genes_fasta_smallset.faa

# Input file with candidate compound SMILES in a column called original_compound
compounds=./tests/full_workflow_test/random_smiles.csv

# Other parameters for MAGI that may be useful
cpu_count=8
min_diameter=12 # this is the Retro Rules diameter, see https://retrorules.org/doc for details
output_directory=./tests/$name/output_$name
## Log files
logfile_name=./tests/$name/log_magi_run_$name.txt
error_log_name=./tests/$name/error_log_magi_run_$name.txt

#####################################################################################
source activate magi_2 #if this does not work, use conda activate magi

echo "Starting MAGI at $(date)"
python $magi_path/workflow_2/compound_to_reaction.py \
    --compounds $compounds \
    --fasta $fasta \
    --diameter $min_diameter \
    --cpu_count $cpu_count \
    --use_precomputed_reactions True \
    --output $output_directory > $logfile_name 2> $error_log_name
if [ $? -eq 0 ]; then # Check if previous MAGI run did not fail
python $magi_path/workflow_2/gene_to_reaction.py --not_first_script --output $output_directory >> $logfile_name 2>> $error_log_name
else exit 1; fi
if [ $? -eq 0 ]; then
python $magi_path/workflow_2/reaction_to_gene.py --not_first_script --output $output_directory >> $logfile_name 2>> $error_log_name
else exit 1; fi
if [ $? -eq 0 ]; then
python $magi_path/workflow_2/scoring.py --not_first_script --output $output_directory >> $logfile_name 2>> $error_log_name
else exit 1; fi
echo "Done"
