#!/bin/bash

# Full path to location where magi is installed
magi_path=/path_to_magi/magi

# Input fasta file with protein sequences and unique identifiers in the header
# Example of a header:
# >gene_1 some description
fasta_file=../input_files/example_fasta_name.fa

# Input file with either candidate compounds in a column called original_compound
# Or with m/z values in the column original_compound
compounds_file=../input_files/example_compounds_name.csv

# Other parameters for MAGI that may be useful
cpu_count=1
output_directory=./output
logfile_name=log_magi_run.txt
error_log_name=error_log_magi_run.txt

#####################################################################################
source activate magi #if this does not work, use conda activate magi

# Run MAGI
echo "Starting MAGI at $(date)"
python $magi_path/workflow/magi_workflow_gene_to_reaction.py --fasta $fasta_file --compounds $compounds_file --output $output_directory --cpu_count $cpu_count > $logfile_name 2> $error_log_name
if [ $? -eq 0 ]; then # Check if previous MAGI run did not fail
python $magi_path/workflow/magi_workflow_accurate_mass_search.py --not_first_script --output $output_directory >> $logfile_name 2>> $error_log_name
else exit 1; fi
if [ $? -eq 0 ]; then
python $magi_path/workflow/magi_workflow_compound_to_reaction.py --not_first_script --output $output_directory >> $logfile_name 2>> $error_log_name
else exit 1; fi
if [ $? -eq 0 ]; then
python $magi_path/workflow/magi_workflow_reaction_to_gene.py --not_first_script --output $output_directory >> $logfile_name 2>> $error_log_name
else exit 1; fi
if [ $? -eq 0 ]; then
python $magi_path/workflow/magi_workflow_scoring.py --not_first_script --output $output_directory >> $logfile_name 2>> $error_log_name
else exit 1; fi
echo "Done"
