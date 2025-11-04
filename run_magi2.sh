#!/bin/bash

if [ "$#" -ne 5 ]; then
  echo "Usage: script.sh name_of_run fasta_sequence compound_table output_dir magi_source_dir" >&2
  exit 1
fi

# Name of the run
name=$1
sequences=$2
compounds=$3
output_dir=$4
magi_source_dir=$5
if [ ! -e $output_dir/$name ]; then
	mkdir $output_dir/$name
fi

# Full path to location where magi is installed
magi_path=$magi_source_dir

# Input fasta file with protein sequences and unique identifiers in the header
# Example of a header:
# >gene_1 some description
fasta=$sequences

# Input file with candidate compound SMILES in a column called original_compound
compounds=$compounds

# Other parameters for MAGI that may be useful
cpu_count=8
min_diameter=12 # this is the Retro Rules diameter, see https://retrorules.org/doc for details
output_directory=$output_dir/$name/output_$name
## Log files
logfile_name=$output_dir/$name/log_magi_run_$name.txt
error_log_name=$output_dir/$name/error_log_magi_run_$name.txt

#####################################################################################
module load conda
conda activate magi_2 #if this does not work, use conda activate magi

echo "Starting MAGI at $(date)"

python $magi_path/workflow_2/compound_to_reaction.py \
    --compounds $compounds \
    --fasta $fasta \
    --diameter $min_diameter \
    --cpu_count $cpu_count \
    --use_precomputed_reactions True \
    --output $output_directory \
    > $logfile_name 2> $error_log_name

if [ $? -eq 0 ]; then # Check if previous MAGI run did not fail
python $magi_path/workflow_2/gene_to_reaction.py \
    --not_first_script \
    --output $output_directory \
    >> $logfile_name 2>> $error_log_name
else exit 1; fi

if [ $? -eq 0 ]; then
python $magi_path/workflow_2/reaction_to_gene.py \
    --not_first_script \
    --output $output_directory \
    >> $logfile_name 2>> $error_log_name
else exit 1; fi

if [ $? -eq 0 ]; then
python $magi_path/workflow_2/scoring.py \
    --not_first_script \
    --output $output_directory \
    >> $logfile_name 2>> $error_log_name
else exit 1; fi

echo "Done"
