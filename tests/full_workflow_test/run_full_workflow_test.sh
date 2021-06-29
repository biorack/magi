#!/bin/bash
set -e # To exit when error occurs
output_directory=./test_output_files

# Run MAGI
echo "Starting MAGI at $(date)"
python ../../workflow/magi_workflow_gene_to_reaction.py \
    --fasta ./s_coelicolor_genes_fasta_tinyset.faa \
    --compounds ./s_coelicolor_pactolus_data_smallset.csv \
    --output $output_directory \
    --annotations ./s_coelicolor_genes_annotations_smallset.tab \
    --cpu_count 4 --mute
echo "Gene to reaction workflow finished without error."

python ../../workflow/magi_workflow_accurate_mass_search.py --not_first_script --output $output_directory
echo "Accurate mass search workflow finished without error."

python ../../workflow/magi_workflow_compound_to_reaction.py --not_first_script --output $output_directory
echo "Compound to reaction workflow finished without error."

python ../../workflow/magi_workflow_reaction_to_gene.py --not_first_script --output $output_directory
echo "Reaction to gene workflow finished without error."

python ../../workflow/magi_workflow_scoring.py --not_first_script --output $output_directory
echo "Scoring workflow finished without error."
echo "Done"
