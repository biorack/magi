#!/bin/bash

time python /project/projectdirs/metatlas/projects/metatlas_reactions/workflow/magi_workflow.py \
--fasta ./s_coelicolor_genes_fasta_smallset.faa \
--compounds ./s_coelicolor_pactolus_data_smallset.csv \
--annotations ./s_coelicolor_genes_annotations_smallset.tab \
--output ./test_output_files \
--cpu_count 10 --final_weights 1 2 3 4 --blast_filter 100 --reciprocal_closeness 90 --chemnet_penalty 2 --mute
