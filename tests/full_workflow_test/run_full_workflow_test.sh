#!/bin/bash

time python /project/projectdirs/metatlas/projects/metatlas_reactions/workflow/magi_workflow_20170519.py \
--fasta ./s_coelicolor_genes_fasta_smallset.faa \
--compounds ./s_coelicolor_pactolus_data_smallset.csv \
--annotations ./s_coelicolor_genes_annotations_smallset.tab \
--output ./test_output_files \
--cpu_count 10 --mute
