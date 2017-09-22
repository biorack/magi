#!/bin/bash

time python ../../workflow/magi_workflow_20170519.py \
--fasta ./s_coelicolor_genes_fasta_smallset.faa \
--output ./test_output_files \
--cpu_count 4 --mute
