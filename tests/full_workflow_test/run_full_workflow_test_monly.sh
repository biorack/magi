#!/bin/bash

time python ../../workflow/magi_workflow.py \
--compounds ./s_coelicolor_pactolus_data_smallset.csv \
--output ./test_output_files \
--cpu_count 4 --mute
