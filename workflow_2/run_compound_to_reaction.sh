#!/bin/bash

# Run new compound to reaction
python compound_to_reaction.py \
--compounds /global/cscratch1/sd/leegwate/magi_2/EMA_HILICz-v11_Annotations20180824_with_SMILES.csv \
--output ./output \
--retro_rules /global/cscratch1/sd/leegwate/magi_2/retrorules_rr02_flat_all.tsv