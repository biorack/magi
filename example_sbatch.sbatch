#!/bin/bash -l
## This sbatch header submits to CORI
##
## To submit to Denovo, the following parameters are different:
## --constraint is removed
## --partition is different, or removed in lieu for specific resource requests
## --account is removed, or changed
## the module load line is either removed or changed as appropriate.

#SBATCH --account=m1541
#SBATCH --job-name=magi
#SBATCH --time=03:50:00
#SBATCH --output=<PATH>/log_out.txt
#SBATCH --error=<PATH>/log_err.txt
#SBATCH --partition=realtime
#SBATCH --constraint=haswell
#SBATCH --license=project
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=ALL

module load python/2.7-anaconda-4.4

python <PATH>/magi/workflow/magi_workflow_20170519.py \
--fasta <PATH>/fasta_file.faa \          # required
--compounds <PATH>/metabolite_file.csv \ # required
--level 2 \                       # optional
--final_weights 1.0 1.0 1.0 1.0 \ # optional
--blast_filter 85 \               # optional
--reciprocal_closeness 75 \       # optional
--chemnet_penalty 4.0 \           # optional
--output <PATH> --mute            # 'output' is required, 'mute' is recommended
