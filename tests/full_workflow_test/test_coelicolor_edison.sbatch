#!/bin/bash -l
#SBATCH --account=m2650
#SBATCH --job-name=magi_test
#SBATCH --time=0:10:00
#SBATCH --output=/project/projectdirs/metatlas/projects/metatlas_reactions/tests/full_workflow_test/edisontestlog_out.txt
#SBATCH --error=/project/projectdirs/metatlas/projects/metatlas_reactions/tests/full_workflow_test/edisontestlog_err.txt
#SBATCH --partition=debug
## uncomment constraint flag for jobs on Cori
## #SBATCH --constraint=haswell
#SBATCH --license=project
#SBATCH --mail-user=oerbilgin@lbl.gov
#SBATCH --mail-type=BEGIN,END,FAIL

module load python/2.7-anaconda

time python /project/projectdirs/metatlas/projects/metatlas_reactions/workflow/magi_workflow.py \
--fasta /project/projectdirs/metatlas/projects/metatlas_reactions/tests/full_workflow_test/s_coelicolor_genes_fasta_smallset.faa \
--compounds /project/projectdirs/metatlas/projects/metatlas_reactions/tests/full_workflow_test/s_coelicolor_pactolus_data_smallset.csv \
--annotations /project/projectdirs/metatlas/projects/metatlas_reactions/tests/full_workflow_test/s_coelicolor_genes_annotations_smallset.tab \
--cpu_count 24 \
--final_weights 1 2 3 4 \
--blast_filter 100 \
--reciprocal_closeness 100 \
--chemnet_penalty 2 \
--output /project/projectdirs/metatlas/projects/metatlas_reactions/tests/full_workflow_test/test_output_files --mute
