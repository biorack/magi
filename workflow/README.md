# The MAGI Workflow
* The entire process is wrapped in `MAGI_NERSCworkflow<date>.py`
* Helper functions called in the workflow `.py` file are stored in `helpertools/`
* Important database files are stored in `database/`
* BLAST binary files are stored in `blastbin/`
* `TEMPLATE_QSUB.qsub` is a read-only template with which to build Genepool submission scripts
* `test_coelicolor.qsub` is a Genepool submission script that runs the entire workflow on *Streptomyces coelicolor* C18-neg data, but only for the first 24 compounds
    * Running the job (`qsub test_coelicolor.qsub`) should take less than 5 minutes, but make sure to look through the header and change things like email, working directory, etc
    * Running the job will store output data here: `/global/project/projectdirs/openmsi/projects/temp_chem_net_data/MAGI_data/s_coelicolor_genes_fasta_<TODAYS DATE>` 
    * `coelicolor_test_data/` stores the test-run input files
    
## The workflow.py script has many available arguments
You can see a summary by just running the .py file with -h as an argument, but here's a list anyway

* REQUIRED
    * `-f, --fasta`, path to fasta file of genes in sample
    * `-g, --gene_info`, path to gene info table downloaded from IMG that corresponds to fasta
    * `-p, --pactolus`, path to pactolus results file
* OPTIONAL
    * `-n, --cpu_count`, number of cpus to use for multiprocessing, type=int, default=24
    * `-o, --output`, path to a custom output
    * `-t, --test`, run MAGI only on the first # of pactolus compounds, type=int
    * `-l, --level`, how many levels deep to search the chemical network, type=int, choices=[1,2,3], default=2)
    * `--mute`, mutes pandas warnings just calling this flag turns it on; no other input needed
* JUMP-STARTERS (if you want to start the script after a major computation, *e.g.* your run crashed)
    * `--level0_sdf`, magi_level0.pkl file if already computed; must also provide --level0_mdf file!)
    * `--level0_mdf`, magi_level0_metadata.pkl file if already computed; must also provide --level0_sdf file!
    * `--neighbor`, magi_level#_neighbors.pkl file if already computed
    * `--neighbor_blast`, magi_level#_multimerge.pkl file if already computed
