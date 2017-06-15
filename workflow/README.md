# The MAGI Workflow
* The entire process is wrapped in `magi_workflow_<date>.py`
* Helper functions called in the workflow `.py` file are stored in `helpertools/`
* Important database files are stored in `database/`
* BLAST binary files are stored in `blastbin/`
* `TEMPLATE_QSUB.qsub` is a read-only template with which to build Genepool submission scripts
* `test_coelicolor.qsub` is a Genepool submission script that runs the entire workflow on *Streptomyces coelicolor* C18-neg data, but only for the first 24 compounds
    * Running the job (`qsub test_coelicolor.qsub`) should take less than 5 minutes, but make sure to look through the header and change things like email, working directory, etc
    * Running the job will store output data here: `/global/project/projectdirs/openmsi/projects/temp_chem_net_data/MAGI_data/s_coelicolor_genes_fasta_<TODAYS DATE>` 
    * `coelicolor_test_data/` stores the test-run input files

# Running the workflow
1. Ensure that your local settings setup is correct and points to all the appropriate data files
2. Only two input files are needed: the fasta file and the compounds file. Ensure that they are in the correct format (see below)
3. There will be several output files
    1. `gene_blast.pkl` is the table of *all* gene-to-reaction blast hits; the top hits have not been filtered yet
    2. `gene_to_reaction.pkl` is the table where the gene-to-reaction blast hits have been connected to reactions, and only the top hits were kept. Reactions are the index of the reaction database.
    3. `compound_to_reaction.pkl` is the table describing all the reactions that the input compounds are associated with. Reactions are the index of the reaction database.
    4. `reaction_to_gene.pkl` is the table describing the blast results of reactions from compound_to_reaction to the input FASTA sequences
    5. `merged_before_score.h5` is the merge of compound_to_reaction and reaction_to_gene on the reactions, and to gene_to_reaction on the genes, but before the MAGI score was calculated
    6. `magi_results.csv` is the full table of results that includes all the gene to reaction/compound results
    7. `magi_compound_results.csv` is a table describing only the best MAGI score for each compound suggestion. This table is merged back with the input compound table on `original_compound` and `compound_score`

## The workflow.py script has many available arguments
You can see a summary by just running the .py file with -h as an argument, but here's a list anyway

* REQUIRED
    * `-f, --fasta`, path to fasta file of genes in sample. The FASTA headers cannot have any `>` characters inside the header. The header should have the following format: `> UNIQUE_GENE_ID OTHER_INFORMATION` noting the space that delimits the unique gene ID from other information.
    * `-c, --compounds`, path to the compounds table. This table must have one column named `original_compound` that has the InChI keys for the observed compounds. There may also be a column named `compound_score` that scores each compound. If this column does not exist, all compounds will automatically be given a score of 1.0.
* OPTIONAL
    * `-a, --annotations`, path to a gene annotation table. This table must have one column named `Gene_ID` that corresponds to the gene identifier in the FASTA header for each gene.
    * `-n, --cpu_count`, number of cpus to use for multiprocessing, type=int, default=0, where 0 uses the maximum number of CPUs on the machine
    * `-o, --output`, path to a custom output
    * `-l, --level`, how many levels deep to search the chemical network, type=int, choices=[1,2,3], default=2)
    * `-t, --tautomer`, whether or not to include tautomers in the search, choices=[True, False], default=True
    * `--mute`, mutes pandas warnings just calling this flag turns it on; no other input needed
    * `--pactolus`, flag that tells the script the compound input is a pactolus file. Runs a custom converter to rename the pactolus input.
* JUMP-STARTERS (if you want to start the script after a major computation, *e.g.* your run crashed)
    * `--gene_to_reaction`, path to gene_to_reaction file, must be in pickle format)
    * `--compound_to_reaction`, path to compound_to_reaction file, must be in pickle format
    * `--reaction_to_gene`, path to reaction_to_gene file, must be in pickle format
    * `--merged_before_score`, path to merged_before_score table, must be in hdf5 format,\
    with the key "merged_before_score"
