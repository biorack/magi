Metabolite Atlas Reaction Database
==================================

Metagenomics and single-cell sequencing have enabled, for the first time, glimpses into the vast metabolic potential of Earth’s collective biological systems.  Yet, for the most part we can’t accurately predict nor identify the products of most biosynthetic pathways. Most of what we know of microbial biochemistry is based on characterization of a few model microorganisms, and these findings have been extended through sequence correlations to the rest of sequence space. Unfortunately, these extrapolations have questionable validity for the vast majority of environmental microbes and therefore requires fundamentally different approaches for directly linking novel sequences to their biochemical functions.

Our vision is to systematically explore dark-biochemistry, using state-of-the-art workflows that integrate large-scale DNA synthesis with metabolomics, high-performance computing, and chemoinformatics.  Bioinformatics mining of the over 500 billion unique genes catalogued by the DOE Joint Genome Institute can be used to prioritize high-novelty candidate biosynthesis clusters. Through synthetic biology approaches candidate clusters can be refactored and expressed in model organisms for characterization of the resulting biochemical activities and products with mass spectrometry. When integrated with novel chemoinformatic algorithms, this creates a closed-loop cycle of design, build, test, and learn for systematically mapping biochemical space.  


Features
--------
- A database of all publicly available reactions
- Compounds in those reactions
- Atom mapping of products and reactants
- Comparison of all reactions
- API for accessing this information programmaticallly

Local Installation
------------------

.. code-block:: bash

    $ git clone https://<username>@bitbucket.org/bpbowen/metatlas_reactions.git
    $ cd metatlas_reactions
    $ pip install .

Local Settings
---------------
Local settings folder should have 3 files:
- local_settings.py
- user_settings.py
- __init__.py

local_settings.py should just have one line in it describing the name of the user_settings.py file:
.. code-block:: python
	$ SETTINGS_FILE = 'user_settings'

user_settings.py should have the following paths and variables defined:
.. code-block:: python
	blastbin = '' # path to BLAST binary

 refseq_path = 		'' # path to reaction reference sequence library
 refseq_db = 		'' # path to BLAST database for reference sequence library
 mrs_reaction_path = '' # path to metabolite-reaction-refseq database
 compounds_df = 		'' # path to compounds database
 mst_path = 			'' # path to chemical similarity network graph
 chemnet_pickle = 	'' # path to chemical similarity network descriptions

 magi_results_storage = 	'' # path to where to store magi outputs and blast databases
 repo_location = 		'' # path to repo location
	# The next 2 lines are only required if you are interfacing with magi_web repo 
 magiwebsuperuser = '' # admin username for magi_web
 magiwebsuperuserpass = '' # admin password for magi_web
