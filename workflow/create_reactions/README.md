# Creation of the reference reactions, sequences, and compounds from original sources

First, clone the repo: https://github.com/biorack/md2st and create a standardized_molecules.csv table that has both original structures and a standardized representation.  The original structure is what you would find in a chemical reaction and the standardized representation is necessary for creating a service to link out to other databases.  Due to tautomers and charged forms of molecules, each thing we think of as a metabolite, can exist in 100s of discrete structural representations.

Next, the following codes will create a reaction json file for each reaction in metacyc and rhea.

1. create_rhea_sequences.py
2. create_biocyc_sequences.py
3. create_rhea_reactions.py
4. create_biocyc_reactions.py

