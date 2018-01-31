# Creation of the reference reactions, sequences, and compounds from original sources

First, clone the repo: https://github.com/biorack/md2st and create a standardized_molecules.csv table that has both original structures and a standardized representation.  The original structure is what you would find in a chemical reaction and the standardized representation is necessary for creating a service to link out to other databases.  Due to tautomers and charged forms of molecules, each thing we think of as a metabolite, can exist in 100s of discrete structural representations.

Next, the following codes will create a reaction json file for each reaction in metacyc and rhea.

1. rhea_sequences.py
2. biocyc_sequences.py
3. rhea_reactions.py
4. biocyc_reactions.py

These codes are not meant to be run very often.  Care should be taken at each step to ensure that the necessary sequences and compounds are found.

For BioCyc reactions, stoichiometry is calculated using the ChemPy.  For about 1/3 of the reactions, the values are higher than you would expect.  Its likely that for some of these a lower stoichiometry exists that satisfies the balanced reaction.
