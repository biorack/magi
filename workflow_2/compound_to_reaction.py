### Test for new compound to reaction search
import os
import sys
import argparse
import datetime
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

import workflow_helpers as mg

def parse_arguments():
    """
    This is the MAGI argument parser that is used in all workflows. 
    It checks if all required arguments are passed, if numbers fall within MAGI-approved ranges and if files exist.

    Outputs
    -------
    An argparse.args object containing all arguments. 
    """
    try:
        """parse arguments"""
        parser = argparse.ArgumentParser()
        # required arguments
        parser.add_argument('-c', '--compounds', type=mg.is_existing_file,
            help='path to observed compounds file')
        parser.add_argument('-o', '--output', 
            help='path to a custom output', 
            type=str)
        parser.add_argument('--retro_rules', 
            help='path to retro rules database', # TODO: read this path from local settings file 
            type=str)
        parser.add_argument('--diameter', 
            help="diameter to use for retro rules reactions", type = int) # TODO: add check to be in range and even
        parser.add_argument('--intermediate_files',
            help='What directory within --output to store intermediate files',
            type=str, default='intermediate_files')
        args = parser.parse_args()
    except argparse.ArgumentTypeError as ex:
        sys.exit(ex.message)
    return args

### Read and prepare data
def read_retro_rules(Retro_rules_db_path, diameter):
    if not os.path.exists(Retro_rules_db_path):
        sys.exit("retro rules path is wrong")
    retro_rules = pd.read_csv(Retro_rules_db_path,sep='\t')
    retro_rules = retro_rules[retro_rules["Diameter"] == diameter]
    retro_rules["Reaction"] = retro_rules["Rule_SMARTS"].apply(AllChem.ReactionFromSmarts)
    retro_rules.reset_index(inplace=True,drop=True)
    return retro_rules

def read_compounds_data(compounds_path):
    if not os.path.exists(compounds_path):
        sys.exit("EMA_adjusted_path is wrong")
    compounds_data = pd.read_csv(compounds_path)
    return compounds_data

def mol_from_smiles(smiles):
    """
    Create Rdkit.Chem.Mol object from SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # TODO: exit or move mol to log unsearched compounds?
        sys.exit("{} cannot be converted to rdkit.Mol object".format(smiles))
    return mol

def calculate_fingerprint_similarity(mol1, mol2, fingerprint_radius = 3):
    """
    This function calculates the Morgan molecular fingerprint similarity for two molecules
    For full documentation, see https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
    
    Inputs
    -------
    mol1:   A rdkit.Chem.Mol object
    mol2:   A rdkit.Chem.Mol object
    fingerprint_radius: Fingerprint radius for generating Morgan fingerprints

    Outputs
    -------
    dice_similarity: The similarity between the two molecular fingerprints.
    TODO: check if useFeatures=True) needs to be added to getting the fingerprints.
    """

    fingerprint1 = AllChem.GetMorganFingerprint(mol1, fingerprint_radius)
    fingerprint2 = AllChem.GetMorganFingerprint(mol2, fingerprint_radius)
    dice_similarity = DataStructs.DiceSimilarity(fingerprint1, fingerprint2)
    return dice_similarity

### Run C2R
def compound_to_reaction(molecule_smiles, rules_to_use, c2r_output_file):
    molecule = mol_from_smiles(molecule_smiles)
    compound_to_reaction = []
    for index, (reaction, substrate_smiles, reaction_id) in rules_to_use[["Reaction", "Substrate_SMILES", "Reaction_ID"]].iterrows():
        results = [Chem.MolToSmiles(x[0]) for x in reaction.RunReactants((molecule,))]
        if len(results) > 0:
            compound_to_reaction.append([molecule_smiles, substrate_smiles,reaction_id])
    if len(compound_to_reaction) > 0:
        compound_to_reaction = pd.DataFrame(compound_to_reaction)
        compound_to_reaction.columns = ["Molecule_SMILES", "Substrate_SMILES", "Reaction_ID"]
        # Store data
        # TODO: find way to also store header
        compound_to_reaction.to_csv(c2r_output_file, mode = 'a', header=False)
        return compound_to_reaction
    else:
        print("No reactions found for {}".format(molecule_smiles))
        return None

def main():
    # Parse arguments and read input
    print("Starting compound to reaction search at {}".format(str(datetime.datetime.now())))
    args = parse_arguments()
    print("loading retro rules database")
    sys.stdout.flush()
    retro_rules = read_retro_rules(args.retro_rules, args.diameter)
    print("loading compounds")
    compounds_data = read_compounds_data(args.compounds)
    sys.stdout.flush()
    # Make output dir
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.exists(os.path.join(args.output, args.intermediate_files)):
        os.mkdir(os.path.join(args.output, args.intermediate_files))
    # Run compound to reaction search
    compound_to_reaction_total = []
    print("starting c2r searches")
    sys.stdout.flush()
    for molecule_smiles in compounds_data["SMILES"]:
        print("Starting with {} at {}".format(molecule_smiles, str(datetime.datetime.now())))
        sys.stdout.flush()
        c2r = compound_to_reaction(molecule_smiles, retro_rules, os.path.join(args.output, "compound_to_reaction.csv"))
        compound_to_reaction_total.append(c2r)
    compound_to_reaction_total = pd.concat(compound_to_reaction_total)
    # Store data
    print("storing final data")
    compound_to_reaction_total.to_csv(os.path.join(args.output, "compound_to_reaction_total.csv"))
    print("Done with compound to reaction search at {}".format(str(datetime.datetime.now())))

if __name__ == "__main__":
    main()