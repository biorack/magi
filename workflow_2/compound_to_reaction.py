### Test for new compound to reaction search
import os
import sys
import argparse
import datetime
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import SaltRemover
from molvs import tautomer

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
            help="Minimum diameter to use for retro rules reactions", type = int, default = 10) # TODO: add check to be in range and even
        parser.add_argument('--fingerprint', 
            help="fingerprint radius for Morgan molecular fingerprint", type = int, default = 3)
        parser.add_argument('--similarity_cutoff', 
            help="Minimum similarity cutoff", type = float, default = 0.6)
        parser.add_argument('--intermediate_files',
            help='What directory within --output to store intermediate files',
            type=str, default='intermediate_files')
        args = parser.parse_args()
    except argparse.ArgumentTypeError as ex:
        sys.exit(ex.message)
    return args

### Read and prepare data
def read_retro_rules(Retro_rules_db_path, min_diameter):
    """
    Read retro rules database and:
    - pre-convert reactions and substrates to rdkit objects
    - remove all entries below the minimum diameter
    - group by reaction ID
    """
    if not os.path.exists(Retro_rules_db_path):
        #TODO: move to argparse?
        sys.exit("retro rules path is wrong")
    retro_rules = pd.read_csv(Retro_rules_db_path,sep='\t')
    retro_rules = retro_rules[retro_rules["Diameter"] >= min_diameter]
    retro_rules["Reaction"] = retro_rules["Rule_SMARTS"].apply(AllChem.ReactionFromSmarts)
    #TODO: Do this only once and never again
    #retro_rules["Substrate_SMILES"] = retro_rules["Substrate_SMILES"].apply(desalt_canonicalize_and_neutralize_smiles)
    retro_rules["Substrate"] = retro_rules["Substrate_SMILES"].apply(mol_from_smiles)
    retro_rules = retro_rules.groupby(by=["Reaction_ID", "Substrate_SMILES"])
    return retro_rules

""" contribution from Hans de Winter on https://www.rdkit.org/docs/Cookbook.html"""
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y)) for x,y in patts]

_reactions=None
def NeutraliseCharges(mol, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return mol

def canonicalize_tautomer(mol):
    """Returns the canonicalized form of a tautomer"""
    canon = tautomer.TautomerCanonicalizer()
    mol = canon.canonicalize(mol)
    return mol

def desalt_canonicalize_and_neutralize_smiles(smiles):
    """
    Removes salts from molecule
    Turns tautomers in canonical form
    Neutralizes charged molecules
    """
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(Chem.MolFromSmiles(smiles), dontRemoveEverything=True)
    mol = canonicalize_tautomer(mol)
    smiles = Chem.MolToSmiles(NeutraliseCharges(mol))
    return smiles

def read_compounds_data(compounds_path):
    if not os.path.exists(compounds_path):
        sys.exit("Compounds_path is wrong")
    compounds_data = pd.read_csv(compounds_path)
    compounds_data["SMILES"] = compounds_data["original_compound"].apply(desalt_canonicalize_and_neutralize_smiles)
    # TODO: allow other input?
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

def mol_matches_reaction(molecule, reaction):
    """
    Return true if the molecule can be used in the reaction.
    molecule: Rdkit.Chem.Mol object
    reaction: Rdkit.Chem.Reaction object
    """
    results = [Chem.MolToSmiles(x[0]) for x in reaction.RunReactants((molecule,))]
    if len(results) > 0:
        return True
    else:
        return False

### Run C2R
def compound_to_reaction(molecule_smiles, rules_to_use, min_diameter, c2r_output_file, fingerprint_radius=3, similarity_cutoff=0.6):
    """
    Find reactions in which a molecule can be used as a substrate
    - For the lowest retro rules diameters, find all reactions in which the compound can be used.
    - For those reactions, check if the reaction also uses the molecule with a higher diameter.
    - Write compound to reaction to intermediate file
    - Return compound to reactions
    """
    molecule = mol_from_smiles(molecule_smiles)
    compound_to_reaction = []
    
    # Find reactions per group of reaction-substrates
    for (reaction_id, substrate_smiles), rules_df in rules_to_use:
        #calculate fingerprint similarity for the substrate of the reaction and the molecule of interest
        substrate = rules_df.iloc[0]["Substrate"]
        similarity = calculate_fingerprint_similarity(substrate, molecule, fingerprint_radius)
        
        if similarity >= similarity_cutoff:
            diameters_to_check = sorted(rules_df["Diameter"].unique())
            # Assume that > 0.99 similarity is the same molecule, so substrate matches molecule
            if similarity > 0.99:
                reaction_matched = True
                diameter_to_store = diameters_to_check[-1]
            # else, check if compound matches lowest reaction
            else:
                reaction = rules_df[rules_df["Diameter"] == diameters_to_check[0]]["Reaction"].values[0]
                reaction_matched = mol_matches_reaction(molecule, reaction)
                # if yes, check consecutive reactions
                if reaction_matched:
                    diameter_to_store = diameters_to_check[0]
                    for diameter in diameters_to_check[-1:0:-1]:
                        if diameter > diameter_to_store:
                            reaction = rules_df[rules_df["Diameter"] == diameter]["Reaction"].values[0]
                            reaction_matched = mol_matches_reaction(molecule, reaction)
                            if reaction_matched:
                                diameter_to_store = diameter
            # store reaction with highest matching diameter
            if reaction_matched:          
                compound_to_reaction.append([molecule_smiles, substrate_smiles,reaction_id, similarity,diameter_to_store])

    # Store results and return data as dataframe
    if len(compound_to_reaction) > 0:
        # TODO: filter for highest diameter per reaction
        compound_to_reaction = pd.DataFrame(compound_to_reaction)
        compound_to_reaction.columns = ["Molecule_SMILES", "Substrate_SMILES", "Reaction_ID", "Similarity", "Diameter"]
        # Store data
        # TODO: find way to also store header
        compound_to_reaction.to_csv(c2r_output_file, mode = 'a', header=False, index=False)
        return compound_to_reaction
    else:
        print("No reactions found for {}".format(molecule_smiles))
        return None

def main():
    # Parse arguments and read input
    start_time = datetime.datetime.now()
    print("Starting compound to reaction search at {}".format(str(start_time)))
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
        c2r = compound_to_reaction(molecule_smiles = molecule_smiles, 
                                    rules_to_use = retro_rules, 
                                    min_diameter = args.diameter, 
                                    c2r_output_file = os.path.join(args.output, args.intermediate_files, "compound_to_reaction.csv"), 
                                    fingerprint_radius = args.fingerprint, 
                                    similarity_cutoff = args.similarity_cutoff)
        compound_to_reaction_total.append(c2r)
    compound_to_reaction_total = pd.concat(compound_to_reaction_total)
    # Store data
    print("storing final data")
    compound_to_reaction_total.to_csv(os.path.join(args.output, "compound_to_reaction_total.csv"), index=False)
    print("Finished in {}".format(str(datetime.datetime.now() - start_time)))
    print("Done with compound to reaction search at {}".format(str(datetime.datetime.now())))

if __name__ == "__main__":
    main()