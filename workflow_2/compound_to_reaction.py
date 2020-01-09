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
import json

sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import workflow_helpers as mg

## Global parameters. None if not needed
precomputed_reactions = {}

### Read and prepare data
def read_retro_rules(min_diameter=0):
    """
    Read retro rules database and:
    - pre-convert reactions and substrates to rdkit objects
    - remove all entries below the minimum diameter
    - group by reaction ID
    """
    # TODO: get this from local settings file
    Retro_rules_db_path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "database/compound_to_reaction.csv")
    if not os.path.exists(Retro_rules_db_path):
        #TODO: move to argparse?
        sys.exit("retro rules path is wrong")
    retro_rules = pd.read_csv(Retro_rules_db_path)
    retro_rules = retro_rules[retro_rules["Diameter"] >= min_diameter]
    ## Get all reaction objects
    rule_smarts_all = list(set(retro_rules["Rule_SMARTS"]))
    rule_smarts_dict = {}
    for smarts in rule_smarts_all:
        reaction = AllChem.ReactionFromSmarts(smarts)
        rule_smarts_dict[smarts] = reaction
    def reaction_lookup(smarts):
        return rule_smarts_dict[smarts]
    retro_rules["Reaction"] = retro_rules["Rule_SMARTS"].apply(reaction_lookup)

    # Get all substrate objects
    substrate_smiles_all = list(set(retro_rules["Substrate_SMILES"]))
    substrate_smiles_dict = {}
    for smiles in substrate_smiles_all:
        mol = mol_from_smiles(smiles)
        substrate_smiles_dict[smiles] = mol
    def substrate_lookup(smiles):
        return substrate_smiles_dict[smiles]
    retro_rules["Substrate"] = retro_rules["Substrate_SMILES"].apply(substrate_lookup)
    retro_rules = retro_rules.groupby(by=["Reaction_ID", "Substrate_SMILES"])
    return retro_rules

def load_precomputed_reactions(precomputed_reaction_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"database/precomputed_compound_to_reactions.txt")):
    """Load the precomputed compound to reaction dictionary
    Input is the path to the precomputed reactions file
    Object is stored in global variable precomputed_reactions"""
    if not os.path.exists(precomputed_reaction_path):
        sys.exit("precomputed reaction path is wrong")
    with open(precomputed_reaction_path, "r") as json_file:
        global precomputed_reactions 
        precomputed_reactions = json.load(json_file)

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

def prepare_smiles(smiles):
    """
    Removes salts from molecule
    Turns tautomers in canonical form
    Neutralizes charged molecules
    Add hydrogens
    """
    saltremover = SaltRemover.SaltRemover()
    mol = saltremover.StripMol(Chem.MolFromSmiles(smiles), dontRemoveEverything=True)
    mol = canonicalize_tautomer(mol)
    mol = NeutraliseCharges(mol)
    mol = Chem.AddHs(Chem.RemoveHs(mol))
    smiles = Chem.MolToSmiles(mol)
    return smiles

def read_compounds_data(compounds_path):
    if not os.path.exists(compounds_path):
        sys.exit("Compounds_path is wrong")
    compounds_data = pd.read_csv(compounds_path)
    compounds_data["SMILES"] = compounds_data["original_compound"].apply(prepare_smiles)
    # TODO: allow other input?
    return compounds_data

def mol_from_smiles(smiles):
    """
    Create Rdkit.Chem.Mol object from SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(Chem.RemoveHs(mol))
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
    results = reaction.RunReactants((molecule,))
    if len(results) > 0:
        return True
    else:
        return False

def lookup_precomputed_reactions(molecule_smiles, min_diameter, similarity_cutoff=0.6):
    """This function returns precomputed reactions for a compound SMILES.
    Compounds were precomputed with a minimum similarity cutoff of 0.6 and no minimum diameter.
    The fingerprint radius was 3. 

    Input is a MAGI-prepared SMILES, the minimum diameter for a reaction and a minimum similarity cutoff.
    Results with a lower cutoff than 0.6 cannot be obtained, because they were not precomputed.
    Output is a dataframe with five columns:
    ["Molecule_SMILES", "Substrate_SMILES", "Reaction_ID", "Similarity", "Diameter"]
    """
    if molecule_smiles in precomputed_reactions.keys():
        c2r_data = precomputed_reactions[molecule_smiles]
        c2r_data = pd.DataFrame(c2r_data, columns = ["Molecule_SMILES", "Substrate_SMILES", "Reaction_ID", "Similarity", "Diameter"])
        # Filter on min diameter and similarity cutoff
        c2r_data = c2r_data[c2r_data["Diameter"] >= min_diameter]
        c2r_data = c2r_data[c2r_data["Similarity"] >= similarity_cutoff]
        return c2r_data
    else:
        return None

### Run C2R
def compound_to_reaction(molecule_smiles, rules_to_use, min_diameter, c2r_output_file=None, fingerprint_radius=3, similarity_cutoff=0.6, use_precomputed = True):
    """
    Find reactions in which a molecule can be used as a substrate
    - For the lowest retro rules diameters, find all reactions in which the compound can be used.
    - For those reactions, check if the reaction also uses the molecule with a higher diameter.
    - Write compound to reaction to intermediate file
    - Return compound to reactions
    """
    # Lookup molecule in precomputed compounds
    if use_precomputed:
        compound_to_reaction = lookup_precomputed_reactions(molecule_smiles = molecule_smiles, 
                                                            min_diameter = min_diameter, 
                                                            similarity_cutoff=similarity_cutoff)
        if compound_to_reaction is not None:
            #TODO: store this in intermediate file too?
            return compound_to_reaction
    # If molecule is not known, perform compound to reaction search. 
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
    # Parse arguments and prepare MAGI run
    start_time = datetime.datetime.now()
    print("Starting compound to reaction search at {}".format(str(start_time)))
    magi_parameters = mg.general_magi_preparation()

    # Read data
    unknown_compounds_present = True
    if unknown_compounds_present:
        print("loading retro rules database")
        sys.stdout.flush()
        retro_rules = read_retro_rules(magi_parameters["diameter"])
    print("loading compounds")
    compounds_data = read_compounds_data(magi_parameters["compounds"])
    sys.stdout.flush()
    print("loading precomputed reactions")
    load_precomputed_reactions()
    sys.stdout.flush()

    # Run compound to reaction search
    compound_to_reaction_total = []
    print("starting c2r searches")
    c2r_output_file = os.path.join(magi_parameters["intermediate_files_dir"], "compound_to_reaction.csv")
    pd.DataFrame(data = None, columns = ["Molecule_SMILES", "Substrate_SMILES", "Reaction_ID", "Similarity", "Diameter"]).to_csv(c2r_output_file, index = False)  
    sys.stdout.flush()
    for index, row in compounds_data[["SMILES", "original_compound"]].iterrows():
        molecule_smiles = row["SMILES"]
        original_compound = row["original_compound"]
        print("Starting with {} at {}".format(original_compound, str(datetime.datetime.now())))
        sys.stdout.flush()
        c2r = compound_to_reaction(molecule_smiles = molecule_smiles, 
                                    original_compound = original_compound,
                                    rules_to_use = retro_rules, 
                                    min_diameter = magi_parameters["diameter"], 
                                    c2r_output_file = c2r_output_file, 
                                    fingerprint_radius = magi_parameters["fingerprint"], 
                                    similarity_cutoff = magi_parameters["similarity_cutoff"], 
                                    use_precomputed=magi_parameters["use_precomputed_reactions"])
        compound_to_reaction_total.append(c2r)
    compound_to_reaction_total = pd.concat(compound_to_reaction_total)
    # Store data
    print("storing final data")
    #TODO: also store compound scores?
    compound_to_reaction_total.to_csv(os.path.join(magi_parameters["output"], "compound_to_reaction_total.csv"), index=False)
    print("Finished in {}".format(str(datetime.datetime.now() - start_time)))
    print("Done with compound to reaction search at {}".format(str(datetime.datetime.now())))

if __name__ == "__main__":
    main()