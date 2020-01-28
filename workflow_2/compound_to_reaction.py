### Test for new compound to reaction search
import os
import sys
import argparse
import datetime
import pandas as pd
import multiprocessing as mp

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import inchi
from rdkit.Chem import SaltRemover
from molvs import tautomer

import sqlite3
import json

sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import workflow_helpers as mg

## Global parameters.
precomputed_reactions = {}
unknown_compounds_present = False # Find a way to only open the retro rules database if there are unknown compounds present
retro_rules = "Load the retro rules database with the read_retro_rules function" #TODO: think if I want this to be global
Retro_rules_reactions = "Load the retro rules database with the read_retro_rules_reactions_from_db function"
Retro_rules_substrates = "Load the retro rules database with the read_retro_rules_substrates_from_db function"
compounds_metadata = None # To store metadata from the compounds. (TODO: Write and read to file to conserve memory?)
useHs = True #Temp, to fix problems with explicit and implicit hydrogens in rdkit

def reaction_from_smarts(smarts):
    """
    This function creates a Rdkit reaction object from a SMARTS string
    """
    reaction = AllChem.ReactionFromSmarts(smarts)
    if reaction is None:
        # TODO: exit or move mol to log unsearched compounds?
        sys.exit("{} cannot be converted to rdkit reaction object".format(smarts))
    return reaction
def reaction_from_binary(binary_object):
    """Convert rdkit reaction binary to rdkit reaction"""
    reaction = Chem.rdChemReactions.ChemicalReaction(binary_object)
    return reaction
def mol_from_binary(binary_object):
    """Convert rdkit reaction mol to rdkit mol"""
    mol = Chem.Mol(binary_object)
    return mol

def lookup_similar_substrates(molecule, retro_rules_substrates=None, similarity_cutoff = 0.6, fingerprint_radius = 3):
    """Compare the molecule of interest to the database with retro rules substrates and 
    return a dictionary of substrate_IDs of similar substrates. 
    The key is a substrate ID and the value is the similarity score."""
    matching_substrate_ids = {}
    if type(retro_rules_substrates) == str:
        sys.exit(retro_rules_substrates)
    if type(retro_rules_substrates) == None:
        retro_rules_substrates = Retro_rules_substrates

    # Check for each substrate if it is similar to the molecule of interest
    for index, substrate in retro_rules_substrates.iterrows():
        similarity = calculate_fingerprint_similarity(substrate["substrate_rdkit_object"], molecule, fingerprint_radius)
        if similarity > similarity_cutoff:
            matching_substrate_ids[substrate["substrate_ID"]] = similarity
    return matching_substrate_ids

def read_retro_rules_from_db(path_to_database="/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/database/C2R_test_database.db", min_diameter=0, use_mp = False):
    """
    Read retro rules database and:
    - pre-convert reactions and substrates to rdkit objects
    - remove all entries below the minimum diameter
    - group by reaction ID
    """
    with sqlite3.connect(path_to_database) as connection:
        # Read retro rules from database
        global Retro_rules_reactions
        global retro_rules
        query = "SELECT * FROM Retro_rules_reactions \
                WHERE Retro_rules_reactions.diameter >= {}".format(min_diameter)
        Retro_rules_reactions = pd.read_sql_query(query, connection)
        # Turn reactions and substrates into mols
        if use_mp:
            cpu_count = mp.cpu_count() // 2 # use half of available CPUs
            print("using {} cpus".format(cpu_count))
            sys.stdout.flush()
            try: 
                p = mp.Pool(cpu_count)
                Retro_rules_reactions["reaction_rdkit_object"] = p.map(reaction_from_binary, Retro_rules_reactions["reaction_rdkit_object"])
            finally:
                p.close()
                p.terminate()
        else:
            Retro_rules_reactions["reaction_rdkit_object"] = Retro_rules_reactions["reaction_rdkit_object"].apply(reaction_from_binary)
        # Group by reaction ID and Canonical_SMILES
        retro_rules = Retro_rules_reactions
        Retro_rules_reactions = Retro_rules_reactions.groupby(by="retro_rules_ID")

### Read and prepare data
def read_retro_rules_substrates_from_db(path_to_database="/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/database/C2R_test_database.db"):
    """
    Read retro rules substrates database and
    pre-convert substrates to rdkit objects
    """
    # Read retro rules from database
    with sqlite3.connect(path_to_database) as connection:
        global Retro_rules_substrates
        connection = sqlite3.connect(path_to_database)
        query = "SELECT * FROM Retro_rules_substrates"
        Retro_rules_substrates = pd.read_sql_query(query, connection)
        # Turn reactions and substrates into mols
        Retro_rules_substrates["substrate_rdkit_object"] = Retro_rules_substrates["substrate_rdkit_object"].apply(mol_from_binary)

def read_retro_rules(min_diameter=0, useHs = useHs):
    """
    Read retro rules database and:
    - pre-convert reactions and substrates to rdkit objects
    - remove all entries below the minimum diameter
    - group by reaction ID
    """
    global retro_rules
    # TODO: get this from local settings file
    if useHs:
        Retro_rules_db_path="/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/retro_rules_with_canonical_and_inchi.csv"
    else:
        Retro_rules_db_path="/Users/northenlab/Desktop/h_leegwater/internship/magi_2_database/retro_rules_with_canonical_and_inchi.csv"
    if not os.path.exists(Retro_rules_db_path):
        #TODO: move to argparse?
        sys.exit("retro rules path is wrong")
    retro_rules = pd.read_csv(Retro_rules_db_path)
    retro_rules = retro_rules[retro_rules["Diameter"] >= min_diameter]
    ## Get all reaction objects
    rule_smarts_all = list(set(retro_rules["Rule_SMARTS"]))
    rule_smarts_dict = {}
    for smarts in rule_smarts_all:
        reaction = reaction_from_smarts(smarts)
        rule_smarts_dict[smarts] = reaction
    def reaction_lookup(smarts):
        return rule_smarts_dict[smarts]
    retro_rules["Reaction"] = retro_rules["Rule_SMARTS"].apply(reaction_lookup)

    # Get all substrate objects
    substrate_smiles_all = list(set(retro_rules["Canonical_SMILES"]))
    substrate_smiles_dict = {}
    for smiles in substrate_smiles_all:
        mol = mol_from_smiles(smiles)
        substrate_smiles_dict[smiles] = mol
    def substrate_lookup(smiles):
        return substrate_smiles_dict[smiles]
    retro_rules["Substrate"] = retro_rules["Canonical_SMILES"].apply(substrate_lookup)
    retro_rules = retro_rules.groupby(by=["Reaction_ID", "Canonical_SMILES"])

def load_precomputed_reactions(precomputed_reaction_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"database/precomputed_c2r_new.txt")):
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
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

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

def remove_stereochemistry(mol):
    """
    Removes stereochemistry centers from a molecule"""
    Chem.RemoveStereochemistry(mol)
    return mol

def prepare_smiles(smiles, useHs = useHs):
    """
    Remove stereochemistry
    Removes salts from molecule
    Turns tautomers in canonical form
    Neutralizes charged molecules
    Add hydrogens
    """
    saltremover = SaltRemover.SaltRemover()
    mol = saltremover.StripMol(Chem.MolFromSmiles(smiles), dontRemoveEverything=True)
    mol = NeutraliseCharges(mol)
    mol = remove_stereochemistry(mol)
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol)) # To fix kekulization bug
    mol = canonicalize_tautomer(mol)
    if useHs:
        mol = Chem.AddHs(Chem.RemoveHs(mol))
    else:
        mol = Chem.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol)
    return smiles

def read_compounds_data(compounds_path):
    """
    This function will read an input file with candidate compounds
    Output is a pandas DataFrame with a column MAGI-prepared SMILES and a compound score
    """
    global compounds_metadata
    if not os.path.exists(compounds_path):
        sys.exit("Compounds_path is wrong")
    compounds_data = mg.load_compound_results(compounds_path)
    compounds_metadata = compounds_data.copy()
    compounds_data["SMILES"] = compounds_data["original_compound"].apply(prepare_smiles)
    compounds_data = compounds_data[["original_compound", "SMILES"]]
    return compounds_data

def mol_from_smiles(smiles, useHs = useHs):
    """
    Create Rdkit.Chem.Mol object from SMILES and add hydrogens
    """
    mol = Chem.MolFromSmiles(smiles)
    if useHs:
        mol = Chem.AddHs(Chem.RemoveHs(mol))
    else:
        mol = Chem.RemoveHs(mol)
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

def lookup_precomputed_reactions(molecule_inchikey, min_diameter, similarity_cutoff=0.6):
    """
    This function returns precomputed reactions for a compound SMILES.
    Compounds were precomputed with a minimum similarity cutoff of 0.6 and no minimum diameter.
    The fingerprint radius was 3. 

    Input is a MAGI-prepared SMILES, the minimum diameter for a reaction and a minimum similarity cutoff.
    Results with a lower cutoff than 0.6 cannot be obtained, because they were not precomputed.
    Output is a dataframe with five columns:
    TODO: do this
    ["Molecule_SMILES", "Substrate_SMILES", "Reaction_ID", "Similarity", "Diameter"]
    """
    if molecule_inchikey in precomputed_reactions.keys():
        c2r_data = precomputed_reactions[molecule_inchikey]
        c2r_data = pd.DataFrame(c2r_data, columns = ["Molecule_SMILES", "Substrate_SMILES", "Reaction_ID", "Similarity", "Diameter"])
        # Filter on min diameter and similarity cutoff
        c2r_data = c2r_data[c2r_data["Diameter"] >= min_diameter]
        c2r_data = c2r_data[c2r_data["Similarity"] >= similarity_cutoff]
        return c2r_data
    else:
        return None

def get_inchi_key(molecule):
    return inchi.MolToInchiKey(molecule)
### Run C2R
def compound_to_reaction(molecule_smiles, original_compound, rules_to_use=Retro_rules_reactions, substrates_to_use = Retro_rules_substrates, min_diameter=0, c2r_output_file=None, fingerprint_radius=3, similarity_cutoff=0.6, use_precomputed = True):
    """
    Find reactions in which a molecule can be used as a substrate
    - For the lowest retro rules diameters, find all reactions in which the compound can be used.
    - For those reactions, check if the reaction also uses the molecule with a higher diameter.
    - Write compound to reaction to intermediate file
    - Return compound to reactions
    """
    # Generate inchi key for the molecule
    molecule = mol_from_smiles(molecule_smiles)
    molecule_inchikey = get_inchi_key(molecule)
    # Lookup molecule in precomputed compounds
    if use_precomputed:
        compound_to_reaction = lookup_precomputed_reactions(molecule_inchikey = molecule_inchikey, 
                                                            min_diameter = min_diameter, 
                                                            similarity_cutoff=similarity_cutoff)
        if compound_to_reaction is not None:
            #TODO: store this in intermediate file too?
            return compound_to_reaction
    # If molecule is not known, perform compound to reaction search on the user's original compound, not on the MAGI-cleaned molecule smiles. 
    # Risk: problems with tautomers. The Substrate object in retro rules is a canonicalized version of the original compound 
    compound_to_reaction = []
    molecule = mol_from_smiles(original_compound)

    # Find substrates to which the molecule of interest matches
    matching_substrates = lookup_similar_substrates(molecule, retro_rules_substrates=substrates_to_use, similarity_cutoff = similarity_cutoff, fingerprint_radius = fingerprint_radius)
    # Find reactions per group of reaction-substrates
    for (reaction_id, rules_df) in rules_to_use:
        #calculate fingerprint similarity for the substrate of the reaction and the molecule of interest
        #substrate = rules_df.iloc[0]["Substrate"]
        substrate_ID = rules_df.iloc[0]["substrate_ID"]
        if substrate_ID in matching_substrates.keys():
            similarity = matching_substrates[substrate_ID]
        else:
            similarity = 0
        
        if similarity >= similarity_cutoff:
            rules_df = rules_df[rules_df["substrate_ID"]==substrate_ID]
            diameters_to_check = sorted(rules_df["diameter"].unique())
            # Assume that > 0.99 similarity is the same molecule, so substrate matches molecule
            if similarity > 0.99:
                reaction_matched = True
                diameter_to_store = max(diameters_to_check)
            # else, check if compound matches lowest reaction
            else:
                reaction = rules_df[rules_df["diameter"] == diameters_to_check[0]]["reaction_rdkit_object"].values[0]
                reaction_matched = mol_matches_reaction(molecule, reaction)
                # if yes, check consecutive reactions
                if reaction_matched:
                    diameter_to_store = diameters_to_check[0]
                    for diameter in diameters_to_check[-1:0:-1]:
                        if diameter > diameter_to_store:
                            reaction = rules_df[rules_df["diameter"] == diameter]["reaction_rdkit_object"].values[0]
                            reaction_matched = mol_matches_reaction(molecule, reaction)
                            if reaction_matched:
                                diameter_to_store = diameter
            # store reaction with highest matching diameter
            if reaction_matched:          
                compound_to_reaction.append([original_compound, rules_df.iloc[0]["reaction_ID"], similarity])

    # Store results and return data as dataframe
    if len(compound_to_reaction) > 0:
        # TODO: filter for highest diameter per reaction
        compound_to_reaction = pd.DataFrame(compound_to_reaction)
        compound_to_reaction.columns = ["original_compound", "reaction_ID", "Similarity"]
        # Store data
        # TODO: find way to also store header
        compound_to_reaction.to_csv(c2r_output_file, mode = 'a', header=False, index=False)
        return compound_to_reaction
    else:
        print("No reactions found for {}".format(molecule_smiles))
        return None

def compound_to_reaction_scoring(compound, reaction, compound_score, diameter, fingerprint_similarity):
    """This function uses the fingerprint similarity and the reaction diameter to score 
    the likelyhood of a reaction occurring."""
    diameter_max = 16.0 # maximum diameter in the retro rules database
    diameter_score = float(diameter) / diameter_max
    fingerprint_similarity_score = fingerprint_similarity # already between 0 and 1
    compound_score = compound_score # original compound score
    c2r_score = (diameter_score + fingerprint_similarity_score + compound_score) / 3

def format_output(compound_to_reaction_total):
    """
    Merge metadata to the compound to reaction table and return full data frame

    Inputs
    ------
    A data frame with 
    - a compound smiles
    - the reaction ID to which it is matched
    - the similarity score between the compound and the original substrate of the reaction

    Outputs
    ------
    Additional info:
    - the retro rules substrate and reaction IDs
    - the diameter of the reaction
    - input compound metadata?
    """
    ## Merge with reactions
    compound_to_reaction_total = compound_to_reaction_total.merge(
        retro_rules[["reaction_ID","substrate_ID","retro_rules_ID", "diameter"]],
                                how="left", on="reaction_ID")
    compound_to_reaction_total.rename({"retro_rules_ID":"retro_rules_reaction_ID"}, axis=1, inplace=True) 
    ## Merge with substrates
    compound_to_reaction_total = compound_to_reaction_total.merge(
        Retro_rules_substrates[['substrate_ID', 'retro_rules_ID']],
                                how="left", on="substrate_ID")
    compound_to_reaction_total.rename({"retro_rules_ID":"retro_rules_substrate_ID", "inchi_key":"substrate_inchi_key"}, axis=1, inplace=True) 
    ## Merge with initial metadata
    compound_to_reaction_total = compound_to_reaction_total.merge(
            compounds_metadata, how="left", on="original_compound")
    return compound_to_reaction_total

def main():
    global retro_rules
    # Parse arguments and prepare MAGI run
    start_time = datetime.datetime.now()
    print("Starting compound to reaction search at {}".format(str(start_time)))
    magi_parameters = mg.general_magi_preparation()

    # Read data
    unknown_compounds_present = True #TODO: build a check for this? Or move it? Only load full RR database when needed?
    if unknown_compounds_present:
        print("!!! Loading retro rules database")
        sys.stdout.flush()
        read_retro_rules_from_db(min_diameter = magi_parameters["diameter"])
        read_retro_rules_substrates_from_db()
    print("!!! Loading compounds")
    compounds_data = read_compounds_data(magi_parameters["compounds"])
    sys.stdout.flush()
    if magi_parameters["use_precomputed_reactions"]:
        print("!!! Loading precomputed reactions")
        load_precomputed_reactions()
    else:
        print("!!! Not using precomputed reactions.")
    sys.stdout.flush()

    # Run compound to reaction search
    compound_to_reaction_total = []
    print("!!! Starting c2r searches")
    c2r_output_file = os.path.join(magi_parameters["intermediate_files_dir"], "compound_to_reaction.csv")
    pd.DataFrame(data = None, columns = ["Molecule_SMILES", "reaction_ID", "Similarity"]).to_csv(c2r_output_file, index = False)  
    sys.stdout.flush()
    for index, row in compounds_data[["SMILES", "original_compound"]].iterrows():
        molecule_smiles = row["SMILES"]
        original_compound = row["original_compound"]
        print("!!! Starting with {} at {}".format(original_compound, str(datetime.datetime.now())))
        sys.stdout.flush()
        c2r = compound_to_reaction(molecule_smiles = molecule_smiles, 
                                    original_compound = original_compound,
                                    rules_to_use = Retro_rules_reactions, 
                                    substrates_to_use = Retro_rules_substrates,
                                    min_diameter = magi_parameters["diameter"], 
                                    c2r_output_file = c2r_output_file, 
                                    fingerprint_radius = magi_parameters["fingerprint"], 
                                    similarity_cutoff = magi_parameters["similarity_cutoff"], 
                                    use_precomputed=magi_parameters["use_precomputed_reactions"])
        compound_to_reaction_total.append(c2r)
    compound_to_reaction_total = pd.concat(compound_to_reaction_total)

    # Maybe do some cool scoring here?

    # Format output
    compound_to_reaction_total = format_output(compound_to_reaction_total)
    # Store data
    outputfilename = os.path.join(magi_parameters["output"], "compound_to_reaction_total.csv")
    print("!!! Saving compound to reaction results to {}".format(outputfilename))
    #TODO: also store compound scores?
    compound_to_reaction_total.to_csv(outputfilename, index=False)
    print("!!! Finished in {}".format(str(datetime.datetime.now() - start_time)))
    print("!@# Done with compound to reaction search at {}".format(str(datetime.datetime.now())))

if __name__ == "__main__":
    main()