### Test for new compound to reaction search
import os
import sys
import argparse
import datetime
import pandas as pd
import multiprocessing as mp
from functools import partial

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
my_settings = mg.get_settings()
unknown_compounds_present = False # Find a way to only open the retro rules database if there are unknown compounds present
retro_rules = "Load the retro rules database with the read_retro_rules function" #TODO: think if I want this to be global
Retro_rules_reactions = "Load the retro rules database with the read_retro_rules_reactions_from_db function"
Retro_rules_substrates = "Load the retro rules database with the read_retro_rules_substrates_from_db function"
compounds_metadata = None # To store metadata from the compounds. (TODO: Write and read to file to conserve memory?)
useHs = False #Temp, to fix problems with explicit and implicit hydrogens in rdkit

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

def read_retro_rules_from_db(path_to_database=my_settings.magi_database, min_diameter=0, use_mp = False):
    """
    Read retro rules database and:
    - pre-convert reactions and substrates to rdkit objects
    - remove all entries below the minimum diameter
    - group by reaction ID
    """
    global Retro_rules_reactions
    global retro_rules

    with sqlite3.connect(path_to_database) as connection:
        # Read retro rules from database
        query = "SELECT * FROM Retro_rules_reactions \
                WHERE Retro_rules_reactions.diameter >= {}".format(min_diameter)
        retro_rules = pd.read_sql_query(query, connection)
        # Turn reactions and substrates into mols
        if use_mp:
            cpu_count = mp.cpu_count() // 2 # use half of available CPUs
            print("using {} cpus".format(cpu_count))
            sys.stdout.flush()
            try: 
                p = mp.Pool(cpu_count)
                retro_rules["retro_rules_smarts"] = p.map(reaction_from_smarts, retro_rules["retro_rules_smarts"])
            finally:
                p.close()
                p.terminate()
        else:
            retro_rules["retro_rules_smarts"] = retro_rules["retro_rules_smarts"].apply(reaction_from_smarts)
        retro_rules.rename(columns = {"retro_rules_smarts":"reaction_rdkit_object"}, inplace = True)
        # Group by reaction ID and Canonical_SMILES
        Retro_rules_reactions = retro_rules.groupby(by=["retro_rules_ID", "substrate_ID"])

### Read and prepare data
def read_retro_rules_substrates_from_db(path_to_database=my_settings.magi_database):
    """
    Read retro rules substrates database and
    pre-convert substrates to rdkit objects
    """
    global Retro_rules_substrates
    # Read retro rules from database
    with sqlite3.connect(path_to_database) as connection:
        connection = sqlite3.connect(path_to_database)
        query = "SELECT * FROM Retro_rules_substrates"
        Retro_rules_substrates = pd.read_sql_query(query, connection)
    
    # Turn reactions and substrates into mols
    Retro_rules_substrates["retro_rules_smiles"] = Retro_rules_substrates["retro_rules_smiles"].apply(mol_from_smiles)
    Retro_rules_substrates.rename(columns = {"retro_rules_smiles":"substrate_rdkit_object"}, inplace = True)

def _InitialiseNeutralisationReactions():
    """ contribution from Hans de Winter on https://www.rdkit.org/docs/Cookbook.html"""
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
    if mol is None:
        print("Warning! Molecule could not be processed: {}".format(smiles))
        return smiles, None
    mol = canonicalize_tautomer(mol)
    if useHs:
        mol = Chem.AddHs(Chem.RemoveHs(mol))
    else:
        mol = Chem.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol)
    return smiles, mol

def preprocess_compounds_data(compounds_path, cpu_count):
    """
    This function will read an input file with candidate compounds
    Output is a pandas DataFrame with :
    - the original compound SMILES
    - a MAGI-prepared SMILES 
    - an rdkit.mol object
    - an InChI Key of the MAGI-prepared SMILES 
    - a compound score
    """
    # Open compounds data and store metadata in data frame
    global compounds_metadata
    if not os.path.exists(compounds_path):
        sys.exit("Compounds_path is wrong")
    compounds_data = mg.load_compound_results(compounds_path)
    compounds_metadata = compounds_data.copy()

    # Prepare a MAGI SMILES
    # TODO: fix multiprocessing
    #if cpu_count > 1:
    #    try: 
    #        p = mp.Pool(cpu_count)
    #        compounds_data["SMILES", "Mol"] = p.map(prepare_smiles, compounds_data["original_compound"])
    #    finally:
    #        p.close()
    #        p.terminate()
    #else:
    compounds_data["SMILES"], compounds_data["Mol"] = zip(*compounds_data["original_compound"].apply(prepare_smiles))
    # Drop compounds that could not be prepared and have None as Mol object
    compounds_data = compounds_data[compounds_data["Mol"].notnull()]

    # Prepare a InChI Key
    compounds_data["inchi_key"] = compounds_data["Mol"].apply(get_inchi_key)

    compounds_data = compounds_data[["original_compound", "SMILES", "compound_score", "Mol", "inchi_key"]]
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
    """
    fingerprint1 = AllChem.GetMorganFingerprint(mol1, fingerprint_radius, useFeatures=True)
    fingerprint2 = AllChem.GetMorganFingerprint(mol2, fingerprint_radius, useFeatures=True)
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

def lookup_matching_reactions_for_one_compound(original_compound, compound_score, inchi_key, similarity_cutoff = 0.6, path_to_database=my_settings.magi_database):
    """
    This function looks up if reactions were precomputed for a given inchi key.
    Input
    original_compound and compound_score are metadata that will be merged to each reaction found
    inchi_key: the inchi key to look up
    similarity_cutoff: reactions for which the compound of interest and the original substrate are less similar than this cutoff, will not be considered

    Output:
    A data frame with reactions. Columns are "original_compound","reaction_ID","similarity" and "compound_score"
    None if no precomputed reactions are found. 
    """

    # Get molecule ID
    with sqlite3.connect(path_to_database) as connection:
        query = "SELECT molecule_ID FROM Precomputed_molecules WHERE inchi_key = '{}'".format(inchi_key)
        precomputed_molecule = pd.read_sql_query(query, connection)
        if precomputed_molecule.shape[0] > 1:
            # This should not happen, this means that the database itself has duplicates.
            print("WARNING, inchi key {} found multiple times in the Precomputed_molecules database.".format(inchi_key))
        if precomputed_molecule.shape[0] == 0:
            # No precomputed reactions found
            print("No precomputed reactions found for {}".format(original_compound))
            return None
        
        # Get reactions and format output data frame
        precomputed_molecule = precomputed_molecule.iloc[0,0]
        query = "SELECT reaction_ID, similarity FROM Precomputed_reactions \
        WHERE molecule_ID = '{}' AND \
        similarity > {}".format(precomputed_molecule, similarity_cutoff)
        reactions= pd.read_sql_query(query, connection)
    reactions["original_compound"] = original_compound
    reactions["compound_score"] = compound_score
    reactions = reactions[["original_compound","reaction_ID","similarity","compound_score"]]

    return reactions

def get_inchi_key(molecule):
    return inchi.MolToInchiKey(molecule)
### Run C2R
def find_precomputed_reactions(compounds_data, c2r_output_file, min_diameter = 12, similarity_cutoff =0.6):
    """
    This function will look up reactions from previous MAGI runs for all compounds.
    Input: 
    compounds_data: a compounds data frame with at least columns original_compound, compound_score and inchi_key
    min_diameter: the minimum diameter for reactions
    similarity_cutoff: a similarity cutoff used to find matching reactions

    Output:
    A pandas dataframe with precomputed reactions. Headers are "original_compound","reaction_ID","similarity" and "compound_score"
    A pandas dataframe with not precomputed reactions. This is a subset of the original data frame with compounds for which no reactions were precomputed.
    """
    precomputed_reactions = []
    not_precomputed_compounds = []
    # find reactions
    for row in compounds_data.iterrows():
        ix = row[0]; row = row[1]
        c2r_subset = lookup_matching_reactions_for_one_compound(
                        original_compound = row["original_compound"], 
                        compound_score = row["compound_score"], 
                        inchi_key = row["inchi_key"], 
                        similarity_cutoff = similarity_cutoff)
        if c2r_subset is not None:
            precomputed_reactions.append(c2r_subset)
            # Store data
            c2r_subset.to_csv(c2r_output_file, mode = 'a', header=False, index=False)
        else:
            not_precomputed_compounds.append(ix)
    not_precomputed_compounds = compounds_data.iloc[not_precomputed_compounds]
    # Change data type to single pandas dataframe
    if len(precomputed_reactions) > 0:
        precomputed_reactions = pd.concat(precomputed_reactions)
    elif len(precomputed_reactions) == 1:
        precomputed_reactions = precomputed_reactions[0]
   # if len(not_precomputed_compounds) > 0:
   #     not_precomputed_compounds = pd.concat(not_precomputed_compounds)
   # elif len(not_precomputed_compounds) == 1:
   #     not_precomputed_compounds = not_precomputed_compounds[0]
    return precomputed_reactions, not_precomputed_compounds

def find_new_reactions(compounds_data, c2r_output_file, min_diameter = 12, cpu_count = 1, fingerprint_radius = 3, similarity_cutoff = 0.6):
    """
    This function will calculate to which reactions a compound matches based on the retro rules database

    Input is a dataframe for which new reactions should be found and a MAGI parameters dict
    Output is a data frame with the original compound, reactions and similarity scores
    """
    # Load retro rules database
    print("!!! Loading retro rules database")
    sys.stdout.flush()
    read_retro_rules_from_db(min_diameter = min_diameter)
    read_retro_rules_substrates_from_db()

    # Set up multiprocessing
    if cpu_count > 1:
        try:
            p = mp.Pool(cpu_count)
            c2r = p.map(
                partial(compound_to_reaction, 
                        rules_to_use = Retro_rules_reactions, 
                        substrates_to_use = Retro_rules_substrates,
                        min_diameter = min_diameter, 
                        c2r_output_file = c2r_output_file, 
                        fingerprint_radius = fingerprint_radius, 
                        similarity_cutoff = similarity_cutoff), 
                        compounds_data[["SMILES", "original_compound", "compound_score"]].iterrows())
        finally:
            # This should also close the multiprocessing pool in case of an exception
            p.close()
            p.terminate()
    else:
        c2r = []
        for row in compounds_data[["SMILES", "original_compound", "compound_score"]].iterrows():
            c2r_subset = compound_to_reaction(row,  
                            rules_to_use = Retro_rules_reactions, 
                            substrates_to_use = Retro_rules_substrates,
                            min_diameter = min_diameter, 
                            c2r_output_file = c2r_output_file, 
                            fingerprint_radius = fingerprint_radius, 
                            similarity_cutoff = similarity_cutoff)
            c2r.append(c2r_subset)
    
    if all(df is None for df in c2r):
    #    sys.exit("No reactions found for any of the compounds.")
        print("No new reactions found for any of the compounds.")   
        c2r=[] #TODO what to do when no reactions are found for any compound?
    else:
        c2r = pd.concat(c2r) 
    return c2r

def compound_to_reaction(canonical_and_original_smiles_and_compound_score, rules_to_use=Retro_rules_reactions, substrates_to_use = Retro_rules_substrates, min_diameter=0, c2r_output_file=None, fingerprint_radius=3, similarity_cutoff=0.6, use_precomputed = True):
    """
    Find reactions in which a molecule can be used as a substrate
    - For the lowest retro rules diameters, find all reactions in which the compound can be used.
    - For those reactions, check if the reaction also uses the molecule with a higher diameter.
    - Write compound to reaction to intermediate file
    - Return compound to reactions
    """
    molecule_smiles, original_compound, compound_score = canonical_and_original_smiles_and_compound_score[1]
    sys.stdout.flush()
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
        substrate_ID = rules_df.iloc[0]["substrate_ID"]
        if substrate_ID in matching_substrates.keys():
            similarity = matching_substrates[substrate_ID]
        else:
            similarity = 0
        
        if similarity >= similarity_cutoff:
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
                            higher_reaction_matched = mol_matches_reaction(molecule, reaction)
                            if higher_reaction_matched:
                                diameter_to_store = diameter
            # store reaction with highest matching diameter
            if reaction_matched:          
                compound_to_reaction.append([original_compound, rules_df[rules_df["diameter"] == diameter_to_store]["reaction_ID"], similarity, compound_score])
    # Store results and return data as dataframe
    if len(compound_to_reaction) > 0:
        # TODO: filter for highest diameter per reaction
        compound_to_reaction = pd.DataFrame(compound_to_reaction)
        compound_to_reaction.columns = ["original_compound", "reaction_ID", "similarity", "compound_score"]
        try:
            compound_to_reaction["reaction_ID"] = compound_to_reaction["reaction_ID"].astype(int)
        except ValueError:
            print("ERROR: Something wrong with reaction IDs in MAGI database. Exiting...")
        # Store data
        compound_to_reaction.to_csv(c2r_output_file, mode = 'a', header=False, index=False)
        return compound_to_reaction
    else:
        print("No reactions found for {}".format(original_compound))
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
        Retro_rules_substrates[['substrate_ID', 'retro_rules_ID', 'retro_rules_smiles']],
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
    magi_parameters["cpu_count"] = 1 # Quick fix for multiprocessing bug
    
    # Load compounds data and add SMILES, MOLs and InChI Keys
    print("!!! Loading compounds")
    sys.stdout.flush()
    compounds_data = preprocess_compounds_data(magi_parameters["compounds"], magi_parameters["cpu_count"])

    # prepare for output
    c2r_output_file = os.path.join(magi_parameters["intermediate_files_dir"], "compound_to_reaction.csv")
    mg.write_intermediate_file_path(magi_parameters["output_dir"], "compound_to_reaction_path", c2r_output_file)
    print("!!! Saving compound to reaction results to {}".format(c2r_output_file))
    sys.stdout.flush()
    # Write header to output file
    pd.DataFrame(data = None, columns = ["original_compound", "reaction_ID", "similarity", "compound_score"]).to_csv(c2r_output_file, index = False)  
    
    # Check if all compounds have precomputed reactions
    if magi_parameters["use_precomputed_reactions"]:
        print("!!! Finding precomputed reactions for input compounds")
        reactions, not_precomputed_compounds = find_precomputed_reactions(
            compounds_data = compounds_data, 
            c2r_output_file = c2r_output_file,
            min_diameter = magi_parameters["diameter"], 
            similarity_cutoff = magi_parameters["similarity_cutoff"])
        print(str(len(reactions))+" precomputed reactions found")
    else:
        print("!!! Not using precomputed reactions.")
        not_precomputed_compounds = compounds_data
        reactions = []
        # set reactions as empty df? --> yes
    sys.stdout.flush()

    # Run compound to reaction search for all non-precomputed compounds
    if len(not_precomputed_compounds) > 0:
   
        # Run compound to reaction search
        print("!!! Starting new compound to reaction searches for "+str(len(not_precomputed_compounds))+" compounds"); 
        sys.stdout.flush()
        new_reactions = find_new_reactions(not_precomputed_compounds, c2r_output_file)
        if len(new_reactions) > 0:
            print(str(len(new_reactions))+" new reactions found");
            if len(reactions) > 0:
                reactions = pd.concat([reactions, new_reactions])
            else:
                reactions = new_reactions
                del new_reactions
    if len(reactions) == 0:
        sys.exit("No reactions found for any of the compounds.")
        
    # Maybe do some cool scoring here?

    # Format output
    if magi_parameters["compound_to_reaction_only"]:
        reactions = format_output(reactions)
        # Store data
        compound_to_reaction_path = os.path.join(magi_parameters["intermediate_files_dir"], "compound_to_reaction_total.csv")
        print("!!! Saving final compound to reaction results to {}".format(compound_to_reaction_path))
        #TODO: also store compound scores?
        reactions.to_csv(compound_to_reaction_path, index=False)
    print("!!! Finished in {}".format(str(datetime.datetime.now() - start_time)))
    print("!@# Done with compound to reaction search at {}".format(str(datetime.datetime.now())))

if __name__ == "__main__":
    main()
