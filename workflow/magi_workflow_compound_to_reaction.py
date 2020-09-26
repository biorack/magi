"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.1 compound to reaction workflow script

Lines beginning with "@@@" are input job parameters
Lines beginning with "!!!" are verbose log info
Lines beginning with "!@#" are checkpoints/announcements

Required inputs is a Compounds file.
The Compounds file should be in some standard table format (CSV, tab-
delimited, etc.), and is required to have a column named 
"original_compound" (case-sensitive). In this column should be standard
InChI Keys representing a compound structure.

If you have only m/z values, run the accurate mass search workflow before 
running this script.

The Compounds table may also have a column named "compound_score" 
(case-sensitive), where the user can provide a score for each 
compound-row. If this column name does not exist, one will be created
and populated with 1.0
"""
import sys
import os
import argparse
import multiprocessing as mp
import pandas as pd
import time
import pickle
import re
from rdkit import Chem
from molvs.standardize import enumerate_tautomers_smiles
import networkx as nx
import workflow_helpers_new as mg
from functools import partial

def find_reactions_of_compound(inchikey, rxn_db, 
                               compound_col='allcpd_ikeys'):
    """
    Given an inchikey input, find reactions that contain the compound

    Inputs
    ------
    inchikey: can be one, two, or three block inchi
    rxn_db: reaction database table
    compound_col: column in rxn_db that contains a list of inchikeys of
                  all compounds involved in that reaction. Elements in
                  this column must be a string for now (e.g. not a list)

    Outputs
    -------
    A list of reaction indices, or None if the inchikey was not found in
    any reactions.
    """
    # TODO: use join tables instead of .str.contains() method

    reactions = rxn_db[rxn_db[compound_col].str.contains(
        inchikey)]
    if reactions.shape[0] != 0:
        reaction_idx_list = reactions.index.tolist()
        return reaction_idx_list
    else:
        return None

def neighbor_finder(inchikey, cpd_group_lookup, chemical_network, cpd_group=None, level=2):
    """
    finds neighbors of a given compound to the given level in the
    chemical network

    Inputs
    ------
    inchikey:           Full or partial inchikey for the compound to be
                        searched
    chemical_network:   MST chemical network
    cpd_group:          if the compound group of inchikey is already
                        known, this will override the cpd_group
                        neighbor_finder
    level:              to what level to search the network

    Outputs
    -------
    returns a list of tuples describing the network level and inchikeys:
    (level, [list of inchikeys at that level])
    """

    if cpd_group is None:
        # find the compound group/node number
        # flatten the inchikey
        inchikey_query = inchikey.split('-')[0]

        # find the compound group of the inchikey
        # there should only be one!!! If there are multiple, this will
        # only take the first group found
        for idx, istring in enumerate(cpd_group_lookup):
            if inchikey_query in istring:
                cpd_group = idx
                break
    neighbor_groups = []
    if cpd_group is not None:
        # compound group number is also the node ID; find the compound
        # groups that are neighbors in the MST
        neighbor_node_dict = nx.single_source_shortest_path_length(
            chemical_network, cpd_group, level)
        transformed_neighbor_node_dict = {}
        for k in list(set(neighbor_node_dict.values())):
            transformed_neighbor_node_dict[k] = []
        for k, v in neighbor_node_dict.items():
            if v != 0:
                transformed_neighbor_node_dict[v].append(k)
        for level, node_list in transformed_neighbor_node_dict.items():
            neighbor_inchikey_list = []
            for neighbor_cpd_group in node_list:
                neighbor_inchikey_str = cpd_group_lookup[neighbor_cpd_group]
                for neighbor_inchikey in neighbor_inchikey_str.split('///'):
                    neighbor_inchikey_list.append(neighbor_inchikey)
            neighbor_groups.append((level, neighbor_inchikey_list))

    else:
        print( '695 WARNING: Could not find "%s" in the chemical network' \
              % (inchikey))

    return neighbor_groups

def enumerate_compound_results(original_compound, compound_results,
                               reaction_idx_list, level=0,
                               neighbor='', note=''):
    """
    Inputs
    ------
    original_compound: initial compound that began the search
    compound_results:   dictionary with the following keys:
                        'original_compound', 'level', 'neighbor',
                        'reaction_id', 'note' where each key-value is a
                        list
    reaction_idx_list:  the ouptut of find_reaction_of_compound()
    level:              The chemical network search level
    neighbor:           The chemical network neighbor of
                        original_compound

    Outputs:
    --------
    An updated compound_results dictionary
    """

    if reaction_idx_list is not None:
        for rxn_idx in reaction_idx_list:
            compound_results['original_compound'].append(original_compound)
            compound_results['level'].append(level)
            compound_results['neighbor'].append(neighbor)
            compound_results['reaction_id'].append(rxn_idx)
            compound_results['note'].append(note)
    else:
        compound_results['original_compound'].append(original_compound)
        compound_results['level'].append(level)
        compound_results['neighbor'].append(neighbor)
        compound_results['reaction_id'].append('')
        compound_results['note'].append(note)

    return compound_results    

def mol_from_inchikey(inchikey, reference_compounds):
    """
    Returns an RdKit Mol object from an inchikey input via the compound
    DataFrame

    Inputs
    ------
    inchikey: a standard InChI key

    Outputs
    -------
    Rdkit Mol object or None, if the inchikey input was not found in the
    compound DataFrame
    """

    matched_inchis = reference_compounds[reference_compounds['inchi_key'].str.contains(
        inchikey, regex=False)]['inchi'].values
    if matched_inchis.any():
        # in case there are duplicates, just take the first one
        inchi = str(matched_inchis[0])
        compound_mol = Chem.MolFromInchi(inchi)
        return compound_mol
    else:
        return None
    
def tautomer_finder(compound_mol, result='split', raise_errors=False):
    """
    enumerates tautomers of a given compound

    MolVS tautomer enumerator flattens compounds, so the default only
    returns the first block of the inchikey

    Inputs
    ------
    compound_mol: RdKit Mol object of the compound to tautomerize
    result: "split" returns the first block of an inchikey
            "full" returns the full inchikey
            "smiles" returns the smiles string
            "inchi" returns the InChI string 
            "mol" returns RdKit Mol object
    raise_errors: when False, it still prints the error as a warning

    Outputs
    -------
    tautomer_list: list of unique tautomers, in the output format
                   specified by result argument
    """

    if not isinstance(compound_mol, type(Chem.Mol())):
        raise RuntimeError('The input is not an rdkit Mol object!')

    if not isinstance(result, str):
        raise RuntimeError('"result" arg must be a string!')
    if result.lower() not in ['split', 'full', 'smiles', 'inchi', 'mol']:
        raise RuntimeError('"result" arg must be either "split", "full", \
            "smiles", "inchi", or "mol"')

    compound_smiles = Chem.MolToSmiles(compound_mol, isomericSmiles=True)

    # some compounds break the tautomerizer
    try:
        enumerated_tautomers = enumerate_tautomers_smiles(compound_smiles)
    except TypeError as e:
        if e.message == 'tuple indices must be integers, not NoneType':
            enumerated_tautomers = []
    except Exception as e:
        if raise_errors is False:
            print( 'WARNING: %s could not be tautomerized; %s' \
                    % (Chem.InchiToInchiKey(Chem.MolToInchi(compound_mol)),
                        e.args))
            enumerated_tautomers = [compound_smiles]
        else:
            raise

    # convert smiles into inchikeys
    tautomer_list = []
    for ts in enumerated_tautomers:
        if result.lower() in ['split', 'full']:
            i = Chem.InchiToInchiKey(Chem.MolToInchi(Chem.MolFromSmiles(ts)))
            if result.lower() == 'split':
                tautomer_list.append(i.split('-')[0])
            elif result.lower() == 'full':
                tautomer_list.append(i)
            else:
                raise RuntimeError(
                    'could not determine what to do with %s' % (result))
        elif result.lower() == 'smiles':
            tautomer_list.append(ts)
        elif result.lower() == 'inchi':
            tautomer_list.append(Chem.MolToInchi(Chem.MolFromSmiles(ts)))
        elif result.lower() == 'mol':
            tautomer_list.append(Chem.MolFromSmiles(ts))
        else:
            raise RuntimeError(
                'could not determine what to do with %s' % (result))

    return list(set(tautomer_list))

def find_direct_reactions(search_inchikey, inchikey, reference_compounds, c2r, mrs_reaction, tautomer_legacy):
    """
    This function finds reactions in which an inchi key is involved.

    Inputs:
    ------
    search_inchikey: The first two parts of the inchi key of interest.
    inchikey:        The full inchi key of interest.
    tautomer_legacy: Bool that says if tautomer legacy should be used or not.
    """
    # if tautomer flag, do it the legacy way (don't use precomputed c2r)
    # useful when precomputing a new chemical database and/or chemical network
    if tautomer_legacy:
        # find any direct matches
        direct_reaction_idx_list = find_reactions_of_compound(search_inchikey, rxn_db=mrs_reaction)

        # initialize the results dict
        compound_results = {
            'original_compound': [],
            'level': [],
            'neighbor': [],
            'reaction_id': [],
            'note': []
            }
        compound_results = enumerate_compound_results(
            inchikey, compound_results, direct_reaction_idx_list,
            level=0, neighbor='', note='direct')
        # make an rdkit mol of the compound
        compound_mol = mol_from_inchikey(inchikey, reference_compounds)
        if compound_mol is None:
            print( 'WARNING: Could not find "%s" in the compound database; \
                skipping tautomer and neighbor searching' % (inchikey))
            return pd.DataFrame(compound_results)

        # get tautomers of the compound
        tautomer_list = tautomer_finder(compound_mol)

        # find reactions for the tautomers
        if len(tautomer_list) > 0:
            tautomer_search_pattern = '|'.join(tautomer_list)
            tautomer_reaction_idx_list = find_reactions_of_compound(tautomer_search_pattern, rxn_db=mrs_reaction)
            compound_results = enumerate_compound_results(
                inchikey, compound_results, tautomer_reaction_idx_list, level=0,
                neighbor='', note='flat tautomer')
            compound_results_list = None
    else:
        # get level zero results from precomputed results
        compound_results = c2r[search_inchikey].copy()
        compound_results['level'] = [0]*len(compound_results['original_compound'])
        compound_results['neighbor'] = ['']*len(compound_results['original_compound'])
        compound_results_list = [pd.DataFrame(compound_results)]
        compound_results = None
    return compound_results_list, compound_results

def find_neighbor_reactions(compound_results_list, compound_results, inchikey, reference_compounds, c2r, mrs_reaction, chemical_network, cpd_group_lookup, tautomer_legacy, neighbor_level):
    """
    This function finds reactions in which neighbors of inchi keys are involved.
    """    
    neighbor_groups = neighbor_finder(inchikey, cpd_group_lookup, chemical_network=chemical_network, cpd_group=None, level=neighbor_level)

    # next look for reaction matches to neighbors
    for level, neighbor_compound_list in neighbor_groups:
        for neighbor_inchikey in neighbor_compound_list:
            if tautomer_legacy:
                # first get the direct matches
                reaction_idx_list = find_reactions_of_compound(
                    neighbor_inchikey, rxn_db=mrs_reaction)
                compound_results = enumerate_compound_results(
                    inchikey, compound_results, reaction_idx_list,
                    level=level, neighbor=neighbor_inchikey, note='direct')
                # then the flat tautomer searches
                neighbor_mol = mol_from_inchikey(neighbor_inchikey, reference_compounds)
                if neighbor_mol is not None:
                    tautomer_list = tautomer_finder(neighbor_mol)
                    if len(tautomer_list) > 0:
                        tautomer_search_pattern = '|'.join(tautomer_list)
                        tautomer_reaction_idx_list = \
                            find_reactions_of_compound(tautomer_search_pattern, rxn_db=mrs_reaction)
                        compound_results = enumerate_compound_results(
                            inchikey, compound_results,
                            tautomer_reaction_idx_list, level=level,
                            neighbor=neighbor_inchikey, note='flat tautomer')
                compound_results_list = None
            else:
                search_inchikey = '-'.join(neighbor_inchikey.split('-')[:2])
                try:
                    tmp_compound_results = c2r[search_inchikey].copy()
                except KeyError:
                    continue
                tmp_compound_results['level'] = [level]*len(tmp_compound_results['original_compound'])
                tmp_compound_results['neighbor'] = tmp_compound_results['original_compound']
                tmp_compound_results['original_compound'] = [inchikey]*len(tmp_compound_results['original_compound'])
                compound_results_list.append(pd.DataFrame(tmp_compound_results))
                compound_results = None

    return compound_results_list, compound_results

def connect_compound_to_reaction(inchikey, reference_compounds, c2r, mrs_reaction, chemical_network, cpd_group_lookup, tautomer_legacy=False, neighbor_level=2):
    """
    Connects a compound, its neighbors, and optionally its and its
    neighbors' tautomers to a reaction

    replaces part of old searcher_thorough() function and incorporates
    neighbor_searcher

    Inputs
    ------
    inchikey:       compound's inchikey
    mrs_reaction:   reaction database dataframe
    tautomer_legacy:       boolean; whether or not to find tautomers of input
                    compound (and it's neighbors)
    neighbor_level: to what level to search the chemical network
                    0: does not search the network
                    1: finds the immediate neighbors
                    2: finds the neighbor's neighbors
                    etc.

    Outputs
    -------
    a dataframe describing the compound and what reactions it connected
    through, through which neighbors if applicable. The final DataFrame
    is sorted by direct/flat tautomer, and then
    original_compound-reaction duplicates are dropped.

    """
    # regex expression to test inchikey input
    # Character 9 of the second block must be "S" (standard inchi)
    # Character 10 of the second block must be "A" (version 1)
    inchikey_test = '[A-Z]{14}-[A-Z]{8}SA-[A-Z]{1}'
    # test validity of input
    if not re.match(inchikey_test, inchikey):
        raise RuntimeError('"%s" is not a valid inchikey or is not the correct\
         version of inchi : xxxxxxxxxxxxxx-xxxxxxxxSA-X' % (inchikey))

    # convert inchikey into two-block
    search_inchikey = '-'.join(inchikey.split('-')[:2])

    # Find direct reactions for the compound
    compound_results_list, compound_results = find_direct_reactions(search_inchikey, inchikey, reference_compounds, c2r, mrs_reaction, tautomer_legacy)

    # find neighbors
    if neighbor_level != 0:
        compound_results_list, compound_results = find_neighbor_reactions(compound_results_list, compound_results, inchikey, reference_compounds, c2r, mrs_reaction, chemical_network, cpd_group_lookup, tautomer_legacy, neighbor_level)

    if tautomer_legacy:
        compound_reaction_result_df = pd.DataFrame(compound_results)
    else:
        compound_reaction_result_df = pd.concat(compound_results_list)

    # now perform cleanup on the df
    # first sort the df so that all the direct hits come first
    compound_reaction_result_df.sort_values(
        ['note', 'level'], ascending=[True, True], inplace=True)
    # then drop any duplicated compound-reaction associations, keeping
    # only the direct hits
    compound_reaction_result_df.drop_duplicates(
        ['original_compound', 'reaction_id'], inplace=True)

    return compound_reaction_result_df

def load_objects(use_tautomer_legacy = False):
    """
    Load to memory:
    If the tautomer legacy is false:
    c2r is a compounds to reaction dictionary with the inchikey as key and 
        a dictionary with note, original compound and reaction id as values. 

    If the tautomer legacy is true:
    mrs_reaction is a Pandas DataFrame with information on reactions. 
        Only loaded to memory if tautomer_legacy is False.
    reference_compounds is a Pandas DataFrame with MAGI-compatible Inchi keys, 
        inchis, mono-isotopic molecular weight and the compound group in which it belongs.
        Only loaded to memory if tautomer_legacy is False.

    cpd_group_lookup is a list of all compound groups. 
        The row number is the ID of the compound group in the chemical network.
    chemical network is a graph with compound groups as nodes and similarities as edges.
    
    """
    my_settings = mg.get_settings()
    # Either c2r or reference compounds is needed. This saves some memory
    if use_tautomer_legacy:
        print( '!!! loading compound table')
        # Loading reference_compounds
        reference_compounds = mg.load_dataframe(my_settings.compounds_df)
        # with additional metacyc reactions manually added
        # (those that had compounds with R-groups)
        mrs_reaction_path = my_settings.mrs_reaction_path    
        print( '!!! MRS-Reaction: {}'.format(mrs_reaction_path))
        mrs_reaction = mg.load_dataframe(mrs_reaction_path)
        c2r = None
    else:
        print("!!! loading compound-reactions table")
        # Loading c2r
        c2r = pd.read_pickle(my_settings.c2r)
        reference_compounds = None
        mrs_reaction = None

    # load the MST chemical network
    with open(my_settings.mst_path, 'r') as f:
        chemical_network = nx.read_graphml(f, node_type=int)
    #chemnet files
    print( '!!! loading chemnet files')
    cpd_group_lookup = pd.read_pickle(my_settings.chemnet_pickle)

    return c2r, mrs_reaction, reference_compounds, cpd_group_lookup, chemical_network

def workflow(compounds_to_search, tautomer_legacy, neighbor_level, cpu_count, intermediate_files_dir):
    """
    This workflow looks for compounds and finds the reactions

    Inputs: 
    ------
    compounds_to_search: pandas dataframe with compounds of interest in original_compound colum
    tautomer_legacy: whether or not to use pre-compouted tautomers of compounds
    neighbor_level: level to search in the chemical similarity network
    cpu_count: number of CPUs to use
    intermediate_files_dir: path to folder to store intermediate files

    Outputs:
    ------
    compound_to_reaction_path: path to the compound_to_reaction pickle file containing the dataframe.
    """
    c2r, mrs_reaction, reference_compounds, cpd_group_lookup, chemical_network = load_objects(use_tautomer_legacy = tautomer_legacy)

    print( '\n!@# Conducting compound to reaction search | TLOG %s' % (time.time()))
    sys.stdout.flush()
    start = time.time()

    input_compounds = compounds_to_search['original_compound'].unique()
    if cpu_count == 1: #Don't set up multiprocessing if only one CPU is used.
        # Get a DataFrame with reaction info for each input_compound
        print("!!! No multiprocessing used for compound to reaction search")
        out = map(partial(connect_compound_to_reaction, 
                            reference_compounds=reference_compounds, 
                            c2r=c2r, 
                            mrs_reaction = mrs_reaction,
                            chemical_network=chemical_network, 
                            cpd_group_lookup=cpd_group_lookup,
                            tautomer_legacy=tautomer_legacy, 
                            neighbor_level=neighbor_level), input_compounds)
    else:
        try:
            p = mp.Pool(cpu_count)
            out = p.map(partial(connect_compound_to_reaction, 
                                reference_compounds=reference_compounds, 
                                c2r=c2r, 
                                mrs_reaction = mrs_reaction,
                                chemical_network=chemical_network, 
                                cpd_group_lookup=cpd_group_lookup,
                                tautomer_legacy=tautomer_legacy, 
                                neighbor_level=neighbor_level), input_compounds) 
        finally:
            # This should also close the multiprocessing pool in case of an exception
            p.close()
            p.terminate()
    
    compound_to_reaction = pd.concat(out)
    del out
    compound_to_reaction.reset_index(inplace=True, drop=True)

    # connect the compound score
    compound_to_reaction = pd.merge(compounds_to_search, compound_to_reaction, 
                                    on='original_compound', how='inner')
    print( '!@# compound_to_reaction table done in %s minutes'\
                %((time.time()-start)/60))

    compound_to_reaction_path = os.path.join(intermediate_files_dir, 
                                                'compound_to_reaction.pkl')
    compound_to_reaction.to_pickle(compound_to_reaction_path)

    print( '!!! compound_reaction table saved to %s'\
            % (compound_to_reaction_path))

    return compound_to_reaction_path

def format_output(compound_to_reaction_path, output_dir, intermediate_files_dir):
    """
    This function reads the compound to reaction pickle and 
    formats it to be a standalone output file with useful information.
    It writes the file magi_compound_results.csv to the specified output directory.

    Inputs
    ------
    compound_to_reaction_path: path to the compound_to_reaction.pkl pickle dataframe.
    intfile_path: path to the temporary storage place of MAGI files.
    output_dir: path to the location where the magi_compound_results.csv file needs to be saved.
    """
    compound_to_reaction = pd.read_pickle(compound_to_reaction_path)
    # fill NA with empty strings and clean up reaction ID column 
    compound_to_reaction.reaction_id = compound_to_reaction.reaction_id.fillna(-1)
    compound_to_reaction.reaction_id = compound_to_reaction.reaction_id.replace('',-1)
    compound_to_reaction.reaction_id = compound_to_reaction.reaction_id.astype(int)
    compound_to_reaction.reaction_id = compound_to_reaction.reaction_id.replace(-1,'')
    compound_to_reaction.fillna('', inplace = True)
    
    #TODO: Merge some useful information to the data frame.
    compound_to_reaction.to_csv(os.path.join(output_dir, 
                                            'magi_compound_results.csv'))
    print( '!!! compound_reaction table saved to %s'\
            % (os.path.join(output_dir, 'magi_compound_results.csv')))
    with open(os.path.join(intermediate_files_dir, "timer.txt"),"r") as timerfile:
            main_start = float(timerfile.read())
    print( '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60))

def main():
    # Parse arguments and prepare for compounds to reaction workflow
    magi_parameters = mg.general_magi_preparation()

    if magi_parameters["gene_to_reaction_only"]:
        print("Not performing MAGI compound to reaction workflow")
    else:
        if magi_parameters["is_mass_search"]:
            magi_parameters["compounds"] = os.path.splitext(magi_parameters["compounds"])[0] + '_mass_searched.csv'

        # Open file with compounds that need to be searched and prepare the data frame
        compounds_to_search = mg.load_compound_results(
                compounds_file = magi_parameters["compounds"], 
                pactolus = magi_parameters["pactolus"], 
                output_dir = magi_parameters["output_dir"], 
                intermediate_files_dir = magi_parameters["intermediate_files_dir"])

        # Run compound to reaction workflow 
        compound_to_reaction_path = workflow(
                compounds_to_search = compounds_to_search, 
                tautomer_legacy = magi_parameters["legacy"], 
                neighbor_level = magi_parameters["level"], 
                cpu_count = magi_parameters["cpu_count"], 
                intermediate_files_dir = magi_parameters["intermediate_files_dir"])    
        mg.write_intermediate_file_path(magi_parameters["output_dir"], "compound_to_reaction_path", compound_to_reaction_path)

        # Format output if this is the last step of the MAGI run
        if magi_parameters["compound_to_reaction_only"]:
            format_output(
                compound_to_reaction_path = compound_to_reaction_path, 
                output_dir = magi_parameters["output_dir"], 
                intermediate_files_dir = magi_parameters["intermediate_files_dir"])

if __name__ == "__main__":
    main()