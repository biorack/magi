"Run compound to reaction search in separate script"

import sys
import os
import multiprocessing as mp
import pandas as pd
import time
import pickle
import re
from rdkit import Chem
from molvs.standardize import enumerate_tautomers_smiles
import networkx as nx
import workflow_helpers_new as mg
    
def load_compound_results(compounds_file, pactolus, output_dir): #TODO: this to workflow helpers
    """ load compound results"""
    print( '\n!!! LOADING COMPOUNDS')
    compounds = mg.load_dataframe(compounds_file)
    # auto-rename pactolus columns
    if pactolus:
        compounds = mg.reformat_pactolus(compounds)
    # remove any missing compounds
    compounds = compounds[~pd.isnull(compounds['original_compound'])]
    compounds.fillna('', inplace=True)

    if 'original_compound' not in compounds.columns:
        raise RuntimeError('Could not find "original_compound" as a column, please\
            rename the column corresponding to inchi keys for the compounds')
    
    # remove compounds not in the database or network
    print( '!!! Scrubbing compounds')
    compounds['adj'] = compounds['original_compound'].apply(
        lambda x: '-'.join(x.split('-')[:2]))
    
    my_settings = mg.get_settings()
    reference_compounds = mg.load_dataframe(my_settings.compounds_df)
    reference_compounds['adj'] = reference_compounds['inchi_key'].apply(
            lambda x: '-'.join(x.split('-')[:2]))

    filt = compounds.merge(reference_compounds, on='adj', how='left', suffixes=('', '_db'))
    # categorize their reason for not being searched
    not_in_db = filt[pd.isnull(filt['cpd_group'])]
    not_in_db['not_searched_reason'] = 'Not in metabolite database'
    not_in_net = filt[filt['cpd_group'] < 0]
    not_in_net['not_searched_reason'] = 'Not in similarity network yet'
    # combine into one table
    not_searched = pd.concat([not_in_db, not_in_net])
    # make the columns same as user input
    cols = compounds.columns[~compounds.columns.str.contains('adj')].tolist()
    cols.append('not_searched_reason')
    not_searched = not_searched[cols]
    # inform the user and save file
    if not_searched.shape[0] > 0:
        print( 'WARNING: some input compounds were not found in the metabolite database or chemical network; please report these compounds! (see log_unsearched_compounds.csv)')
        print( '!@# {} Compounds not being searched; see log_unsearched_compounds.csv'.format(not_searched['original_compound'].unique().shape[0]))
        not_searched.to_csv(os.path.join(output_dir,
            'log_unsearched_compounds.csv'), index=False)

    to_search = filt[filt['cpd_group'] > 0]['original_compound'].unique()
    compounds = compounds[compounds['original_compound'].isin(to_search)]

    u_cpds = compounds['original_compound'].unique()
    print( '!@# {} total input compounds to search\n'.format(len(u_cpds)))

    if 'compound_score' not in compounds.columns:
        print( 'WARNING: "compound_score" not found as a column; assuming that\
            there is no score for compounds, and setting the compound scores \
            to 1.0')
        compounds['compound_score'] = 1.0
    else:
        compounds['compound_score'] = compounds['compound_score'].apply(float)
    return compounds

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
        for k, v in neighbor_node_dict.iteritems():
            if v != 0:
                transformed_neighbor_node_dict[v].append(k)
        for level, node_list in transformed_neighbor_node_dict.iteritems():
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

def mol_from_inchikey(inchikey, compounds):
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

    matched_inchis = compounds[compounds['inchi_key'].str.contains(
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
    except TypeError, e:
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

def connect_compound_to_reaction(inchikey, reference_compounds, c2r, net, cpd_group_lookup, tautomer=False, neighbor_level=2):
    """
    Connects a compound, its neighbors, and optionally its and its
    neighbors' tautomers to a reaction

    replaces part of old searcher_thorough() function and incorporates
    neighbor_searcher

    Inputs
    ------
    inchikey:       compound's inchikey
    mrs_reaction:   reaction database dataframe
    tautomer:       boolean; whether or not to find tautomers of input
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

    # if tautomer flag, do it the legacy way (don't use precomputed c2r)
    # useful when precomputing a new chemical database and/or chemical network
    if tautomer:
        # find any direct matches
        direct_reaction_idx_list = find_reactions_of_compound(search_inchikey)

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
            tautomer_reaction_idx_list = find_reactions_of_compound(
                tautomer_search_pattern)
            compound_results = enumerate_compound_results(
                inchikey, compound_results, tautomer_reaction_idx_list, level=0,
                neighbor='', note='flat tautomer')
    else:
        # get level zero results from precomputed results
        compound_results = c2r[search_inchikey].copy()
        compound_results['level'] = [0]*len(compound_results['original_compound'])
        compound_results['neighbor'] = ['']*len(compound_results['original_compound'])
        compound_results_list = [pd.DataFrame(compound_results)]

    # find neighbors
    if neighbor_level != 0:
        neighbor_groups = neighbor_finder(inchikey, cpd_group_lookup, chemical_network=net, cpd_group=None, level=neighbor_level)

        # next look for reaction matches to neighbors
        for level, neighbor_compound_list in neighbor_groups:
            for neighbor_inchikey in neighbor_compound_list:
                if tautomer:
                    # first get the direct matches
                    reaction_idx_list = find_reactions_of_compound(
                        neighbor_inchikey)
                    compound_results = enumerate_compound_results(
                        inchikey, compound_results, reaction_idx_list,
                        level=level, neighbor=neighbor_inchikey, note='direct')
                    if tautomer:
                        # then the flat tautomer searches
                        neighbor_mol = mol_from_inchikey(neighbor_inchikey)
                        if neighbor_mol is not None:
                            tautomer_list = tautomer_finder(neighbor_mol)
                            if len(tautomer_list) > 0:
                                tautomer_search_pattern = '|'.join(tautomer_list)
                                tautomer_reaction_idx_list = \
                                    find_reactions_of_compound(tautomer_search_pattern)
                                compound_results = enumerate_compound_results(
                                    inchikey, compound_results,
                                    tautomer_reaction_idx_list, level=level,
                                    neighbor=neighbor_inchikey, note='flat tautomer')
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
    if tautomer:
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
    compound_reaction_result_df.fillna('', inplace=True)
    return compound_reaction_result_df

def load_objects():
    my_settings = mg.get_settings()
    with open(my_settings.c2r, 'r') as fid:
        c2r = pickle.load(fid)
    # with additional metacyc reactions manually added
    # (those that had compounds with R-groups)
    mrs_reaction_path = my_settings.mrs_reaction_path
    
    print( '!!! MRS-Reaction: {}'.format(mrs_reaction_path))
    mrs_reaction = mg.load_dataframe(mrs_reaction_path)
    print( '!!! loading compound table')
    reference_compounds = mg.load_dataframe(my_settings.compounds_df)
    # load the MST chemical network
    with open(my_settings.mst_path, 'r') as f:
        chemical_network = nx.read_graphml(f, node_type=int)
    #chemnet files
    print( '!!! loading chemnet files')
    with open(my_settings.chemnet_pickle, 'r') as fid:
        cpd_group_data = pickle.load(fid)
    cpd_group_lookup = cpd_group_data
    return c2r, mrs_reaction, reference_compounds, cpd_group_lookup, net

def workflow(legacy, level, compound_to_reaction, compounds_file, main_start, cpu_count, fasta_file, output_dir, intermediate_files_dir, pactolus):
    compounds = load_compound_results(compounds_file, pactolus, output_dir)
    """Perform compound to reaction search"""
    c2r, mrs_reaction, reference_compounds, cpd_group_lookup, net = load_objects()
    
    def connect_compound_to_reaction_mp_helper(inchikey, 
                                            tautomer=legacy, 
                                            neighbor_level=level, 
                                            reference_compounds = reference_compounds,
                                            cpd_group_lookup =cpd_group_lookup,
                                            c2r=c2r,
                                            net=net):
        try:
            out = connect_compound_to_reaction(inchikey,
                                               reference_compounds=reference_compounds, 
                                               c2r=c2r, 
                                               net=net, 
                                               cpd_group_lookup=cpd_group_lookup,
                                            tautomer=tautomer, 
                                            neighbor_level=neighbor_level)
        except Exception as e:
            print( inchikey )
            sys.stdout.flush()
            raise RuntimeError('offending inchikey: %s; error message: %s' \
                                %(inchikey, e.args))
        return out
    
    if compound_to_reaction is None:
        print( '\n!@# Conducting compound to reaction search | TLOG %s' % (time.time()))
        sys.stdout.flush()
        start = time.time()

        input_compounds = compounds['original_compound'].unique()
        if cpu_count == 1: #Quick fix for multiprocessing problems
            #To do: fix this
            print("!!! No multiprocessing used for compound to reaction search")
            out = map(connect_compound_to_reaction_mp_helper, input_compounds)
        else:
            p = mp.Pool(cpu_count)
            out = p.map(connect_compound_to_reaction_mp_helper, input_compounds)
            p.close()
            p.terminate()
        
        compound_to_reaction = pd.concat(out)
        del out
        compound_to_reaction.reset_index(inplace=True, drop=True)
    
        # connect the compound score
        compound_to_reaction = pd.merge(compounds, compound_to_reaction, 
                                        on='original_compound', how='inner')
        print( '!@# compound_to_reaction table done in %s minutes'\
                    %((time.time()-start)/60))
        # if no fasta file then just save these results and quit
        if fasta_file is None:
            compound_to_reaction.to_csv(os.path.join(output_dir, 
                                                    'magi_compound_results.csv'))
            print( '!!! compound_reaction table saved to %s'\
                    % (os.path.join(output_dir, 'magi_compound_results.csv')))
            print( '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60))
            sys.exit()
        else:
            compound_to_reaction.to_pickle(os.path.join(intermediate_files_dir, 
                                                    'compound_to_reaction.pkl'))
    
            print( '!!! compound_reaction table saved to %s'\
                    % (os.path.join(intermediate_files_dir, 'compound_to_reaction.pkl')))
    else:
        compound_to_reaction = pd.read_pickle(compound_to_reaction)
        print( '\n!@# compound_to_reaction successfully loaded')
    return compounds, compound_to_reaction


#def parse_arguments():
#    return "arguments are parsed"
#
#def format_output():
#    return "output is formatted"
#
##def main():
##    parse_stuff()
##    workflow()
##    format_output()
#
#if __name__ == "__main__":
#    #main()
#    print("Are you sure? This is not ready yet")