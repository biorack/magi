import pandas as pd
import re as re
import requests
from bs4 import BeautifulSoup as soup
import os
import sys
import imp
import pickle
import hashlib

from io import open

from chempy import balance_stoichiometry

metatlas_rxn_path = '/Users/bpb/repos/build_magi_ref_data/'
sys.path.append(metatlas_rxn_path)

# import refdata_tools.biopax as bp
import refdata_tools.data_loading_3 as metdata
import imp
# from local_settings import local_settings as settings_loc
# my_settings = getattr(__import__('local_settings', fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)
# sys.path.append(my_settings.path_to_magi)
# import metatlas_reactions.reaction_objects as mr
import ast
import json


import pickle
rhea_pickle_file = '../../data/rhea/biopax/rhea-biopax_full.pkl'
if not os.path.isfile(rhea_pickle_file):
    reactions, metabolites, references = bp.read_biopax_level2('../../data/rhea/biopax/rhea-biopax_full.owl')
    with open(rhea_pickle_file,'wb') as fid:
        pickle.dump({'reactions':reactions,
                    'metabolites':metabolites,
                    'references':references},
                    fid,
                   protocol=pickle.HIGHEST_PROTOCOL)
else:
    with open(rhea_pickle_file,'rb') as fid:
        rhea_pickle = pickle.load(fid)
    reactions = rhea_pickle['reactions']
    metabolites = rhea_pickle['metabolites']
    references = rhea_pickle['references']


cpd_df = pd.read_csv('../../data/standardized_molecules_2.csv')
# cpd_df = pd.read_csv('../../data/chebi_and_metacyc_molecules.csv')
cpd_df.set_index('original_id',drop=True,inplace=True)

# %%time
metdata = imp.reload(metdata)
all_reactions = []
for r in reactions:
    """
       
    """
    rxn = metdata.Reaction()
    rxn.source = 'RHEA:' + metdata.get_rhea_id_from_rhea_reaction(r)
#     rxn.rhea_url = 'http://www.rhea-db.org/reaction?id='+metdata.get_rhea_id_from_rhea_reaction(r)
    direction = metdata.get_direction_from_rhea_reaction(r)
    left_key = 'left'
    right_key = 'right'
    direction = direction.lower()
#      Rhea reactions are labeled as 'bidirectional', 'left to right', 'right to left', and 'undefined'
    if 'left to right' in direction:
        direction = 'left-to-right'
    elif 'right to left' in direction:
        direction = 'left-to-right'
        left_key = 'right'
        right_key = 'left'
    elif 'bidirectional' in direction:
        direction = 'reversible'
    else:
        direction = 'undefined'
    rxn.direction = direction
    
    #Do sequences
    prots = metdata.get_uniprot_id_from_rhea_reaction(r)
    if len(prots) > 0:
        seqs = []
        for protein_id in prots:
            # it is necessary that all prots be made and stored prior to this step.
            # if your protein isn't found then extract sequences and store them using
            # a different workflow.
            path='../../data/rhea/uniprot_sequences/'
            filename = os.path.join(path,protein_id+'.pkl')
            with open(filename,'rb') as fid:
                fasta = pickle.load(fid)
            seq = metdata.ProteinSequence()
            seq.md5=hashlib.md5(fasta['sequence'].encode('utf-8')).hexdigest()
            seq.sequence=fasta['sequence'].encode('utf-8')
            seq.source=protein_id
            seqs.append(seq)
        rxn.sequences = seqs            
    
    rxn.description = metdata.get_rhea_reaction_text_from_rhea_reaction(r)
    rxn.name = metdata.get_rhea_reaction_text_from_rhea_reaction(r)
    
    for rxn_cpd in metdata.get_compounds_from_rhea_reaction(r,metabolites,side=right_key):
        cpd = metdata.Compound()
        cpd.common_name = rxn_cpd['name']
        coefficient = metdata.StoichiometricCoefficient()
        coefficient.coefficient=rxn_cpd['stoichiometry']

        if rxn_cpd['chebi_id'] is not None:
            cpd.chebi_id = rxn_cpd['chebi_id']
            if rxn_cpd['chebi_id'] in cpd_df.index:
                cpd.isomeric_smiles = cpd_df.loc[rxn_cpd['chebi_id'],'original_smiles']
                cpd.formula = cpd_df.loc[rxn_cpd['chebi_id'],'formula']
                cpd.standardized_inchikey = cpd_df.loc[rxn_cpd['chebi_id'],'standardized_inchikey']
                cpd.standardized_smiles = cpd_df.loc[rxn_cpd['chebi_id'],'standardized_smiles']

        coefficient.compound=cpd
        rxn.right_compounds.append(coefficient) #initially set stoichiometry to 1
        
    for rxn_cpd in metdata.get_compounds_from_rhea_reaction(r,metabolites,side=left_key):
        cpd = metdata.Compound()
        cpd.common_name = rxn_cpd['name']
        coefficient = metdata.StoichiometricCoefficient()
        coefficient.coefficient=rxn_cpd['stoichiometry']

        if rxn_cpd['chebi_id'] is not None:
            cpd.chebi_id = rxn_cpd['chebi_id']
            if rxn_cpd['chebi_id'] in cpd_df.index:
                cpd.isomeric_smiles = cpd_df.loc[rxn_cpd['chebi_id'],'original_smiles']
                cpd.formula = cpd_df.loc[rxn_cpd['chebi_id'],'formula']
                cpd.standardized_inchikey = cpd_df.loc[rxn_cpd['chebi_id'],'standardized_inchikey']
                cpd.standardized_smiles = cpd_df.loc[rxn_cpd['chebi_id'],'standardized_smiles']

        coefficient.compound=cpd
        rxn.left_compounds.append(coefficient) #initially set stoichiometry to 1
    
    rxn.balanced = True #Assume Rhea has it right
    
    all_reactions.append(rxn)



# Save each reaction as a JSON file

import ast
import json
for rxn in all_reactions:
    filename=os.path.join('../../data/reaction_json4/',rxn.source.replace(':','_')+'.json')
    d = ast.literal_eval(str(rxn).replace("nan","None"))
    json_str = json.dumps(d, indent=4, sort_keys=True)#,default=default)
    with open(filename,'w') as fid:
        fid.write(unicode(json_str))
#     my_str = item.to_json(filename=os.path.join('../../data/reaction_json/',item.unique_id+'.json'))