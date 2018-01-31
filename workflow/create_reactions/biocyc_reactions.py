from __future__ import print_function

import sys
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

import pandas as pd
import numpy as np
from sqlite3 import dbapi2 as sqlite
import re
import copy
import requests
import glob
import hashlib

from chempy import balance_stoichiometry

from io import open

import data_loading as metdata

url_pattern = "href=([^>]+)"
import imp

pd.set_option('display.max_rows', 500)
pd.set_option("display.max_colwidth", 1000000)



"""
#Specify locations of flatfile databases
# Using Metacyc 21.5 flat files

Get them here: http://bioinformatics.ai.sri.com/ecocyc/dist/flatfiles-52983746/

Username: biocyc-flatfiles

Password: ####-##### 

Search your email for actual password.  Subject = Pathway Tools and BioCyc distribution available
"""
biocyc_paths = ['../../data/tier1/metacyc_21.5/data/','../../data/tier1/Eco_21.5/data/','../../data/tier1/Sco_17.5/data/']

# location of table with original compounds and standardized forms
standardized_molecule_path = '../../data/standardized_molecules_2.csv'


monomer_path='/users/bpb/repos/build_magi_ref_data/data/tier1/monomer_sequences'



# #Make a compounds dataframe (name and ID only)
# my_files = [os.path.join(path,'compounds.dat') for path in biocyc_paths]
# cpd_name_df = metdata.parse_biocyc_flat_file(my_files,attributes =  ['UNIQUE-ID','COMMON-NAME'])
# cpd_name_df['UNIQUE-ID'] = cpd_name_df['UNIQUE-ID'].apply(lambda x: x[0] if len(x)>0 else '')
# cpd_name_df['COMMON-NAME'] = cpd_name_df['COMMON-NAME'].apply(lambda x: x[0] if len(x)>0 else '')
# cpd_name_df.drop_duplicates(subset='UNIQUE-ID',inplace=True)
# cpd_name_df.set_index('UNIQUE-ID',drop=True,inplace=True)

#Make a reactions DataFrame
my_files = [os.path.join(path,'reactions.dat') for path in biocyc_paths]
rxn_df = metdata.parse_biocyc_flat_file(my_files,attributes = ['UNIQUE-ID','ENZYMATIC-REACTION','LEFT','REACTION-DIRECTION','RIGHT'])
rxn_df['REACTION-DIRECTION'] = rxn_df['REACTION-DIRECTION'].apply(lambda x: x[0] if len(x)>0 else '')
rxn_df['UNIQUE-ID'] = rxn_df['UNIQUE-ID'].apply(lambda x: x[0] if len(x)>0 else '')
rxn_df.drop_duplicates(subset='UNIQUE-ID',inplace=True)

# Check all reactions to ensure all compounds are accounted for

cpd_df = pd.read_csv(standardized_molecule_path)
cpd_df.set_index('original_id',drop=True,inplace=True)

# get the few sequences provided in flatfile download and store in prot_df

fasta_info = []
my_files = [os.path.join(path,'protseq.fsa') for path in biocyc_paths]
fasta_info = []
all_prots_dfs = []
for my_file in my_files[1:2]:
    with open(my_file,'r',encoding='utf-8', errors='ignore') as fid:
        fasta_info = fid.read()
    prots = [tt.strip().split('\n') for tt in fasta_info.split('\n>') if len(tt)>0]
    for i,p in enumerate(prots):
        prots[i][0] = p[0].split('|')[-1].split(' ')[0]
        prots[i][1] = ''.join(p[1:])
        prots[i] = prots[i][:2]
    all_prots_dfs.append(pd.DataFrame(prots,columns=['UNIQUE-ID','sequence']))
prot_df = pd.concat(all_prots_dfs)
    


# get enzyme complex composition and merge sequences provided in flatfile distribution

all_enzymes_dfs = []
for path in biocyc_paths:
    temp = pd.read_csv(os.path.join(path,'enzymes.col'),sep='\t',comment='#',encoding='ISO-8859-1')
    all_enzymes_dfs.append(temp[['UNIQUE-ID','NAME','SUBUNIT-COMPOSITION']])
enzyme_df = pd.concat(all_enzymes_dfs)


all_subunits = []
for index,row in enzyme_df.iterrows():
    if not pd.isnull(row['SUBUNIT-COMPOSITION']):
        subunits = [re.sub('\d+\*','',s.strip()) for s in row['SUBUNIT-COMPOSITION'].split(',')]
        for s in subunits:
            all_subunits.append(s)

all_subunits = pd.unique(all_subunits)

all_subunits_df = pd.DataFrame(all_subunits,columns=['protein_id'])
all_subunits_df = all_subunits_df.merge(prot_df,how='outer',left_on='protein_id',right_on='UNIQUE-ID')
all_subunits_df['sequence'] = all_subunits_df['sequence'].replace(pd.np.nan,'',regex=True)
all_subunits_df.drop('UNIQUE-ID',axis=1,inplace=True)

# populate dataframe with sequences in files

for i,row in all_subunits_df.iterrows():
    monomer_file = os.path.join(monomer_path,'%s.txt'%row.protein_id)
    if os.path.isfile(monomer_file):
        with open(monomer_file,'r') as fid:
            seq_str = fid.read()
            all_subunits_df.loc[i,'sequence'] = "".join(seq_str.split())
all_subunits_df['sequence'] = all_subunits_df['sequence'].replace(pd.np.nan,'',regex=True)


# scrape the majority of other sequences using web APIs

# for i,row in all_subunits_df.iterrows():
#     metdata.scrape_sequences_not_in_flatfile(row['protein_id'],monomer_path=monomer_path,do_print=True)

# print(len(pool_resuts),all_subunits_df.shape)
# # all_subunits_df.to_csv('../../data/tier1_sequences_protein-subunits_4.csv.gz',compression='gzip')


"""
Some get screwed up.  The value in sequence looks like this: <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 T...

"""
idx = all_subunits_df.sequence.str.contains('<!DOCTYPE html PUBLIC ')
print('This many are screwed up',sum(idx))

"""
Save the subunit sequences each time you do something. If something major changes, dump them to monoer sequence files
"""
all_subunits_df.to_csv('../../data/tier1_sequences_protein-subunits_4.csv.gz',compression='gzip',encoding='utf-8')
# all_subunits_df = pd.read_csv('../../data/tier1_sequences_protein-subunits_4.csv.gz')
# monomer_path='/users/bpb/repos/build_magi_ref_data/data/tier1/monomer_sequences'
# for i,row in all_subunits_df.iterrows():
#     if (not pd.isnull(row.sequence)) and (len(row.sequence)>0):
#         monomer_file = os.path.join(monomer_path,'%s.txt'%row.protein_id)
#         with open(monomer_file,'w') as fid:
#             fid.write(row.sequence)

# make reactions with sequences, stoichiometry, and compounds

all_reactions = []
balanced = 0
undefined = 0
illogical = 0
identical = 0
for i,row in rxn_df.iterrows():
    rxn = metdata.Reaction()
    rxn.balanced=False
    # Create Descriptive Information
    rxn.source = 'BIOCYC:%s'%row['UNIQUE-ID'] #store the reaction id here.  makes it easy to associate enzyme with reaction later

    # TODO ADD DIRECTIONS
    left_key = 'LEFT'
    right_key = 'RIGHT'
    direction = row['REACTION-DIRECTION']
    direction = direction.lower()
    if 'left-to-right' in direction:
            direction = 'left-to-right'
    elif 'right-to-left' in direction:
            direction = 'left-to-right'
            left_key = 'RIGHT'
            right_key = 'LEFT'
    else:
        direction = 'reversible'
    rxn.direction = direction

    for rxn_cpd in row[right_key]:
        cpd = metdata.Compound()
        if rxn_cpd in cpd_df.index:
            cpd.isomeric_smiles = cpd_df.loc[rxn_cpd,'original_smiles']
            cpd.common_name = cpd_df.loc[rxn_cpd,'name']
            cpd.formula = cpd_df.loc[rxn_cpd,'formula']
            cpd.standardized_inchikey = cpd_df.loc[rxn_cpd,'standardized_inchikey']
            cpd.standardized_smiles = cpd_df.loc[rxn_cpd,'standardized_smiles']
            cpd.metacyc_id = rxn_cpd
        else:
            cpd.common_name = rxn_cpd
            cpd.metacyc_id = rxn_cpd
        coefficient = metdata.StoichiometricCoefficient()
        coefficient.coefficient=1
        coefficient.compound=cpd
        rxn.right_compounds.append(coefficient) #initially set stoichiometry to 1

    for rxn_cpd in row[left_key]:
        cpd = metdata.Compound()
        if rxn_cpd in cpd_df.index:
            cpd.isomeric_smiles = cpd_df.loc[rxn_cpd,'original_smiles']
            cpd.common_name = cpd_df.loc[rxn_cpd,'name']
            cpd.formula = cpd_df.loc[rxn_cpd,'formula']
            cpd.standardized_inchikey = cpd_df.loc[rxn_cpd,'standardized_inchikey']
            cpd.standardized_smiles = cpd_df.loc[rxn_cpd,'standardized_smiles']
            cpd.metacyc_id = rxn_cpd
        else:
            cpd.common_name = rxn_cpd
            cpd.metacyc_id = rxn_cpd
        coefficient = metdata.StoichiometricCoefficient()
        coefficient.coefficient=1
        coefficient.compound=cpd
        rxn.left_compounds.append(coefficient) #initially set stoichiometry to 1
    
    # Get all enzymes and get all sequences
    seqs = []
    for enzyme_id in rxn_df.loc[i,'ENZYMATIC-REACTION']:
        # For each reaction, get the enzyme ID
        subunit_text = enzyme_df.loc[enzyme_df['UNIQUE-ID'] == enzyme_id,'SUBUNIT-COMPOSITION']
        try:
            subunit_text = subunit_text.values[0]
        except:
            subunit_text = pd.np.nan
        #For each enzyme ID get the subunit text string
        if not pd.isnull(subunit_text):
            subunits = [re.sub('\d+\*','',s.strip()) for s in subunit_text.split(',')]
            for s in subunits:
                seq = all_subunits_df.loc[all_subunits_df.protein_id == s,'sequence'].values[0]
                md5=hashlib.md5(seq.encode('utf-8')).hexdigest()
                metseq = metdata.ProteinSequence()
                metseq.sequence=seq
                metseq.source=s
                metseq.md5=md5
                seqs.append(metseq)
    rxn.sequences = seqs
    
    lhs = []
    rhs = []   
    for cpds in rxn.left_compounds:
        lhs.append(cpds.compound.formula)
    for cpds in rxn.right_compounds:
        rhs.append(cpds.compound.formula)
    
    try:
        # an incompletely defined molecule will have empty string for formula
        if min([len(f) for f in lhs+rhs]) > 0:
            if set(lhs) == set(rhs):
                # isomerization reactions etc are identical
                for cpds in rxn.left_compounds:
                    cpds.coefficient = 1
                for cpds in rxn.right_compounds:
                    cpds.coefficient = 1
                rxn.balanced = True
                identical += 1
            else:
                b_lhs,b_rhs = balance_stoichiometry(set(lhs),set(rhs),underdetermined=1)                        
                for k,v in b_lhs.items():
                    if v>1:
                        for cpds in rxn.left_compounds:
                            if k == cpds.compound.formula:
                                cpds.coefficient = v
                for k,v in b_rhs.items():
                    if v>1:
                        for cpds in rxn.right_compounds:
                            if k == cpds.compound.formula:
                                cpds.coefficient = v
                rxn.balanced = True
                balanced += 1
        else:
            # There are undefined molecules in this equiation
            rxn.balanced = False
            undefined += 1
    except:
        # The system can not be simplified into integer coefficients
        rxn.balanced = False
        illogical += 1        
    
    all_reactions.append(rxn)

print(balanced,undefined,illogical,identical)

# from chempy import balance_stoichiometry
# balance_stoichiometry({u'C7H5O3-', u'O2', u'C21H27N7O14P2-2', u'H+'},{u'C7H5O4-', u'C21H26N7O14P2-', u'H2O'},underdetermined=1)

# balance_stoichiometry({'C22H16O8', u'O2', u'C21H26N7O17P3-4', u'H+'},
#                       {'C22H18O11', u'H2O', u'C21H25N7O17P3-3'},
#                       underdetermined=True)

# Save each reaction object as a json file

import ast
import json
outdir = '../../data/reaction_json4/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
for rxn in all_reactions:
    filename=os.path.join(outdir,rxn.source.replace(':','_')+'.json')
    d = ast.literal_eval(str(rxn).replace("nan","None")) #some smiles are nan and can't be evaled
    json_str = json.dumps(d, indent=4, sort_keys=True)#,default=default)
    with open(filename,'w') as fid:
        fid.write(unicode(json_str))
#     my_str = item.to_json(filename=os.path.join('../../data/reaction_json/',item.unique_id+'.json'))