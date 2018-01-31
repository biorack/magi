import pandas as pd
import re as re
import os
import sys
import imp
import json
import glob
import itertools
import multiprocessing as mp

from chempy import balance_stoichiometry  # Main reaction in NASA's booster rockets:


metatlas_rxn_path = '/Users/bpb/repos/build_magi_ref_data/'
sys.path.append(metatlas_rxn_path)

import refdata_tools.biopax as bp
import refdata_tools.data_loading as metdata

from local_settings import local_settings as settings_loc
my_settings = getattr(__import__('local_settings', fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)
sys.path.append(my_settings.path_to_magi)
import metatlas_reactions.reaction_objects as mr

from chempy import balance_stoichiometry

from cobra import Model, Reaction, Metabolite, Gene
from cobra.io import write_sbml_model
from cobra.io import save_json_model

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

json_path = '../../data/reaction_json/'
json_files = glob.glob(os.path.join(json_path,'*.json'))

my_file = json_files[100]
print(my_file)
with open(my_file,'r') as fid:
    d = json.load(fid)



def make_compound(dd):
    """
    makes a cobra compound object from dict
    """
    compound = Metabolite(id=dd['unique_id'],name=dd['name'],formula=dd['formula'],charge=dd['charge'])        
    compound.notes = dd
    if (compound.formula == '') & (compound.notes['inchi'] != ''):
        mol = Chem.MolFromInchi(compound.notes['inchi'])
        formula = rdMolDescriptors.CalcMolFormula(mol)
        compound.formula = formula
    return compound


def balance_compounds(mets):
    lhs = set()
    rhs = set()
    lhs_names = set()
    rhs_names = set()
    has_formula = True
    balanced = False
    for k,v in mets.items():
        if v<1:
            lhs.add(k.formula)
            lhs_names.add(k.name)
            if len(k.formula) == 0:
                has_formula = False
        else:
            rhs.add(k.formula)
            rhs_names.add(k.name)
            if len(k.formula) == 0:
                has_formula = False

    if has_formula:
        try:
            b_lhs,b_rhs = balance_stoichiometry(lhs,rhs)
            for k,v in mets.items():
                if v<1:
                    mets[k] = b_lhs[k.formula] * -1
                else:
                    mets[k] = b_rhs[k.formula]
            balanced = True
            changed_stoichiometry = False
            for k,v in mets.items():
                if abs(v) != 1:
                    changed_stoichiometry = True
            if changed_stoichiometry:
                print('############ Stoichiometry Updated #############')
                print(b_lhs,b_rhs)
        except ValueError:
            print('############ Unsatisfiable Reaction #############')
    else:
        print('############ Missing Chemical Formula #############')
    print(lhs,rhs)
    print(lhs_names,rhs_names)
    return mets,balanced




model = Model('magi reference data')

for my_file in json_files[:100]:
    print(my_file)
    with open(my_file,'r') as fid:
        d = json.load(fid)

    rxn = Reaction(id=d['unique_id'],
                   name=d['name'])
    
    rxn.notes = {k:v for (k,v) in d.items() if (not 'compounds' in k) & (not 'sequences' in k)}
    print(rxn.notes['metacyc_url'],rxn.notes['rhea_url'])
    mets = {}
    for cpds in d['left_compounds']:
        mets[make_compound(cpds['compound'])] = -1*cpds['coefficient']
    for cpds in d['right_compounds']:
        mets[make_compound(cpds['compound'])] = cpds['coefficient']
        
    mets,balanced = balance_compounds(mets)
    rxn.notes['balanced'] = balanced
    
    rxn.add_metabolites(mets)

    if len(d['sequences']) > 0:
        gene_str = '('
        gene_str = ' or '.join([dd['unique_id'] for dd in d['sequences']])
        gene_str = '( %s )'%gene_str
        rxn.gene_reaction_rule = gene_str
        for gene in rxn.genes:
            for magi_gene in d['sequences']:
                if gene.id == magi_gene['unique_id']:
                    gene.name = magi_gene['source']
                    gene.notes = magi_gene

    model.add_reactions([rxn])
    
    print(' ')




# Now there are things in the model
print('%i reaction' % len(model.reactions))
print('%i metabolites' % len(model.metabolites))
print('%i genes' % len(model.genes))    


#This doesn't export notes
# write_sbml_model(model, "/Users/bpb/Downloads/magi_ref_data.xml")

save_json_model(model, "/Users/bpb/Downloads/magi_ref_data.json")

with open('/Users/bpb/Downloads/magi_ref_data.json','r') as fid:
    d = json.load(fid)
# d

# # Best practise: SBML compliant IDs
# model = Model('example_model')

# reaction = Reaction('3OAS140')
# reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
# reaction.subsystem = 'Cell Envelope Biosynthesis'
# reaction.lower_bound = 0.  # This is the default
# reaction.upper_bound = 1000.  # This is the default

# ACP_c = Metabolite(
#     'ACP_c',
#     formula='C11H21N2O7PRS',
#     name='acyl-carrier-protein',
#     compartment='c')
# omrsACP_c = Metabolite(
#     '3omrsACP_c',
#     formula='C25H45N2O9PRS',
#     name='3-Oxotetradecanoyl-acyl-carrier-protein',
#     compartment='c')
# co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
# malACP_c = Metabolite(
#     'malACP_c',
#     formula='C14H22N2O10PRS',
#     name='Malonyl-acyl-carrier-protein',
#     compartment='c')
# h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
# ddcaACP_c = Metabolite(
#     'ddcaACP_c',
#     formula='C23H43N2O8PRS',
#     name='Dodecanoyl-ACP-n-C120ACP',
#     compartment='c')

# h_c.notes() #contains dict



# reaction.add_metabolites({
#     malACP_c: -1.0,
#     h_c: -1.0,
#     ddcaACP_c: -1.0,
#     co2_c: 1.0,
#     ACP_c: 1.0,
#     omrsACP_c: 1.0
# })

# reaction.reaction  # This gives a string representation of the reaction

# g = Gene()
# g.notes()

# reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
# list(reaction.genes)[0].name = 'asdf'
# for g in reaction.genes:
#     print(g.id,g.name)

# model.add_reactions([reaction])

# # Now there are things in the model
# print('%i reaction' % len(model.reactions))
# print('%i metabolites' % len(model.metabolites))
# print('%i genes' % len(model.genes))

# # Iterate through the the objects in the model
# print("Reactions")
# print("---------")
# for x in model.reactions:
#     print("%s : %s" % (x.id, x.reaction))

# print("")
# print("Metabolites")
# print("-----------")
# for x in model.metabolites:
#     print('%9s : %s' % (x.id, x.formula))

# print("")
# print("Genes")
# print("-----")
# for x in model.genes:
#     associated_ids = (i.id for i in x.reactions)
#     print("%s is associated with reactions: %s" %
#           (x.id, "{" + ", ".join(associated_ids) + "}"))

# model.add_reactions([reaction])

# # Now there are things in the model
# print('%i reaction' % len(model.reactions))
# print('%i metabolites' % len(model.metabolites))
# print('%i genes' % len(model.genes))

# # for i in range(100):
# #     with open(json_files[i],'r') as fid:
# #         d = json.load(fid)
# #     print(d['sequences'])


# r.traits()



# type(o)

# isinstance(o,mr.MetList)

# import time
# import datetime
# s = "2017-08-08T08:03:10"
# time.mktime(datetime.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S").timetuple())


# d['left_compounds'][0]['compound']

# my_file = json_files[0]
# print(my_file)
# with open(my_file,'r') as fid:
#     d = json.load(fid)
# # def reaction_dict_to_object(d):
# r = mr.Reaction()
# for k,o in r.traits().items():
#     if k in d.items():
#         if isinstance(o,mr.MetList):# we will do the lists last
#             pass
#         elif isinstance(o,mr.MetInt) & isinstance(d[k],str):
#             # its a date string that has to be converted to timestamp
#             setattr(r,k,time.mktime(datetime.datetime.strptime(d[k], "%Y-%m-%dT%H:%M:%S").timetuple()))
#         else:
#             setattr(r,k,d[k])
#     for reaction_compound in d['left_compounds']:
#         compound = mr.Compound(reaction_compound['compound'])
#         r.left_compounds.append(mr.StoichiometricCoefficient(coefficient=1,compound=compound))
# # return r
# # r2 = reaction_dict_to_object(d)
# print(r)


# # %%time
# mr = imp.reload(mr)
# my_file = json_files[0]
# print(my_file)

# r = mr.reaction_from_json(json_files[0])

# reac, prod = balance_stoichiometry({'NH4ClO4', 'Al'}, {'Al2O3', 'HCl', 'H2O', 'N2'},)

# def make_compound(dd):
#     """
#     makes a metatlas compound object from dict
#     """
#     compound = mr.Compound()
#     for k,o in compound.traits().items():
#         if k in dd:
#             if isinstance(o,mr.MetList):# compounds don't have any lists
#                 pass
#             elif isinstance(o,mr.MetInt) & isinstance(dd[k],str):
#                 # its a date string that has to be converted to timestamp
#                 setattr(compound,k,time.mktime(datetime.datetime.strptime(dd[k], "%Y-%m-%dT%H:%M:%S").timetuple()))
#             else:
#                 setattr(compound,k,dd[k])
#     return compound
# make_compound(d['left_compounds'][0]['compound'])

# %%time
# pool_size = 6
# p = mp.Pool(pool_size)

# pool_results = p.map(mr.reaction_from_json, json_files[:30000])
# p.close()
# p.join()

# rxn_ik=[]
# for r in pool_results:
#     rxn_ik.append(metdata.reaction_to_inchikey(r))

# print(len(rxn_ik),len(pd.unique(rxn_ik)))

# temp_ik = [k.replace('-I','-U').replace('-R','-U')[1:] for k in rxn_ik]
# print(len(temp_ik),len(pd.unique(temp_ik)))

# metdata = imp.reload(metdata)
# metdata.reaction_to_inchikey(r)

# for c in itertools.chain(r.left_compounds,r.right_compounds):
#     print(c.inchi)
#     print('')

