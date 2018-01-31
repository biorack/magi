import pandas as pd
import re as re
import requests
from bs4 import BeautifulSoup as soup
import os
import sys
import hashlib
import imp
import pickle

metatlas_rxn_path = '/Users/bpb/repos/build_magi_ref_data/'
sys.path.append(metatlas_rxn_path)

import refdata_tools.data_loading as metdata

# from local_settings import local_settings as settings_loc
# my_settings = getattr(__import__('local_settings', fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)
# sys.path.append(my_settings.path_to_magi)
# import metatlas_reactions.reaction_objects as mr





xref_df = pd.read_csv('../../data/rhea/rhea2xrefs.tsv',sep='\t')

rhea_uniprot_ids = xref_df[xref_df.DB=='UNIPROT'].ID.unique()
print(len(rhea_uniprot_ids))

# path='../../data/rhea/uniprot_sequences/'
# """
# Redo-pickle everything to python2.  Must be run in python3.

# """
# for r in rhea_uniprot_ids:
#     filename = os.path.join(path,r+'.pkl')
#     if os.path.isfile(filename):
#         with open(filename,'rb') as fid:
#             fasta = pickle.load(fid)
#         with open(filename,'wb') as fid:
#             pickle.dump({'header':fasta['header'],'sequence':fasta['sequence']},
#                         fid,
#                         protocol=2)

"""
These have 'javascript' in them:

../../data/rhea/uniprot_sequences/B1JGA5.pkl
../../data/rhea/uniprot_sequences/Q9K8E7.pkl
../../data/rhea/uniprot_sequences/A1RJ17.pkl
../../data/rhea/uniprot_sequences/C1AK32.pkl
../../data/rhea/uniprot_sequences/P48848.pkl
../../data/rhea/uniprot_sequences/A7X2S3.pkl
../../data/rhea/uniprot_sequences/Q6WI70.pkl
../../data/rhea/uniprot_sequences/B7J921.pkl
../../data/rhea/uniprot_sequences/Q977Z3.pkl


These have 'html' in them:
../../data/rhea/uniprot_sequences/B1JGA5.pkl
../../data/rhea/uniprot_sequences/Q9K8E7.pkl
../../data/rhea/uniprot_sequences/A1RJ17.pkl
../../data/rhea/uniprot_sequences/C1AK32.pkl
../../data/rhea/uniprot_sequences/P48848.pkl
../../data/rhea/uniprot_sequences/A7X2S3.pkl
../../data/rhea/uniprot_sequences/Q6WI70.pkl
../../data/rhea/uniprot_sequences/B7J921.pkl
../../data/rhea/uniprot_sequences/Q977Z3.pkl

"""

files = ['../../data/rhea/uniprot_sequences/B1JGA5.pkl',
'../../data/rhea/uniprot_sequences/Q9K8E7.pkl',
'../../data/rhea/uniprot_sequences/A1RJ17.pkl',
'../../data/rhea/uniprot_sequences/C1AK32.pkl',
'../../data/rhea/uniprot_sequences/P48848.pkl',
'../../data/rhea/uniprot_sequences/A7X2S3.pkl',
'../../data/rhea/uniprot_sequences/Q6WI70.pkl',
'../../data/rhea/uniprot_sequences/B7J921.pkl',
'../../data/rhea/uniprot_sequences/Q977Z3.pkl']
for filename in files:
    uniprot = os.path.basename(filename).split('.')[0]
    print(uniprot)
    fasta_header,fasta_seq = metdata.get_uniparc_from_uniprot(uniprot)
    print(fasta_header,fasta_seq)
    with open(filename,'wb') as fid:
        pickle.dump({'header':fasta_header,'sequence':fasta_seq},
                            fid,
                            protocol=2)

path='../../data/rhea/uniprot_sequences/'
# """
# Check everyfile for scraping error

# """
for r in rhea_uniprot_ids:
    filename = os.path.join(path,r+'.pkl')
    if os.path.isfile(filename):
        with open(filename,'rb') as fid:
            fasta = pickle.load(fid)
        if 'html' in fasta['sequence']:
            print filename

# """
# Test that all the sequences are readable.  This is something to be run very carefully
# """

# def get_sequence(uniprot,path='../../data/rhea/uniprot_sequences/'):
#     filename = os.path.join(path,uniprot+'.pkl')
#     try:
#         with open(filename,'rb') as fid:
#             fasta = pickle.load(fid)
# #             print('worked',fasta['header'])
#     except:
#         try:
#             fasta_header,fasta_seq = metdata.get_uniparc_from_uniprot(uniprot)
#             with open(filename,'wb') as fid:
#                 pickle.dump({'header':fasta_header,'sequence':fasta_seq},
#                             fid,
#                             protocol=2)
#                 print('got new one')
#         except:
#             print(filename,'failed')

# # for r in rhea_uniprot_ids[:1000]:
# #     get_sequence(r)

# import multiprocessing as mp
# pool_size = 10
# p = mp.Pool(pool_size)
# p.map(get_sequence, rhea_uniprot_ids)
# p.close()
# p.join()


mr.database.executable

mr.store([mr.ProteinSequence()])

prots = []
for uniprot in rhea_uniprot_ids:
    path='../../data/rhea/uniprot_sequences/'
    filename = os.path.join(path,uniprot+'.pkl')
    with open(filename,'rb') as fid:
        fasta = pickle.load(fid)
    premade_prot = mr.retrieve('proteinsequences',md5=hashlib.md5(fasta['sequence'].encode('utf-8')).hexdigest())
    if len(premade_prot) == 0:
        prots.append(mr.ProteinSequence(sequence=fasta['sequence'],
                                        source=fasta['header'][1:],
                                        md5=hashlib.md5(fasta['sequence'].encode('utf-8')).hexdigest()))

mr.store(prots)