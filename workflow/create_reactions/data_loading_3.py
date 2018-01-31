from __future__ import print_function

import pandas as pd
import re as re
import requests
import time
import sys
import os
import copy

from bs4 import BeautifulSoup
from rdkit import Chem

from collections import defaultdict

import getpass
import uuid
import pprint
import json

import six

import functools
# import dataset as dataset

# Use python 3 style file open
from io import open

REACTION_DIRECTION = ('undefined','left-to-right','reversible')

class RxnObject():
    def __init__(self, **kwargs):
        """Set the default attributes."""
        super(MetatlasObject, self).__init__(**kwargs)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        names = sorted(dir(self))
        names = [n for n in names if not n.startswith('_')]
        list_of_objects = []
        for n in names:
            t = type(getattr(self,n))
            list_of_objects.append((n,getattr(self,n)))
        
        state = dict(list_of_objects)
        return pprint.pformat(state,width=10000)


class Compound(RxnObject):
    """

    """
    def __init__(self, **kwargs):
        self.isomeric_smiles = None
        self.standardized_inchikey = None
        self.standardized_smiles = None
        self.chebi_id = None
        self.metacyc_id = None
        self.formula = None
        self.name = None
 
class ProteinSequence(RxnObject):
    """
    Amino acid sequence.  
    """
    def __init__(self, **kwargs):
        self.sequence = None
        self.source = None
        self.md5 = None


class StoichiometricCoefficient(RxnObject):
    """
    An integer coefficient for reaction stoichiometry
    associate
    """
    def __init__(self, **kwargs):
        self.coefficient = None
        self.compound = Compound


class Reaction(RxnObject):
    """
    A Reaction ties together one or more reference sequences to compounds.
    """
    def __init__(self, **kwargs):
        self.sequences = []#[ProteinSequence) #All sequences that are part of the potential enzyme complex for this reaction from any species
        self.left_compounds = []
        self.right_compounds = []
        self.direction = None#MetEnum(REACTION_DIRECTION,'undefined',help='directionality of reaction')
        self.source = None#MetUnicode(help='')
        self.balanced = False#MetBool(False,help='is reaction stoichiometry balanced')
        self.name = None
        self.description = None


def get_uniprot_id_from_rhea_reaction(rxn):
    """
    
    """
    prots = []
    for x in rxn['xref']:
        if 'UNIPROT:' in x:
            prots.append(x.split(':')[-1])
    return prots

def get_direction_from_rhea_reaction(rxn):
    """
    Rhea reactions are labeled as 'bidirectional', 'left to right', 'right to left', and 'undefined'
    
    Possibly, a good idea to get polymerization here too, but determine later
    # comment {'RHEA:Direction=bidirectional', 'RHEA:Formuled=true', 'RHEA:Polymerization=false', 'RHEA:Transport=false', 'RHEA:Class of reactions=false', 'RHEA:Chemically balanced=true', 'RHEA:Mapped=true', 'RHEA:Status=approved'}
    """
    return [d.split('=')[-1] for d in rxn['comment'] if 'RHEA:Direction' in d][-1]

def get_EC_from_rhea_reaction(rxn):
    """
    
    """
    return rxn['EC']

def get_rhea_id_from_rhea_reaction(rxn):
    """
    grabs a rhea id
    """
    return rxn['rdf_name'].split('/')[-1]

def get_rhea_reaction_text_from_rhea_reaction(rxn):
    """
    this grabs a simple text of the reaction with common names
    """
    return rxn['name']

def get_compounds_from_rhea_reaction(rxn,metabolites,side='left'):
    """
    Inputs:
    rxn:
        a parsed xml from the biopax file
    metabolites:
        parsed metabolites from the xml biopax file
     side:
        'left' or 'right'
    
    
    
    """
    out_compounds = []
    for k,v in rxn[side].items():
        compound = {}
        chebi_id = metabolites[k]['xref'][0].replace('_',':').split('/')[-1]
        compound['stoichiometry'] = v
        compound['name'] = metabolites[k]['name']
        compound['isomeric_smiles'] = None
        compound['formula'] = None
        if 'CHEBI' in chebi_id:
            compound['chebi_id'] = chebi_id
        else:
            compound['chebi_id'] = None

        out_compounds.append(compound)
    return out_compounds



def parse_biocyc_flat_file(my_files,attributes =  ['UNIQUE-ID','COMMON-NAME']):
    """
     We will only use the following

     ['UNIQUE-ID',
     'COMMON-NAME']

    """
    header_lines = []
    all_lines = []
    for my_file in my_files:
        with open(my_file,'r',encoding='utf-8', errors='ignore') as fid:
            for line in fid.readlines():
                if line.startswith('#'):
                    header_lines.append(line)
                else:
                    all_lines.append(line)
    header_lines = '\n'.join(header_lines)
    all_entities = '\n'.join(all_lines).split('//\n')

    # attributes = [a.strip() for a in header_lines.split('Attributes:\n')[-1].replace('#','').split('\n')]
    # attributes = [a for a in attributes if len(a)>0 ]
    
    empty_dict = {}
    for a in attributes:
        empty_dict[a] = ''

    all_dicts = []
    for a in all_entities:
        new_entity = copy.deepcopy(empty_dict)
        for a_attr in a.split('\n'):
            a_attr_key_value = a_attr.split(' - ')
            if a_attr_key_value[0] in new_entity:
                try:
                    new_entity[a_attr_key_value[0]].append(a_attr_key_value[1])
                except:
                    new_entity[a_attr_key_value[0]] = [a_attr_key_value[1]]
        all_dicts.append(new_entity)
    df = pd.DataFrame(all_dicts)
    return df

def scrape_sequences_not_in_flatfile(monomer_id,
                                     monomer_path='/users/bpb/repos/build_magi_ref_data/data/tier1/monomer_sequences',
                                    do_print=False):

    """
    The majority of protein subunits do not contain sequences in the flat file download.
    You must scrape them from the website in at least two steps
    1) Get dblinks from biocyc API
    2) Get sequence from uniprot ID
    3) Use other IDs contained in dblinks (two ncbi options are implemented here)

    Get missing ones by:
    https://biocyc.org/gene?orgid=META&id=MONOMER-20281
    and search for things like this:

    <dblink>
    <dblink-db>UNIPROT</dblink-db>
    <dblink-oid>I3ZNU9</dblink-oid>
    <dblink-relationship>unification</dblink-relationship>
    <dblink-URL>http://www.uniprot.org/uniprot/I3ZNU9</dblink-URL>
    </dblink>


    another scenario is that dblinks will look like this:
    http://www.ncbi.nlm.nih.gov/protein/AAY37623

    Seems about 1 in 20

    These fasta sequences can be obtained slightly differently

    https://www.ncbi.nlm.nih.gov/protein/AAY37623?report=fasta&format=text

    This doesn't work directly.  It opens a link to redirect to a new URL:

    https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=AAY37623&db=protein&report=fasta&extrafeat=0&fmt_mask=0&maxplex=1&sendto=t&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000

    Seems about 1 in 50 missing ones can be obtained from NCBI:

        https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?val=NP_416397&db=protein&report=fasta&extrafeat=0&fmt_mask=0&maxplex=1&sendto=t&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000

    Get missing ones by:
    https://biocyc.org/gene?orgid=META&id=PHOSPHO-CHEB
    and search for things like this:

    <dblink-URL>
        http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=NP_416397
    </dblink-URL>

    """
    
    if not os.path.exists(monomer_path):
        os.makedirs(monomer_path)
    
    monomer_file = os.path.join(monomer_path,'%s.txt'%monomer_id)
    
    if not os.path.isfile(monomer_file):
        fasta_seq = ''
        r = requests.get('https://websvc.biocyc.org/getxml?META:%s'%monomer_id)
        uniprot_str = 'www.uniprot.org/uniprot'
        uniprot_url_re = re.findall(uniprot_str + '[^<]+',r.text)
        entrez_str = 'www.ncbi.nlm.nih.gov/protein/'
        entrez_url_re = re.findall(entrez_str + '[^<]+',r.text)
        entrez_str_2 = 'www.ncbi.nlm.nih.gov/entrez/viewer.fcgi'
        entrez_url_2_re = re.findall(entrez_str + '[^<]+',r.text)
        r.close()

        if len(uniprot_url_re)>0:
            url = 'https://'+uniprot_url_re[0] + '.fasta'
            r2 = requests.get(url)
            fasta_header,fasta_seq = fasta_text_to_seq(r.text)
            r2.close()
        elif len(entrez_url_re)>0:
            url_params = '&db=protein&report=fasta&extrafeat=0&fmt_mask=0&maxplex=1&sendto=t&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000'
            url = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id='+entrez_url_re[0].strip().split('/')[-1] + url_params
            r2 = requests.get(url)
            fasta_header,fasta_seq = fasta_text_to_seq(r2.text)
            r2.close()
        elif len(entrez_url_2_re)>0:
            url_params = '&db=protein&report=fasta&extrafeat=0&fmt_mask=0&maxplex=1&sendto=t&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000'
            url = 'https://'+entrez_url_2_re[0].strip() + url_params
            r2 = requests.get(url)
            fasta_header,fasta_seq = fasta_text_to_seq(r2.text)
            r2.close()  
        if do_print:
            print(fasta_seq)
        with open(monomer_file,'w') as fid:
            fid.write(fasta_seq)

def fasta_text_to_seq(my_str):
    """

    """
    split_str = my_str.strip().split('\n')
    fasta_header = split_str[0]
    fasta_seq = ''.join(split_str[1:])
    return fasta_header,fasta_seq

def get_uniparc_from_uniprot(uniprot_id):
    """

    """
    r = requests.get('http://www.uniprot.org/uniparc/?query='+uniprot_id)
    bs = BeautifulSoup(r.content,"lxml")
    result = bs.find_all("tr", {"class": " entry selected-row"})
    uniparc_ids = [s.replace('id="','') for s in re.findall(r'id="[^"]+',unicode.join(u'\n',map(unicode,result)))]
    r = requests.get('http://www.uniprot.org/uniparc/' + uniparc_ids[0] + '.fasta')
    fasta_header,fasta_seq = fasta_text_to_seq(r.text)
    return fasta_header,fasta_seq
