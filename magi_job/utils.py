"""
Utility functions for preprocessing magi jobs from webserver.
"""

from __future__ import print_function
import pandas as pd
import os
import subprocess
import json
import requests
import sys
# load local settings
sys.path.insert(
    0,
    '/project/projectdirs/metatlas/projects/metatlas_reactions/')
from local_settings import local_settings as settings_loc
my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)

def retrieve_jobs(
		username=my_settings.magiwebsuperuser,
		password=my_settings.magiwebsuperuserpass,
		base_url='https://magi-dev.nersc.gov/',
		admin_suffix='admin/login'
	):
	"""
	Retrieves jobs json from magi web.
	Output is a python list of dictionaries representing the json of all
	jobs.

    TODO: filter json by jobs that need to be run
	"""

	auth_url = base_url + admin_suffix

	client = requests.session()
	# Retrieve the CSRF token first
	client.get(auth_url,verify=False)  # sets cookie
	csrftoken = client.cookies['csrftoken']
	login_data = {'username': username,
	              'password': password,
	              'csrfmiddlewaretoken': csrftoken}
	client.post(auth_url, data=login_data,
		headers=dict(Referer=auth_url)
		)

	r = client.get(os.path.join(base_url,'admin/ids/?filter=all&json=True'))
	# Keep jobs that don't have a job script yet
	# load jobs json
	all_jobs = json.loads(r.text)
	return all_jobs
def mirror_inputs(all_jobs, verbose=False):
    for job in all_jobs:
        base_url='https://magi-dev.nersc.gov/files/input'
        dir_root='/project/projectdirs/metatlas/projects/magi_tasks'

        try:
            fasta_file = job['fields']['fasta_file'].split('input/')[1]
            metabolite_file = job['fields']['metabolite_file'].split(
                'input/')[1]
        except IndexError:
            print('WARNING: job %s input path not properly formatted; skipping' 
                % (job['pk']))
            continue
        job_path = '/'.join(fasta_file.split('/')[:-1])
        job_path = os.path.join(dir_root, job_path)

        if os.path.isfile(os.path.join(dir_root, fasta_file)) and \
            os.path.isfile(os.path.join(dir_root, metabolite_file)):
            if verbose:
                print('job inputs exist for %s' % (job['pk']))
            continue

        if not os.path.isdir(job_path):
            if verbose:
                print('making %s' % (job_path))
            os.makedirs(job_path)

        if not os.path.isfile(os.path.join(dir_root, fasta_file)):
            if verbose:
                print('getting %s' % os.path.join(base_url, fasta_file))
            fasta_data = requests.get(os.path.join(base_url, fasta_file))
            if verbose:
                print('writing %s' % os.path.join(dir_root, fasta_file))
            with open(os.path.join(dir_root, fasta_file), 'w') as f:
                f.write(fasta_data.text)
        if not os.path.isfile(os.path.join(dir_root, metabolite_file)):
            if verbose:
                print('getting %s' % (os.path.join(base_url, metabolite_file)))
            metabolite_data = requests.get(os.path.join(base_url,
                                           metabolite_file))
            if verbose:
                print('writing %s' % (os.path.join(dir_root, metabolite_file)))
            with open(os.path.join(dir_root, metabolite_file), 'w') as f:
                f.write(metabolite_data.text)
    return None

def jobs_to_script(
		all_jobs,
		dir_root='/project/projectdirs/metatlas/projects/magi_tasks'
	):
	"""
	Determine which jobs need an SBATCH script.
	
	Inputs
	all_jobs: a list of dictionaries that represents the json of magi
			  web jobs.
	dir_root: the full path prefix to where the input files are located
			  at NERSC.

	Output is a list of python dictionaries representing the json of all
	jobs that do not have an SBATCH script already written for them, and
	a boolean describing whether or not there are any jobs that require
	accurate mass searching (so that main program can load up the
	compound dataframe for accurate mass searching)
	"""

	to_script = []
	mass_search = []
	for i, job in enumerate(all_jobs):
	    year = job['fields']['uploaded_at'].split('-')[0]
	    month = job['fields']['uploaded_at'].split('-')[1]
	    pk = job['pk']
	    job_path = os.path.join(dir_root, year, month, pk)
	    if not os.path.isfile(os.path.join(job_path, 'job_script.sbatch')):
	        # adjust paths
	        try:
	            job['fields']['fasta_file'] = os.path.join(dir_root, job['fields']['fasta_file'].split('input/')[1])
	            job['fields']['metabolite_file'] = os.path.join(dir_root, job['fields']['metabolite_file'].split('input/')[1])
	        except IndexError:
	            print('WARNING: could not find input file(s) for job %s index %s; not making a job script' % (pk, i))
	            continue
	        
	        # check to make sure the files exist
	        if not os.path.isfile(job['fields']['fasta_file']):
	            print('WARNING: fasta file for job %s does not exist; not making a job script' % (pk))
	            continue
	        if not os.path.isfile(job['fields']['metabolite_file']):
	            print('WARNING: fasta file for job %s does not exist; not making a job script' % (pk))
	            continue
	        
	        # add job json to todo list
	        mass_search.append(job['fields']['is_mass_search'])
	        to_script.append(job)
	mass_search = any(mass_search)
	return to_script, mass_search

def protein_translate(seq, warnings=False):
    """
    Translates an input DNA or RNA sequence to amino acid

    Inputs
    ------
    seq: DNA or RNA sequence
    warnings: print out warnings or not
    na: nucleic acid language('dna' or 'rna')
    
    Outputs
    -------
    protein_seq: Translated protein sequence
    """
    codons = {
		'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T',
		'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
		'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H',
		'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
		'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L',
		'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
		'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G',
		'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
		'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S',
		'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
		'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
        }
    # check input sequence
    if len(seq) % 3 != 0:
        raise RuntimeError(
        	'The nucelotide sequence is not a multiple of 3: %s' % (seq)
        	)
    
    # replace all Uracil with Thymine
    seq = seq.upper()
    seq = seq.replace('U', 'T')

    # make index tuples list
    cidx = range(0, len(seq)+3, 3)
    cidx = zip(cidx[:-1], cidx[1:])
    
    # translate
    protein_seq = ''
    for i, j in cidx:
        codon = seq[i:j]
        aa = codons[codon]
        protein_seq += aa
    
    # check translation for suspicious activity
    if warnings:
        if protein_seq[0] != 'M':
            print('WARNING: sequence does not start with Methionine')
        if protein_seq[-1] != '*':
            print('WARNING: sequence does not end with a STOP codon')
        if '*' in protein_seq[:-1]:
            counter = []
            for i, aa in enumerate(protein_seq[:-1]):
                if aa == '*':
                    counter.append(i * 3)
            print(
                'WARNING: sequence contains internal STOP codons at DNA positions: %s'
                 % (counter))
    
    # return the translation
    return protein_seq

def determine_fasta_language(job_data, translate=True):
    """
    Assuming this is a valid fasta file,
    determines if the fasta is protein or nucleic acid
    if translate = True, translates the file to protein if it is dna/rna
    
    Input is dict of json of one job.
    
    alters the fasta_file field in job_data json if translation occured.
    
    returns job_data
    """
    file_path = job_data['fields']['fasta_file']
    # read fasta file
    with open(file_path, 'r') as f:
        file_data = f.read()
        
    # convert newlines
    for newline in ['\r\n', '\n\r']:
        if newline in file_data:
            file_data.replace(newline, '\n')
    if '\r' in file_data:
        file_data.replace('\r', '\n')
    
    # parse gene sequences into one long string and convert to DNA if RNA
    genes = file_data.split('>')[1:]
    seqs = ''.join([''.join(g.split('\n')[1:]).replace('-', '') for g in genes]).upper().replace('U', 'T')
    
    # determine what letters are in the gene sequences
    letters = pd.Series(list(set(seqs)))
    # first check for DNA
    if len(letters) <= 4 and letters.str.contains('[ACTG]').all():
        answer = 'dna'
    # then amino acids 
    elif letters.str.contains('[ACDEFGHIKLMNPQRSTVWY*]').all():
        answer = 'protein'
    else:
        no = letters[~letters.str.contains('[ACDEFGHIKLMNPQRSTVWY*]')].values
        raise RuntimeError('Could not determine if FASTA is nucleotide or protein. Offending character(s): %s' % (no))
    
    # translate if desired
    if translate and answer == 'dna':       
        # translate the genes
        new_data = ''
        gene_list = file_data.split('>')[1:]
        for gene in gene_list:
            header = gene.split('\n')[0]
            seq = gene.split('\n')[1:]
            seq = ''.join([i for i in seq if i != ''])
            protein = protein_translate(seq)
            new_data += '>' + header + '\n'
            new_data += protein + '\n\n'
        
        # save the new file
        new_filename = file_path.split('/')[-1].split('.')[0] + '_translated.faa'
        new_filepath = '/'.join(file_path.split('/')[:-1]) + '/%s' %(new_filename)
        with open(new_filepath, 'w') as f:
            f.write(new_data)
        # change the field in job_data
        job_data['fields']['fasta_file'] = new_filepath
        return job_data
    else:
        return job_data

def job_script(job_data):
    """
    uses job data json to create a magi job submission script for cori
    """
    account_id = 'm2650' # metatlas
    partition = 'debug'
    
    # where to write the job script to
    out_path = '/'.join(job_data['fields']['fasta_file'].split('/')[:-1])
    
    ###################
    # TEMPORARY BLOCK #
    ###################
    job_data['fields']['score_weights'] = [1, 2, 3, 4]
    job_data['fields']['blast_cutoff'] = 100
    job_data['fields']['reciprocal_cutoff'] = 75
    job_data['fields']['chemnet_penalty'] = 3
    ###################
    # TEMPORARY BLOCK #
    ###################

    # need to convert this into string so we can join it later
    job_data['fields']['score_weights'] = [str(i) for i in job_data['fields']['score_weights']]
    
    # need to interpret tautomer flag
    if  job_data['fields']['is_tautomers']:
        tautomer = '--tautomer'
    else:
        tautomer = ''

    job_lines = [
        '#!/bin/bash -l',
        '#SBATCH --account=%s' % (account_id),
        '#SBATCH --job-name=%s' % (job_data['pk'].split('-')[0]),
        '#SBATCH --time=0:10:00',
        '#SBATCH --output=%s/log_out.txt' % (out_path),
        '#SBATCH --error=%s/log_err.txt' % (out_path),
        '#SBATCH --partition=%s' % (partition),
        '#SBATCH --constraint=haswell',
        '#SBATCH --license=project',
        '#SBATCH --mail-user=%s' %(job_data['fields']['email']),
        '#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT',
        '',
        'module load python/2.7-anaconda',
        '',
        'time python /project/projectdirs/metatlas/projects/metatlas_reactions/workflow/magi_workflow_20170519.py \\',
        '--fasta %s \\' % (job_data['fields']['fasta_file']),
        '--compounds %s \\' % (job_data['fields']['metabolite_file']),
        '--level %s \\' % (job_data['fields']['network_level']),
        # not sure if this line will break anything at nersc
        # if it does, put it at the end of the previous line
        '%s \\' % (tautomer),
        '--final_weights %s \\' % (' '.join(job_data['fields']['score_weights'])),
        '--blast_filter %s \\' % (job_data['fields']['blast_cutoff']),
        '--reciprocal_closeness %s \\' % (job_data['fields']['reciprocal_cutoff']),
        '--chemnet_penalty %s \\' % (job_data['fields']['chemnet_penalty']),
        '--output %s/output_files --mute' % (out_path),
    ]
    
    job = '\n'.join(job_lines)

    # change umask temporarily; don't want job script to be world-read
    old_mask = os.umask(007)
    # write job
    with open(os.path.join(out_path, 'job_script.sbatch'), 'w') as f:
        f.write(job)
    os.umask(old_mask)
    return None

def ppm_window(mass, ppm=5, result='bounds'):
    """
    Given a mass and a ppm error, returns lower and upper bounds 
    corresponding to that ppm window. 
    Inputs
    ------
    mass:   monoisotopic mass
    ppm:    the ppm error to construct the window
    result: "bounds" returns a list of lower and upper bounds
            "error" returns the amu value of the error, if you wanted to
            add and/or subtract the error as you wish
    Outputs
    -------
    Either a list of lower and upper mass bounds, or a single value
    corresponding to the amu error
    """
    error = ppm/1e6 * mass
    lower_bound = mass - error
    upper_bound = mass + error
    if result.lower() == 'bounds':
        return [lower_bound, upper_bound]
    elif result.lower() == 'error':
        return error
    else:
        raise RuntimeError(
            '%s is not a valid result argument' %(result)
            )

def accurate_mass_match(mass, compound_df=None, ppm=5, extract='inchi_key'):
    """
    Accurate mass searching against a compound database.
    Inputs
    ------
    mass:           An accurate monoisotopic mass
    compound_df:    A dataframe of compounds. Must have a column named
                    "mono_isotopic_molecular_weight"
    ppm:            ppm error to allow the mass matcher
    extract:        What compound information to return. Must correspond
                    to a valid column in compound_df
    
    Outputs
    -------
    cpd:            List of compounds that were matched
    """

    err = ppm_window(mass, ppm=ppm, result='error')

    potential_compounds = compound_df[
                            abs(compound_df['mono_isotopic_molecular_weight'] \
                            - mass) <= err]
    cpds = potential_compounds[extract].values.tolist()
    if len(cpds) == 0:
        cpds = None
    return cpds

def accurate_mass_search_wrapper(job_data, reference_compounds):
    """
    performs accurate mass search using unique_compounds table
    """

    # math table for all adducts
    adduct_table = {
        'p1' : -1.007276,
        'p2' : -18.033823,
        'p3' : -22.989218,
        'p4' : -38.963158,
        'p8' : - 42.033823,
    }

    search_ppm = job_data['fields']['ppm']

    # get adducts according to polarity
    if job_data['fields']['polarity'] == 'pos':
        adducts = job_data['fields']['adducts_pos'].split(',')
    elif job_data['fields']['polarity'] == 'pos':
        adducts = job_data['fields']['adducts_neg'].split(',')
    else:
        raise RuntimeError('Could not understand polarity')

    # load compound table (should be masses in original_compounds)
    compounds = pd.read_csv(job_data['fields']['metabolite_file'])

    # rename original_compounds column
    columns = compounds.columns.values
    columns[columns == 'original_compound'] = 'original_mz'
    compounds.columns = columns

    # set up data container
    data = {
        'original_compound': [],
        'searched_adduct': [],
        'original_mz' : []
    }
    # accurate mass search and store results
    for mz in compounds['original_mz'].unique():
        for adduct in adducts:
            neutral_mass = mz + adduct_table[adduct]
            found_compounds = accurate_mass_match(neutral_mass,
                                                  compound_df=reference_compounds,
                                                  ppm=search_ppm
                                                 )
            if found_compounds is not None:
                for cpd in found_compounds:
                    data['original_compound'].append(cpd)
                    data['searched_adduct'].append(adduct)
                    data['original_mz'].append(mz)

    # merge with user input and save
    df = pd.DataFrame(data)
    compounds = compounds.merge(df, on='original_mz', how='left')
    
    # save the new table
    new_path = job_data['fields']['metabolite_file'].split('.')[0] + '_mass_searched.csv'
    job_data['fields']['metabolite_file'] = new_path
    compounds.to_csv(job_data['fields']['metabolite_file'])

    return job_data

