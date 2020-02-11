"""
Utility functions for preprocessing magi jobs from webserver.
"""

from __future__ import print_function
import pandas as pd
import os
import subprocess
import json
import requests
import smtplib
from email.mime.text import MIMEText
import sys

# load local settings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from local_settings import local_settings as settings_loc
my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)

MAGI_EMAIL = my_settings.admin_email

def retrieve_jobs(
        username=my_settings.magiwebsuperuser,
        password=my_settings.magiwebsuperuserpass,
        base_url=my_settings.magiweburl,
        admin_logon='admin/login/',
        sift=[('filter', 'all')]
    ):
    """
    Retrieves jobs json from magi web.
    Output is a python list of dictionaries representing the json of all
    jobs.

    sift: list of tuples to construct url query filter.
          first element of tuple is the variable
          second element is the value
          example [('filter', 'all')]:
              admin/ids/?filter=all&json=True
          example [('pk', 'b823fd89-ed34-48d6-a657-258c90542088')]:
              admin/ids/?pk=b823fd89-ed34-48d6-a657-258c90542088&json=True
          example [('email', 'magi_web@lbl.gov'), ('year_lte', 2015)]:
              admin/ids/?email=magi_web@lbl.gov&uploaded_at__year__lte=2015&json=True

    """
    if not isinstance(sift, list):
        raise ValueError('sift argument must be a list of tuples or lists')
    for f in sift:
        if not (isinstance(sift, list) or isinstance(sift, tuple)):
            raise ValueError('sift argument must be a list of tuples or lists')

    auth_url = base_url + admin_logon

    client = requests.session()
    # Retrieve the CSRF token first
    client.get(auth_url)#,verify=False)  # sets cookie
    csrftoken = client.cookies['csrftoken']
    login_data = {'username': username,
                  'password': password,
                  'csrfmiddlewaretoken': csrftoken}
    r = client.post(auth_url, data=login_data,
        headers=dict(Referer=auth_url)
        )
    if r.status_code not in [200, 404]:
        raise RuntimeError(
            'Could not authenticate; status code %s' % (r.status_code))
    
    # build filter string
    filter_string = ''
    for f in sift:
        filter_string += '%s=%s&' %(f[0], f[1])

    get_url = os.path.join(base_url,'admin/ids/?%sjson=True' % (filter_string))
    r = client.get(get_url)
    if r.status_code not in [200]:
        raise RuntimeError(
            'GET request failed; status code %s; url: %s'
            % (r.status_code, get_url))
    # Keep jobs that don't have a job script yet
    # load jobs json
    try:
        all_jobs = json.loads(r.text)
    except ValueError as e:
        if 'No JSON object could be decoded' in e.args:
            return None
        else:
            raise e
    if len(all_jobs) == 0:
        return None
    return all_jobs

def change_params(
        job_id, field, value,
        username=my_settings.magiwebsuperuser,
        password=my_settings.magiwebsuperuserpass,
        base_url=my_settings.magiweburl,
        admin_logon='admin/login/',
        post_suffix='admin/',
        verify=True
        ):
    """
    changes a parameter in the data model at server via POST.

    Inputs
    ------
    job_id: uuid of the magi job
    field: field in the data model to change, e.g. "runflag"
    value: value to change the field to, e.g. "True"
    verify: retrieves the job and verifies that the field was changed
    """

    auth_url = base_url + admin_logon

    client = requests.session()
    # authenticate
    # Retrieve the CSRF token first
    client.get(auth_url,verify=False)  # sets cookie
    csrftoken = client.cookies['csrftoken']
    login_data = {'username': username,
                  'password': password,
                  'csrfmiddlewaretoken': csrftoken}
    r = client.post(auth_url, data=login_data,
        headers=dict(Referer=auth_url)
        )
    if r.status_code not in [200, 404]:
        raise RuntimeError(
            'Could not authenticate; status code %s' % (r.status_code))

    # change parameters
    # Retrieve new CSRF token
    client.get(auth_url,verify=False)  # sets cookie
    csrftoken = client.cookies['csrftoken']
    post_url = auth_url = base_url + 'admin/'
    payload = {
        'job_id': job_id,
        'field': field,
        'value': value,
        'csrfmiddlewaretoken': csrftoken
    }
    r = client.post(post_url, data=payload,
        headers=dict(Referer=auth_url)
        )

    if r.status_code not in [200]:
        raise RuntimeError(
    'Change job POST request failed with status code %s. URL: %s; Payload: %s'
    % (r.status_code, post_url, payload)
    )

    if verify:
        # verify that the field changed
        # get the one job
        job = retrieve_jobs(
            username=username,
            password=password,
            base_url=base_url,
            admin_logon=admin_logon,
            sift=[('pk', job_id)]
            )
        if job is None:
            return False
        if job[0]['fields'][field] == value:
            return True
        else:
            return False
    else:
        return None

def get_job_dir(job_json):
    """
    uses upload time to get what the job directory path should be:
    year/month/primary_key
    """
    uptime = job_json['fields']['uploaded_at']
    y = uptime.split('-')[0]
    m = uptime.split('-')[1]
    pk = job_json['pk']
    job_path = os.path.join(y, m, pk)
    return job_path

def mirror_inputs(all_jobs,
    base_url=my_settings.magiweburl,
    dir_root=my_settings.magi_task_path,
    verbose=False):
    """
    dir_root is where to mirror the files to
    """
    base_url = os.path.join(base_url, 'files', 'input')

    for job in all_jobs:
        fasta_exist = False
        met_exist = False
        if job['fields']['fasta_file'] != '':
            fasta_exist = True
        if job['fields']['metabolite_file'] != '':
            met_exist = True
        if not fasta_exist and not met_exist:
            if verbose:
                print('no input files for %s' % (job['pk']))
            continue
        try:
            if fasta_exist:
                fasta_file = job['fields']['fasta_file'].split('input/')[1]
            if met_exist:
                metabolite_file = job['fields']['metabolite_file'].split(
                'input/')[1]
        except IndexError:
            print('WARNING: job %s input path not properly formatted; skipping' 
                % (job['pk']))
            continue
        job_path = get_job_dir(job)
        job_path = os.path.join(dir_root, job_path)

        # make directory
        if not os.path.isdir(job_path):
            if verbose:
                print('making %s' % (job_path))
            os.makedirs(os.path.join(job_path, 'admin'))

        # copy job's json
        with open(os.path.join(job_path, 'admin', '%s.json' % (job['pk'])), 'w') as f:
            json.dump(job, f)

        # make file if does not exist
        if fasta_exist and not os.path.isfile(os.path.join(dir_root, fasta_file)):
            if verbose:
                print('getting %s' % os.path.join(base_url, fasta_file))
            fasta_data = requests.get(os.path.join(base_url, fasta_file))
            if verbose:
                print('writing %s' % os.path.join(dir_root, fasta_file))
            with open(os.path.join(dir_root, fasta_file), 'w') as f:
                f.write(fasta_data.text)
        if met_exist and not os.path.isfile(os.path.join(dir_root, metabolite_file)):
            if verbose:
                print('getting %s' % (os.path.join(base_url, metabolite_file)))
            metabolite_data = requests.get(os.path.join(base_url,
                                           metabolite_file))
            if verbose:
                print('writing %s' % (os.path.join(dir_root, metabolite_file)))
            with open(os.path.join(dir_root, metabolite_file), 'w') as f:
                f.write(metabolite_data.text)
    return None

def adjust_file_paths(
    all_jobs,
    dir_root=my_settings.magi_task_path):
    """
    creates appropriate full path for input files after being mirrored to disk

    removes jobs that don't have a file they are supposed to
    """
    to_pop = []
    for i, job in enumerate(all_jobs):
        job_path = get_job_dir(job)
        newroot = os.path.join(dir_root, job_path)
        if job['fields']['fasta_file'] != '':
            job['fields']['fasta_file'] = os.path.join(newroot, job['fields']['fasta_file'].split('/')[-1])
            if not os.path.isfile(job['fields']['fasta_file']):
                print('ERROR: fasta file for job %s does not exist' % (job['pk']))
                to_pop.append(i)
                continue
        if job['fields']['metabolite_file'] != '':
            job['fields']['metabolite_file'] = os.path.join(newroot, job['fields']['metabolite_file'].split('/')[-1])
            if not os.path.isfile(job['fields']['metabolite_file']):
                print('ERROR: metabolite file for job %s does not exist' % (job['pk']))
                to_pop.append(i)
                continue
    # remove the jobs that didn't have input files
    if len(to_pop) > 0:
        for i in to_pop:
            all_jobs.pop(i)
    return all_jobs

def jobs_to_script(
        all_jobs,
        dir_root=my_settings.magi_task_path
    ):
    """
    Determine which jobs need a job script.
    
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
        job_path = get_job_dir(job)
        script_path = os.path.join(dir_root, job_path, 'admin')
        if not (os.path.isfile(os.path.join(script_path, 'job_script.sbatch')) or os.path.isfile(os.path.join(script_path, 'job_script.qsub'))):            
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
            'The nucelotide sequence is not a multiple of 3'
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
            file_data = file_data.replace(newline, '\n')
    if '\r' in file_data:
        file_data = file_data.replace('\r', '\n')
    
    # parse gene sequences into one long string and convert to DNA if RNA
    genes = file_data.split('>')[1:]
    seqs = ''.join([''.join(g.split('\n')[1:]).replace('-', '') for g in genes]).upper().replace('U', 'T')
    
    # determine what letters are in the gene sequences
    letters = pd.Series(list(set(seqs)))
    letter_ratios = pd.Series(list(seqs)).value_counts() / float(len(seqs))
    base_sum = letter_ratios[['C', 'G', 'A', 'T']].sum()
    # first check for DNA: if >90% of the letters are CGAT, it's probably DNA
    if base_sum > 0.9:
        answer = 'dna'
    # then amino acids 
    elif letters.str.contains('[ACDEFGHIKLMNPQRSTVWYX*]').all():
        answer = 'protein'
    else:
        no = letters[~letters.str.contains('[ACDEFGHIKLMNPQRSTVWYX*]')].values
        raise RuntimeError('Could not determine if FASTA is nucleotide or protein. Offending character(s): %s; file: %s' % (no, job_data['fields']['fasta_file']))
    
    # translate if desired
    if translate and answer == 'dna':      
        # translate the genes 
        new_data = ''
        gene_list = file_data.split('>')[1:]
        for gene in gene_list:
            header = gene.split('\n')[0]
            seq = gene.split('\n')[1:]
            seq = ''.join([i for i in seq if i != ''])
            try:
                protein = protein_translate(seq)
            except RuntimeError as e:
                change_params(job_data['pk'], 'runflag', 'True')
                msg = 'There was an error translating the DNA sequence with the header:\n'
                msg += '"%s"\n\n' % (header)
                msg += '\nThe error message was:\n%s\n\n' %(e.args)
                msg += 'Please check all your DNA sequences and resubmit your job.\n'
                msg += 'Also note that MAGI does not do any gene calling; do not submit scaffolds or contigs.' 
                msg += '\n\nIf you believe this was in error, please reply to this email.'
                subj = 'MAGI DNA translation error'
                email_user(job_data['fields']['email'], subj, msg)
                email_user(MAGI_EMAIL, subj, msg)
                return None
            new_data += '>' + header + '\n'
            new_data += protein + '\n\n'
        
        # save the new file
        new_filename = os.path.splitext(os.path.basename(file_path))[0]+'_translated.faa'
        new_filepath = os.path.join(os.path.dirname(file_path), new_filename)
        with open(new_filepath, 'w') as f:
            f.write(new_data)
        # change the field in job_data
        job_data['fields']['fasta_file'] = new_filepath
        return job_data
    else:
        return job_data

def job_script(job_data, n_cpd=None):
    """
    uses job data json to create a magi job submission script for cori
    """

    # account_id = 'm2650' # metatlas
    account_id = 'm1541' # openmsi - for realtime
    
    # where to write the job script to
    if job_data['fields']['fasta_file'] != '':
        out_path = os.path.dirname(job_data['fields']['fasta_file'])
    else:
        out_path = os.path.dirname(job_data['fields']['metabolite_file'])

    script_path = os.path.join(out_path, 'admin')
    # prepare score weights
    score_weights = [
        job_data['fields']['score_weight_compound'],
        job_data['fields']['score_weight_reciprocal'],
        job_data['fields']['score_weight_homology'],
        job_data['fields']['score_weight_rxnconnect'],
    ]
    # need to convert this into string so we can join it later
    score_weights = [str(i) for i in score_weights]

    # estimate timing:
    if n_cpd <= 100:
        t_limit = '00:10:00'
        partition = 'debug'
        filetype = 'sbatch'
    else:
        t_limit = '02:00:00'
        partition = 'realtime'
        filetype = 'sbatch'
        # partition = 'genepool'
        # filetype = 'qsub'

    if partition == 'realtime':
        header_lines = [
            '#!/bin/bash -l',
            '#SBATCH --account=%s' % (account_id),
            '#SBATCH --job-name=%s' % (job_data['pk'].split('-')[0]),
            '#SBATCH --time=%s' % (t_limit),
            '#SBATCH --nodes=1',
            '#SBATCH -c 64',
            '#SBATCH --output=%s/log_out.txt' % (out_path),
            '#SBATCH --error=%s/log_err.txt' % (out_path),
            '#SBATCH --partition=%s' % (partition),
            '#SBATCH --constraint=haswell',
            '#SBATCH --license=project',
            '#SBATCH --mail-user=%s' %(MAGI_EMAIL),
            '#SBATCH --mail-type=FAIL,TIME_LIMIT',
            '',
            'source /global/common/software/m2650/python-cori/bin/activate',
            #'module load python/2.7-anaconda-4.4',
            ''
        ]
    elif partition == 'debug':
        header_lines = [
            '#!/bin/bash -l',
            '#SBATCH --account=%s' % (account_id),
            '#SBATCH --job-name=%s' % (job_data['pk'].split('-')[0]),
            '#SBATCH --time=%s' % (t_limit),
            '#SBATCH --output=%s/log_out.txt' % (out_path),
            '#SBATCH --error=%s/log_err.txt' % (out_path),
            '#SBATCH --partition=%s' % (partition),
            '#SBATCH --constraint=haswell',
            '#SBATCH --license=project',
            '#SBATCH --mail-user=%s' %(MAGI_EMAIL),
            '#SBATCH --mail-type=FAIL,TIME_LIMIT',
            '',
            'source /global/common/software/m2650/python-cori/bin/activate',
           # 'module load python/2.7-anaconda-4.4',
            ''
        ]
    elif partition == 'genepool':
        header_lines = [
            '#!/bin/bash',
            '#$ -M %s' % (MAGI_EMAIL),
            '#$ -m a',
            '#$ -l h_rt=%s' % (t_limit),
            '#$ -pe pe_32 32',
            '#$ -l ram.c=7.5G,h_vmem=7.5G',
            '#$ -q exclusive.c',
            '#$ -wd %s' % (out_path),
            '#$ -o %s/log_out.txt' % (out_path),
            '#$ -e %s/log_err.txt' % (out_path),
            '',
            'source /global/common/software/m2650/python-cori/bin/activate',
           # 'module switch python/2.7.4 python/2.7-anaconda_4.3.0',
            ''
            ]

    else:
        raise RuntimeError('%s is an unknown partition' % (partition))
    if job_data['fields']['fasta_file'] != '':
        fasta_file_line = '--fasta %s \\' % (job_data['fields']['fasta_file'])
    else:
        fasta_file_line = '\\'
    if job_data['fields']['metabolite_file'] != '':
        met_file_line = '--compounds %s \\' % (job_data['fields']['metabolite_file'])
    else:
        met_file_line = '\\'
    job_lines = [
        'date -u > %s/start_time.txt' % (os.path.join(out_path, 'admin')),
        '',
        'umask 002',
        '',
        'python /project/projectdirs/metatlas/projects/metatlas_reactions/workflow/helpertools/nersc_memmonitor.py > %s &' % (os.path.join(script_path, 'memory.txt')),
        '',
        'magi_path=/global/homes/p/pasteur/repos/magi',
        'time python $magi_path/workflow/magi_workflow_gene_to_reaction.py \\',
        '%s' % (fasta_file_line),
        '%s' % (met_file_line),
        '--level %s \\' % (job_data['fields']['network_level']),
        '--final_weights %s \\' % (' '.join(score_weights)),
        '--blast_filter %s \\' % (job_data['fields']['blast_cutoff']),
        '--reciprocal_closeness %s \\' % (job_data['fields']['reciprocal_cutoff']),
        '--chemnet_penalty %s \\' % (job_data['fields']['chemnet_penalty']),
        '--output %s --mute' % (out_path),
        '',
        'if [ $? -eq 0 ] && [ ! -f %s/incomplete ]; then' % (os.path.join(out_path, 'admin')),
        '  python $magi_path/workflow/magi_workflow_compound_to_reaction.py --not_first_script --output %s' % (out_path), 
        'else touch %s/incomplete; fi' % (os.path.join(out_path, 'admin')),
        'if [ $? -eq 0 ] && [ ! -f %s/incomplete ]; then' % (os.path.join(out_path, 'admin')),
        '  python $magi_path/workflow/magi_workflow_reaction_to_gene.py --notf_first_script --output %s' % (out_path), 
        'else touch %s/incomplete; fi' % (os.path.join(out_path, 'admin')),
        '  if [ $? -eq 0 ] && [ ! -f %s/incomplete ]; then' % (os.path.join(out_path, 'admin')),
        'python $magi_path/workflow/magi_workflow_scoring.py --not_first_script --output %s' % (out_path), 
        '  else touch %s/incomplete; fi' % (os.path.join(out_path, 'admin')),
        'if [ $? -eq 0 ] && [ ! -f %s/incomplete ]; then' % (os.path.join(out_path, 'admin')),
        '  date -u > %s/end_time.txt' % (os.path.join(out_path, 'admin')),
        'else touch %s/incomplete; fi' % (os.path.join(out_path, 'admin'))
    ]
    
    job = '\n'.join(header_lines) + '\n' + '\n'.join(job_lines) + '\n'

    # change umask temporarily; don't want job script to be world-read
    old_mask = os.umask(007)
    # write job
    with open(os.path.join(script_path, 'job_script.%s') % (filetype), 'w') as f:
        f.write(job)
    os.umask(old_mask)
    return None

def ppm_error(mass, theoretical_mass):
    """
    Returns the ppm error of a given observed mass, and theoretical mass
    """
    try:
        mass = float(mass)
        theoretical_mass = float(theoretical_mass)
    except Exception as e:
        print('one of the arguments is not a number!')
        raise e
    ppm = (mass - theoretical_mass) / theoretical_mass * 1e6
    return abs(ppm)

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

    theoretical = potential_compounds['mono_isotopic_molecular_weight']
    ppm_error = (theoretical - mass) / theoretical * 1e6
    ppm_error = abs(ppm_error)
    potential_compounds['ppm_error'] = ppm_error

    cpds = potential_compounds[[extract, 'ppm_error']].values.tolist()
    if len(cpds) == 0:
        cpds = None
    return cpds

def mz_neutral_transform(val, adduct, transform='neutralize'):
    """
    val: m/z or neutral mass
    adduct: adduct to consider
    transform: 'neutralize' or 'ionize'
      if neutralize, neutralizes 'val' by subtracting adduct
      if ionize, ionizes 'val' by adding adduct
    list of acceptible adducts:
      M+, M+H, M+NH4, M+Na, M+CH3OH+H, M+K, M+ACN+H, M+2Na-H,
      M+IsoProp+H, M+ACN+Na, M+2K-H, M+DMSO+H, M+2ACN+H,
      M+IsoProp+Na+H, 2M+H, 2M+NH4, 2M+Na, 2M+K, 2M+ACN+H,
      2M+ACN+Na, M+3H, M+2H+Na, M+H+2Na, M+3Na, M+2H, M+H+NH4,
      M+H+Na, M+H+K, M+ACN+2H, M+2Na, M+2ACN+2H, M+3ACN+2H, M-H,
      M+Cl, M+FA-H, M+Hac-H, 2M-H, 2M+FA-H, 2M+Hac-H, 3M-H, M-3H,
      M-2H, M-H2O-H, M+Na-2H, M+K-2H, M+Br, M+TFA-H
    """
    acceptible_adducts = [
    'M+','M+H', 'M+NH4', 'M+Na', 'M+CH3OH+H', 'M+K', 'M+ACN+H', 'M+2Na-H',
    'M+IsoProp+H', 'M+ACN+Na', 'M+2K-H', 'M+DMSO+H', 'M+2ACN+H',
    'M+IsoProp+Na+H', '2M+H', '2M+NH4', '2M+Na', '2M+K', '2M+ACN+H',
    '2M+ACN+Na', 'M+3H', 'M+2H+Na', 'M+H+2Na', 'M+3Na', 'M+2H', 'M+H+NH4',
    'M+H+Na', 'M+H+K', 'M+ACN+2H', 'M+2Na', 'M+2ACN+2H', 'M+3ACN+2H', 'M-H',
    'M+Cl', 'M+FA-H', 'M+Hac-H', '2M-H', '2M+FA-H', '2M+Hac-H', '3M-H',
    'M-3H', 'M-2H', 'M-H2O-H', 'M+Na-2H', 'M+K-2H', 'M+Br', 'M+TFA-H'
    ]

    # M + N
    simple = {
      'M+': 0.0000,
      'M+H': 1.007276,
      'M+NH4': 18.033823,
      'M+Na': 22.989218,
      'M+CH3OH+H': 33.033489,
      'M+K': 38.963158,
      'M+ACN+H': 42.033823,
      'M+2Na-H': 44.971160,
      'M+IsoProp+H': 61.06534,
      'M+ACN+Na': 64.015765,
      'M+2K-H': 76.919040,
      'M+DMSO+H': 79.02122,
      'M+2ACN+H': 83.060370,
      'M+IsoProp+Na+H': 84.05511,
      'M-H': -1.007276,
      'M+Cl': 34.969402,
      'M+FA-H': 44.998201,
      'M+Hac-H': 59.013851,
      'M-H2O-H': -19.01839,
      'M+Na-2H': 20.974666,
      'M+K-2H': 36.948606,
      'M+Br': 78.918885,
      'M+TFA-H': 112.985586,
    }
    # 2M + N
    two_M = {
      '2M+H': 1.007276,
      '2M+NH4': 18.033823,
      '2M+Na': 22.989218,
      '2M+K': 38.963158,
      '2M+ACN+H': 42.033823,
      '2M+ACN+Na': 64.015765,
      '2M-H': -1.007276,
      '2M+FA-H': 44.998201,
      '2M+Hac-H': 59.013851,
    }
    # 3M + N
    three_M = {
      '3M-H': -1.007276,
    }
    # M/2 + N
    two_charge = {
      'M+2H': 1.007276,
      'M+H+NH4': 9.520550,
      'M+H+Na': 11.998247,
      'M+H+K': 19.985217,
      'M+ACN+2H': 21.520550,
      'M+2Na': 22.989218,
      'M+2ACN+2H': 42.033823,
      'M+3ACN+2H': 62.547097,
      'M-2H': -1.007276,
    }
    # M/3 + N
    three_charge = {
      'M+3H': 1.007276,
      'M+2H+Na': 8.334590,
      'M+H+2Na': 15.7661904,
      'M+3Na': 22.989218,
      'M-3H': -1.007276,
    }
    transform = transform.lower()
    if transform not in ['neutralize', 'ionize']:
      raise RuntimeError('%s is not an acceptible transformaion;\
           please use "ionize" or "neutralize"' % (transform))
    if adduct not in acceptible_adducts:
      raise RuntimeError('%s not in the list of acceptible adducts'
           % (adduct))
    x = None
    if adduct in simple.keys():
      if transform == 'neutralize':
           x = val - simple[adduct]
      elif transform == 'ionize':
           x = val + simple[adduct]

    if adduct in two_M.keys():
      if transform == 'neutralize':
           x = (val - two_M[adduct]) / 2
      elif transform == 'ionize':
           x = 2 * val + two_M[adduct]

    if adduct in three_M.keys():
      if transform == 'neutralize':
           x = (val - three_M[adduct]) / 3
      elif transform == 'ionize':
           x = 3 * val + three_M[adduct]

    if adduct in two_charge.keys():
      if transform == 'neutralize':
           x = (val - two_charge[adduct]) * 2
      elif transform == 'ionize':
           x = (val / 2) + two_charge[adduct]

    if adduct in three_charge.keys():
      if transform == 'neutralize':
           x = (val - three_charge[adduct]) * 3
      elif transform == 'ionize':
           x = (val / 3) + three_charge[adduct]
    return x

def accurate_mass_search_wrapper(job_data, reference_compounds, max_compounds=25000):
    """
    performs accurate mass search using unique_compounds table
    """
    # allow administrative override of compound limit
    note_path = os.path.join(os.path.dirname(job_data['fields']['metabolite_file']), 'admin')
    if 'cpd_override' in os.listdir(note_path):
        max_compounds = 1e6
    search_ppm = job_data['fields']['ppm']

    # get adducts according to polarity
    if job_data['fields']['polarity'] == 'pos':
        adducts = job_data['fields']['adducts_pos'].split(',')
    elif job_data['fields']['polarity'] == 'neg':
        adducts = job_data['fields']['adducts_neg'].split(',')
    elif job_data['fields']['polarity'] == 'neut':
        adducts = ['']
    else:
        raise RuntimeError('Could not understand polarity %s in %s' % (job_data['fields']['polarity'], job_data['pk']))

    # load compound table (should be masses in original_compounds)
    compounds = pd.read_csv(job_data['fields']['metabolite_file'])

    # rename original_compounds column
    columns = compounds.columns.values
    if 'original_compound' not in columns:
        raise RuntimeError('no original_compound')
    columns[columns == 'original_compound'] = 'original_mz'
    compounds.columns = columns

    # set up data container
    data = {
        'original_compound': [],
        'searched_adduct': [],
        'original_mz': [],
        'ppm_error': [],
        'compound_score': []
    }
    # accurate mass search and store results
    for mz in compounds['original_mz'].unique():
        for adduct in adducts:
            if adduct != '':
                neutral_mass = mz_neutral_transform(mz, adduct)
            else:
                neutral_mass = mz
            found_compounds = accurate_mass_match(neutral_mass,
                                                  compound_df=reference_compounds,
                                                  ppm=search_ppm
                                                 )
            if found_compounds is not None:
                for cpd in found_compounds:
                    data['original_compound'].append(cpd[0])
                    data['ppm_error'].append(cpd[1])
                    data['compound_score'].append(search_ppm + 1 - cpd[1])
                    data['searched_adduct'].append(adduct)
                    data['original_mz'].append(mz)

    # merge with user input and save
    df = pd.DataFrame(data)
    compounds = compounds.merge(df, on='original_mz', how='left')
    
    # save the new table
    new_path = os.path.splitext(job_data['fields']['metabolite_file'])[0] + '_mass_searched.csv'
    job_data['fields']['metabolite_file'] = new_path
    compounds.to_csv(job_data['fields']['metabolite_file'])

    if compounds['original_compound'].drop_duplicates().shape[0] > max_compounds:
        raise RuntimeError('too many compounds')

    return job_data

def email_user(email, subject, text, subtype='plain'):
    """
    emails a MAGI user a specific message
    """
    msg = MIMEText(text, subtype)
    msg['Subject'] = subject
    msg['From'] = 'magi_web@lbl.gov'
    msg['To'] = email

    s = smtplib.SMTP('localhost')
    s.sendmail('magi_web@lbl.gov', email, msg.as_string())

def save_job_params(job_data, fname='too_many_compounds'):
    """
    saves a job's json and input file info for later comparison
    """
    compound_df = pd.read_csv(job_data['fields']['metabolite_file'])
    job_data['n_mz'] = compound_df.shape[0]

    job_path = os.path.join(os.path.dirname(job_data['fields']['metabolite_file']), 'admin')
    with open(os.path.join(job_path, fname), 'w') as f:
        f.write(json.dumps(job_data))

def accurate_mass_checkpoint(job_data, fname='too_many_compounds'):
    """
    performs various checks to determine if accurate mass searching
    should proceed
    1. is the file 'too_many_compounds' present in the job dir?
    2. are the current params different from old params?
    3. is the metabolite table length the same(same number of cpds?)

    returns True if accurate mass searching should proceed
    returns False if the issues were not fixed

    """
    # check if the email file exists
    job_path = os.path.join(os.path.dirname(job_data['fields']['metabolite_file']), 'admin')
    if not os.path.isfile(os.path.join(job_path, fname)):
        return True
    
    # are the relevant params different?
    with open(os.path.join(job_path, fname), 'r') as f:
        data = f.read()
    old_job_data = json.loads(data)

    relevant_keys = ['polarity', 'adducts_pos', 'adducts_neg', 'ppm']
    for k in relevant_keys:
        if job_data['fields'][k] != old_job_data['fields'][k]:
            return True

    # is the number of mz the same?
    compound_df = pd.read_csv(job_data['fields']['metabolite_file'])
    if old_job_data['n_mz'] != compound_df.shape[0]:
        return True

    # if all the checkpoints failed, don't run the accurate mass search
    return False
