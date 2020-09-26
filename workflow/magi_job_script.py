"""
This script contains functions to make a MAGI job script that can be run from the command line.
"""
import os
import sys
def create_job_script(
                    location_to_store_script, 
                    output_path,
                    cpu_count = 1,
                    fasta_file = None,
                    compounds_file = None, 
                    magi_location = os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    level = 2,
                    final_weights = [1.0, 1.0, 1.0, 1.0],
                    blast_filter = 85,
                    reciprocal_closeness = 75,
                    chemnet_penalty = 4
                    ):
    """
    Create a job file for MAGI and store it in location_to_store_script.
    """
    if fasta_file is not None and os.path.exists(fasta_file):
        fasta_file = os.path.abspath(fasta_file)
    if compounds_file is not None and os.path.exists(compounds_file):
        compounds_file = os.path.abspath(compounds_file)
    # Write fasta or compounds file or both to script
    if fasta_file == None and compounds_file == None:
        raise RuntimeError("Specify compounds file or genes file.")
    elif fasta_file == None:
        fasta_compounds_settings = "--compounds {}".format(compounds_file)
    elif compounds_file == None:
        fasta_compounds_settings = "--fasta {}".format(fasta_file)
    else:
        fasta_compounds_settings = "--fasta {} \\\n--compounds {} \\".format(fasta_file, compounds_file)
    
    # This is the script
    job_lines = [
        '#!/bin/bash',
        '',
        'magi_path={}'.format(magi_location),
        'source activate magi || conda activate magi', 
        '',
        'mkdir {}'.format(output_path),
        'date -u > {}/start_time.txt'.format(output_path),
        'python $magi_path/workflow/magi_workflow_gene_to_reaction.py \\',
        fasta_compounds_settings,
        '--cpu_count {} \\'.format(cpu_count),
        '--level {} \\'.format(level),
        '--final_weights {} \\'.format(' '.join([str(weight) for weight in final_weights])),
        '--blast_filter {} \\'.format(blast_filter),
        '--reciprocal_closeness {} \\'.format(reciprocal_closeness),
        '--chemnet_penalty {} \\'.format(chemnet_penalty),
        '--output {} --mute'.format(output_path),
        '',
        'if [ $? -eq 0 ] && [ ! -f {}/incomplete.txt ]; then'.format(output_path),
        '  python $magi_path/workflow/magi_workflow_compound_to_reaction.py --not_first_script --output %s' % (output_path), 
        'else touch {}/incomplete.txt; fi'.format(output_path),
        'if [ $? -eq 0 ] && [ ! -f {}/incomplete.txt ]; then'.format(output_path),
        '  python $magi_path/workflow/magi_workflow_reaction_to_gene.py --not_first_script --output %s' % (output_path), 
        'else touch {}/incomplete.txt; fi'.format(output_path),
        'if [ $? -eq 0 ] && [ ! -f {}/incomplete.txt ]; then'.format(output_path),
        'python $magi_path/workflow/magi_workflow_scoring.py --not_first_script --output %s' % (output_path), 
        'else touch {}/incomplete.txt; fi'.format(output_path),
        'if [ $? -eq 0 ] && [ ! -f {}/incomplete.txt ]; then'.format(output_path),
        '  date -u > {}/end_time.txt'.format(output_path),
        'else touch {}/incomplete.txt; fi'.format(output_path)]
    
    ## Join job lines
    job = '\n'.join(job_lines) + '\n'

    # write job
    with open(os.path.join(location_to_store_script, 'job_script.sh'), 'w') as jobfile:
        jobfile.write(job)
    return None
