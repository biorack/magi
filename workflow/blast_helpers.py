import os
import sys
import pandas as pd
import numpy as np
import subprocess

import workflow_helpers_new as mg

## Load data needed for BLASTing
my_settings = mg.get_settings()
mrs_reaction = mg.load_mrs_reaction()

refseq_dbpath = my_settings.refseq_db
# set up joiner tables
a = mrs_reaction.index.values
b = mrs_reaction['refseq_id'].values
to_blowup = zip(a,b)
data = []
for pair in to_blowup:
    rxn_id = pair[0]
    refseqs = pair[1]
    for r in refseqs.split('|'):
        if r != '':
            data.append([rxn_id, r])
rxn_refseq_join = pd.DataFrame(data, columns=['reaction_id', 'refseq_id'])

def partition_indexes(totalsize, numberofpartitions, offset=0):
    """
    Used to split up an iterable into multiple chunks
    This function returns the indices to create chunks of roughly the
    same size

    From:
    http://stackoverflow.com/questions/33878939/
        get-indices-of-roughly-equal-sized-chunks

    Inputs
    ------
    totalsize: length of the iterable to partition
    numberofpartitions: number of partitions to split the iterable into
    offset: number to add to the partition indices. Used when you're
            partitioning a slice of something, for example if you want
            to parititon array[100:200] into 4 partitions, this function
            would return [[0, 25], [25, 50], [50, 75], [75, 100]] if
            offset were kept default. So you would set offset to 100,
            which would return
            [[100, 125], [125, 150], [150, 175], [175, 200]]
    """

    chunk_indices = []
    # Compute the chunk size (using floor division (rounding down))
    chunksize = totalsize // numberofpartitions
    # How many chunks need an extra 1 added to the size?
    remainder = totalsize - chunksize * numberofpartitions
    a = 0
    for i in range(numberofpartitions):
        b = a + chunksize + (i < remainder)
        # Yield the inclusive-inclusive range
#         chunk_indices.append([a, b - 1])
        chunk_indices.append([a, b])
        a = b

    offset_chunk_indices = []
    for indices in chunk_indices:
        i = indices[0] + offset
        j = indices[1] + offset
        offset_chunk_indices.append([i, j])

    return offset_chunk_indices

def refseq_to_reactions(blast_results, refseq_col):
    """
    enumerates the reaction(s) that a refseq is associated with

    Inputs
    ------
    blast_results: output of multi_blast()
    refseq_col: column name that has the reaction reference sequences
                either "query acc." or "subject acc."
                "query acc." for reaction_to_gene
                "subject acc." for gene_to_reaction

    Outputs
    -------
    result_table: pandas dataframe that connects BLAST results to a
                  reaction ID in the reaction database
    """

    if refseq_col not in blast_results.columns:
        raise RuntimeError('%s is not a column in the dataframe!'
                           % (refseq_col))

    result_table = pd.merge(
        blast_results, rxn_refseq_join, 
        left_on=refseq_col, right_on='refseq_id', how='left')
    result_table.drop('refseq_id', axis=1, inplace=True)

    # keep only the best hit to each reaction, for each query
    result_table.sort_values(
        ['query acc.', 'e_score'],
        ascending=[True, False]).drop_duplicates(['query acc.', 'reaction_id'])

    return result_table

def keep_top_blast(escore_series, filt=0.85):
    """
    Function to keep only the top scores of a blast result.
    Meant to be an "apply" function for a pandas groupby object.
    Can also be used as a regular "apply" function.

    Ex. Usage:
    gene_groups = df.groupby('query acc.')
    multidx = gene_groups['e_score'].apply(keep_top_blast).index
    idx = multidx.levels[1]
    df.loc[idx] # this is now the filtered df

    Inputs
    ------
    escore_series: Pandas Series
    filt: lower cutoff of max value to keep

    Outputs
    -------
    final: Pandas Series
    """
    max_score = escore_series.max()
    final = escore_series[escore_series >= (max_score * filt)]
    return final

def tabulate_blast(results):
    """
    Convert a blast results text output into a dataframe

    Inputs
    ------
    results: BLAST results output when run with the following parameter:
             -outfmt \'10 qacc sacc qcovs length ppos evalue bitscore\'

    Outputs
    -------
    df: Pandas dataframe summarizing the BLAST results
    """

    cols = ['query acc.', 'subject acc.', '% query coverage per subject',
            'alignment length', '% positives', 'evalue', 'bit score']
    result_list = results.split('\n')[:-1]
    data = []
    for r in result_list:
        elements = r.split(',')
        data.append(elements)
    df = pd.DataFrame(data, columns=cols)
    return df

def multi_blast(query_list, query_full_table, database_path, result_path,
                cpu=2, raise_blast_error=True):
    """
    Large-scale BLASTP with multiple processes.
    Partitions the query list into N lists, then uses a shell script to
    open several BLASTP processes, saves them to disk and then collects
    them.

    Inputs
    ------
    query_list: list of indices to query_full_table that you want blasted
    query_full_table:   pandas DataFrame representing protein sequences
                        Must have a column named "sequence"
    database: path to blast database to be used
    result_path: where you want results saved
    cpu: number of processes to open
    scriptpath: path to shell script that opens multiple blast processes
    raise_blast_error: boolean, False will still print( BLAST errors/warnings,
                        but will not raise an error

    Outputs
    -------
    blast_results: pandas DataFrame that stores all blast results.
    """
    my_settings = mg.get_settings()
    blastbin = my_settings.blastbin
    db_name = os.path.basename(database_path)

    # don't open more processes than necessary
    if cpu > len(query_list):
        cpu = len(query_list)

    # set up the blast search strings
    idxes = partition_indexes(len(query_list), cpu)
    mplist = []
    for pair in idxes:
        i, j = pair
        seq = ''
        for gid in query_list[i:j]:
            if gid in query_full_table.index:
                seq += '>' + gid + '\n'
                seq += query_full_table.loc[gid, 'sequence']
                seq += '\n'
            else:
                print( '%s not in sequence database!' % (gid))
        mplist.append(seq)

    # call a shell script that opens n blasts,
    # instead of opening pool of workers
    # first need to save the input seq as a temporary file
    cwd = os.path.join(result_path, 'multi_blast_files')
    if not os.path.isdir(cwd):
        os.makedirs(cwd)
    for i, seq in enumerate(mplist):
        with open('%s/tmp_seq_%s.faa' % (cwd, i), 'w') as f:
            f.write(seq)

    # then call the job. TODO: fix prints for windows. They are not accurate
    scriptpath = os.path.join(blastbin, 'recip_blaster.sh')
    print( '!!! blast script: {}'.format(scriptpath))
    print( '!!! # processes to open: {}'.format(cpu))
    print( '!!! database path: {}'.format(database_path))
    print( '!!! results stored: {}'.format(cwd))
    sys.stdout.flush()
    blaster_file = "{}__{}.txt".format(os.path.join(result_path, 'blasterr'), db_name)
    if os.name == 'nt': #Check if the operating system is windows or linux/mac
        # This is for windows
        subprocess.call('{0} -query  {1} -db {2} -outfmt "10 qacc sacc qcovs length ppos evalue bitscore" -evalue 1 -max_target_seqs 10 > {3} &'.format(
        os.path.join(blastbin, "blastp.exe"),
        os.path.join(cwd, "tmp_seq_0.faa"),
        database_path,
        os.path.join(cwd, "tmp_out_blasted_0.txt")), shell=True)
    else:
        subprocess.call(
        '%s %s %s %s %s 2> %s'
        % (scriptpath, cpu, database_path, cwd, my_settings.repo_location,blaster_file),
        shell=True)
        with open(blaster_file, 'r') as f:
            msg = f.read()
        if msg != '':
            print( '!@# WARNING: BLAST script error message:')
            print( msg)
            print( '-'*80)
            if raise_blast_error:
                raise RuntimeError('blast had an error message')
    sys.stdout.flush()

    # collect the results
    results = ''
    for i, seq in enumerate(mplist):
        with open(os.path.join(cwd, 'tmp_out_blasted_%s.txt' % i), 'r') as f:
            results += f.read()
    results_table = tabulate_blast(results)

    # convert e-value to blast score
    results_table['e_score'] = results_table['evalue'].apply(float).apply(
        np.log10) * -1
    # convert perfect blast scores to 200
    idx = results_table[results_table['e_score'] == np.inf].index
    results_table.loc[idx, 'e_score'] = 200.

    # clean up temp folders
    # don't do this anymore. makes debugging really hard.
    # maybe do a cleanup step every once in a while via cronjob
    # rmtree(cwd)

    return results_table