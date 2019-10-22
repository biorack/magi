import os
import pandas as pd
import numpy as np
import sys
import subprocess
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from local_settings import local_settings as settings_loc

def load_dataframe(fname, filetype=None, key=None):
    """
    Uses the appropriate pandas function to load a file based on the
    file extension:
    .pkl or .pickle: pickle file
    .xls or .xlsx: excel file
    .csv: comma separated file
    .txt, .tab, or .tsv: tab separated file
    .h5 or .hdf5: HDF5 formatted files. To load these, you must also
                  pass a key argument!

    Inputs
    ------
    fname: path to the file to be loaded as a pandas dataframe
    filetype: Used to override the autodetection of a filetype based on
              its file extension. Accepts file extensions listed above,
              but without the preceding "." (e.g. pass it "pkl")
    key: They key for the table to be loaded in the HDF5 file. Only
         required/used when loading HDF5 tables. 

    Outputs
    -------
    df: the loaded pandas dataframe

    WARNING: for binary file types (pickle and hdf5), be wary of
    saving and loading from different versions of pandas, as this will
    very likely break the loader.
    """

    if filetype is None:
        file_ext = os.path.splitext(fname)[1][1:]
    else:
        file_ext = filetype
    if file_ext in ['pkl', 'pickle']:
        df = pd.read_pickle(fname)
    elif file_ext in ['xls', 'xlsx']:
        df = pd.read_excel(fname)
    elif file_ext in ['csv']:
        df = pd.read_csv(fname)
    elif file_ext in ['txt', 'tab', 'tsv']:
        df = pd.read_csv(fname, sep='\t')
    elif file_ext in ['h5', 'hdf5']:
        if key is None:
            raise IOError('"key" argument must be used when loading\
                an HDF5 file')
        else:
            df = pd.read_hdf(fname, key)
    else:
        raise IOError('could not infer what type of file %s is... please \
            make it a csv, tab, or pickle file')

    # remove rows that are empty
    df = df[~pd.isnull(df).all(axis=1)]
    return df

def get_settings():
    my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)
    return my_settings

def reformat_pactolus(df, original_compound=None, compound_score=None):
    """
    Reformats pactolus output table to be direcly portable into MAGI.
    1.  changes column "inchi_key_y" or "inchi_key" to
        "original_compound"
    2.  changes column "score" to "compound_score"

    Inputs
    ------
    df: Pactolus output table as Pandas DataFrame
    original_compound: column name corresponding to the original
                       compound inchi key. If None, will look for
                       matches to both "inchi_key_y" or "inchi_key".
    compound_score: column name corresponding to the pactolus score. If
                    None, will look for match to "score"

    Outputs:
    df: reformatted Pactolus table
    """

    cols = df.columns.values
    if original_compound is None:
        old_cols = pd.Series(['inchi_key_y', 'inchi_key'])
        tmp = old_cols[old_cols.isin(cols)].values
        if len(tmp) == 0:
            raise RuntimeError('no columns named "inchi_key_y" or "inchi_key" \
                in Pactolus table to rename')
        elif len(tmp) > 1:
            raise RuntimeError('both "inchi_key_y" and "inchi_key" found in \
                Pactolus table columns')
        else:
            original_compound = tmp[0]

    i = np.argwhere(cols == original_compound)
    if len(i) == 0:
        raise RuntimeError('There is no column named "%s", please double \
            check your input.' % (original_compound))
    else:
        i = i[0][0]
    cols[i] = 'original_compound'
    
    if compound_score is None:
        compound_score = 'score'
    
    i = np.argwhere(cols == compound_score)
    if len(i) == 0:
        raise RuntimeError('There is no column named "%s", please double \
            check your input.' % (original_compound))
    else:
        i = i[0][0]
    cols[i] = 'compound_score'
    df.columns = cols
    return df

def ec_parse(x):
    """
    Cleans up IMG gene table's Enzyme column by pulling out and
    separating EC numbers. Can be used as a pandas.apply() function.

    Inputs
    ------
    x: a string

    Outputs
    -------
    out: a string of EC numbers separated by a "|" character, or only the
         "|" character if the input is null or empty
    """

    if pd.isnull(x):
        return '|'

    out = '|'
    if x != '':
        if '<<>>' in x:
            one = x.split('<<>>')
            for two in one:
                try:
                    ec = two.split('EC:')[1].split('=')[0]
                except:
                    ec = ''
                out += ec+'|'
        else:
            try:
                ec = x.split('EC:')[1].split('=')[0]
            except:
                ec = ''
            out += ec+'|'
        return out
    else:
        return out
    
def load_genome(fasta, intfile_path, annotation_file=None):
    """
    This function will take some standard fasta a input and convert it
    to an appropriate genome dataframe.
    This function will also construct a blast database for the genome.

    Inputs
    ------
    fasta: path to standard fasta file of all the genes.
           FASTA sequence header format should be:
           >UNIQUE_GENE_IDENTIFIER OTHER_INFORMATION
           Note the space between the unique identifier and other info

    intfile_path: path to the temporary storage place of the database. If None, the
               BLAST database is not made, and the gene table is not
               stored.

    annotation_file: path to an optional .TAB file that contains
                     annotations for the genes in the fasta file. One
                     column must be named "Gene_ID" and correspond to
                     the UNIQUE_GENE_IDENTIFIER in the fasta headers.

    Outputs
    -------
    genome: gene sequence table that is merged with the annotation table
            if one is provided. This table is saved as a pickle file to
            intfile_path/gene_fastas/filename.pkl
    db_path: path to the genome's BLAST database
    """
    # TODO: find a way to handle windows text files (\n\r for new lines)
    # TODO: find a way to handle additional '>' characters in the header

    # process the fasta file
    with open(fasta, 'r') as f:
        genes = f.read()
    if genes.strip()[0] != '>':
        raise RuntimeError('%s does not appear to be a FASTA file' % (fasta))

    data = []
    for entry in genes.split('>')[1:]:
        header = '>' + entry.splitlines()[0]
        gene_id = header.split(' ')[0][1:]
        sequence = ''.join(entry.splitlines()[1:])
        data.append([gene_id, header, sequence])
    genome = pd.DataFrame(data, columns=['Gene_ID', 'header', 'sequence'])
    if genome['Gene_ID'].duplicated().any():
        first_dup = genome[genome['Gene_ID'].duplicated()].head(1)
        first_header = first_dup.iloc[0, 1]
        first_identifier = first_dup.iloc[0, 0]
        raise RuntimeError('There are duplicated Gene_ID fields! please check\
            your unique gene identifiers in the fasta headers. For the first\
            FASTA sequence of %s, I parsed the unique gene identifier to be \
            %s' % (repr(first_header), repr(first_identifier)))

    if annotation_file is not None:
        # Make gene info table
        annotation_table = load_dataframe(annotation_file)
        # turn gene IDs into strings to keep things consistent
        annotation_table['Gene_ID'] = annotation_table['Gene_ID'].apply(str)

        # clean up the nans
        annotation_table.fillna(value='', inplace=True)

        # Really only if the annotation file was an IMG-derived table
        try:
            annotation_table['EC'] = annotation_table['Enzyme'].apply(ec_parse)

            # move the new EC column to front for ease of visualization
            cols = annotation_table.columns.tolist()
            newcols = cols[-1:] + cols[:-1]
            annotation_table = annotation_table[newcols]
        except KeyError:
            print( 'Could not find a column corresponding to EC annotations, \
                    skipping EC parsing')

        genome = pd.merge(genome, annotation_table, on='Gene_ID', how='left')

    print( '!@# FASTA file loaded with {} genes'.format(len(genome)))

    genome.set_index('Gene_ID', inplace=True, drop=True)
    genome_name = os.path.splitext(os.path.basename(fasta))[0]
    if intfile_path is not None:
        gene_fasta_path = os.path.join(intfile_path, 'gene_fastas')
        if not os.path.isdir(gene_fasta_path):
            os.makedirs(gene_fasta_path)
        gene_seq_path = os.path.join(gene_fasta_path,
            '%s_sequences.pkl' % (genome_name))
        # gene_seq_path = '%s/gene_fastas/%s_sequences.pkl' \
        #                 % (intfile_path, genome_name)
        genome.to_pickle(gene_seq_path)
        print( '!!! saved gene_sequence table here: {}'.format(gene_seq_path))

        # Make the blast database of the fasta
        # this command makes the blast database for a given fasta file.
        blastbin = get_settings(); blastbin = blastbin.blastbin
        makeblastdb_path = os.path.join(blastbin, 'makeblastdb')
        fasta_path = fasta
        db_path = os.path.join(intfile_path, 'BLAST_dbs',(os.path.splitext(os.path.basename(fasta_path))[0]+'.db'))
        print( '!!! blast database stored here: {}'.format(db_path))
        if os.name == 'nt': #Check if the operating system is windows or linux/mac
            make_db_command = '%s -in %s -out %s -dbtype prot' \
                % (makeblastdb_path+'.exe', fasta_path, db_path)
        else:
            make_db_command = '%s -in %s -out %s -dbtype prot' \
                % (makeblastdb_path, fasta_path, db_path)
        blastp = subprocess.Popen(
            make_db_command,
            shell=True, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if blastp.stdout is not None:
            print( blastp.stdout.read())
        if blastp.stderr is not None:
            print( blastp.stderr.read())
    else:
        db_path = None
    return genome, db_path



def reciprocal_agreement(df, forward_name='reaction_id_r2g',
                         reverse_name='reaction_id_g2r',
                         closeness_threshold=0.75):
    """
    Finds and scores reciprocal agreement converging on a reaction
    between forward and reverse homology searches.
    
    Inputs
    ------
    df: dataframe containing forward and reverse homology searches
    forward_name: column name in df corresponding to the reaction ID for
                  the forward search
    reverse_name: column name in df corresponding to the reaction ID for
                  the reverse search
    closeness_threshold: Determines cutoff to call a reciprocal search
                         as "close" if the reaction IDs between forward
                         and reverse searches do not agree.
                         For example, if there is no agreement but the
                         forward score is 100 and reverse score is 90,
                         this would be scored as "close" if the
                         closeness_threshold was less than or equal to
                         0.9

    Outputs
    -------
    df: input df but with a new column named "reciprocal_score"
        corresponding to scored reciprocal agreement:
        2.0: agreement between forward and reverse searches
        1.0: forward and reverse searches disagree on reaction, but the
             homology scores are close as judged by closeness_threshold
        0.1: either a forward or reverse search does not exist (homology
             score below threshold)
        0.01: forward and reverse searches disagree on reaction and are
              not close
    """
    # find agreement
    agree_idx = df[df['reaction_id_r2g'] == df['reaction_id_g2r']].index
    df.loc[agree_idx, 'reciprocal_score'] = 2.
    # disagreement
    disagree = df[df['reaction_id_r2g'] != df['reaction_id_g2r']].index
    slc = df.loc[disagree]
    # close disagreements get a medium score
    close = (
        slc[['e_score_r2g', 'e_score_g2r']].min(axis=1)
        >= (slc[['e_score_r2g', 'e_score_g2r']].max(axis=1) * closeness_threshold))
    close_idx = slc.loc[close].index
    df.loc[close_idx, 'reciprocal_score'] = 1
    # very different disagreements get a low score
    wrong_idx = df[pd.isnull(df['reciprocal_score'])].index
    df.loc[wrong_idx, 'reciprocal_score'] = 0.01
    # if one direction did not get a blast score, 
    # change reciprocal score to 0.1 - not wrong, but not close either.
    incomparable_idx = df.loc[pd.isnull(df[['e_score_r2g',
        'e_score_g2r']]).any(axis=1)].index
    df.loc[incomparable_idx, 'reciprocal_score'] = 0.1

    return df


def homology_score(df, forward_name='e_score_r2g', reverse_name='e_score_g2r'):
    """
    Calculates homology score for a compound-gene association:
    Takes the sum of forward and reverse scores, subtracts the difference
    between them: (forward + reverse) - abs(forward-reverse)

    deprecated:
    average(forward_score, reverse_score) - abs(forward_score - reverse_score)

    Inputs
    ------
    df: dataframe containing forward and reverse homology scores
    forward_name: column name in df corresponding to the forward score
    reverse_name: column name in df corresponding to the reverse score

    Output
    ------
    array of combined homology scores
    """

    forward = df[forward_name].values.astype(float)
    reverse = df[reverse_name].values.astype(float)

    score = (forward + reverse) - abs(forward-reverse)

    return score

def magi_score(s, w=np.asarray([np.nan, np.nan])):
    """
    Calculates the geometric mean of an array of numbers.
    Used to calculate one score even if sub-scores are scaled differently.

    Inputs
    ------
    s: 1D array-like list of scores
    weights: 1D array-like list of weights

    Outputs
    -------
    The (weighted) geometric mean
    """
    
    if not isinstance(s, np.ndarray):
        s = np.asarray(s)
    s = s.astype(float)

    if not isinstance(w, np.ndarray):
        w = np.asarray(w)
    # if no weights provided, make them all ones
    if np.isnan(w).all():
        w = np.ones(s.shape)
    w = w.astype(float)

    # infer summing axis
    a = len(s.shape) - 1
    if a > 1:
        raise RuntimeError(
            'Scores array has too many dimensions (%s)' % (s.shape)
            )
    if w.shape != s.shape:
        raise RuntimeError(
            'Weights array does not have same dimensions as Scores array: %s vs %s' % (w.shape, s.shape)
            )

    return np.exp(np.sum(np.multiply(w, np.log(s)), axis=a) / np.sum(w, axis=a))
