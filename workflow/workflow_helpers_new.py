import os
import pandas as pd
import numpy as np
import sys
import subprocess
import datetime
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from local_settings import local_settings as settings_loc

def make_output_dirs(output_dir=None, fasta_file=None, compounds_file=None, intermediate_files='intermediate_files'):
    """set up where the results will be stored"""
    if output_dir is None:
        # autoname the directory based on fasta, or compound file
        # this will change eventually
        if fasta_file is not None:
            experiment_name = os.path.splitext(os.path.basename(fasta_file))[0]
            experiment_dir = os.path.abspath(os.path.dirname(fasta_file))
        else:
            experiment_name = os.path.splitext(os.path.basename(compounds_file))[0]
            experiment_dir = os.path.abspath(os.path.dirname(compounds_file))
        today = datetime.datetime.now()
        experiment_name += today.strftime('_%Y%m%d')
        experiment_path = os.path.join(experiment_dir, experiment_name)
    else:
        experiment_path = output_dir
    
    output_dir = os.path.abspath(experiment_path)
    intermediate_files_dir = os.path.join(output_dir, intermediate_files)
    
    print( '!!! Saving all results here: {}'.format(output_dir))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(intermediate_files_dir):
        os.makedirs(intermediate_files_dir)
    
    return output_dir, intermediate_files_dir

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
            make it a csv, tab, or pickle file'%fname)

    # remove rows that are empty
    df = df[~pd.isnull(df).all(axis=1)]
    return df

def load_compound_results(compounds_file, pactolus, output_dir, intermediate_files_dir): 
    """ load compound results"""
    print( '\n!!! LOADING COMPOUNDS')
    compounds = load_dataframe(compounds_file)
    # auto-rename pactolus columns
    if pactolus:
        compounds = reformat_pactolus(compounds)
    # remove any missing compounds
    compounds = compounds[~pd.isnull(compounds['original_compound'])]
    compounds.fillna('', inplace=True)

    if 'original_compound' not in compounds.columns:
        raise RuntimeError('Could not find "original_compound" as a column, please\
            rename the column corresponding to inchi keys for the compounds')
    
    # remove compounds not in the database or network
    print( '!!! Scrubbing compounds')
    compounds['adj'] = compounds['original_compound'].apply(
        lambda x: '-'.join(x.split('-')[:2]))
    
    my_settings = get_settings()
    reference_compounds = load_dataframe(my_settings.compounds_df)
    reference_compounds['adj'] = reference_compounds['inchi_key'].apply(
            lambda x: '-'.join(x.split('-')[:2]))

    filt = compounds.merge(reference_compounds, on='adj', how='left', suffixes=('', '_db'))
    # categorize their reason for not being searched
    not_in_db = filt[pd.isnull(filt['cpd_group'])]
    not_in_db['not_searched_reason'] = 'Not in metabolite database'
    not_in_net = filt[filt['cpd_group'] < 0]
    not_in_net['not_searched_reason'] = 'Not in similarity network yet'
    # combine into one table
    not_searched = pd.concat([not_in_db, not_in_net])
    # make the columns same as user input
    cols = compounds.columns[~compounds.columns.str.contains('adj')].tolist()
    cols.append('not_searched_reason')
    not_searched = not_searched[cols]
    # inform the user and save file
    if not_searched.shape[0] > 0:
        print( 'WARNING: some input compounds were not found in the metabolite database or chemical network; please report these compounds! (see log_unsearched_compounds.csv)')
        print( '!@# {} Compounds not being searched; see log_unsearched_compounds.csv'.format(not_searched['original_compound'].unique().shape[0]))
        not_searched.to_csv(os.path.join(output_dir,
            'log_unsearched_compounds.csv'), index=False)

    to_search = filt[filt['cpd_group'] > 0]['original_compound'].unique()
    compounds = compounds[compounds['original_compound'].isin(to_search)]

    u_cpds = compounds['original_compound'].unique()
    print( '!@# {} total input compounds to search\n'.format(len(u_cpds)))

    if 'compound_score' not in compounds.columns:
        print( 'WARNING: "compound_score" not found as a column; assuming that\
            there is no score for compounds, and setting the compound scores \
            to 1.0')
        compounds['compound_score'] = 1.0
    else:
        compounds['compound_score'] = compounds['compound_score'].apply(float)

    compounds.to_pickle(os.path.join(intermediate_files_dir, 'scrubbed_compounds.pkl'))
    return compounds

def get_settings():
    my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)
    return my_settings

def load_mrs_reaction():
    my_settings = get_settings()
    mrs_reaction_path = my_settings.mrs_reaction_path
    print( '!!! MRS-Reaction: {}'.format(mrs_reaction_path))
    mrs_reaction = load_dataframe(mrs_reaction_path)
    return mrs_reaction


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
