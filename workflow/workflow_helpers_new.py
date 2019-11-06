import os
import pandas as pd
import numpy as np
import sys
import subprocess
import warnings
import datetime
import argparse
from multiprocessing import cpu_count as counting_cpus

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from local_settings import local_settings as settings_loc

def parse_arguments():
    def is_existing_file(filepath):
        """Checks if a file exists and return absolute path if it exists"""
        if not os.path.exists(filepath):
            msg = "{0} does not exist".format(filepath)
            raise argparse.ArgumentTypeError(msg)
        else:
            return os.path.abspath(filepath)
    def percentage_values_to_decimal(percentage):
        """Turns the blast filter and reciprocal closeness percentages 
        into decimal numbers"""
        try:
            percentage = int(percentage)
        except:
            msg = "Please enter an integer value"
            raise argparse.ArgumentTypeError(msg)        
        if percentage > 100:
            msg = "Max value is 100"
            raise argparse.ArgumentTypeError(msg)
        elif percentage < 0:
            msg = "Value cannot be negative"
            raise argparse.ArgumentTypeError(msg)
        else:
            decimal = percentage/100.
        return decimal
    
    def positive_number(number):
        """Checks if none of the number/numbers are negative"""
        try:
            number = float(number)
        except:
            msg = "Please enter a numeric value"
            raise argparse.ArgumentTypeError(msg)        
        if number < 0:
            msg = "Value cannot be negative"
            raise argparse.ArgumentTypeError(msg)
        else:
            return number
    
    def set_cpu_count(cpu_count):
        max_cpu = counting_cpus()  
        if cpu_count == 0:
            cpu_count = max_cpu    
        if cpu_count > max_cpu:
            msg = "ERROR: You have exceeded the cpus on this machine ({})".format(max_cpu)
            raise argparse.ArgumentTypeError(msg)
        return cpu_count
    try:
        """parse arguments"""
        parser = argparse.ArgumentParser()
        # required arguments
        required_args = parser.add_argument_group('Required arguments')
        required_args.add_argument('-f', '--fasta', type=is_existing_file,
            help='path to fasta file of genes in sample')
        required_args.add_argument('-c', '--compounds', type=is_existing_file,
            help='path to observed compounds file')
        
        # jump-start the script after certain computations
        start_halfway_args = parser.add_argument_group('Arguments to jump-start the script after certain computations')
        start_halfway_args.add_argument('--gene_to_reaction', type=is_existing_file,
            help='path to gene_to_reaction file, must be in pickle format')
        start_halfway_args.add_argument('--compound_to_reaction', type=is_existing_file,
            help='path to compound_to_reaction file, must be in pickle format')
        start_halfway_args.add_argument('--reaction_to_gene', type=is_existing_file,
            help='path to reaction_to_gene file, must be in pickle format')
        start_halfway_args.add_argument('--merged_before_score', type=is_existing_file,
            help='path to merged_before_score table, must be in hdf5 format,\
            with the key "merged_before_score"')
        # Use this if only a part of the workflow should be run
        stop_halfway_args = parser.add_argument_group('Arguments to run a part of the script')
        stop_halfway_args.add_argument('--gene_to_reaction_only',
            help="Use this parameter if you are only interested in the gene to reaction search", 
            action='store_true', default=False)
        
        
        # optional runtime variables
        optional_args = parser.add_argument_group("Optional runtime variables")
        optional_args.add_argument('-a', '--annotations', type=is_existing_file,
            help='path to annotation file for genes in sample', 
            default=None)
        optional_args.add_argument('-n', '--cpu_count', 
            help='number of cpus to use for multiprocessing. Default is to use max!', 
            type=int, default=0)
        optional_args.add_argument('-o', '--output', 
            help='path to a custom output', 
            type=str)
        optional_args.add_argument('-l', '--level', 
            help='how many levels deep to search the chemical network', 
            type=int, choices=[0,1,2,3], default=2)
        optional_args.add_argument('--legacy', dest='legacy', action='store_true',
            help='use legacy tautomer searching; default is no')
        optional_args.add_argument('--no-legacy', dest='legacy', action='store_false',
            help='use precomputed compound-to-reaction; default is yes')
        optional_args.set_defaults(legacy=False)
        optional_args.add_argument('--mute', 
            help='mutes pandas warnings', 
            action='store_true')
        optional_args.add_argument('--pactolus', 
            help='Flag to tell MAGI that the compounds input is a pactolus file', 
            action='store_true')
        optional_args.add_argument('--test', 
            help='TBD: run MAGI only on the first # of pactolus compounds', 
            type=int)
        optional_args.add_argument('--debug', 
            help='TBD: prints a lot of info', 
            action='store_true')
        optional_args.add_argument('--blast_filter', 
            help='How stringent to filter the top BLAST results, as percent;\
            default is 85 meaning that only BLAST results within 85%% of the top\
            result will be taken.', 
            type=percentage_values_to_decimal, default=0.85)
        optional_args.add_argument('--reciprocal_closeness', 
            help='Cutoff to call a reciprocal disagreement as "close", as percent;\
            default is 75 meaning that a reciprocal disagreement will be classified\
            as "close" if the lower blast score (e score) is within 75%% of the higher\
            score', 
            type=percentage_values_to_decimal, default=0.75)
        optional_args.add_argument('--final_weights', 
            help='Defined weights to weight the final scoring for the scores:\
            compound_score reciprocal_score homology_score reaction_connection', 
            type=positive_number, nargs=4, default=None)
        optional_args.add_argument('--chemnet_penalty', 
            help='Base factor in the chemical network search level penalty', 
            type=positive_number, default=4)
        optional_args.add_argument('--intermediate_files',
            help='What directory within --output to store intermediate files',
            type=str, default='intermediate_files')
        
        # Parameters for accurate mass search
        mass_search_args = parser.add_argument_group("Arguments for the optional accurate mass search")
        mass_search_args.add_argument('--accurate_mass_search',
                            type=str, choices=['pos','neg','neut'], default=None,
                            help = "Perform accurate mass search on m/z values in original_compound column in the input compounds file. \
                            Specify if the masses are measured in negative mode, positive mode or if they have been transformed to neutral masses."
                            )
        mass_search_args.add_argument('--accurate_mass_search_only',
                            action='store_true', default=False,
                            help = "If this parameter is used, MAGI will stop after performing the accurate mass search."
                            )
        mass_search_args.add_argument('--adduct_file',
                            type=str, default=None,
                            help="optionally specify which adducts to investigate. If not specified, it will search for M+,M+H,M+NH4,M+Na in positive mode or M-H,M+Cl,M+FA-H,M+Hac-H in negative mode."
                            )
        mass_search_args.add_argument('--ppm_cutoff', 
            help='The ppm cutoff for the accurate mass search. Default is 10 ppm.', 
            type=int, default=10)
        
        args = parser.parse_args()
        
        # Check parameters and set number of required CPUs
        if args.fasta is None and args.compounds is None:
            raise argparse.ArgumentTypeError('ERROR: either FASTA or metabolites file is required')
        args.cpu_count = set_cpu_count(args.cpu_count)
    except argparse.ArgumentTypeError as ex:
        print(ex.message)
        sys.exit(1)
    return args

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

def print_version_info():
    """Print versions of modules that may be troublesome for magi."""
    print('!!! Python version:'+ sys.version)
    print('!!! numpy version: '+ np.__version__)
    print('!!! pandas version:'+ pd.__version__)
    #print('!!! pickle version:'+ pickle.__version__)
    print('#'*80)

def print_parameters(args):
    print('~~~~~PARAMETERS~~~~~~')

    # print your paths in stdout for logging
    print('@@@ FASTA file input: %s' %(args.fasta))
    print('@@@ Compound input: %s' %(args.compounds))
    if args.annotations is not None:
        print( '@@@ annotations input: %s' %(args.annotations))
    # print parameters
    if args.fasta is not None: 
        print( '@@@ BLAST filter: %s' % (args.blast_filter))
    if args.compounds is not None:
        print( '@@@ Using precomputed compound results: %s' % (not args.legacy))
        print( '@@@ Chemnet search to level %s' % (args.level))
        print( '@@@ Reciprocal closeness: %s' % (args.reciprocal_closeness))
        print( '@@@ Chemnet base penalty: %s' % (args.chemnet_penalty))
    print( '@@@ MAGI score weights: %s' % (args.final_weights))
    print( '@@@ Using %s CPUs' % (args.cpu_count))
    
    if args.gene_to_reaction is not None:
        print( '@@@ gene_to_reaction input: %s' %(args.gene_to_reaction))
    
    if args.compound_to_reaction is not None:
        print( '@@@ compound_to_reaction input: %s' %(args.compound_to_reaction))
    
    if args.reaction_to_gene is not None:
        print( '@@@ reaction_to_gene input: %s' %(args.reaction_to_gene))
    if args.mute:
        print( '!!! Warnings are muted')
        warnings.filterwarnings('ignore')
        
def general_magi_preparation():
    """
    This function prepares for a MAGI run. It:
        - parses arguments
        - makes an output file directory,
        - prints versions of possibly troublesome modules 
        - prints input parameters.
    It returns a dictionary with parameters for the MAGI run.
    """
    args = parse_arguments()
    print_version_info()
    print_parameters(args)
    output_dir, intermediate_files_dir = make_output_dirs(output_dir=args.output, 
                                                          fasta_file=args.fasta, 
                                                          compounds_file=args.compounds, 
                                                          intermediate_files=args.intermediate_files)
    magi_parameters = vars(args)
    magi_parameters["output_dir"] = output_dir
    magi_parameters["intermediate_files_dir"] = intermediate_files_dir
    return magi_parameters
    
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
