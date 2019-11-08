"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.1 scoring workflow script

Lines beginning with "@@@" are input job parameters
Lines beginning with "!!!" are verbose log info
Lines beginning with "!@#" are checkpoints/announcements

This workflow merges information from various magi sub-workflows and scores
compounds-reaction-gene combinations. The MAGI score is the geometric mean of 
the homology score, reciprocal score, reaction connection score and compound score. 

Required inputs are gene_to_reaction, compound_to_reaction
and reaction_to_gene intermediate files, or a merged_before_score file
or a directory in which all intermediate files are stored.
"""

import os
import sys
import argparse
import time
import pandas as pd
import numpy as np
import workflow_helpers_new as mg

def parse_arguments():
    def is_existing_file(filepath):
        """Checks if a file exists and return absolute path if it exists"""
        if not os.path.exists(filepath):
            msg = "{0} does not exist".format(filepath)
            raise argparse.ArgumentTypeError(msg)
        else:
            return os.path.abspath(filepath)
    
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

    try:
        """parse arguments"""
        parser = argparse.ArgumentParser()
        # required arguments
        parser.add_argument('--gene_to_reaction', type=is_existing_file, 
            help='path to gene_to_reaction file, must be in pickle format')
        parser.add_argument('--compound_to_reaction', type=is_existing_file,
            help='path to compound_to_reaction file, must be in pickle format')
        parser.add_argument('--reaction_to_gene', type=is_existing_file, 
            help='path to reaction_to_gene file, must be in pickle format')
        
        parser.add_argument('--merged_before_score', type=is_existing_file,
            help='path to merged_before_score table, must be in hdf5 format,\
            with the key "merged_before_score"')
        
        # optional runtime variables
        parser.add_argument('-o', '--output', 
            help='path to a custom output', required = False,
            type=str)

        parser.add_argument('--final_weights', 
            help='Defined weights to weight the final scoring for the scores:\
            compound_score reciprocal_score homology_score reaction_connection', 
            type=positive_number, nargs=4, default=None)
        parser.add_argument('--chemnet_penalty', 
            help='Base factor in the chemical network search level penalty', 
            type=positive_number, default=4)
        parser.add_argument('--reciprocal_closeness', 
            help='Cutoff to call a reciprocal disagreement as "close", as percent;\
            default is 75 meaning that a reciprocal disagreement will be classified\
            as "close" if the lower blast score (e score) is within 75%% of the higher\
            score', 
            type=percentage_values_to_decimal, default=0.75)
        
        parser.add_argument('--intermediate_files',
            help='What directory within --output to store intermediate files',
            type=str, default='intermediate_files')
        args = parser.parse_args()
    except argparse.ArgumentTypeError as ex:
        print(ex.message)
        sys.exit(1)
    return args

def merge_g2r_and_r2g_searches(compound_to_reaction, reaction_to_gene, gene_to_reaction, intermediate_files_dir):
    """
    This part of the workflow merges tables
    """
    # Load data frames
    gene_to_reaction_df = pd.read_pickle(gene_to_reaction)
    compound_to_reaction_df = pd.read_pickle(compound_to_reaction)
    compound_to_reaction_df["reaction_id"] = compound_to_reaction_df["reaction_id"].astype(float)
    reaction_to_gene_df = pd.read_pickle(reaction_to_gene)
    reaction_to_gene_df["reaction_id"] = reaction_to_gene_df["reaction_id"].astype(float)
    
    # Start merging
    print( '\n!@# Merging final table | TLOG %s' % (time.time()))
    sys.stdout.flush()
    start = time.time()
    
    compound_to_reaction_df = compound_to_reaction_df[['reaction_id',
                                'compound_score',
                                'original_compound', 'level', 'neighbor',
                                'note']]
    reaction_to_gene_df = reaction_to_gene_df[['subject acc.', 'reaction_id',
                                'e_score']]
    
    merged_dataframe = pd.merge(compound_to_reaction_df, reaction_to_gene_df, 
                                on='reaction_id', how='left')
    del reaction_to_gene_df
    del compound_to_reaction_df
    
    # Drop duplicate rows
    merged_dataframe.drop_duplicates(inplace=True)

    gene_to_reaction_df = gene_to_reaction_df[['query acc.', 'reaction_id',
                                                    'e_score']]

    gene_to_reaction_df.drop_duplicates(inplace=True)

    # Make an integrated dataframe, joining on the gene
    merged_dataframe = pd.merge(merged_dataframe, gene_to_reaction_df, 
        left_on='subject acc.', right_on='query acc.', 
        suffixes=('_r2g', '_g2r'), how='outer')
    del gene_to_reaction_df
    
    merged_dataframe.reset_index(inplace=True, drop=True)
    merged_dataframe.drop_duplicates(inplace=True)

    # Clean up NaNs in string columns and replace it with empty strings ('')
    def check_str(x):
        if isinstance(x, str):
            return True
        else:
            return False

    for column in merged_dataframe.columns:
        if len(merged_dataframe[column].apply(type).unique()) > 1:
            string_checked = merged_dataframe[column].apply(check_str)
            if string_checked.any():
                merged_dataframe[column].fillna('', inplace=True)

    # Clean up neighbor column and reaction to gene id column
    merged_dataframe['neighbor'] = merged_dataframe['neighbor'].astype(str)

    # Write to file
    merged_dataframe.to_hdf(os.path.join(intermediate_files_dir, 'merged_before_score.h5'),
        'merged_before_score', mode='w', format='table',
        complib='blosc', complevel=9)

    print( '!@# Final Merged table done in %s minutes'\
        %((time.time() - start) / 60))
    print( '!!! Final Merged table saved to %s'\
            % (os.path.join(intermediate_files_dir, 'merged_before_score.h5')))

    return merged_dataframe

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

def magi_score(scores, weights=np.asarray([np.nan, np.nan])):
    """
    Calculates the geometric mean of an array of numbers.
    Used to calculate one score even if sub-scores are scaled differently.

    Inputs
    ------
    scores: 1D array-like list of scores
    weights: 1D array-like list of weights

    Outputs
    -------
    The (weighted) geometric mean
    """
    
    if not isinstance(scores, np.ndarray):
        scores = np.asarray(scores)
    scores = scores.astype(float)

    if not isinstance(weights, np.ndarray):
        weights = np.asarray(weights)
    # if no weights provided, make them all 1
    if np.isnan(weights).all():
        weights = np.ones(scores.shape)
    weights = weights.astype(float)

    # Count dimensions of the score table
    a = len(scores.shape) - 1
    if a > 1:
        raise RuntimeError(
            'Scores array has too many dimensions (%s)' % (scores.shape)
            )
    if weights.shape != scores.shape:
        raise RuntimeError(
            'Weights array does not have same dimensions as Scores array: %s vs %s' % (weights.shape, scores.shape)
            )
        
    # Multiply the scores by their respective weights and sum the result
    # Divide this by the total sum of the weights
    # Take the exponential (e^x) of this value.
    geometric_mean = np.exp(np.sum(np.multiply(weights, np.log(scores)), axis=a) / np.sum(weights, axis=a))
    
    return geometric_mean

def calculate_scores(df, reciprocal_closeness, final_weights, chemnet_penalty, start_time):
    """
    This part of the script calculates:
        - the reciprocal agreement score, 
        - the homology score, 
        - the reaction connection score 
        - the final MAGI integrated score.
    
    Inputs
    ------
    df: The data frame with g2r and c2g merged
    reciprocal_closeness:  The closeness treshold to use to score 
                           reciprocal agreement
    final_weights:         The weights for compound_score, reciprocal_score,
                    homology_score and reaction_connection in that order
    chemnet_penalty:        penalty for finding compound in the 
                chemical similarity network instead of the compound itself
    start_time:             Time value to print the time

    Outputs
    -------
    Data fame with a homology_score, reaction_connection and MAGI_score column added.
    """
    print( '\n!@# Calculating final scores | TLOG %s' % (time.time()))
    sys.stdout.flush()
    
    # score reciprocal agreement
    df = reciprocal_agreement(df, closeness_threshold=reciprocal_closeness)
    
    # calculate homology score
    score = homology_score(df)
    # the nulls get a really low score
    score[pd.isnull(score)] = 1
    df['homology_score'] = score
    
    # reaction connection score says if the compound got connected to any
    # reaction in the database. Can't have zero because that messes up
    # geometric mean, so added a small number.
    df['reaction_connection'] = df[['reaction_id_r2g', 'reaction_id_g2r']]\
                                    .apply(pd.notnull).sum(axis=1) + 0.01
    
    # calculate final MAGI integrated score
    scoring_data = ['compound_score', 'reciprocal_score', \
                    'homology_score', 'reaction_connection']
    
    scores = df[scoring_data].values
    if final_weights is not None:
        weights = np.asarray([final_weights] * scores.shape[0])
        scores = magi_score(scores, weights)
    else:
        scores = magi_score(scores)
    df['MAGI_score'] = scores / (chemnet_penalty ** df['level'].values)
    print( '!@# Scoring done in %s minutes' %((time.time() - start_time) / 60))
    return df
    
    
    
    
def format_table(df):
    """
    This function will:
        - Write gene IDs to strings if they were floats
    """
    print( '\n!@# Formatting final table | TLOG %s' % (time.time()))
    start = time.time()
    sys.stdout.flush()
          
    # find gene ids that are floats, convert those to strings, without the decimal
    float_entries = df['subject acc.'].apply(lambda x: isinstance(x, float))
    df.loc[float_entries, 'subject acc.'] = df.loc[float_entries, \
                    'subject acc.'].apply(lambda x: "{:.0f}".format(x))
       
    # sort the final table and drop key duplicates
    df = df.sort_values(
        ['original_compound', 'MAGI_score'], 
        ascending=[True, False]
        ).drop_duplicates(
            ['original_compound', 'level', 'neighbor', 'compound_score',
             'reciprocal_score', 'query acc.', 'reaction_id_r2g',
             'reaction_id_g2r']
             )
    
    # Add database IDs for the g2r and r2g searches
    mrs_reaction_path = mg.get_settings(); mrs_reaction_path = mrs_reaction_path.mrs_reaction_path
    mrs_reaction = mg.load_dataframe(mrs_reaction_path)
    df = df.merge(mrs_reaction[['database_id']],
        left_on='reaction_id_r2g', right_index=True, how='left')
    df = df.merge(mrs_reaction[['database_id']],
        left_on='reaction_id_g2r', right_index=True, how='left',
        suffixes=('_r2g', '_g2r'))
    
    # Rename query acc column to gene_id column
    cols = df.columns.values
    idx = pd.np.argwhere(cols == 'query acc.')[0][0]
    cols[idx] = 'gene_id'
    df.columns = cols
    
    # Shuffle column order and remove reaction ID and subject acc columns
    df = df[['MAGI_score','gene_id', 'original_compound', 'neighbor',
        'note', 'compound_score','level','homology_score','reciprocal_score',
        'reaction_connection', 'e_score_r2g','database_id_r2g', 'e_score_g2r',
        'database_id_g2r']]
    return df, start

def save_outputs(df, start, output_dir, intermediate_files_dir):
    # save the full dataframe
    df.to_csv(os.path.join(output_dir, 'magi_results.csv'), index=False)
    print( 'full results saved to {}'.format(os.path.join(output_dir, 'magi_results.csv')))
    
    # save a compound-centric dataframe, where only the best row for each
    # original_compound was chosen (this is only for compound scoring, do
    # not use this for any kind of gene function analysis!)
    
    compound_centric = df[pd.notnull(df['original_compound'])]\
                         .sort_values('MAGI_score', ascending=False)\
                         .drop_duplicates(['original_compound', 'compound_score'])
    input_compounds = pd.read_pickle(os.path.join(intermediate_files_dir, "scrubbed_compounds.pkl"))
    compound_centric = pd.merge(
        compound_centric, input_compounds,
        on=['original_compound', 'compound_score'],
        how='right')
    compound_centric.to_csv(os.path.join(output_dir, 
        'magi_compound_results.csv'), index=False)
    
    # Save a gene-centric dataframe.
    gene_centric = df.sort_values(['MAGI_score', 'e_score_g2r'], 
        ascending=[False, False])\
        .drop_duplicates(['gene_id', 'database_id_g2r'])
    gene_centric.to_csv(os.path.join(output_dir,
        'magi_gene_results.csv'), index=False)
    
    print( '!@# MAGI Scoring done in %s minutes' %((time.time() - start) / 60))
    with open(os.path.join(intermediate_files_dir, "timer.txt"),"r") as timerfile:
        main_start = float(timerfile.read())
    print( '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60))
    print( '!!! final results stored to %s' \
            %(os.path.join(output_dir, 'magi_results.csv')))
    
    
def main():
    # Parse arguments and prepare for reaction to gene workflow
    magi_parameters = mg.general_magi_preparation()

    # Get paths to intermediate files and check if all files exist
    if magi_parameters["gene_to_reaction"] is not None:
        gene_to_reaction_path = mg.is_existing_file(magi_parameters["gene_to_reaction"])
    else:
        gene_to_reaction_path = mg.get_intermediate_file_path(magi_parameters["intermediate_files_dir"], "gene_to_reaction_path")
    if magi_parameters["compound_to_reaction"] is not None:
        compound_to_reaction_path = mg.is_existing_file(magi_parameters["compound_to_reaction"])
    else:
        compound_to_reaction_path = mg.get_intermediate_file_path(magi_parameters["intermediate_files_dir"], "compound_to_reaction_path")
    if magi_parameters["reaction_to_gene"] is not None:
        reaction_to_gene_path = mg.is_existing_file(magi_parameters["reaction_to_gene"])
    else:
        reaction_to_gene_path = mg.get_intermediate_file_path(magi_parameters["intermediate_files_dir"], "reaction_to_gene_path")

    # Merge the gene to reaction and reaction to gene tables
    if magi_parameters["merged_before_score"] is None:
        merged_dataframe = merge_g2r_and_r2g_searches(gene_to_reaction = gene_to_reaction_path,
                                   compound_to_reaction = compound_to_reaction_path,
                                   reaction_to_gene = reaction_to_gene_path,
                                   intermediate_files_dir = magi_parameters["intermediate_files_dir"])
    else:
        merged_dataframe = pd.read_hdf(magi_parameters["merged_before_score"])
        print( '\n!@# merged_before_score successfully loaded')
    
    # Calculate MAGI scores
    start_time = time.time()
    merged_dataframe = calculate_scores(merged_dataframe, 
                                       magi_parameters["reciprocal_closeness"], 
                                       magi_parameters["final_weights"], 
                                       magi_parameters["chemnet_penalty"], start_time)
    # Merge more info to data frame and shuffle column order
    merged_dataframe, time_start = format_table(merged_dataframe)
    # Save final output
    save_outputs(merged_dataframe, time_start, magi_parameters["output_dir"], magi_parameters["intermediate_files_dir"])

if __name__ == "__main__":
    main()    