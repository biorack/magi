## Scoring version 2.0

import os
import sys
import argparse
import time
import pandas as pd
import numpy as np
import workflow_helpers as mg
import sqlite3

my_settings = mg.get_settings()

def read_compounds_metadata(compounds_metadata_path):
    """ Read everything in the input csv file"""
    compounds_metadata = mg.load_dataframe(compounds_metadata_path)
    return compounds_metadata

def read_compound_to_reaction(compound_to_reaction_path, path_to_database=my_settings.magi_database):
    ''' 
    Get reaction diameter, similarity score and original compound score
    Final data frame has the original_compound, rhea_ID, similarity_score,
    reaction_diameter and compound_score.
    '''
    compound_to_reaction = pd.read_csv(compound_to_reaction_path)

    ## Get reaction diameter from reaction ID
    with sqlite3.connect(path_to_database) as connection:
        query = "SELECT Retro_rules_reactions.reaction_ID, Retro_rules_reactions.diameter, Retro_rules_to_rhea_reactions.rhea_ID FROM \
            Retro_rules_reactions INNER JOIN Retro_rules_to_rhea_reactions ON Retro_rules_reactions.retro_rules_ID = Retro_rules_to_rhea_reactions.retro_rules_ID"
        database_info = pd.read_sql_query(query, connection)
    compound_to_reaction = compound_to_reaction.merge(database_info, how = "left", on="reaction_ID")
    compound_to_reaction.drop("reaction_ID", axis = 1, inplace = True)

    ## Fix problems with unavailable reactions and add penalty for compounds not matched to any reaction
    compound_to_reaction.drop_duplicates(inplace = True)
    compound_to_reaction["rhea_ID"] = compound_to_reaction["rhea_ID"].astype("str")
    return compound_to_reaction

def read_gene_to_reaction(gene_to_reaction_path):
    """
    Columns in reaction to gene are:
    query acc.	subject acc.	% query coverage per subject	alignment length	% positives	evalue	bit score	e_score	rhea_ID
    """
    gene_to_reaction = pd.read_csv(gene_to_reaction_path)
    gene_to_reaction = gene_to_reaction[["query acc.", "subject acc.", "e_score", "rhea_ID"]]
    gene_to_reaction.columns = ["gene_ID", "g2r_reference_protein", "e_score", "rhea_ID"]
    gene_to_reaction.sort_values(by="e_score", ascending = False, inplace = True)
    gene_to_reaction.drop_duplicates(inplace = True)

    gene_to_reaction["rhea_ID"] = gene_to_reaction["rhea_ID"].astype("str")
    return gene_to_reaction

def read_reaction_to_gene(reaction_to_gene_path):
    """
    Columns in reaction to gene are:
    query acc.	subject acc.	% query coverage per subject	alignment length	% positives	evalue	bit score	e_score	rhea_ID
    Columns to keep are: r2g_reference_protein gene_ID e_score rhea_ID
    """
    reaction_to_gene = pd.read_csv(reaction_to_gene_path)
    reaction_to_gene = reaction_to_gene[["query acc.", "subject acc.", "e_score", "rhea_ID"]]
    reaction_to_gene.columns = ["r2g_reference_protein", "gene_ID", "e_score", "rhea_ID"]
    reaction_to_gene.sort_values(by="e_score", ascending = False, inplace = True)
    reaction_to_gene.drop_duplicates(inplace = True)

    reaction_to_gene["rhea_ID"] = reaction_to_gene["rhea_ID"].astype("str")
    return reaction_to_gene

def merge_g2r_and_r2g_searches(compound_to_reaction, reaction_to_gene, gene_to_reaction, intermediate_files_dir):
    """
    This part of the workflow merges tables
    """
    # Start merging
    print( '\n!@# Merging final table | TLOG %s' % (time.time()))
    sys.stdout.flush()
    start = time.time()
    
    # Merge compound to reaction and reaction to gene info
    merged_dataframe = pd.merge(compound_to_reaction, reaction_to_gene, 
                                on='rhea_ID', how='left')
    
    # Make an integrated dataframe, joining on the gene
    merged_dataframe = pd.merge(merged_dataframe, gene_to_reaction, 
        on="gene_ID", suffixes=('_r2g', '_g2r'), how='outer')
    
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
    agree_idx = df[df['rhea_ID_r2g'] == df['rhea_ID_g2r']].index
    df.loc[agree_idx, 'reciprocal_score'] = 2.
    # disagreement
    disagree = df[df['rhea_ID_r2g'] != df['rhea_ID_g2r']].index
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

def calculate_scores(merged_dataframe, reciprocal_closeness, final_weights, start_time):
    """
    This part of the script calculates:
        - the homology score, which consists of two homology components and a reciprocal agreement score.
        - the reaction connection score, which represents how well a compound is matched to a reaction.
            It consists of the similarity to the original substrate of the reaction and the Retro Rules reaction diameter.
        - the final MAGI integrated score.
    It also uses the original user's compound score
    
    Inputs
    ------
    merged_dataframe: The data frame with g2r and c2g merged
    reciprocal_closeness:  The closeness treshold to use to score 
                           reciprocal agreement
    final_weights:         The weights for compound_score, compound_similarity_score, 
                            reaction_diameter_score, reciprocal_score,
                            homology_score and reaction_connection in that order

    start_time:             Time value to print the time

    Outputs
    -------
    Data fame with a homology_score, reaction_connection and MAGI_score column added.
    """
    print( '\n!@# Calculating final scores | TLOG %s' % (time.time()))
    sys.stdout.flush()
    
    # score reciprocal agreement
    merged_dataframe = reciprocal_agreement(merged_dataframe, closeness_threshold=reciprocal_closeness)
    
    # calculate homology score
    homology_scores = homology_score(merged_dataframe)
    # the nulls get a really low score
    homology_scores[pd.isnull(homology_scores)] = 1
    merged_dataframe['homology_score'] = homology_scores
    
    # reaction connection score says if the compound got connected to any
    # reaction in the database. Can't have zero because that messes up
    # geometric mean, so added a small number.
    merged_dataframe['reaction_connection'] = merged_dataframe[['rhea_ID_r2g', 'rhea_ID_g2r']]\
                                    .apply(pd.notnull).sum(axis=1) + 0.01
    
    # calculate final MAGI integrated score
    scoring_data = ['compound_score', 'similarity', 
                    'diameter', 'reciprocal_score', 
                    'homology_score', 'reaction_connection']
    
    scores = merged_dataframe[scoring_data].values
    if final_weights is not None:
        weights = np.asarray([final_weights] * scores.shape[0])
        merged_dataframe['MAGI_score'] = magi_score(scores, weights)
    else:
        merged_dataframe['MAGI_score'] = magi_score(scores)
    print( '!@# Scoring done in %s minutes' %((time.time() - start_time) / 60))
    return merged_dataframe

def format_table(df, compounds_metadata):
    """
    This function will:
        - Write gene IDs to strings if they were floats
        - Sort the dataframe based on the MAGI score and drop duplicates
    """
    print( '\n!@# Formatting final table | TLOG %s' % (time.time()))
    start = time.time()
    sys.stdout.flush()
          
    # find gene ids that are floats, convert those to strings, without the decimal
    float_entries = df['gene_ID'].apply(lambda x: isinstance(x, float))
    df.loc[float_entries, 'gene_ID'] = df.loc[float_entries, \
                    'gene_ID'].apply(lambda x: "{:.0f}".format(x))
       
    # sort the final table and drop matches with more than one reference protein
    df = df.sort_values('MAGI_score', ascending=False).drop_duplicates(
            ['original_compound', 'compound_score',
             'reciprocal_score', 'gene_ID', 'rhea_ID_r2g',
             'rhea_ID_g2r']
             )
    
    # Merge compound input metadata to the data frame
    if "compound_score" in compounds_metadata.columns:
        df = df.merge(compounds_metadata, how = "right", on = ["original_compound", "compound_score"]).sort_values('MAGI_score', ascending=False)
    else:
        df = df.merge(compounds_metadata, how = "right", on = "original_compound").sort_values('MAGI_score', ascending=False)

    df['MAGI_score'] = df['MAGI_score'].fillna(0)
    return df, start

def save_outputs(df, start, output_dir, intermediate_files_dir):
    # save the full dataframe
    df.to_csv(os.path.join(output_dir, 'magi_results.csv'), index=False)
    print( '!!! Full results saved to {}'.format(os.path.join(output_dir, 'magi_results.csv')))
    
    # save a compound-centric dataframe, where only the best row for each
    # original_compound was chosen (this is only for compound scoring, do
    # not use this for any kind of gene function analysis!)
    
    compound_centric = df[pd.notnull(df['original_compound'])]\
                         .sort_values('MAGI_score', ascending=False)\
                         .drop_duplicates(['original_compound', 'compound_score'])

    compound_centric.to_csv(os.path.join(output_dir, 
        'magi_compound_results.csv'), index=False)
    print( '!!! Compound results saved to {}'.format(os.path.join(output_dir, 'magi_compound_results.csv')))

    # Save a gene-centric dataframe.
    gene_centric = df.sort_values(['MAGI_score', 'e_score_g2r'], 
        ascending=[False, False])\
        .drop_duplicates(['gene_ID', 'rhea_ID_g2r'])
    gene_centric.to_csv(os.path.join(output_dir,
        'magi_gene_results.csv'), index=False)
    print( '!!! Gene results saved to {}'.format(os.path.join(output_dir, 'magi_gene_results.csv')))
    
    print( '!@# MAGI Scoring done in %s minutes' %((time.time() - start) / 60))
    with open(os.path.join(intermediate_files_dir, "timer.txt"),"r") as timerfile:
        main_start = float(timerfile.read())
    print( '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60))

def main():
    # Parse arguments and prepare for reaction to gene workflow
    magi_parameters = mg.general_magi_preparation()
    
    # Only run scoring if gene_to_reaction_only and compound_to_reaction_only are false
    if magi_parameters["gene_to_reaction_only"] or magi_parameters["compound_to_reaction_only"]:
        print("Not performing MAGI scoring workflow")
    else:
        if magi_parameters["merged_before_score"] is None:
            # Get paths to intermediate files and check if all files exist
            if magi_parameters["gene_to_reaction"] is not None:
                gene_to_reaction_path = mg.is_existing_file(magi_parameters["gene_to_reaction"])
            else:
                gene_to_reaction_path = mg.get_intermediate_file_path(magi_parameters["output_dir"], "gene_to_reaction_path")
            if magi_parameters["compound_to_reaction"] is not None:
                compound_to_reaction_path = mg.is_existing_file(magi_parameters["compound_to_reaction"])
            else:
                compound_to_reaction_path = mg.get_intermediate_file_path(magi_parameters["output_dir"], "compound_to_reaction_path")
            if magi_parameters["reaction_to_gene"] is not None:
                reaction_to_gene_path = mg.is_existing_file(magi_parameters["reaction_to_gene"])
            else:
                reaction_to_gene_path = mg.get_intermediate_file_path(magi_parameters["output_dir"], "reaction_to_gene_path")

            # Read intermediate files
            gene_to_reaction = read_gene_to_reaction(gene_to_reaction_path)
            compound_to_reaction = read_compound_to_reaction(compound_to_reaction_path)
            reaction_to_gene = read_reaction_to_gene(reaction_to_gene_path)
            # Merge the gene to reaction and reaction to gene tables
            merged_dataframe = merge_g2r_and_r2g_searches(gene_to_reaction = gene_to_reaction,
                                    compound_to_reaction = compound_to_reaction,
                                    reaction_to_gene = reaction_to_gene,
                                    intermediate_files_dir = magi_parameters["intermediate_files_dir"])
        else:
            merged_dataframe = pd.read_hdf(magi_parameters["merged_before_score"])
            print( '\n!@# merged_before_score successfully loaded')
        
        # Calculate MAGI scores
        start_time = time.time()
        merged_dataframe = calculate_scores(merged_dataframe, 
                                        magi_parameters["reciprocal_closeness"], 
                                        magi_parameters["final_weights"], start_time)
        # Merge more info to data frame and shuffle column order
        if magi_parameters["compounds"] is not None:
            compounds_metadata_path = mg.is_existing_file(magi_parameters["compounds"])
        else:
            compounds_metadata_path = mg.get_intermediate_file_path(magi_parameters["output_dir"], "compounds")
        compounds_metadata = read_compounds_metadata(compounds_metadata_path)
        merged_dataframe, time_start = format_table(merged_dataframe, compounds_metadata)
        # Save final output
        save_outputs(merged_dataframe, time_start, magi_parameters["output_dir"], magi_parameters["intermediate_files_dir"])

if __name__ == "__main__":
    main()