"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.1 accurate mass search workflow script

Lines beginning with "@@@" are input job parameters
Lines beginning with "!!!" are verbose log info
Lines beginning with "!@#" are checkpoints/announcements

Required inputs is a Compounds file.
The Compounds file should be in some standard table format (CSV, tab-
delimited, etc.), and is required to have a column named 
"original_compound" (case-sensitive). In this column should be m/z values which need to be mass-searched.

The output will be a compounds.csv file with the suffix _mass_searched. This file has inchi keys in
a column called original_compound, and the m/z values can now be found in a column called original_mz.
"""
import sys
import os
import pandas as pd
import workflow_helpers_new as mg

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

def accurate_mass_search(features_to_search, polarity, adducts, ppm_cutoff, reference_compounds):
    """
    Given an input filename, finds all compounds in the reference_compounds database matching within the given ppm error.
    
    Inputs
    ------
    features_to_search: Dataframe with m/z values. m/z values need to be in a column called 'original_compound'
    polarity: either pos, neg or neut
    adducts: list of adducts to search.
    ppm_cutoff: ppm error cuttof for accurate mass search
    reference_compounds: default is the unique_compound_groups_magi.pkl file in local settings (compounds_df).
    
    Outputs
    -------
    mass_searched_filename: filename of dataframe with columns query_mass, target_mass, ppm,
        original_compound, compound_score corresponding to the mass_list
        masses, found masses, ppm difference, inchikeys for those
        compounds, and a function of the ppm error as a compound score

    compound_score = 1
    """
    # rename original_compounds column
    columns = features_to_search.columns.values
    if 'original_compound' not in columns:
        raise RuntimeError('no original_compound')
    columns[columns == 'original_compound'] = 'original_mz'
    features_to_search.columns = columns

    # set up data container
    data = {
        'original_compound': [],
        'searched_adduct': [],
        'original_mz': [],
        'ppm_error': [],
        'compound_score': []
    }
    # accurate mass search and store results
    for mz in features_to_search['original_mz'].unique():
        for adduct in adducts:
            if adduct != '':
                neutral_mass = mz_neutral_transform(mz, adduct)
            else:
                neutral_mass = mz
            found_compounds = accurate_mass_match(neutral_mass,
                                                  compound_df=reference_compounds,
                                                  ppm=ppm_cutoff
                                                 )
            if found_compounds is not None:
                for cpd in found_compounds:
                    data['original_compound'].append(cpd[0])
                    data['ppm_error'].append(cpd[1])
                    data['compound_score'].append(1)
                    data['searched_adduct'].append(adduct)
                    data['original_mz'].append(mz)

    # merge with user input 
    df = pd.DataFrame(data)
    features_to_search = features_to_search.merge(df, on='original_mz', how='left')    
    return features_to_search

def prepare_adducts(polarity, adducts_pos = 'M+,M+H,M+NH4,M+Na', adducts_neg = 'M-H,M+Cl,M+FA-H,M+Hac-H'):
    """
    Make list of adducts to search for

    Inputs:
    ------
    adducts_pos:   string to read positive adducts from. 
    adducts_neg:   string to read positive adducts from. 
    polarity:       pos, neg or neut
    
    Outputs:
    ------
    adducts:        list of adducts to use. Empty list if the polarity is neutral.
    """
    # get adducts according to polarity
    if polarity == 'pos':
        adducts = adducts_pos.split(',')
    elif polarity == 'neg':
        adducts = adducts_neg.split(',')
    elif polarity == 'neut':
        adducts = ['']
    else:
        # Error when polarity is not in pos, neg or neut
        msg = 'Could not understand polarity %s' % (polarity)
        sys.exit(msg)
    return adducts

def workflow(compounds_file, adducts_pos, adducts_neg, polarity, ppm_cutoff, reference_compounds):
    """
    Run an accurate mass search on m/z values and return candidate compounds

    Inputs
    ------
    compounds_file: absolute path to the file with m/z values that need to be searched
    adduct_file:    list of adducts to use.
    polarity:       polarity of the search. Either positive, negative or neutral
    ppm_cutoff:     ppm cutoff to use for the mass search
    """

    print("\n!!! Performing accurate mass search for {}".format(compounds_file))
    adducts = prepare_adducts(adducts_pos = adducts_pos, adducts_neg = adducts_neg, polarity = polarity)

    # Read compounds file
    features_to_search = pd.read_csv(compounds_file)
    mass_searched_compounds = accurate_mass_search(features_to_search, polarity, adducts, ppm_cutoff, reference_compounds)
    return mass_searched_compounds

def format_output(mz_filename, output_dir, mass_searched_compounds):
    """
    Write mass_searched_compounds dataframe to file
    """
    # save the new table
    mass_searched_filename = os.path.splitext(mz_filename)[0] + '_mass_searched.csv'
    mass_searched_compounds.to_csv(mass_searched_filename, index = False)
    print("\n!!! Accurate mass search done. Mass-searched file stored in {}".format(mass_searched_filename))    

def main():
    # Parse arguments and prepare for compounds to reaction workflow
    magi_parameters = mg.general_magi_preparation()
    if magi_parameters["is_mass_search"]:
        # Check mass-search-specific parameter
        if magi_parameters["polarity"] is None:
            sys.exit("Please specify polarity")
        if magi_parameters["adducts_pos"] is None and magi_parameters["adducts_neg"] is None:
            sys.exit("Please specify adducts")

        # Load reference compounds list with reference masses
        my_settings = mg.get_settings()
        reference_compounds = mg.load_dataframe(my_settings.compounds_df)

        # Run workflow
        mass_searched_compounds = workflow(compounds_file = magi_parameters["compounds"], 
                    adducts_pos = magi_parameters["adducts_pos"], 
                    adducts_neg = magi_parameters["adducts_neg"],
                    polarity = magi_parameters["polarity"], 
                    ppm_cutoff = magi_parameters["ppm"], 
                    reference_compounds = reference_compounds)
        # Write to file
        format_output(mz_filename = magi_parameters["compounds"], 
                    output_dir = magi_parameters["output_dir"],
                    mass_searched_compounds = mass_searched_compounds)

if __name__ == "__main__":
    main()
