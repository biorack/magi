"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.0a workflow script

Lines beginning with "@@@" are input job parameters
Lines beginning with "!!!" are verbose log info
Lines beginning with "!@#" are checkpoints/announcements

Required inputs are a FASTA file and a Compounds file.
The FASTA file should be in standard FASTA format, E.g.:
'''fasta
> UNIQUE_GENE_ID OTHER_INFORMATION
AMINO_ACID_SEQUENCE

> UNIQUE_GENE_ID OTHER_INFORMATION
AMINO_ACID_SEQUENCE
'''
*Please note the space in between UNIQUE_GENE_ID and OTHER_INFORMATION*

The FASTA file will be used to make a BLAST database, and to make a
gene table with these columns:
Gene_ID        | header                           | sequence
-----------------------------------------------------------------------
UNIQUE_GENE_ID | UNIQUE_GENE_ID OTHER_INFORMATION | AMINO_ACID_SEQUENCE

The Compounds file should be in some standard table format (CSV, tab-
delimited, etc.), and is required to have a column named 
"original_compound" (case-sensitive). In this column should be standard
InChI Keys representing a compound structure.

The Compounds table may also have a column named "compound_score" 
(case-sensitive), where the user can provide a score for each 
compound-row. If this column name does not exist, one will be created
and populated with 1.0
"""
# TODO: make fixed typing in columns, particularly reaction_id and neighor
# TODO: find the best way to level adjust the MAGI score. 10^n is too strong

import sys
import argparse
import os
import warnings
import multiprocessing as mp
import pandas as pd
import numpy as np
import time
import pickle
import datetime

# load local settings
sys.path.insert(
    0,
    '/global/homes/e/erbilgin/repos/magi/')
from local_settings import local_settings as settings_loc
my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)

# print versions of troublesome modules
print '!!! Python version:', sys.version
print '!!! numpy version: ', np.__version__
print '!!! pandas version:', pd.__version__
print '!!! pickle version:', pickle.__version__
print '#'*80

# parse arguments
parser = argparse.ArgumentParser()
# required
parser.add_argument('-f', '--fasta', 
	help='path to fasta file of genes in sample', 
	required=True)
parser.add_argument('-c', '--compounds', 
	help='path to observed compounds file', 
	required=True)

# jump-start the script after certain computations
parser.add_argument('--gene_to_reaction', 
	help='path to gene_to_reaction file, must be in pickle format')
parser.add_argument('--compound_to_reaction', 
	help='path to compound_to_reaction file, must be in pickle format')
parser.add_argument('--reaction_to_gene', 
	help='path to reaction_to_gene file, must be in pickle format')
parser.add_argument('--merged_before_score', 
	help='path to merged_before_score table, must be in hdf5 format,\
	with the key "merged_before_score"')

# optional runtime variables
parser.add_argument('-a', '--annotations', 
	help='path to annotation file for genes in sample', 
	default=None)
parser.add_argument('-n', '--cpu_count', 
	help='number of cpus to use for multiprocessing. Default is to use max!', 
	type=int, default=0)
parser.add_argument('-o', '--output', 
	help='path to a custom output', 
	type=str)
parser.add_argument('-l', '--level', 
	help='how many levels deep to search the chemical network', 
	type=int, choices=[1,2,3], default=2)
parser.add_argument('-t', '--tautomer', 
	help='include tautomers in search; default is true', 
	choices=[True, False], default=True)
parser.add_argument('--mute', 
	help='mutes pandas warnings', 
	action='store_true')
parser.add_argument('--pactolus', 
	help='Flag to tell MAGI that the compounds input is a pactolus file', 
	action='store_true')
parser.add_argument('--test', 
	help='TBD: run MAGI only on the first # of pactolus compounds', 
	type=int)
parser.add_argument('--debug', 
	help='TBD: prints a lot of info', 
	action='store_true')
parser.add_argument('--blast_filter', 
	help='How stringent to filter the top BLAST results, as percent;\
	default is 85 meaning that only BLAST results within 85% of the top\
	result will be taken.', 
	type=int, choices=range(0, 101), default=85)
parser.add_argument('--reciprocal_closeness', 
	help='Cutoff to call a reciprocal disagreement as "close", as percent;\
	default is 75 meaning that a reciprocal disagreement will be classified\
	as "close" if the lower blast score (e score) is within 75% of the higher\
	score', 
	type=int, choices=range(0, 101), default=75)
parser.add_argument('--final_weights', 
	help='Defined weights to weight the final scoring for the scores:\
	compound_score reciprocal_score homology_score reaction_connection', 
	type=float, nargs=4, default=None)
parser.add_argument('--chemnet_penalty', 
	help='Base factor in the chemical network search level penalty', 
	type=float, default=4)

args = parser.parse_args()

# Convert paths to absolute paths
args.fasta = os.path.abspath(args.fasta)
args.compounds = os.path.abspath(args.compounds)

print '~~~~~PARAMETERS~~~~~~'

# print your paths in stdout for logging
print '@@@ FASTA file input: %s' %(args.fasta)
print '@@@ compound input: %s' %(args.compounds)
if args.annotations is not None:
	args.annotations =  os.path.abspath(args.annotations)
	print '@@@ annotations input: %s' %(args.annotations)

# convert parameters
args.blast_filter = args.blast_filter / 100.
args.reciprocal_closeness = args.reciprocal_closeness / 100.

# check parameters
if args.blast_filter < 0:
	raise RuntimeError('argument blast_filter cannot be negative!')
if args.reciprocal_closeness < 0:
	raise RuntimeError('argument reciprocal_closeness cannot be negative!')
if args.final_weights is not None:
	for w in args.final_weights:
		if w < 0:
			raise RuntimeError('argument final_weights cannot be negative!')
if args.chemnet_penalty < 0:
	raise RuntimeError('argument chemnet_penalty cannot be negative!')

# print parameters
print '@@@ Searching Tautomers: %s' % (args.tautomer)
print '@@@ Chemnet search to level %s' % (args.level)
print '@@@ BLAST filter: %s' % (args.blast_filter)
print '@@@ Reciprocal closeness: %s' % (args.reciprocal_closeness)
print '@@@ MAGI score weights: %s' % (args.final_weights)
print '@@@ chemnet base penalty: %s' % (args.chemnet_penalty)
print '@@@ Using %s CPUs' % (args.cpu_count)

# setup multiprocessing helpers based on input params
def keep_top_blast_helper(x, param=args.blast_filter):
	"""
	x is the normal input
	param is the defined parameter
	"""
	return mg.keep_top_blast(x, filt=param)



max_cpu = mp.cpu_count()

if args.cpu_count > max_cpu:
	raise RuntimeError('You have exceeded the cpus on this machine (%s)!' 
		%(max_cpu))
if args.cpu_count == 0:
	args.cpu_count = max_cpu

if args.gene_to_reaction is not None:
	args.gene_to_reaction =  os.path.abspath(args.gene_to_reaction)
	if not os.path.isfile(args.gene_to_reaction):
		raise IOError('%s does not exist!' %(args.gene_to_reaction))
	else:
		print '@@@ gene_to_reaction input: %s' %(args.gene_to_reaction)

if args.compound_to_reaction is not None:
	args.compound_to_reaction =  os.path.abspath(args.compound_to_reaction)
	if not os.path.isfile(args.compound_to_reaction):
		raise IOError('%s does not exist!' %(args.compound_to_reaction))
	else:
		print '@@@ compound_to_reaction input: %s' %(args.compound_to_reaction)

if args.reaction_to_gene is not None:
	args.reaction_to_gene =  os.path.abspath(args.reaction_to_gene)
	if not os.path.isfile(args.reaction_to_gene):
		raise IOError('%s does not exist!' %(args.reaction_to_gene))
	else:
		print '@@@ reaction_to_gene input: %s' %(args.reaction_to_gene)

if args.mute:
	print '!!! Warnings are muted'
	warnings.filterwarnings('ignore')

# import workflow helpers after all the argument checking
sys.path.insert(
	0,
	os.path.join(my_settings.repo_location, 'workflow/helpertools'))
print '!!! importing workflow helpers'
import workflow_helpers as mg

# path to MAGI data storage
MAGI_PATH = my_settings.magi_results_storage

# set up where the results will be stored
if args.output is None:
	experiment_name = args.fasta.split('/')[-1].split('.')[0]
	today = datetime.datetime.now()
	experiment_name += today.strftime('_%Y%m%d')
	experiment_path = '%s/%s' %(MAGI_PATH, experiment_name)
else:
	experiment_path = args.output

experiment_path = os.path.abspath(experiment_path)

print '!!! Saving all results here:', experiment_path
if not os.path.isdir(experiment_path):
	os.makedirs(experiment_path)

main_start = time.time() # overall program timer

# load genome
print '\n!!! LOADING GENOME'
genome, genome_db_path = mg.load_genome(args.fasta, MAGI_PATH, 
										annotation_file=args.annotations)

# load pactolus results
print '\n!!! LOADING COMPOUNDS'
compounds = mg.load_dataframe(args.compounds)
# auto-rename pactolus columns
if args.pactolus:
	compounds = mg.reformat_pactolus(compounds)
# remove any missing compounds
compounds = compounds[~pd.isnull(compounds['original_compound'])]
compounds.fillna('', inplace=True)

if 'original_compound' not in compounds.columns:
	raise RuntimeError('Could not find "original_compound" as a column, please\
		rename the column corresponding to inchi keys for the compounds')
u_cpds = compounds['original_compound'].unique()
print '!@#', len(u_cpds), 'total input compounds to search\n'

if 'compound_score' not in compounds.columns:
	print 'WARNING: "compound_score" not found as a column; assuming that\
		there is no score for compounds, and setting the compound scores \
		to 1.0'
	compounds['compound_score'] = 1.0
else:
	compounds['compound_score'] = compounds['compound_score'].apply(float)

# Conduct gene to reaction search
if args.gene_to_reaction is None:
	print '!@# Conducting gene to reaction search'
	start = time.time()
	gene_blast = mg.multi_blast(genome.index, genome, mg.refseq_dbpath, 
		experiment_path, raise_blast_error=False, cpu=args.cpu_count)

	print '!@# Homology searching done in %s minutes' \
			%((time.time() - start) / 60)
	gene_blast.to_pickle(os.path.join(experiment_path, 'gene_blast.pkl'))
	print '!!! g2r blast results saved to %s' \
			%(os.path.join(experiment_path, 'g2r_blast.pkl'))

	start = time.time()
	gene_to_reaction = mg.refseq_to_reactions(gene_blast, 'subject acc.')
	del gene_blast
	gene_groups = gene_to_reaction.groupby('query acc.')
	multidx = gene_groups['e_score'].apply(keep_top_blast_helper).index
	idx = multidx.levels[1]
	gene_to_reaction_top = gene_to_reaction.loc[idx]
	del gene_to_reaction
	print '!@# gene_to_reaction table completed in %s minutes' \
			%((time.time() - start) / 60)

	gene_to_reaction_top.to_pickle(os.path.join(experiment_path, 
											'gene_to_reaction.pkl'))
	print '!!! gene to reaction results saved to %s' \
			%(os.path.join(experiment_path, 'gene_to_reaction.pkl'))
else:
	gene_to_reaction_top = pd.read_pickle(args.gene_to_reaction)
	print '\n!@# gene_to_reaction successfully loaded'
del genome

# compound to reaction search
if args.compound_to_reaction is None:
	print '\n!@# Conducting compound to reaction search'
	sys.stdout.flush()
	start = time.time()

	def connect_compound_to_reaction_mp_helper(inchikey, 
											tautomer=args.tautomer, 
											neighbor_level=args.level):
	    try:
	    	out = mg.connect_compound_to_reaction(inchikey, 
									    	tautomer=tautomer, 
									    	neighbor_level=neighbor_level)
	    except Exception as e:
	    	print inchikey
	    	sys.stdout.flush()
	    	raise RuntimeError('offending inchikey: %s; error message: %s' \
	    						%(inchikey, e.args))
	    return out

	input_compounds = compounds['original_compound'].unique()

	p = mp.Pool(args.cpu_count)
	out = p.map(connect_compound_to_reaction_mp_helper, input_compounds)
	p.close()
	p.terminate()

	compound_to_reaction = pd.concat(out)
	del out
	compound_to_reaction.reset_index(inplace=True, drop=True)

	# connect the compound score
	compound_to_reaction = pd.merge(compounds, compound_to_reaction, 
									on='original_compound', how='inner')

	compound_to_reaction.to_pickle(os.path.join(experiment_path, 
											'compound_to_reaction.pkl'))

	print '!@# compound_to_reaction table done in %s minutes'\
			%((time.time()-start)/60)
	print '!!! compound_reaction table saved to %s'\
			% (os.path.join(experiment_path, 'compound_to_reaction.pkl'))
else:
	compound_to_reaction = pd.read_pickle(args.compound_to_reaction)
	print '\n!@# compound_to_reaction successfully loaded'

# reaction to gene search
if args.reaction_to_gene is None:
	print '\n!@# Conducting reaction to gene search'
	sys.stdout.flush()
	start = time.time()

	# set up a list of reference sequences to blast against the genome
	reactions = compound_to_reaction[compound_to_reaction\
						['reaction_id'] != '']['reaction_id'].tolist()
	reactions_refseqs = mg.mrs_reaction.loc[reactions, 'refseq_id']
	del reactions
	reactions_refseqs = reactions_refseqs[reactions_refseqs != '']
	rseq_list = []
	for reaction in reactions_refseqs:
	    for rseq in reaction.split('|'):
	        if rseq != '':
	            rseq_list.append(rseq)
	rseq_list = list(set(rseq_list))

	# rseq_list is the "query_list" for multi_blast()
	# query_full_table is the refseq table
	# database_path is the path to the genome's blast database
	print '!!!', len(rseq_list), 'reference sequences to search'
	sys.stdout.flush()

	reaction_to_gene_blast = mg.multi_blast(rseq_list, mg.refseq, 
		genome_db_path, experiment_path, cpu=args.cpu_count, 
		raise_blast_error=False)

	reaction_to_gene = mg.refseq_to_reactions(reaction_to_gene_blast,
		'query acc.')
	del reaction_to_gene_blast

	reaction_to_gene.to_pickle(os.path.join(experiment_path,
		'reaction_blast.pkl'))
	print '!!! r2g blast results saved to %s' \
			%(os.path.join(experiment_path, 'r2g_blast.pkl'))

	reaction_groups = reaction_to_gene.groupby('query acc.')
	multidx = reaction_groups['e_score'].apply(keep_top_blast_helper).index
	del reaction_groups
	idx = multidx.levels[1]
	reaction_to_gene_top = reaction_to_gene.loc[idx]
	del reaction_to_gene
	reaction_to_gene_top.to_pickle(os.path.join(experiment_path, 
											'reaction_to_gene.pkl'))
	print '!@# reaction_to_gene table done in %s minutes'\
			%((time.time()-start)/60)
	print '!!! reaction_to_gene table saved to %s'\
			% (os.path.join(experiment_path, 'reaction_to_gene.pkl'))
else:
	reaction_to_gene_top = pd.read_pickle(args.reaction_to_gene)
	print '\n!@# reaction_to_gene successfully loaded'

if args.merged_before_score is None:
	print '\n!@# Merging final table'
	sys.stdout.flush()
	start = time.time()

	compound_to_gene = pd.merge(compound_to_reaction, reaction_to_gene_top, 
								on='reaction_id', how='left')
	del reaction_to_gene_top
	del compound_to_reaction

	compound_to_gene_small = compound_to_gene[['subject acc.', 'reaction_id',
								'e_score', 'compound_score',
								'original_compound', 'level', 'neighbor',
								'note']]
	del compound_to_gene

	# okay to drop duplicates, because i only care about these columns 
	# anyway; if these are duplicated then other information doesn't really 
	# matter or can easily be re-expanded by joining 
	compound_to_gene_small.drop_duplicates(inplace=True)

	gene_to_reaction_small = gene_to_reaction_top[['query acc.', 'reaction_id',
													'e_score']]
	del gene_to_reaction_top
	gene_to_reaction_small.drop_duplicates(inplace=True)

	# Make an integrated dataframe, joining on the gene
	df = pd.merge(compound_to_gene_small, gene_to_reaction_small, 
		left_on='subject acc.', right_on='query acc.', 
		suffixes=('_r2g', '_g2r'), how='outer')

	df.reset_index(inplace=True, drop=True)
	df.drop_duplicates(inplace=True)

	# Clean up reaction_id_r2g column
	idx = df[df['reaction_id_r2g'] == ''].index
	df.loc[idx, 'reaction_id_r2g'] = np.nan
	df['reaction_id_r2g'] = df['reaction_id_r2g'].astype(float)

	# Clean up stupid NaNs in string columns
	def check_str(x):
	    if isinstance(x, str):
	        return True
	    else:
	        return False

	for c in df.columns:
	    if len(df[c].apply(type).unique()) > 1:
	        string_checked = df[c].apply(check_str)
	        if string_checked.any():
	            df[c].fillna('', inplace=True)

	# Clean up neighbor column
	df['neighbor'] = df['neighbor'].astype(str)

	df.to_hdf(os.path.join(experiment_path, 'merged_before_score.h5'),
		'merged_before_score', mode='w', format='table',
		complib='blosc', complevel=9)

	print '!@# Final Merged table done in %s minutes'\
		%((time.time() - start) / 60)
	print '!!! Final Merged table saved to %s'\
			% (os.path.join(experiment_path, 'merged_before_score.pkl'))
else:
	del reaction_to_gene_top
	del compound_to_reaction
	del gene_to_reaction_top
	df = pd.read_pickle(args.merged_prescore, 'merged_before_score')
	print '\n!@# merged_before_score successfully loaded'

print '\n!@# Calculating final scores...'
start = time.time()
sys.stdout.flush()

# score reciprocal agreement
df = mg.reciprocal_agreement(df, closeness_threshold=args.reciprocal_closeness)

# calculate homology score
score = mg.homology_score(df)
# the nulls get a really low score
score[pd.isnull(score)] = 1
df['homology_score'] = score

# reaction connection score says if the compound got connected to any
# reaction in the database. Can't have zero because that messes up
# geometric mean, so added a small number.
df['reaction_connection'] = df[['reaction_id_r2g', 'reaction_id_g2r']]\
								.apply(pd.notnull).sum(axis=1) + 0.01

print '!@# Pre-scoring done in %s minutes' %((time.time() - start) / 60)
print '!@# Calculating final integrated MAGI score'
sys.stdout.flush()
start = time.time()

# calculate final MAGI integrated score
scoring_data = ['compound_score', 'reciprocal_score', \
				'homology_score', 'reaction_connection']
scores = []
to_score = df[scoring_data].values
data = []
for s in to_score:
    data.append(mg.magi_score(s, weights=args.final_weights))
scores.append(data)
df['MAGI_score'] = scores[0] / (args.chemnet_penalty ** df['level'].values)

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

df = df.merge(mg.mrs_reaction[['database_id']],
	left_on='reaction_id_r2g', right_index=True, how='left')
df = df.merge(mg.mrs_reaction[['database_id']],
	left_on='reaction_id_g2r', right_index=True, how='left',
	suffixes=('_r2g', '_g2r'))
cols = df.columns.values
idx = pd.np.argwhere(cols == 'query acc.')[0][0]
cols[idx] = 'gene_id'
df.columns = cols
df = df[['MAGI_score','gene_id', 'original_compound', 'neighbor',
	'note', 'compound_score','level','homology_score','reciprocal_score',
	'reaction_connection', 'e_score_r2g','database_id_r2g', 'e_score_g2r',
	'database_id_g2r']]

# save the full dataframe
df.to_csv(os.path.join(experiment_path, 'magi_results.csv'), index=False)
print 'full results saved to', os.path.join(experiment_path, 'magi_results.csv')
# save a compound-centric dataframe, where only the best row for each
# original_compound was chosen (this is only for compound scoring, do
# not use this for any kind of gene function analysis!)

compound_centric = df[pd.notnull(df['original_compound'])]\
					 .sort_values('MAGI_score', ascending=False)\
					 .drop_duplicates(['original_compound', 'compound_score'])
compound_centric = pd.merge(
	compound_centric, compounds,
	on=['original_compound', 'compound_score'],
	how='right')
compound_centric.to_csv(os.path.join(experiment_path, 
	'magi_compound_results.csv'), index=False)

gene_centric = df.sort_values(['MAGI_score', 'e_score_g2r'], 
	ascending=[False, False])\
	.drop_duplicates(['gene_id', 'database_id_g2r'])
gene_centric.to_csv(os.path.join(experiment_path,
	'magi_gene_results.csv'), index=False)

print '!@# MAGI Scoring done in %s minutes' %((time.time() - start) / 60)
print '\n!@# MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60)
print '!!! final results stored to %s' \
		%(os.path.join(experiment_path, 'magi_results.csv'))
