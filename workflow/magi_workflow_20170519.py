"""
Metabolites Annotations and Genes Integrated (MAGI)

MAGI 1.0a workflow script

"""

import argparse
import os
import warnings
import multiprocessing as mp

# parse arguments
parser = argparse.ArgumentParser()
# required
parser.add_argument('-f', '--fasta', 
	help='path to fasta file of genes in sample', 
	required=True)
parser.add_argument('-a', '--annotations', 
	help='path to annotation file for genes in sample', 
	default=None)
parser.add_argument('-c', '--compounds', 
	help='path to observed compounds file', 
	required=True)
# optional runtime variables
parser.add_argument('-n', '--cpu_count', 
	help='number of cpus to use for multiprocessing. Default is to use max!', 
	type=int, default=0)
parser.add_argument('-o', '--output', 
	help='path to a custom output', 
	type=str)
parser.add_argument('--test', 
	help='run MAGI only on the first # of pactolus compounds', 
	type=int)
parser.add_argument('-l', '--level', 
	help='how many levels deep to search the chemical network', 
	type=int, choices=[1,2,3], default=2)
parser.add_argument('-t', '--tautomer', 
	help='include tautomers in search; devault is true', 
	choices=[True, False], default=True)
parser.add_argument('--mute', 
	help='mutes pandas warnings', 
	action='store_true')
# jump-start the script after certain computations
# after the initial level 0 tautomer and blast search
parser.add_argument('--gene_to_reaction', 
	help='')
parser.add_argument('--compound_to_reaction', 
	help='')
parser.add_argument('--reaction_to_gene', 
	help='')

args = parser.parse_args()

# print your paths in stdout for logging
print '@@@ FASTA file input: %s' %(args.fasta)
print '@@@ compound input: %s' %(args.compounds)
if args.annotations is not None:
	print '@@@ annotations input: %s' %(args.annotations)

max_cpu = mp.cpu_count()

if args.cpu_count > max_cpu:
	raise RuntimeError('You have exceeded the cpus on this machine (%s)!' 
		%(max_cpu))
if args.cpu_count == 0:
	args.cpu_count = max_cpu

if args.gene_to_reaction is not None:
	if not os.path.isfile(args.gene_to_reaction):
		raise IOError('%s does not exist!' %(args.gene_to_reaction))
	else:
		print '@@@ gene_to_reaction input: %s' %(args.gene_to_reaction)

if args.compound_to_reaction is not None:
	if not os.path.isfile(args.compound_to_reaction):
		raise IOError('%s does not exist!' %(args.compound_to_reaction))
	else:
		print '@@@ compound_to_reaction input: %s' %(args.compound_to_reaction)

if args.reaction_to_gene is not None:
	if not os.path.isfile(args.reaction_to_gene):
		raise IOError('%s does not exist!' %(args.reaction_to_gene))
	else:
		print '@@@ reaction_to_gene input: %s' %(args.reaction_to_gene)

if args.mute:
	print '!!! Warnings are muted !!!'
	warnings.filterwarnings('ignore')

# import these after all the argument checking, because it takes so long
import sys
sys.path.insert(0, '/project/projectdirs/metatlas/projects/metatlas_reactions/workflow/helpertools')
import magitools2 as mg
import pandas as pd
import numpy as np
import time
import pickle
import datetime
print '\n'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#local_settings file stuff

# path to MAGI data storage
MAGI_PATH = '/global/project/projectdirs/openmsi/projects/temp_chem_net_data/MAGI_data'

# set up where the results will be stored
if args.output is None:
	experiment_name = args.fasta.split('/')[-1].split('.')[0]
	today = datetime.datetime.now()
	experiment_name += today.strftime('_%Y%m%d')
	experiment_path = '%s/%s' %(MAGI_PATH, experiment_name)
else:
	experiment_path = args.output

print '@@@ Saving all results here:', experiment_path
if not os.path.isdir(experiment_path):
	os.makedirs(experiment_path)

print '@@@ Using %s CPUs' %(args.cpu_count)

main_start = time.time() # overall program timer

#~~~~~~~~~~~~~~~~~~~~~~~~~

# load genome
print '\n*** LOADING GENOME ***'
genome, genome_db_path = mg.load_genome(args.fasta, MAGI_PATH, 
										annotation_file=args.annotations)

# load pactolus results
print '\n*** PACTOLUS RESULTS ***'
compounds = mg.load_dataframe(args.compounds)
compounds.fillna('', inplace=True)

if 'original_compound' not in compounds.columns:
	raise RuntimeError('Could not find "original_compound" as a column, please\
		rename the column corresponding to inchi keys for the compounds')
u_cpds = compounds['original_compound'].unique()
print len(u_cpds), 'total compounds to search'

if 'compound_score' not in compounds.columns:
	print 'WARNING: "compound_score" not found as a column; assuming that\
		there is no score for compounds, and setting the compound scores \
		to 1.0'
	compounds['compound_score'] = 1.0
else:
	compounds['compound_score'] = compounds['compound_score'].apply(float)

# Conduct gene to reaction search
print 'Conducting gene to reaction search'
start = time.time()
gene_blast = mg.multi_blast(genome.index, genome, mg.refseq_dbpath, 
	experiment_path, raise_blast_error=False, cpu=args.cpu_count)
print '... Homology searching done in %s minutes' %((time.time() - start) / 60)
gene_blast.to_csv(os.path.join(experiment_path, 'gene_blast.csv'))
print '... scored blast results saved to %s' \
		%(os.path.join(experiment_path, 'gene_blast.csv'))

gene_to_reaction = mg.refseq_to_reactions(gene_blast, 'subject acc.', 
	cpu=args.cpu_count)
gene_groups = gene_to_reaction.groupby('query acc.')
multidx = gene_groups['e_score'].apply(mg.keep_top_blast).index
idx = multidx.levels[1]
gene_to_reaction_top = gene_to_reaction.loc[idx]
print '... gene_to_reaction table completed in %s minutes' \
		%((time.time() - start) / 60)

gene_to_reaction_top.to_csv(os.path.join(experiment_path, 
										'gene_to_reaction.csv'))

# compound to reaction search
print 'Conducting compound to reaction search'
start = time.time()

def connect_compound_to_reaction_mp_helper(inchikey, 
										tautomer=args.tautomer, 
										neighbor_level=args.level):
    out = mg.connect_compound_to_reaction(inchikey, 
								    	tautomer=tautomer, 
								    	neighbor_level=neighbor_level)
    return out

input_compounds = compounds['original_compound'].unique()

p = mp.Pool(args.cpu_count)
out = p.map(connect_compound_to_reaction_mp_helper, input_compounds)
p.close()
p.terminate()

compound_to_reaction = pd.concat(out)
compound_to_reaction.reset_index(inplace=True, drop=True)

# connect the compound score
compounds['original_compound'] = compounds['original_compound'].apply(
									lambda x: '-'.join(x.split('-')[:2]))
compound_to_reaction = pd.merge(compounds, compound_to_reaction, 
								on='original_compound', how='inner')

compound_to_reaction.to_csv(os.path.join(experiment_path, 
										'compound_to_reaction.csv'))

print '... compound_to_reaction table done in %s minutes and saved to %s'\
		%((time.time()-start)/60, os.path.join(experiment_path, 
											'compound_to_reaction.csv'))

# reaction to gene search
print 'Conducting reaction to gene search'
start = time.time()

# set up a list of reference sequences to blast against the genome
reactions = compound_to_reaction[compound_to_reaction['reaction_id'] != '']\
					['reaction_id'].tolist()
reactions_refseqs = mg.mrs_reaction.loc[reactions, 'refseq_id']
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
print len(rseq_list), 'reference sequences to search'

reaction_to_gene_blast = mg.multi_blast(rseq_list, mg.refseq, genome_db_path, 
	experiment_path, cpu=args.cpu_count, raise_blast_error=False)

reaction_to_gene = mg.refseq_to_reactions(reaction_to_gene_blast, 'query acc.', 
	cpu=args.cpu_count)

reaction_groups = reaction_to_gene.groupby('query acc.')
multidx = reaction_groups['e_score'].apply(mg.keep_top_blast).index
idx = multidx.levels[1]
reaction_to_gene_top = reaction_to_gene.loc[idx]
reaction_to_gene_top.to_csv(os.path.join(experiment_path, 
										'reaction_to_gene.csv'))
print '... reaction_to_gene table done in %s minutes and saved to %s'\
		%((time.time()-start)/60, os.path.join(experiment_path, 
											'reaction_to_gene.csv'))

compound_to_gene = pd.merge(compound_to_reaction, reaction_to_gene_top, 
							on='reaction_id', how='left')

# Make an integrated dataframe, joining on the gene
df = pd.merge(compound_to_gene, gene_to_reaction_top, 
	left_on='subject acc.', right_on='query acc.', 
	suffixes=('_r2g', '_g2r'), how='outer')
df.reset_index(inplace=True, drop=True)

print 'Calculating final scores...'

# score reciprocal agreement
# find agreement
agree_idx = df[df['reaction_id_r2g'] == df['reaction_id_g2r']].index
df.loc[agree_idx, 'reciprocal_score'] = 2.

# disagreement
disagree = df[df['reaction_id_r2g'] != df['reaction_id_g2r']].index
slc = df.loc[disagree]

# close disagreements get a medium score
close = slc[['e_score_r2g', 'e_score_g2r']].min(axis=1) > (slc[['e_score_r2g', 'e_score_g2r']].max(axis=1)*0.75)
close_idx = slc.loc[close].index
df.loc[close_idx, 'reciprocal_score'] = 1

# very different disagreements get a low score
wrong_idx = df[pd.isnull(df['reciprocal_score'])].index
df.loc[wrong_idx, 'reciprocal_score'] = 0.01

# if one direction did not get a blast score, change reciprocal score to 0.1 - not wrong, but not close either.
incomparable_idx = df.loc[pd.isnull(df[['e_score_r2g', 'e_score_g2r']]).any(axis=1)].index
df.loc[incomparable_idx, 'reciprocal_score'] = 0.1

df.sort_values(['original_compound', 'level'])

# calculate homology score
score = mg.homology_score(df)
# the nulls get a really low score
score[pd.isnull(score)] = 1
df['homology_score'] = score

# adjust compound score by leveled search
df['level_adjusted_compound_score'] = df['compound_score'].values / \
											(10 ** df['level'].values)
# calculate final MAGI integrated score
scoring_data = ['level_adjusted_compound_score', 'reciprocal_score', 'homology_score']
scores = []
to_score = df[scoring_data].values
data = []
for s in to_score:
    data.append(mg.magi_score(s, weights=None))
scores.append(data)
df['MAGI_score'] = scores[0] / (10. ** df['level'].values)

# merge with annotation table
# find gene ids that are floats, convert those to strings, without the decimal
float_entries = df['subject acc._r2g'].apply(lambda x: isinstance(x, float))
df.loc[float_entries, 'subject acc._r2g'] = df.loc[float_entries, \
				'subject acc._r2g'].apply(lambda x: "{:.0f}".format(x))
# merge with genome table
df = df.merge(genome, left_on='subject acc._r2g', right_index=True, how='left')

# sort the final table
df.sort_values(['original_compound', 'MAGI_score'], 
				ascending=[True, False], inplace=True)

df.to_csv(os.path.join(experiment_path, 'magi_results.csv'))

print 'MAGI analysis complete in %s minutes' %((time.time() - main_start) / 60)
print 'final results stored to %s' \
		%(os.path.join(experiment_path, 'magi_results.csv'))