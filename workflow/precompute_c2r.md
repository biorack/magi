# Precomputing compound-to-reaction results
This should be performed whenever the metabolite database or the chemical similarity network are updated

This took a little longer than 6 hours on genepool for the first version of the metabolite database and network

## 1. prepare the input list
In a notebook:
```python
import pandas as pd

df = pd.read_pickle('<path>/magi/workflow/database/unique_compounds_groups_magi.pkl')
df = df[['inchi_key']]
df.columns = ['original_compound']
df.to_csv('<path>/all_compounds.csv')
```

## 2. Submit the job
The script command should look like this. the `--legacy` flag is critical here!
```bash
$ python <path>/magi_workflow_20170519.py \
$ --compounds <path>/all_compounds.csv \
$ --level 0 \
$ --output <path>/all_cpds
$ --legacy
```

## 3. Clean up the results file and save it
In a notebook:
```python
import pandas as pd
import pickle

c2rdf = pd.read_csv('<path>/all_cpds/magi_compound_results.csv')
c2rdf = c2rdf[['original_compound', 'note', 'reaction_id']]
c2rdf['index'] = c2rdf['original_compound'].apply(lambda x: '-'.join(x.split('-')[:2]))
c2rdf.set_index('index', inplace=True)
c2rdf.sort_index(inplace=True)

with open('<path>/magi/workflow/database/c2r.pkl', 'w') as f:
    pickle.dump(c2r, f)
```
