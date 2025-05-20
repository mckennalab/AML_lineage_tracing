#===========================================================================

### TREES 
# Generating trees from allele table from 1_LineageProcessing and clones assigned in 2_TagClustering 
# Save trees as newick strings 
# Plot trees with metadata 

#===========================================================================


import pandas as pd
import numpy as np 

from AML_lineage_utilities import *

import seaborn as sns
import matplotlib.pyplot as plt
import re

from glob import glob
import os
import subprocess

import cassiopeia as cas

#===========================================================================

# 1. Make rmWT allele table (remove cells with no editing)

### Input: Allele table 
### Output: Allele table w/ cells without editing removed. 

res_at_clean = res_at[['cellBC', 'exp', 'intBC', 'allele', 'r1', 'r2', 'r3','r4', 'r5', 'r6', 'r7', 'r8', 'UMI', 'readCount', 'founder', 'GEX_bc', 'FB_bc','status', 'cas_clone']]

edits_melted = pd.melt(res_at_clean[ ['cellBC'] + ['r' + str(i) for i in range(1, 9)]], id_vars = 'cellBC')[['cellBC', 'value']].drop_duplicates()
ncells = len(set(edits_melted.cellBC))
print(ncells)

# one line per unique edit + cell combo -- so any cells that are removed when NONE are filtered out means thats the only edit the cell has
ncells_rm_NONE = len(set(edits_melted.loc[edits_melted['value'] != 'NONE'].cellBC))
print(ncells_rm_NONE)

print(f'{ncells - ncells_rm_NONE} cells have NO editing and are exluded from the rmWT newick tree.')

mask = edits_melted.loc[edits_melted['value'] != 'NONE'].cellBC

res_at_rmWT = res_at[res_at['cellBC'].isin(mask)]
res_at_rmWT.to_csv(f'{adir}/241211_prefinal_cas_founder_combined_annotated_allele_table_rmWT.csv')


rs_at_rmWT_clean = res_at_rmWT[['cellBC', 'exp', 'intBC', 'allele', 'r1', 'r2', 'r3','r4', 'r5', 'r6', 'r7', 'r8',
                               'UMI', 'readCount', 'founder', 'GEX_bc','FB_bc','status', 'cas_clone']]
rs_at_rmWT_clean.to_csv(f'{direct}/final_tables_res1/hl60s_start_end_combined_lineage_castags_res1_allele_table_rmWT.csv')

#===========================================================================

# 2. Build trees using the cassiopeia vanilla greedy solver. 

### Input: Allele table 
### Output: Newick files per clone. 

input_at = pd.read_csv(f'{direct}/final_tables_res1/hl60s_start_end_combined_lineage_castags_res1_allele_table_rmWT.csv')
tag_key = 'cas_clone'
indel_priors = cas.pp.compute_empirical_indel_priors(input_at, grouping_variables=['intBC', tag_key])
solver = 'Vanilla'

indat = list(np.unique(input_at[[tag_key]]))
pbar = tqdm(total = len(indat))

tree_outdir = f'{direct}/final_trees_res1'

ensure_directory_exists(tree_outdir)
                       
for cl in indat:
    if input_at.loc[input_at[tag_key] == cl].cellBC.nunique() == 1:
        print(f'Only one cell found. Skipping {cl}')
        continue
    cas_tree = make_vanilla_tree(input_at, cl, 
                                            tag_colid=tag_key, 
                                            rep_thresh=0.95, solve=False)#, indel_priors=indel_priors)
    if cas_tree is None:
        continue
    res = cassiopeia_solve_tree(cas_tree, solver=solver)
    write_newick(f'{tree_outdir}/with_diversity_cutoff/hl60s_start_end_combined_{cl}_{tag_key}.{solver}_tree_rmWT.newick', res.get_newick())

    pbar.update(1)

pbar.close()

#===========================================================================

# 3. Combine into a single megatree 

### Input: List of paths to newick strings to combine into a single tree. 
### Output: Single newick MEGATREE file. 

import re
import os
import networkx as nx 

founders = input_at.founder.unique()
tree_outdir = f'{direct}/final_trees_res1/with_diversity_cutoff'
prefix = 'hl60s_start_end_combined'

pbar = tqdm(total = len(founders))

for founder in founders:
    pattern = fr".*{founder}.\d+_cas_clone.Vanilla_tree_rmWT.newick$"
    newick_trees = [f for f in os.listdir(f'{tree_outdir}')  if re.match(pattern, f)]
    
    if len(newick_trees) == 0:
        print(f'No clonal trees found for founder {founder}. Exiting.')
        continue 
        
    tree_top = []
    included = []
    for i in range(0, len(newick_trees)):
        with open(f'{tree_outdir}/{newick_trees[i]}', 'r') as f:
            tree = f.read()

        res = cas.data.CassiopeiaTree(tree=tree)
        if res is None:
            continue

        clone_id = re.findall(fr'{founder}.\d+', newick_trees[i])[0]
        map_names = dict(zip(res.internal_nodes, [f'{clone_id}' + i for i in res.internal_nodes]))
        G = res.get_tree_topology()
        G = nx.relabel_nodes(G, map_names)

        included.append(clone_id)
        tree_top.append(G)

    tree_tops = join_digraphs(tree_top, founder, included)

    combined_tree = cas.data.CassiopeiaTree(
                            tree=tree_tops, 
                            root_sample_name=founder)

    print(f'{len(included)} included in final megatree for founder {founder}.')

    nwk = combined_tree.get_newick()
    write_newick( f'{tree_outdir}/{prefix}_{founder}_{solver}_rmWT_megatree.newick', nwk)
    pbar.update(1)
pbar.close()

#===========================================================================

# 4. Plotting 

### Input: Tree or list of trees with associated metadata. 
### Output: Circular tree plots. 

## Load metadata 
md = pd.read_csv(f'{direct}/final_tables_res1/hl60s_start_end_combined_cell_table_w_founder_and_casclone_annotation.csv')
md.set_index('cellBC', inplace=True)

### General Colors 
group_colors = ['#4C8BDB', '#FFCE00', '#292930', '#9B9CA0'] ### I think we changed resistant to firebrick3 ... 
exp_groups = ['Control', 'AraC', 'rs29', 'rs30'] 
hex_dict_group = dict(zip(list(range(0, len(exp_groups))), group_colors))
cmap=create_custom_cmap(hex_dict_group) 
value_mapping = dict(zip(exp_groups, hex_dict_group.keys()))
unique_values = exp_groups

## Load tree 
with open(f'{tree_outdir}/{nwk}', 'r') as f:
        tree = f.read()
    
res = cas.data.CassiopeiaTree(tree=tree)
if res is None:
    continue
    
res.cell_meta = md.loc[res.leaves, ]

# Create color map and legend 
cas_cmap=create_custom_cmap(hex_dict_cas_group) 
value_mapping = dict(zip(unique_values, hex_dict_cas_group.keys()))
color_dict = dict(zip(unique_values, relevant_colors))
handles = []
for label, color in color_dict.items():
    # Create a custom legend handle (filled circle)
    handle = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10)
    handles.append(handle)
        
fig, ax = cas.pl.plot_matplotlib(res, meta_data=['cas_clone', 'founder', 'sample'] , 
                                    internal_node_kwargs={'s':0}, leaf_kwargs={'s':0.1},  
                                    branch_kwargs={'linewidth':0.5, 'c':'black'},
                                    value_mapping=value_mapping,
                                    categorical_cmap=cas_cmap, 
                                     ontinuous_cmap='viridis')
ax = plt.gca()
ax.set_aspect(1.0)
ax.legend(handles=handles, labels=unique_values, loc='upper center',
            bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=4, fontsize=8)
    
plt.title(f'AML HL60s | {clone_id}')
plt.savefig(f'{tree_outdir}/{nwk[:-7]}_by_sample.pdf')
plt.close()