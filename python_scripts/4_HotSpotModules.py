from datetime import datetime
current_date = datetime.now()
savedate = current_date.strftime("%y%m%d")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D
import scanpy as sc

import re
from statsmodels.stats.multitest import multipletests
from cnmf import cNMF
import random
from scipy import stats
import os
from itertools import combinations 
# clean up 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import hotspot 
from ete3 import Tree

from tqdm import tqdm

#===========================================================================

# Load Adata 
adata = sc.read_h5ad(counts_file)
adata.var['geneMask'] = ~adata.var_names.str.startswith(("MT-", "RPS", "RPL"))
keep_genes = adata.var['geneMask']  # Extract boolean mask
adata_subset = adata[:, keep_genes]

### iterate through trees of interest 
root = '/dartfs/rc/lab/M/McKennaLab/projects/hannah/aml/analysis/human_cell_lines/hl60/lineage/final_trees_res1/with_diversity_cutoff/megatrees_newicks'
trees = [os.path.join(root, name) for name in  os.listdir(root)]

rds_label = 'sfc_hl60_full_joinged_121324_rmMTRB'

for tree_path in trees:
    tree = Tree(tree_path, format=1)
    leaves = set()
    for tn in tree.traverse('postorder'):
        if tn.is_leaf():
            leaves.add(tn.name)
    #treeid = re.findall('F\\d+.\\d+|F\\d+_.*megatree', os.path.basename(tree_path))
    treeid = re.findall('H\d+|H_.*META', os.path.basename(tree_path))
    if len(treeid) == 0:
        print(f'Skipping {tree_path}.')
        continue
    else:
        treeid = treeid[0]
    
    if 'META' in treeid:
        treeid = treeid.replace('_Vanilla_rmWT', '')

    print(f'tree {treeid} has {len(leaves)} leaves.')

    ### Get indices in adata to match the newick string. 
    is_valid = [x in leaves for x in adata.obs_names]
    is_valid_indices = np.nonzero(is_valid)[0]
    valid_barcodes = [adata.obs_names[i] for i in is_valid_indices]
    ldata = adata[valid_barcodes]
    sc.pp.filter_genes(ldata, min_cells=10)

    hs = hotspot.Hotspot(ldata, model='normal', tree=tree)
    try:
        hs.create_knn_graph(weighted_graph=False, n_neighbors=10,)
        hs_results = hs.compute_autocorrelations(jobs=1)

    except:
        print(f'An error occured for {treeid}. Moving on') 
        continue 
    
    hs_genes = hs_results.index[hs_results.FDR < 0.05]

    if hs_genes.shape[0] <= 1:
        print(f'No genes found for {treeid}. Moving on.')
        continue
        
    lcz = hs.compute_local_correlations(hs_genes, jobs=1)

    ### 2. Call modules 
    modules = hs.create_modules(min_gene_threshold=10, core_only=True, fdr_threshold=0.05)
    modules.value_counts()

    if len(modules.value_counts()) == 1:
        print(f'No modules found for {treeid}. Moving on.')
        continue 
    else:
        print(f'Found {len(modules.value_counts()) - 1} for {treeid}.')

    hs.plot_local_correlations()
    plt.savefig(f'{outdir}/{rds_label}_{treeid}_hotspot_local_correlations_mingenethresh_10.pdf')

    results = hs.results.join(hs.modules)
    results.to_csv(f'{outdir}/{rds_label}_{treeid}_hotspot_local_correlations_mingenethresh_10_results.csv')

    module_scores = hs.calculate_module_scores()
    ldata.raw = ldata
    sc.pp.scale(ldata)
    sc.tl.pca(ldata)
    sc.pp.neighbors(ldata)
    sc.tl.umap(ldata)

    module_cols = []
    for c in module_scores.columns:
        key = f"Module {c}"
        ldata.obs[key] = module_scores[c]
        module_cols.append(key)

    fig = sc.pl.umap(ldata, color=module_cols + ['group'], frameon=False, vmin=-1, vmax=1, return_fig = True,  ncols = 2)
    plt.savefig(f'{outdir}/{rds_label}_{treeid}_hotspot_local_correlations_mingenethresh_10_UMAP.pdf')