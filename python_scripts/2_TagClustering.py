
#===========================================================================

### TAG CLUSTERING 
## Cluster founder tags with additional filtering to identify founders in HL60 dataset
## Cluster cas-tags cleaned from 1_LineagePreprocessing 
## Clean cas-tag / founders --> cas-tags should uniquely map to a single founder. 

## Final dataset generation. 

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

#===========================================================================

# Assign founders to cells in absense of founder whitelist 

### Input: Allele table from 1_LineagePreprocessing.py
### Output: Dataframe with cells assigned to founders based on their founder-tags (intBCs). 

# 1. Clean founder tags 
lin_combined = barcodes_pre_final_df ### final allele table from 1_LineagePreprocessing.py
lin_combined['propUMI'] = lin_combined['UMI'] / lin_combined.groupby('cellBC').UMI.transform(sum)
lin_combined['propreadCount'] = lin_combined['readCount'] / lin_combined.groupby('cellBC').readCount.transform(sum)

# Plot histogram of proportion reads to decide background cutoff 
sns.histplot(lin_combined, x = 'propreadCount', log_scale=True)
plt.axvline(x=0.015, c='red')

filt_by_reads = bc_final_resolved.loc[bc_final_resolved['propreadCount'] > 0.015]

# aggregate info per intBC
agglin_ = filt_by_reads.groupby(['intBC']).agg({'UMI':'sum', 'readCount':'sum', 'cellBC':'count'}).reset_index()
agglin_['rc_per_cell'] = agglin_['readCount'] / agglin_['cellBC']
agglin_['umi_per_cell'] = agglin_['UMI'] / agglin_['cellBC']

# filter by number of cells per intBC and number of UMIs per cell 
# to get robust set of intBCs to identify founder groups 
ncells_cutoff, umi_cutoff = 5, 3
agglin_ft = agglin_.loc[(agglin_['cellBC'] > ncells_cutoff) & (agglin_['umi_per_cell'] > umi_cutoff)]

filt_by_reads_sf = filt_by_reads.loc[filt_by_reads['intBC'].isin(set(agglin_ft['intBC']))]
intBC_to_founders = filt_by_reads_sf[['cellBC','experiment', 'intBC', 'readCount']].drop_duplicates()

# 2. Louvain clustering 
res = do_clusters(intBC_to_Founders, tag='intBC', res=1, algorithm='louvain')

# 3. Visualize cluster assignments -- by proportion of cells within a founder that have a tag. 
expanded_founders = expand_dataframe(res, 'intBC')
expanded_founders['nCells'] = expanded_founders.groupby('louvain').cellBC.transform('nunique') # number of cells in the founder 
expanded_founders_count = expanded_founders.groupby(['intBC', 'nCells', 'louvain']).cellBC.count().reset_index() # number of cells with the tag
expanded_founders_count['normCount'] = expanded_founders_count['cellBC'] / expanded_founders_count['nCells']          
expanded_piv = expanded_founders_count.pivot(columns = 'intBC', index = 'louvain', values = 'normCount').fillna(0)
sns.clustermap(expanded_piv, cmap='rocket_r', figsize=(6, 6))

sns.histplot(expanded_founders_count, x = 'normCount', log_scale=True)
plt.axvline(x = 0.1, c='red')

# 4. Discretize founder tags to founders by filtering out background tags that are found low numbers of cells within a founder. 
expanded_founders['total_cells_in_founder'] = expanded_founders.groupby([ 'louvain']).cellBC.transform('nunique')
counts_by_founder = expanded_founders.groupby([ 'louvain','total_cells_in_founder', 'intBC']).cellBC.nunique().reset_index()
counts_by_founder['all_cells_with_ft'] = counts_by_founder.groupby('intBC').cellBC.transform(sum)
counts_by_founder['norm_count_1'] = counts_by_founder['cellBC'] / counts_by_founder['all_cells_with_ft']

### only keep cellBC+founder-tag combos where the number of cells with the tag in a founder is greater than 50% of all cells with the tag. 
### cutoff decided based on distribution of norm_count_1. will be higher with fewer founders. 

cutoff = 0.5 
counts_by_founder_f = counts_by_founder.loc[counts_by_founder["norm_count_1"] > cutoff]
print(f'{counts_by_founder.intBC.nunique() - counts_by_founder_f.intBC.nunique()} ambiguous intBCs were removed.')
print(counts_by_founder_f.intBC.nunique() == counts_by_founder_f.intBC.shape[0])

sns.histplot(counts_by_founder, x = 'norm_count_1', log_scale=True, bins=50)
plt.axvline(x=cutoff, c='red')

# 5. Save result. Cols:  cellBC, intBC, founder assignment, number of cells with tag in founder, total number of cells in founder
res_clean = expanded_founders.merge(counts_by_founder_f[['intBC', 'louvain']], on=['louvain', 'intBC'])
res_clean.to_csv(f'{adir}/241211_hl60s_intBC_to_founders_clustering_res1_CLEAN.csv')

#===========================================================================

# Final Cleaning - make sure cas-tags map to cells belonging to a single founder. 

### Input: qc allele table and qc cas-tags from 1_LineagePreprocessing.py. 
### Output: Unique mapping of cas-tags per founder and cleaned allele table. 

# 1. Filter against lineage table and resolve cas-tags to founders. 
bc_final_resolved = pd.read_csv(path_to_allele_table) # result of 1_LineagePreprocessing.py
ct_final = pd.read_csv(path_to_qc_castags) # result of 1_LineagePreprocessing.py

all_lin_combined = ct_final.merge(bc_final_resolved, on = ['lin_cellID', 'exp'])
all_lin_combined_cln = all_lin_combined[['lin_cellID', 'exp', 'CasTag_x', 'nUMI', 'aggRC', 
        'intBC', 'allele', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'UMI', 'readCount', 'founder', 'GEX_bc',
        'FB_bc']]

all_lin_combined_cln.rename(columns = {'lin_cellID':'cellBC', 'CasTag_x':'CasTag'}, inplace=True)

agg_by_founder = all_lin_combined_cln.groupby(['CasTag', 'founder', 'exp']).agg({'nUMI':'sum', 'aggRC':'sum'}).reset_index()

cutoff = 0.5 # adjust cutoff depending on visualization, expected number of founders, etc. 

# keep cas-tag + founder-tag pairings with highest identity 
agg_by_founder['propUMI'] = agg_by_founder['nUMI'] / agg_by_founder.groupby(['CasTag', 'exp'])['nUMI'].transform(sum)
agg_by_founder['propReads'] =agg_by_founder['aggRC'] / agg_by_founder.groupby(['CasTag', 'exp'])['aggRC'].transform(sum)

plt.hist(agg_by_founder['propUMI'])
plt.axvline(x=cutoff, color='r', linestyle='--', linewidth=2)


cas_x_founder = agg_by_founder.loc[agg_by_founder['propUMI'] > cutoff, ['founder', 'CasTag']] 
clean = all_lin_combined_cln.merge(cas_x_founder, on=['founder', 'CasTag']) # merge with cleaned allele table to get final cleaned set of tags and founders. 

clean[['cellBC', 'CasTag']].to_csv(f'{adir}/241211_hl60s_intBC_to_founders_clustering_res1_CLEAN_castags_and_cells.csv')

#===========================================================================

# 1. Clustering cas-tags 

### Input: cleaned cas-tags from 1_LineagePreprocessing.py that have been discretized to founders and filtered. 
### Output: cells' louvain clusters and sets of cas-tags. 

clean_castags = pd.read_csv(f'{adir}/241211_hl60s_intBC_to_founders_clustering_res1_CLEAN_castags_and_cells.csv', index_col=0)
ct_clust = do_clusters(clean_castags, res=1, algorithm='louvain')

# 2. OPTIONAL cleaning of cas-tag sets to remove loosely assigned cells (i.e. with one tag that is lowly represented in the cluster). 

ct_expanded = expand_dataframe(ct_clust, 'CasTag') # removed singletons (unclustered cells)
ct_expanded['nCells'] = ct_expanded.groupby('louvain').cellBC.transform('nunique')
ct_expanded_count = ct_expanded.groupby(['CasTag', 'nCells', 'louvain']).cellBC.count().reset_index()
ct_expanded_count['normCount'] = ct_expanded_count['cellBC'] / ct_expanded_count['nCells']

ct_expanded_count['total_cells_in_clone'] = ct_expanded_count.groupby([ 'louvain']).cellBC.transform('nunique')
counts_by_clones = ct_expanded_count.groupby([ 'louvain','total_cells_in_clone', 'CasTag']).cellBC.nunique().reset_index()
counts_by_clones['all_cells_with_castag'] = counts_by_clones.groupby('CasTag').cellBC.transform(sum)

counts_by_clones['norm_count_1'] = counts_by_clones['cellBC'] / counts_by_clones['all_cells_with_castag']
counts_by_clones['norm_count_2'] = counts_by_clones['cellBC'] / counts_by_clones['total_cells_in_clone']

cutoff = 0.1 # removes lowly represented castags within clones. 
sns.histplot(counts_by_clones, x = 'norm_count_1', log_scale=True, bins=50)
plt.axvline(x=cutoff, c='red')

counts_by_clones_f = counts_by_clones.loc[counts_by_clones['norm_count_1'] > cutoff]

print(f'{counts_by_clones.CasTag.nunique() - counts_by_clones_f.CasTag.nunique()} ambiguous CasTags were removed.')
print(counts_by_clones_f.CasTag.nunique() == counts_by_clones_f.CasTag.shape[0])

cas_clusters_cleaned = counts_by_clones_f[['louvain', 'CasTag']].merge(ct_expanded) 
cas_clusters_cleaned.to_csv(f'{adir}/241211_hl60s_final_cas_clone_assignments_cleaned.csv')

#===========================================================================

# Final dataset included cellBC + intBC + allele + cas-clone + founder assignments 

### Input: Final cleaned cas-tags with assigned clones and final cleaned founder tags with assigned founder groups. 
### Output: Final annotated allele table with all lineage information. 

# check that the cas-clones are discrete to one founder 
cas_clusters_cleaned = cas_clusters_cleaned.rename(columns={'louvain' : 'cas_clone.1', 'nCells':'cells_in_clone'})
res_clean = res_clean.rename(columns = {'louvain':'founder', 'nCells':'cells_in_founder'})
compare_groups = cas_clusters_cleaned.merge(res_clean, on='cellBC')
compare_groups.groupby('cas_clone.1').founder.nunique().max()

# merge cleaned founder-tags with original allele table and rename founders based on size ranking. 
lin_anno  = res_clean[['cellBC', 'intBC', 'louvain_founder']].merge(lin_combined, on = ['cellBC', 'intBC'])
cell_counts(lin_anno)

lin_grouped = lin_anno.groupby(['louvain_founder']).cellBC.nunique().reset_index()

lin_grouped['ranking'] = lin_grouped['cellBC'].rank(method='first', ascending=False).astype(int)
lin_anno = lin_anno.merge(lin_grouped[['louvain_founder', 'ranking']], on='louvain_founder')
lin_anno['founder'] = ['H' + str(i) for i in lin_anno.ranking]

# merge with lineage on cellBCs and similarly assign cas-clone names based on size ranking within a founder. 
cl_all_ = cas_clusters_cleaned.merge(lin_anno[['cellBC', 'founder']].drop_duplicates())
cl_grouped = cl_all_.groupby(['founder','louvain_cas']).cellBC.nunique().reset_index()

cl_grouped['ranking'] = cl_grouped.groupby('founder')['cellBC'].rank(method='first', ascending=False).astype(int)
cl_all_ = cl_all_.merge(cl_grouped[['louvain_cas', 'ranking']], on='louvain_cas')

# merge to generate final dataset with alleles, founders, and cas-clones
cl_combined_anno = lin_anno[['cellBC', 'exp', 'intBC', 'allele', 'r1',
       'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'UMI', 'readCount', 'founder']].merge(cl_all_[['cellBC', 'louvain_cas', 'ranking']], on='cellBC', how = 'left')

ct_lin = cl_combined_anno.drop_duplicates()
ct_lin['ranking'] = ct_lin.ranking.fillna(0)
ct_lin['ranking'] = ct_lin.ranking.astype(int)
ct_lin['ranking'] = ct_lin.ranking.replace(0, 'unassigned')
ct_lin['cas_clone'] = ct_lin.apply(lambda row: f"{row['founder']}.{str(row['ranking'])}", axis=1)

ct_lin.to_csv(f'{adir}/241211_prefinal_cas_founder_combined_annotated_allele_table.csv')
#===========================================================================


