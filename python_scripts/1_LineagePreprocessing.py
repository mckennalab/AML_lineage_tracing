#===========================================================================

### LINEAGE PREPROCESSING 
# Take output from SingleCellLineage pipeline and do QC with cassiopeia preprocessing pipeline (for lineage) and custom QC for cas-tags. 

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

# Lineage Processing 
## Part 1 - Preprocessing 

### Input: .stats file from SingleCellLineage 
### Output: Allele table with unique cellBC+intBC+indel allele pairings. 

# 1. Initial Cleaning 

# Read in .stats file from SingleCellLineage and extract read info and founder tags. 
barcodes = pd.read_table(barcode_path)
barcodes = extractReadInfo(barcodes)
barcodes['founderTag'] = barcodes.apply(lambda row: getBC(row['mergedReadRef'], row['mergedRead'], r'NNNNNNNNNNNNNN'), axis=1) 

# Remove tags that don't PASS SingleCellLinege filter 
# Remove cells that are not in the 10x whitelist 
barcodes = barcodes.loc[barcodes['keep'] == 'PASS']
barcodes = barcodes.merge(whitelist, left_on='cellBC', right_on='FB_bc')

### OPTIONAL - can also filter against a known cellBC list here. 

# Reformatting, no cells lost here 
barcodes = barcodes[[ 'readName', 'readCount', 'cellBC', 'UMI', 'founderTag',  'GEX_bc', 'FB_bc'] + target_cols]
target_cols = ['target' + str(i) for i in range(1, 9)]

barcodes['allele'] = barcodes[target_cols].agg('_'.join, axis=1)
barcodes['bad_tag'] = barcodes['founderTag'].apply(lambda x: contains_pattern(r'-+', x) | contains_pattern(r'([ACTG])\1{4,}', x))
barcodes = barcodes.loc[barcodes['bad_tag'] == False]

#2. Cassiopeia: Resolve UMI Sequences 
barcodes_step2 = cas.pp.resolve_umi_sequence(
    barcodes, output_directory='./', min_umi_per_cell=2, min_avg_reads_per_umi=2.0, plot=False,
)

# 3. Correct Founder Tags / intBC
barcodes_step3 = starcode_correct(barcodes_step2)

# Reformat column names 
barcodes_step3 = barcodes_step3.rename(columns={'founderTag':'intBC'})
barcodes_step3 = barcodes_step3.rename(columns=dict(zip(['target'+str(i) for i in range(1,9)], ['r'+str(i) for i in range(1,9)])))

# 4. Cassiopeia: Error correct UMIs and filter molecule table 

barcodes_step4 = cas.pp.error_correct_umis(barcodes_step3, max_umi_distance=1)
barcodes_step5 = cas.pp.filter_molecule_table(barcodes_step4, 
                                            output_directory='./', 
                                            min_umi_per_cell = 2,
                                            min_avg_reads_per_umi = 2.0,
                                            min_reads_per_umi = 2, #1, # 99th percentile leaves us with 3. 
                                            intbc_dist_thresh = 0, 
                                            doublet_threshold = None,
                                            allow_allele_conflicts=False)
# Final filter 
barcodes_step6 = utilities.filter_cells(barcodes_step5,min_umi_per_cell=int(3),min_avg_reads_per_umi=2)

#5. Group so one row per cellBC+intBC+allele 

grouping = ["cellBC", "intBC",  "allele"] + ['r'+str(i) for i in range(1,9)]
barcodes_pre_final_df = barcodes_step6.groupby(grouping, as_index=False).agg(
        {"UMI": "count", "readCount": "sum"}
    )

#6. OPTIONAL: Filter intBCs on known Founder Tag whitelist 
refclone_ids = pd.read_csv(path_to_founder_whitelist)
barcodes_final_df = barcodes_pre_final_df.merge(refclone_ids, left_on='intBC', right_on='ID').drop_duplicates()

# In absense of founder whitelist, see 2_TagClustering for additional founder-tag (intBC) clustering to identify founders. 

#===========================================================================

## Part 2 - Post-processing 

### Input: Allele table with unique cellBC+intBC+indel allele pairings. 
### Output: Cleaned allele table with unique cellBC+intBC+indel with final founder assignments.

# Make sure each all reads in a cell map uniquely to one founder population 

barcodes_final_df = pd.read_csv(f'{dir_path}/{prefix}_barcodes_final.csv')

ft = (barcodes_final_df.groupby(['cellBC','clone'])
          .agg({'readCount':'sum', 'UMI':'sum'}))
ft.columns = ['clonal_nReads', 'clonal_nUMI']
ft = ft.reset_index()

if len(ft.cellBC) == len(set(ft.cellBC)):
    print('All remaining cells are unique.')
else:
    print('Non-unique cell-founder pairs remain.')

ft['propUMI'] = ft['clonal_nUMI'] / ft.groupby(['cellBC'])['clonal_nUMI'].transform(sum)
ft['propReads'] = ft['clonal_nReads'] / ft.groupby(['cellBC'])['clonal_nReads'].transform(sum)
sns.histplot(ft, x='propReads')#, hue='clone', multiple='dodge')
plt.axvline(0.95, color='red')

ft['pass_mask'] = ft.apply(lambda x: 'KEEP' if (x['propReads'] >= 0.95) & (x['propUMI'] >= 0.9) else 'FAIL', axis=1)

ft_thresh = ft.loc[(ft['propReads'] >= 0.95) & (ft['propUMI'] >= 0.9)]
print(f'{len(set(ft_thresh.cellBC))} cells survive.')
if len(ft_thresh.cellBC) == len(set(ft_thresh.cellBC)):
    print('All cell-founder pairs are unique.')

bc_final_resolved = barcodes_final_df.merge(ft_thresh[['cellBC', 'clone']], on=['cellBC', 'clone'])

#===========================================================================

# Clone Tag Processing 

### Input: .stats file from SingleCellLineage 
### Output: Cleaned dataframe with cells and cas-tags. 

# 1. Initial Cleaning 

# Read in .stats file from SingleCellLineage and extract read info and clone tags. 
ct_raw = pd.read_table(castag_path)
ct_raw = extractReadInfo(ct_raw)

# Extract raw clone tags 
ct_regex = r'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' # r'NNNNNNNNNNNNNN'
ct_raw['castag_raw'] = ct_raw.apply(lambda row: getBC(row['mergedReadRef'], row['mergedRead'], ct_regex), axis=1) 
ct_raw_cln_cols = ct_raw[['readName','castag_raw','keep', 'matchRate1', 'mergedRead', 'mergedReadRef', 'seq_info', 'readCount', 'cellBC', 'UMI']]

# Filter against 10x whitelist 
ct_raw_cln_cols = ct_raw_cln_cols.loc[ct_raw_cln_cols['cellBC'].isin(whitelist)]

# 2. OPTIONAL: Filter for tags against a pattern
pattern = r'[ACGT-]{3}[AT-][CG-][ACGT-]{3}[AT-][CG-][ACGT-]{3}[AT-][CG-][ACGT-]{3}[AT-][CG-][ACGT-]{3}[AT-][CG-][ACGT-]{3}[AT-][CG-][ACGT-]{3}[AT-][CG-][ACTG-]{3}'

ct_raw_cln_cols['castag_raw'] = ct_raw_cln_cols.castag_raw.apply(lambda x: find_custom(rs29_pattern, x))
ct_cln_final = ct_raw_cln_cols.loc[ct_raw_cln_cols['castag_raw'] != 'FAIL']

# 3. Error correct clone tags with Starcode 
ct_cln_final = ct_cln_final[['readName', 'readCount','exp','cellBC', 'UMI', 'castag_raw']]
ct_cln_final = ct_cln_final.loc[c1498_raw_cmb['castag_raw'].str.count('-') <= 1]  # remove castags with more than one dash 
ct_cln_final['castag_raw'] = ct_cln_final.castag_raw.str.replace('-', '')

ct_cln_corre = starcode_correct(ct_cln_final, tag_id='castag_raw', do_count=True)
ct_cln_corre = ct_cln_corre.rename(columns = {'assignedTag':'CasTag'})

# 4. Aggregate counts per cellBC-CloneTag and filter 

grp_ct_cln_corre = ct_cln_corre.groupby(['cellBC','exp', 'CasTag']).agg({'UMI':'count', 'readCount':'sum'})
grp_ct_cln_corre.columns = ['nUMI', 'aggRC'] 
grp_ct_cln_corre = grp_ct_cln_corre.reset_index()

run_pattern = r'([ACTG])\1{4,}' # this will get rid of nucleotide runs longer than 5 base pairs.  
ct_final = filter_castags(grp_ct_cln_corre, run_pattern = run_pattern)


