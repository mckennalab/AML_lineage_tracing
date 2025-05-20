### Utils functions for AML_Lineage [Saxe et al. in Prep]

import pandas as pd
import numpy as np
import re
import cassiopeia as cas

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import silhouette_score
from tqdm import tqdm 
import datetime
import networkx as nx 

### NOTE: cas-tags == clone-tags 


### Lineage Processing Functions 
def contains_pattern(pattern, string):
    return re.search(pattern, string) is not None

def find_custom(pattern, value):
    res = re.findall(pattern, value)
    if len(res) == 0:
        return 'FAIL'
    else:
        return res[0]

def getBC(x1, x2, search):
    res = re.search(search, x1)
    if res is not None:
        return x2[res.start():res.end()]
    else:
        return '--------------'
    
def extractReadInfo(data, split_cols=['seq_info', 'R1', 'readCount', 'totalReads', 'x2', 'x3', 'x4'], method='10X'):
    split_values = data['readName'].str.split('_', expand=True)
    data = data.reindex(columns=[*data.columns, *split_cols])
    data[split_cols] = split_values

    if method == '10X': 
        data['cellBC'] = data.R1.str[0:16]
        data['UMI'] = data.R1.str[16:28]
        data['readName'] = data['cellBC'] + '_' + data['UMI'] + '_' + data['readCount']
    else:
        data['readName'] =data['R1'].str[:-3] + '_' + data['readCount']
    data['readCount'] = data['readCount'].astype('int')
    
    return data


def make_vanilla_tree(allele_table, group="AllCells", metadata=None, 
                      indel_priors=None, rep_thresh=1.0, tag_colid='assigned_tag', solve=True):
    if group != "AllCells": 
        at_ = allele_table.loc[allele_table[tag_colid] == group]
    else:
        at_ = allele_table
        
    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(at_,  allele_rep_thresh=rep_thresh, 
                                                                                             missing_data_allele='UNKNOWN', 
                                                                                             mutation_priors = indel_priors)
    if character_matrix.shape[0] == 0:
        print('Exiting... No target sites remain after filtering. ')
        return 
    
    cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)
    
    cell_meta = at_.groupby('cellBC').agg({"intBC": 'nunique', 'UMI': 'sum'})

    missing_proportion = (character_matrix == -1).sum(axis=0) / character_matrix.shape[0]
    uncut_proportion = (character_matrix == 0).sum(axis=0) / character_matrix.shape[0]
    n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis=0)

    character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 'uncut_prop', 'n_unique_states']).T

    cas_tree.cell_meta = cell_meta
    cas_tree.character_meta = character_meta
    
    if not solve:
        return cas_tree
    
    vanilla_greedy = cas.solver.VanillaGreedySolver()
    vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)
    
    if metadata is not None:
        cas_tree.cell_meta = cas_tree.cell_meta.merge(metadata, left_index=True, right_index=True)
        
    return cas_tree

def write_newick(path, newick):
    with open(path, 'w') as file:
        file.write(newick)
        
    file.close()
    return

def filter_castags(ct, 
                   rc_thresh = 1, 
                   umi_thresh = 2,
                   prop_thresh = 0.1,
                   run_pattern = r'([ACTG])\1{7,}'):
    '''
    CasTag file with required columns nUMI and aggRC where nUMI is the number of UMIs supporting a castag-cellBC pair 
    and aggRC is the total number of reads supporting that castag-cellBC pair (aggregated over the UMIs). 
    Also required: changing GEX_bc and FB_bc to just cellBC. keeping CasTag required. 
    '''
    ncells = len(set(ct.cellBC))
    print(f'Number of Cells, pre-filtering || {ncells}')
    ct = ct.loc[ct['aggRC'] > rc_thresh]
    print(f'Read Count per cellBC+Tag || {ncells - len(set(ct.cellBC))} cells filtered out.')
    ncells = len(set(ct.cellBC))

    ct = ct.loc[ct["nUMI"] >= umi_thresh]
    print(f'Number UMI per cellBC+Tag || {ncells - len(set(ct.cellBC))} cells filtered out.')
    ncells = len(set(ct.cellBC))
    
    #### Calc Proportion 
    ct['propReads'] = ct['aggRC'] / ct.groupby(['cellBC']).aggRC.transform(sum)

    ### Optional Step -- no filtering yet, just defining the filter column 
    #1. Evaluate whether the cas-tag contains a run (minimum = 5, this defaulted...) 
    ntags = len(set(ct.CasTag))
    ct['containsRun'] = ct.CasTag.apply(lambda x: bool(re.search(run_pattern, x))) # find tags with repeats 
    ct = ct.loc[ct['containsRun'] == False]
    print(f'Removed {ntags - len(set(ct.CasTag))} tags with nucleotide runs || {ncells - len(set(ct.cellBC))} cells filtered out.')
    ntags = len(set(ct.CasTag))
    ncells =  len(set(ct.cellBC))
    
    ### Final aggressive filtering step
    ct = ct.loc[ct["propReads"] > prop_thresh]
    print(f'Min. Prop Reads per Cell || {ncells - len(set(ct.cellBC))} cells filtered out.')
    ncells = len(set(ct.cellBC))
    print(f'Final cell count || {ncells}')
    return ct

import networkx as nx
def join_digraphs(digraphs, new_root, subroot=None):
    '''
    chatGPT | i have more than one digraph that I want to join together at a new root node, 
    I want to code to work for any number of digraphs such that the new root node an out degree
    of n where n is the number of input digraphs.
    '''
    new_digraph = nx.DiGraph()
    
    # Add the new root node to the new digraph
    new_digraph.add_node(new_root)
    
    # Add edges from the new root node to the roots of the original digraphs
    if subroot is not None: 
        for ix in range(len(digraphs)):
            digraph=digraphs[ix]
            root_nodes = [node for node in digraph.nodes if digraph.in_degree(node) == 0]

            new_digraph.add_edge(new_root, subroot[ix])
            
            for root_node in root_nodes:
                new_digraph.add_edge(subroot[ix], root_node)
    else: 
        for digraph in digraphs:
            root_nodes = [node for node in digraph.nodes if digraph.in_degree(node) == 0]            
            for root_node in root_nodes:
                new_digraph.add_edge(new_root, root_node)
    
    # Add all nodes and edges from the input digraphs to the new digraph
    for digraph in digraphs:
        for node, edges in digraph.adjacency():
            new_digraph.add_node(node)
            for edge in edges:
                new_digraph.add_edge(node, edge)
    
    return new_digraph

### Tag Utils 

# import networkx as nx
from itertools import combinations 
from community import community_louvain
import igraph as ig
import leidenalg as la
import subprocess
import pandas as pd
import os 

def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union != 0 else 0

def generate_graph(cooccurrences):
    I = nx.Graph()

    for i, s in enumerate(cooccurrences):
        I.add_node(i, label=cooccurrences.index[i, ])

    for (i, s1), (j, s2) in combinations(enumerate(cooccurrences), 2):
        weight = jaccard_similarity(s1, s2)
        if weight > 0:  # Only add edges with non-zero similarity
            I.add_edge(i, j, weight=weight)
    
    return I

def get_louvain(I, res):
    comms2 = community_louvain.best_partition(I, random_state=1234, resolution=res)
    unique_coms2 = np.unique(list(comms2.values()))
    print(f'{len(unique_coms2)} communities detected with louvain algorithm.' )

    new_comms = {I.nodes[k]['label']:v for k,v in comms2.items()}
    
    return new_comms

def get_leiden(I, res):

    h = ig.Graph.from_networkx(I)
    partition = la.find_partition(h, la.CPMVertexPartition, resolution_parameter =res, weights='weight')
    unique_coms2 = np.unique(partition.membership)
    print(f'{len(unique_coms2)} communities detected with leiden algorithm.' )
    
    new_comms =  {I.nodes[k]['label']:v for k,v in enumerate(partition.membership)}
    
    return new_comms
    
def do_clusters(data, res=0.8):
    cooccurrences = data.groupby('cellBC')['CasTag'].apply(lambda x: set(x))
    I = generate_graph(cooccurrences)
    new_cooccurrences = cooccurrences.reset_index()

    louvain_map = get_louvain(I, res)
    new_cooccurrences['louvain'] = new_cooccurrences.cellBC.map(louvain_map)
    
    leiden_map = get_leiden(I, res)
    new_cooccurrences['leiden'] = new_cooccurrences.cellBC.map(leiden_map)
   
    return new_cooccurrences

def expand_dataframe(df, set_column):
    rows = []
    for _, row in df.iterrows():
        set_values = row[set_column]
        for value in set_values:
            new_row = row.copy()
            new_row[set_column] = value
            rows.append(new_row)
    return pd.DataFrame(rows)

def cell_counts(df, cellkey='cellBC'):
    n_obs = df.shape[0]
    n_cells = len(set(df[cellkey]))

    print(f'Dataset has {n_obs} observations and {n_cells} cells.')
    
    return n_obs, n_cells

def starcode_correct(df, tag_id='founderTag', do_count=False, ref=None):
    in_fasta = 'tmp.txt'
    
    if do_count:
        print('Running tags as counts')
        df_counts = df.groupby([tag_id]).size().reset_index()
        
        if ref is not None:
            ref_data = pd.read_csv(ref, index_col=0)
            ref_data.columns = ['castag_raw', 0]
            df_counts = pd.concat([ref_data, df_counts])
            
        df_counts.to_csv(in_fasta, index=False, header=False, sep='\t')
    else:
        print('Running tags as raw list') 
        df[tag_id].to_csv(in_fasta, index=False, header=False)
    
    starcode_path = '/dartfs/rc/lab/M/McKennaLab/projects/hannah/software/starcode'
    starcodeCommand = [starcode_path + '/starcode', '-d', '1', '-t', '20', '-i', in_fasta, '-o', 'tmp_out.txt', '--print-clusters']
    subprocess.call(starcodeCommand)
    
    starcode_corrected = pd.read_table('tmp_out.txt', sep='\t', header=None)
    starcode_corrected.columns = ['assignedTag', 'clusterSize', 'originalTag']
    sclust_expanded = starcode_corrected.assign(originalTag=starcode_corrected['originalTag'].str.split(',')).explode('originalTag').reset_index(drop=True)
    
    df1 = df.merge(sclust_expanded, left_on=tag_id, right_on='originalTag')
    
    return df1


def ensure_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")
