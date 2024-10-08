import os
import random
import numpy as np
import pandas as pd
import networkx as nx
import numpy as np
from itertools import chain
from joblib import Parallel, delayed
from tqdm import tqdm

def TOPAS(network, seeds, expansion_steps=2, cores=1):
    if network is None or len(network) < 2:
        return None
    if seeds is None:
        return None
    if expansion_steps is None:
        expansion_steps = 2
    if cores is None:
        cores = max(1, os.cpu_count())
    if not isinstance(cores, int) or cores <= 0:
        return None
    
    ###################################
    # STEP 1 - Full Seed Network
    print('STEP 1: Full Seed Network')
    # Create a networkx graph from network
    network.columns = ['gene1','gene2']
    seeds_graph = nx.from_pandas_edgelist(network,'gene1','gene2', create_using=nx.Graph())
    seeds_graph.remove_edges_from(nx.selfloop_edges(seeds_graph))

    seeds = list(map(str, seeds))
    ###################################
    # STEP 2 - Largest Connected Module
    print(f"STEP 2: Largest Connected Module (using {cores} {'cores' if cores > 1 else 'core'})")
    def _extract_lcc(graph):
        # This function takes a graph as input, computes all connected components,
        # and outputs the one with the largest number of seeds.
        components = nx.connected_components(graph)
        largest_cc = max(components, key=len)
        return graph.subgraph(largest_cc)
    
    # Extract the largest seed connected component
    graph_lcc = _extract_lcc(seeds_graph)
    # Update seeds
    seeds = list(set(seeds) & set(graph_lcc.nodes))
    def _sp_compute(source_v, seeds, graph):
        # This function computes all shortest paths between a seed and all other seeds
        # and returns only those which are shorter than the input "expansion_steps".        
        sp_v = []
        for dv in seeds:
            sp = nx.shortest_path(graph, source=source_v, target=dv)[1:-1]
            if len(sp)<=expansion_steps:
                sp_v += sp
            
        return list(set(sp_v))

    # Proceed with the computation of shortest paths
    progress = tqdm(seeds)
    connectors = list(set(chain.from_iterable(Parallel(n_jobs=cores)(delayed(_sp_compute)(seed, seeds, graph_lcc) for seed in progress))))
    # Extract a subgraph composed of seeds and potential connectors
    seeds_graph = graph_lcc.subgraph(list(set(seeds) | set(connectors)))
    
    # If no edge is retrieved, the output is None
    if seeds_graph.number_of_edges() < 1:
        return None
    
    # If more than one component exists, take the largest seed component
    if nx.number_connected_components(seeds_graph) > 1:
        seeds_graph = _extract_lcc(seeds_graph)
    
    # Update seeds
    seeds = list(set(seeds) & set(graph_lcc.nodes))
    
    ###################################
    # STEP 3 - Pruned Module
    print('\nSTEP 3: Pruned Module')
    network_df = pd.DataFrame(np.zeros(seeds_graph.number_of_nodes()),columns=['seed'],index=seeds_graph.nodes())
    network_df['seed'] = np.where(network_df.index.isin(seeds), 1, 0)
    random.seed(1996)
    stationary_probability = nx.pagerank(seeds_graph,personalization=network_df['seed'].to_dict(), alpha=0.75) 
    network_df['stationary_probability'] = list(stationary_probability.values())
    network_df['gene'] = network_df.index
    # Randomize the alphabetical order first, it adds nothing if all probabilities vary
    network_df = network_df.sample(frac=1, random_state=1996)
    # Order by probability afterwards
    network_df.sort_values('stationary_probability', inplace=True)
    network_df.reset_index(drop=True, inplace=True)
    
    # Iteratively test connectors, from the one with the lowest stationary probability to the highest
    for _, row in network_df[network_df['seed'] == 0].iterrows():
        v = row['gene']
        if v not in seeds_graph.nodes:
            continue
        
        # Temporary delete a candidate connector from the module
        tmp_graph = seeds_graph.copy()
        tmp_graph.remove_node(v)
        
        # If removing a vertex generates more than 1 connected component,
        # then check if seed coverage decreases.
        # Else remove the vertex permanently.
        if nx.number_connected_components(tmp_graph) > 1:
            best_cc = 1
            max_seeds = 0
            tmp_graph_lcc = _extract_lcc(tmp_graph)
            
            for i, component in enumerate(nx.connected_components(tmp_graph_lcc)):
                seeds_cc = sum(1 for node in component if node in seeds)
                if seeds_cc > max_seeds:
                    best_cc = i
                    max_seeds = seeds_cc
            
            if (len(seeds) - max_seeds) == 0:
                tmp_graph_lcc_vertices = [v for v in tmp_graph_lcc.nodes if v in tmp_graph_lcc[best_cc]]
                seeds_graph = tmp_graph.subgraph(tmp_graph_lcc_vertices).copy()
        else:
            seeds_graph = tmp_graph.copy()

    module_df = nx.to_pandas_edgelist(seeds_graph)
    module_df.columns = ['gene1','gene2']
    
    return module_df

if __name__ == '__main__':
    ### NETWORK
    cutoff_ppv = 0.8
    network_f = 'data/FC5.0_H.sapiens_compact.gz'
    network_df = pd.read_csv(network_f, sep='\t', usecols=[0,2,3])
    network_df.columns = ['ppv','gene1','gene2']
    network_df['ppv'] = network_df['ppv'].astype(float)

    # Filter the DataFrame based on the cutoff weight
    network_df = network_df[network_df['ppv'] >= cutoff_ppv]
    network_df = network_df.drop('ppv', axis=1).reset_index(drop=True)

    ### SEED GENES
    seed_genes_f = 'data/adrenal_gland_diseases.txt'
    seed_genes = pd.read_csv(seed_genes_f, sep='\t', header=None)[0].to_list()

    ## MaxLink
    candidate_genes = TOPAS(network_df, seed_genes, expansion_steps=2, cores=4)
