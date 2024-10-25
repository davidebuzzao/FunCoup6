import os
import pandas as pd
import numpy as np
import time
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm
from auxiliary_functions import *
### NETWORK
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from django.conf import settings
import yaml
from yaml.loader import SafeLoader

##### pre-computed anubix
def pre_compute_crosstalk(graph, genesets, pathway_id):
    # Iterate over pathways and measure crosstalk
    pathway_genes = genesets[genesets['pathway_id'] == pathway_id]['protein_id'].tolist()
    # Measure crosstalk using edge_boundary and append to the list
    crosstalk_pathway = []
    for node_id in list(graph.nodes):
        crosstalk_pathway.append(len(list(nx.edge_boundary(graph,{node_id}, pathway_genes))))
    return crosstalk_pathway

def precomputeAnubix():
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)

    print("Initiating precomputation of ANUBIX crosstalk")
    if os.path.exists('configFiles/websiteConfig.yml'):
        with open('configFiles/websiteConfig.yml') as f:
            webConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your websiteConfig.yml does not exist")
        exit()

    cutoff_ppv = float(webConfig['anubixCutoff'])
    networksPath="website/static/website/networks/FunCoup6.0/"
    for networkFile in os.listdir(networksPath):
        if networkFile.endswith("_full.gz") and "TRANSFERRED" not in networkFile and 'SARS-CoV-2' not in networkFile:
            spName = networkFile.replace("FC6.0_","").replace("_full.gz","").split(".")[1]
            species_object = Species.objects.filter(species_name__icontains=spName)
            taxid=species_object[0].tax_id
            print(taxid)

            network_df = pd.read_csv(networksPath+"/"+networkFile, compression='gzip', sep='\t', usecols=[2,3,5])
            network_df.columns = ['proteinA_id','proteinB_id','ppv']
            # Filter the DataFrame based on the cutoff weight
            network_df = network_df[network_df['ppv'] >= cutoff_ppv]
            network_df = network_df.drop('ppv', axis=1).reset_index(drop=True)
            network_genes = list(set(network_df['proteinA_id'].to_list()+network_df['proteinB_id'].to_list()))
            graph = nx.from_pandas_edgelist(network_df, 'proteinA_id','proteinB_id')
      
            genesets = pd.DataFrame.from_records(Pathway.objects.filter(Q(species=taxid)).values('protein_id','pathway_id'))
            genesets = genesets[genesets['protein_id'].isin(network_genes)]

            # Count occurrences of each pathway_id
            pathway_counts = genesets['pathway_id'].value_counts()
            pathway_size_cutoff = 5
            filtered_genesets = genesets[genesets['pathway_id'].isin(pathway_counts[pathway_counts >= pathway_size_cutoff].index)].reset_index(drop=True)
            pathways = filtered_genesets['pathway_id'].unique()

            min_bin_size = 150
            degree_f = 'website/static/website/anubix/%s_degree-%i_ppv%.2f.pkl' %(taxid,min_bin_size,cutoff_ppv)
            if os.path.exists(degree_f)==False:
                print('### Extracting network degrees')
                # Calculate degrees and create a DataFrame
                degrees = pd.DataFrame(list(graph.degree()), columns=['node', 'degree'])
                # Sort the DataFrame by degree and node
                degrees = degrees.sort_values(by=['degree', 'node']).reset_index(drop=True)

                # Calculate the bin for each node based on min_bin_size
                degrees['bin'] = (degrees.index // min_bin_size) + 1
                joblib.dump(degrees, degree_f)

            precomputed_crosstalk_f = 'website/static/website/anubix/%s_xtalk_ppv%.2f.pkl' %(taxid,cutoff_ppv)
            if os.path.exists(precomputed_crosstalk_f)==False:      
                print('### Computing Xtalk')
                ## Pre-compute crosstalk
                time_start = time.perf_counter()
                with tqdm_joblib(tqdm(desc='Pre-compute Xtalk', total=len(pathways))) as progress_bar:
                    precomputed_crosstalk_matrix = parallel(delayed(pre_compute_crosstalk)(graph, filtered_genesets, pathway_id) for pathway_id in pathways)
                time_elapsed = (time.perf_counter() - time_start)
                print('### Pre-computed Xtalk in %5.1f secs' %time_elapsed)

                precomputed_crosstalk_df = pd.DataFrame(precomputed_crosstalk_matrix, index=pathways, columns=graph.nodes).T
                precomputed_crosstalk_df.info()
                joblib.dump(precomputed_crosstalk_df, precomputed_crosstalk_f)

