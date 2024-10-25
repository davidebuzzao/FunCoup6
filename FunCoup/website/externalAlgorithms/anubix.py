import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import betabinom
from scipy.optimize import minimize
import networkx as nx
import joblib
### NETWORK
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from django.conf import settings

def mid_p_value(obv, max_val, alpha, beta):
    # Calculate the exact PMF for the observed value
    exact_prob = betabinom.pmf(obv, max_val, alpha, beta)
    
    # Calculate the cumulative sum for values greater than obv
    cumulative_prob = betabinom.sf(obv, max_val, alpha, beta)  # sf gives the survival function (1 - cdf)
    
    # Calculate the mid-p-value
    mid_pvalue = 0.5 * exact_prob + cumulative_prob
    
    return mid_pvalue

def log_likelihood(params, *args):
    alpha, beta = params
    data, n = args
    ll = np.sum(betabinom.logpmf(data, n, alpha, beta))
    return -ll

def fit_beta_binomial(crosstalk_values,length_query_set,length_pathway_dict,optimizeParams=True):
    ### w or w/o ML optimization of alpha and beta params
    # Ex: crosstalk_values = crosstalk.iloc[:,1]
    observed_crosstalk = crosstalk_values.iloc[-1]
    random_crosstalk = crosstalk_values.iloc[:-1]

    length_pathway = length_pathway_dict[crosstalk_values.name]
    max_value = (length_query_set * length_pathway) - min(length_query_set, length_pathway)
    n = max_value  # This should match the R definition of n

    m_1 = np.mean(random_crosstalk)
    m_2 = np.mean(random_crosstalk**2)

    # Initial alpha and beta guess
    denominator = (n * (m_2 / m_1 - m_1 - 1) + m_1)
    if denominator == 0:
        return np.nan, np.nan

    alpha = (n * m_1 - m_2) / denominator
    beta = (n - m_1) * (n - m_2 / m_1) / denominator

    if np.isnan(alpha) or np.isinf(alpha) or np.isnan(beta) or np.isinf(beta):
        return np.nan, np.nan

    if optimizeParams:
        # Optimize alpha and beta using the log-likelihood function
        # It double the execution time!!!
        alpha_init, beta_init = alpha, beta

        result = minimize(log_likelihood, [alpha_init, beta_init], args=(random_crosstalk, n),
                        bounds=((1e-10, None), (1e-10, None)), method='L-BFGS-B')

        if not result.success:
            return np.nan, np.nan

        alpha, beta = result.x

    # Calculate the p-value for the observed crosstalk
    # p_value = 1 - betabinom.cdf(observed_crosstalk, n, alpha_opt, beta_opt)
    p_value = mid_p_value(observed_crosstalk, max_value, alpha, beta)
    
    ## From formula
    # expected_crosstalk = n * alpha_opt / (alpha_opt + beta_opt)
    ## From R anubix implementation
    expected_crosstalk = np.mean(random_crosstalk)

    return p_value, expected_crosstalk

def constrained_sampling(query_degree_dict,degrees,seed):
    sample_genes = []
    for bin in query_degree_dict:
        nodes = degrees[degrees['bin'] == bin]['node']
        count = query_degree_dict[bin]
        np.random.seed(seed)
        sample_genes.extend(list(np.random.choice(nodes, count, replace=False)))
    return sample_genes

def sample_crosstalk(graph, query_degree_dict, degrees, filtered_genesets, pathways):
    # Sample query set of genes
    query_genes = constrained_sampling(query_degree_dict, degrees)
    # Initialize a list to store crosstalk results for this iteration
    crosstalk_iteration = []
    # Iterate over pathways and measure crosstalk
    for pathway_id in pathways:
        pathway_genes = filtered_genesets[filtered_genesets['pathway_id'] == pathway_id]['protein_id'].tolist()
        # Measure crosstalk using edge_boundary and append to the list
        crosstalk_iteration.append(len(list(nx.edge_boundary(graph, query_genes, pathway_genes))))
    return crosstalk_iteration

def compute_crosstalk(precomputed_crosstalk_df,query_set):
    return precomputed_crosstalk_df.loc[query_set, :].sum()

def ANUBIX(query_set,speciesA,ppv_cutoff=0.95,num_iterations=2000,degree_bin_size=150,optimizeParams=False): 
    degree_f = 'website/static/website/anubix/%s_degree-%i_ppv%.2f.pkl' %(speciesA,degree_bin_size,ppv_cutoff)
    if os.path.exists(degree_f):
        # print('### Reading network degrees')
        degrees = joblib.load(degree_f)
    else:
        degrees = []
        print("Could not find precomputed degrees")

    if len(degrees)>0: query_set = [qv for qv in query_set if qv in degrees['node'].to_list()] # Filter out query genes that are not in the network at cutoff
    else: query_set = []

    precomputed_crosstalk_f = 'website/static/website/anubix/%s_xtalk_ppv%.2f.pkl' %(speciesA,ppv_cutoff)
    precomputed_crosstalk_df = []
    if os.path.exists(precomputed_crosstalk_f):
        # print('### Reading pickled Xtalk')
        precomputed_crosstalk_df = joblib.load(precomputed_crosstalk_f)
        # print("Done getting crosstalk")
    else:        
        print("Could not find precomputed crosstalk")

    length_query_set = len(query_set)
    if len(precomputed_crosstalk_df)==0:
        random_strings = [''.join(np.random.choice(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'), size=4)) for _ in range(50)]
        results_df = pd.DataFrame({
            'pathway_id': random_strings,
            'Xtalk': 0,
            'Pvalue': 1.0,
            'FDR': 1.0,
            'FWER': 1.0,
            'expectedCrosstalk': 0
            })
    elif length_query_set==0:
        pathway_ids = list(precomputed_crosstalk_df.columns)
        results_df = pd.DataFrame({
            'pathway_id': pathway_ids,
            'Xtalk': 0,
            'Pvalue': 1.0,
            'FDR': 1.0,
            'FWER': 1.0,
            'expectedCrosstalk': 0
        })
    else:
        query_degree_dict = {}
        for qv in query_set:
            bin = int(degrees[degrees['node']==qv]['bin'])
            query_degree_dict[bin] = query_degree_dict.get(bin,0)
            query_degree_dict[bin] += 1

        ## Observerd crosstalk, Subset rows based on the index names
        observed_crosstalk = compute_crosstalk(precomputed_crosstalk_df,query_set)
        np.random.seed(int(speciesA))
        seeds = np.random.choice(range(1000000), num_iterations, replace=False)
        crosstalk = []
        for seed in seeds:
            sampled_genes = constrained_sampling(query_degree_dict, degrees, seed)
            crosstalk.append(compute_crosstalk(precomputed_crosstalk_df,sampled_genes))
        
        # Convert the crosstalk matrix to a DataFrame for better presentation
        crosstalk = pd.concat(crosstalk,axis=1)

        ## Stack the observation at the bottom
        crosstalk = pd.concat([crosstalk,observed_crosstalk],axis=1).T.reset_index(drop=True)

        # Apply the function to each column (pathway) in the crosstalk matrix
        # p_values = crosstalk.apply(fit_beta_binomial, axis=0)
        network_genes = precomputed_crosstalk_df.index
        genesets = pd.DataFrame.from_records(Pathway.objects.filter(Q(species=speciesA)).values('protein_id','pathway_id'))
        genesets = genesets[genesets['protein_id'].isin(network_genes)]
        length_pathway_dict = genesets['pathway_id'].value_counts().to_dict()
        # Otherwise use pre-computed
        # length_pathway_f = 'website/static/website/anubix/%s_pathway_length.pkl' %speciesA
        # length_pathway_dict = {}
        # if os.path.exists(length_pathway_f):
        #     length_pathway_dict = joblib.load(length_pathway_f)
        # else:        
        #     print("Could not find precomputed pathway lengths")
        p_values = crosstalk.apply(lambda col: fit_beta_binomial(col, length_query_set, length_pathway_dict, optimizeParams), axis=0)

        p_vals = [float(i) for i in p_values.iloc[0].fillna(1).tolist()]
        expCrosstalk=p_values.iloc[1].fillna(1)

        # Adjust p-values using multiple testing correction
        _, fdr, _, _ = multipletests(p_vals, method='fdr_bh')
        _, fwer, _, _ = multipletests(p_vals, method='bonferroni')

        results_df = pd.DataFrame({
            'pathway_id': crosstalk.columns,
            'Xtalk': observed_crosstalk,
            'expectedCrosstalk':expCrosstalk,
            'Pvalue': p_vals,
            'FDR': fdr,
            'FWER': fwer
        }).reset_index(drop=True)

    return results_df


if __name__ == '__main__':
    speciesA = '9606'
    ppv_cutoff=0.95
    num_iteration=10000
    optimizeParam=True
    # degree_f = 'website/static/website/anubix/%s_degree_ppv%.2f.pkl' %(speciesA,ppv_cutoff)
    # degrees = joblib.load(degree_f)
    # network_genes = degrees['node'].to_list()
    # query_set=list(np.random.choice(network_genes, 10, replace=False)) 
    query_set_df = pd.read_csv('website/static/website/anubix/jun_subnetwork.tsv', sep='\t')
    query_set = query_set_df['protein_id'].to_list()
    enrichment_df = ANUBIX(query_set,speciesA,ppv_cutoff=ppv_cutoff,num_iterations=num_iteration, optimizeParams=optimizeParam)
    enrichment_df = enrichment_df.sort_values(by=['FWER', 'Pvalue'], ascending=[True, True])
    enrichment_df.to_csv('website/static/website/anubix/JUNsubnet-EAanubixPython.tsv', sep='\t', index=None)