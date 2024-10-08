import os
import time
import numpy as np
from scipy.stats import fisher_exact
import pandas as pd
from statsmodels.stats.multitest import multipletests

from joblib import Parallel, delayed
### NETWORK
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from django.conf import settings

def EASE(query_set, genesets, PT=20000):
    """
    Run EASE analysis.

    Parameters:
    A (list or array-like): User genes.
    PT (int, optional): Population total. Default is 18000.
    gs (dict): Dictionary containing gene sets.

    Returns:
    pandas.DataFrame: Enrichment results.
    """
    LT = len(query_set)
    PT = int(PT)  # Convert PT to integer (experimental)
    gs_names = list(genesets.keys())
    enrichment_df = []
    for gs_n in gs_names:
        B = genesets[gs_n]
        PH = len(B)
        LH = len(set(query_set).intersection(B))

        modLH = LH if LH != 0 else 1
        # Implement EASE score from DAVID
        contingency_table = np.array([[modLH-1, PH-LH+1],
                                      [LT-LH, PT-LT-(PH-LH)]])
        _, p_value = fisher_exact(contingency_table, alternative='greater')

        enrichment_df.append({
            'pathway_id': gs_n,
            'Overlap': LH,
            'Pvalue': p_value,
        })
    
    enrichment_df = pd.DataFrame(enrichment_df)
    # Adjust p-values using multiple testing correction
    _, enrichment_df['FDR'], _, _ = multipletests(enrichment_df['Pvalue'], method='fdr_bh')
    _, enrichment_df['FWER'], _, _ = multipletests(enrichment_df['Pvalue'], method='bonferroni')

    return enrichment_df


if __name__ == '__main__':

    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus,prefer="threads")

    speciesA = '83333'
    species_objects = Species.objects.get(Q(tax_id=speciesA))
    speciesA_name = species_objects.species_name
    speciesA_name_compact = '%s.%s' %(speciesA_name.split()[0][0],speciesA_name.split()[1])
    network_f = '/mnt/sdb/FunCoup/pyFunCoup/FunCoup/website/static/website/networks/FunCoup6.0/FC6.0_' + speciesA_name_compact + '_full.gz'
    network_df = pd.read_csv(network_f, compression='gzip', sep='\t', usecols=[2,3,5])
    network_df.columns = ['proteinA_id','proteinB_id','ppv']
    # Filter the DataFrame based on the cutoff weight
    cutoff_ppv = 0.95
    network_df = network_df[network_df['ppv'] >= cutoff_ppv]
    network_df = network_df.drop('ppv', axis=1).reset_index(drop=True)
    network_genes = list(set(network_df['proteinA_id'].to_list()+network_df['proteinB_id'].to_list()))
    
    ## Test
    query_set = list(np.random.choice(network_genes, 10, replace=False))

    genesets = pd.DataFrame.from_records(Pathway.objects.filter(Q(species=speciesA)).values('protein_id','pathway_id'))
    genesets = genesets[genesets['protein_id'].isin(network_genes)]
    ## genesets to dict
    genesets_dict = genesets.groupby('pathway_id')['protein_id'].apply(list).to_dict()

    pathway_id = '01503'
    query_set = genesets[genesets['pathway_id'] == pathway_id]['protein_id'].tolist()

    ## PT --> genome coverage
    genome_coverage = 4000

    ## Pre-compute crosstalk
    time_start = time.perf_counter()
    ease_results = EASE(query_set, genesets_dict, PT=genome_coverage)
    time_elapsed = (time.perf_counter() - time_start)
    print('### EASE in %5.1f secs' %time_elapsed)

    ease_results[ease_results['pathway_id']==pathway_id]