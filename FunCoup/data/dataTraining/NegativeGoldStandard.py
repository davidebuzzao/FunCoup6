import os
import django
from collections import Counter
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
import yaml
from yaml.loader import SafeLoader
import random
import pandas as pd
import numpy as np

from data.dataTraining.ExtractOrthologLinks import *

def sampleNegativeGoldstandard(evidenceLink,tax_id,trainingConfig):
    min_n_sample, max_n_sample = trainingConfig['nSample']

    tot = evidenceLink.shape[0]
    if tot < min_n_sample: return None
    
    sample_size = min(tot,max_n_sample)
    random.seed(int(tax_id))
    sample_evidenceLink = evidenceLink.sample(sample_size).reset_index(drop=True)
    return sample_evidenceLink
    # if max_n_sample > sample_size:
    #     return sample_evidenceLink
    # else:
    #     # Filter out over-studied genes above 90th percentiles
    #     combined_p = sample_evidenceLink['proteinA_id'].to_list() + sample_evidenceLink['proteinB_id'].to_list()
    #     combined_p_counts = Counter(combined_p)
    #     combined_p_vals = [val for p,val in combined_p_counts.items()]
    #     # Calculate the 10th and 90th percentiles
    #     # percentile_10 = np.percentile(combined_p_vals, 5)
    #     percentile_95 = np.percentile(combined_p_vals, 95)
    #     # Set thresholds within the 10th and 90th percentiles
    #     selected_p = [p for p,val in combined_p_counts.items() if val <= percentile_95]
    #     # selected_p = [p for p,val in combined_p_counts.items() if percentile_10 <= val <= percentile_90]

    #     filtered_sample_evidenceLink = sample_evidenceLink[sample_evidenceLink['proteinA_id'].isin(selected_p) & \
    #                                             sample_evidenceLink['proteinB_id'].isin(selected_p)].reset_index(drop=True)
    #     # sample_evidenceLink[~sample_evidenceLink['proteinA_id'].isin(selected_p) & \
    #     #                                         ~sample_evidenceLink['proteinB_id'].isin(selected_p)].reset_index(drop=True)
    #     return filtered_sample_evidenceLink

def filterNegativeGoldStandard(negative_goldStandardLink,spAGoldStandard,col_name,instanceConfig):
    speciesA = spAGoldStandard.species.tax_id

    ## Filter out the positive links of the respective species, and also the species ortholog positive links.
    for speciesB in instanceConfig['instance']['species']:
        spBGoldStandard = GoldStandard.objects.filter(Q(type=spAGoldStandard.type) & Q(species__tax_id=speciesB))
        if not spBGoldStandard.exists(): continue
        spBGoldStandard=spBGoldStandard[0]

        if speciesB == speciesA:
            # Extract positive gold standard links from species A
            positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=spBGoldStandard)).values('proteinA_id','proteinB_id'))
            if (positive_goldStandardLink.shape[0]==0): continue
            else: positive_goldStandardLink.columns = ['proteinA_id','proteinB_id']
        else:
            # Extract positive gold standard links from species B
            positive_goldStandardLink = extractOrthologGoldStandardLinks(spBGoldStandard,speciesA)
            if positive_goldStandardLink is None: continue

        positive_goldStandardLink['proteinA_id'] = positive_goldStandardLink['proteinA_id'].astype(int)
        positive_goldStandardLink['proteinB_id'] = positive_goldStandardLink['proteinB_id'].astype(int)
        
        negative_goldStandardLink = \
            pd.merge(negative_goldStandardLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer')\
                .query('_merge=="left_only"')\
                .drop('_merge', axis=1)

    if len(negative_goldStandardLink)>0:
        negative_goldStandardLink = negative_goldStandardLink[col_name].to_list()
        return negative_goldStandardLink
    else:
        return None

def separateNegativePositiveGoldStandard(evidenceLink,spAGoldStandard,col_name,instanceConfig):
    speciesA = spAGoldStandard.species.tax_id

    ## Filter out the positive links of the respective species, and also the species ortholog positive links.
    speciesA_positive_goldStandardLink = pd.DataFrame()
    positive_goldStandardLink = pd.DataFrame()
    for speciesB in instanceConfig['instance']['species']:
        spBGoldStandard = GoldStandard.objects.filter(Q(type=spAGoldStandard.type) & Q(species__tax_id=speciesB))
        if not spBGoldStandard.exists(): continue
        spBGoldStandard=spBGoldStandard[0]

        if speciesB == speciesA:
            # Extract positive gold standard links from species A
            tmp_positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=spBGoldStandard)).values('proteinA_id','proteinB_id'))
            if (tmp_positive_goldStandardLink.shape[0]==0): continue
            else: tmp_positive_goldStandardLink.columns = ['proteinA_id','proteinB_id']
            if speciesA_positive_goldStandardLink.empty: speciesA_positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
            else: speciesA_positive_goldStandardLink = pd.merge(speciesA_positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')
        else:
            # Extract positive gold standard links from species B
            tmp_positive_goldStandardLink = extractOrthologGoldStandardLinks(spBGoldStandard,speciesA)
            if tmp_positive_goldStandardLink is None: continue
        
        tmp_positive_goldStandardLink['proteinA_id'] = tmp_positive_goldStandardLink['proteinA_id'].astype(int)
        tmp_positive_goldStandardLink['proteinB_id'] = tmp_positive_goldStandardLink['proteinB_id'].astype(int)
        if positive_goldStandardLink.empty: positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
        else: positive_goldStandardLink = pd.merge(positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')

    if len(positive_goldStandardLink)>0:
        negative_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)
        positive_goldStandardLink = pd.merge(evidenceLink,speciesA_positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="both"').drop('_merge', axis=1).reset_index(drop=True)

    if len(positive_goldStandardLink)>0 and len(negative_goldStandardLink)>0:
        positive_goldStandardLink = positive_goldStandardLink[col_name].to_list()
        negative_goldStandardLink = negative_goldStandardLink[col_name].to_list()
        return positive_goldStandardLink,negative_goldStandardLink
    else:
        return None,None

## Not used ##
def separateNegativePositiveGoldStandardAllSpecies(evidenceLink,spAGoldStandard,col_name,instanceConfig,goldStandardConfig):
    speciesA = spAGoldStandard.species.tax_id

    ## Filter out the positive links of the respective species, and also the species ortholog positive links.
    positive_goldStandardLink = pd.DataFrame()
    for speciesB in instanceConfig['instance']['species']:
        spBGoldStandard = GoldStandard.objects.filter(Q(type=spAGoldStandard.type) & Q(species__tax_id=speciesB))
        if not spBGoldStandard.exists(): continue
        spBGoldStandard=spBGoldStandard[0]

        if speciesB == speciesA:
            # Extract positive gold standard links from species A
            tmp_positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=spBGoldStandard)).values('proteinA_id','proteinB_id'))
            if (tmp_positive_goldStandardLink.shape[0]==0): continue
            else: tmp_positive_goldStandardLink.columns = ['proteinA_id','proteinB_id']
        else:
            # Extract positive gold standard links from species B
            tmp_positive_goldStandardLink = extractOrthologGoldStandardLinks(spBGoldStandard,speciesA)
            if tmp_positive_goldStandardLink is None: continue
        
        tmp_positive_goldStandardLink['proteinA_id'] = tmp_positive_goldStandardLink['proteinA_id'].astype(int)
        tmp_positive_goldStandardLink['proteinB_id'] = tmp_positive_goldStandardLink['proteinB_id'].astype(int)
        if positive_goldStandardLink.empty: positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
        else: positive_goldStandardLink = pd.merge(positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')

    if len(positive_goldStandardLink)>0:
        negative_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)
        positive_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="both"').drop('_merge', axis=1).reset_index(drop=True)

    if len(positive_goldStandardLink)>0 and len(negative_goldStandardLink)>0:
        positive_goldStandardLink = positive_goldStandardLink[col_name].to_list()
        negative_goldStandardLink = negative_goldStandardLink[col_name].to_list()
        return positive_goldStandardLink,negative_goldStandardLink
    else:
        return None,None

def separateNegativePositiveGoldStandardAllGSAllSpecies(evidenceLink,spAGoldStandard,col_name,instanceConfig,goldStandardConfig):
    speciesA = spAGoldStandard.species.tax_id
    g_type = spAGoldStandard.type
    ## Filter out the positive links of the respective species, and also the species ortholog positive links.
    positive_goldStandardLink = pd.DataFrame()

    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id__in=instanceConfig['instance']['species'])) \
                    if goldstandard_t.type in instanceConfig['instance']['goldStandard']]
    
    for spBGoldStandard in goldstandard_list:
        speciesB = spBGoldStandard.species.tax_id
        if speciesB == speciesA:
            # Extract positive gold standard links from species A
            tmp_positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=spBGoldStandard)).values('proteinA_id','proteinB_id'))
            if (tmp_positive_goldStandardLink.shape[0]==0): continue
            else: tmp_positive_goldStandardLink.columns = ['proteinA_id','proteinB_id']
        else:
            # Extract positive gold standard links from species B
            tmp_positive_goldStandardLink = extractOrthologGoldStandardLinks(spBGoldStandard,speciesA)
            if tmp_positive_goldStandardLink is None: continue
        
        tmp_positive_goldStandardLink['proteinA_id'] = tmp_positive_goldStandardLink['proteinA_id'].astype(int)
        tmp_positive_goldStandardLink['proteinB_id'] = tmp_positive_goldStandardLink['proteinB_id'].astype(int)
        if positive_goldStandardLink.empty: positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
        else: positive_goldStandardLink = pd.merge(positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')

    if len(positive_goldStandardLink)>0:
        negative_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)
        positive_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="both"').drop('_merge', axis=1).reset_index(drop=True)

    if len(positive_goldStandardLink)>0 and len(negative_goldStandardLink)>0:
        positive_goldStandardLink = positive_goldStandardLink[col_name].to_list()
        negative_goldStandardLink = negative_goldStandardLink[col_name].to_list()
        return positive_goldStandardLink,negative_goldStandardLink
    else:
        return None,None

def separateNegativePositiveGoldStandardNoOverfitting(evidenceLink,spAGoldStandard,col_name,instanceConfig,goldStandardConfig):
    speciesA = spAGoldStandard.species.tax_id
    g_type = spAGoldStandard.type
    ## Filter out the positive links of the respective species, and also the species ortholog positive links.
    speciesA_positive_goldStandardLink = pd.DataFrame()
    positive_goldStandardLink = pd.DataFrame()

    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id__in=instanceConfig['instance']['species'])) \
                    if goldstandard_t.type in instanceConfig['instance']['goldStandard']]
    
    for spBGoldStandard in goldstandard_list:
        speciesB = spBGoldStandard.species.tax_id
        if speciesB == speciesA:
            # Extract positive gold standard links from species A
            tmp_positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=spBGoldStandard)).values('proteinA_id','proteinB_id'))
            if (tmp_positive_goldStandardLink.shape[0]==0): continue
            else: tmp_positive_goldStandardLink.columns = ['proteinA_id','proteinB_id']
            if spBGoldStandard.type == g_type:
                if speciesA_positive_goldStandardLink.empty: speciesA_positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
                else: speciesA_positive_goldStandardLink = pd.merge(speciesA_positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')
        else:
            # Extract positive gold standard links from species B
            tmp_positive_goldStandardLink = extractOrthologGoldStandardLinks(spBGoldStandard,speciesA)
            if tmp_positive_goldStandardLink is None: continue
        
        tmp_positive_goldStandardLink['proteinA_id'] = tmp_positive_goldStandardLink['proteinA_id'].astype(int)
        tmp_positive_goldStandardLink['proteinB_id'] = tmp_positive_goldStandardLink['proteinB_id'].astype(int)
        if positive_goldStandardLink.empty: positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
        else: positive_goldStandardLink = pd.merge(positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')

    if len(positive_goldStandardLink)>0:
        negative_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)
        # negative_goldStandardLink = evidenceLink.query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)
        positive_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="both"').drop('_merge', axis=1).reset_index(drop=True)
        ## To avoid "overfitting", we exclude GS links that were used for training
        positive_goldStandardLink = pd.merge(positive_goldStandardLink,speciesA_positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)

    if len(positive_goldStandardLink)>0 and len(negative_goldStandardLink)>0:
        positive_goldStandardLink = positive_goldStandardLink[col_name].to_list()
        negative_goldStandardLink = negative_goldStandardLink[col_name].to_list()
        return positive_goldStandardLink,negative_goldStandardLink
    else:
        return None,None

def separateNegativePositiveGoldStandard2b(evidenceLink,spAGoldStandard,col_name,instanceConfig,goldStandardConfig):
    speciesA = spAGoldStandard.species.tax_id

    ## Filter out the positive links of the respective species, and also the species ortholog positive links.
    # speciesA_positive_goldStandardLink = pd.DataFrame()
    positive_goldStandardLink = pd.DataFrame()

    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id__in=instanceConfig['instance']['species'])) \
                    if goldstandard_t.type in instanceConfig['instance']['goldStandard']]
    
    for spBGoldStandard in goldstandard_list:
        speciesB = spBGoldStandard.species.tax_id
        if speciesB == speciesA:
            # Extract positive gold standard links from species A
            tmp_positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=spBGoldStandard)).values('proteinA_id','proteinB_id'))
            if (tmp_positive_goldStandardLink.shape[0]==0): continue
            else: tmp_positive_goldStandardLink.columns = ['proteinA_id','proteinB_id']
            # if speciesA_positive_goldStandardLink.empty: speciesA_positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
            # else: speciesA_positive_goldStandardLink = pd.merge(speciesA_positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')
        else:
            # Extract positive gold standard links from species B
            tmp_positive_goldStandardLink = extractOrthologGoldStandardLinks(spBGoldStandard,speciesA)
            if tmp_positive_goldStandardLink is None: continue
        
        tmp_positive_goldStandardLink['proteinA_id'] = tmp_positive_goldStandardLink['proteinA_id'].astype(int)
        tmp_positive_goldStandardLink['proteinB_id'] = tmp_positive_goldStandardLink['proteinB_id'].astype(int)
        if positive_goldStandardLink.empty: positive_goldStandardLink = tmp_positive_goldStandardLink.copy()
        else: positive_goldStandardLink = pd.merge(positive_goldStandardLink,tmp_positive_goldStandardLink, on=['proteinA_id','proteinB_id'],how='outer')

    if len(positive_goldStandardLink)>0:
        negative_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)
        # negative_goldStandardLink = evidenceLink.query('_merge=="left_only"').drop('_merge', axis=1).reset_index(drop=True)
        positive_goldStandardLink = pd.merge(evidenceLink,positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="both"').drop('_merge', axis=1).reset_index(drop=True)
        # positive_goldStandardLink = pd.merge(evidenceLink,speciesA_positive_goldStandardLink, on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="both"').drop('_merge', axis=1).reset_index(drop=True)

    if len(positive_goldStandardLink)>0 and len(negative_goldStandardLink)>0:
        positive_goldStandardLink = positive_goldStandardLink[col_name].to_list()
        negative_goldStandardLink = negative_goldStandardLink[col_name].to_list()
        return positive_goldStandardLink,negative_goldStandardLink
    else:
        return None,None

if __name__ == '__main__':
    if os.path.exists('configFiles/exampleGenerate'):
        with open('configFiles/exampleGenerate/instanceConfig.yml') as f:
            instanceConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerate/goldStandardConfig.yml') as f:
            goldStandardConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerate/proteomeConfig.yml') as f:
            proteomeConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerate/evidenceConfig.yml') as f:
            evidenceConfig = yaml.load(f, Loader=SafeLoader)
    print(instanceConfig['instance']['species'])