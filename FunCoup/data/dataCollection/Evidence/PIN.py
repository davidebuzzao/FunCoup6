import os
# import matplotlib.pyplot as plt
import dask
import dask.dataframe as dd
from collections import Counter
from sklearn.preprocessing import MinMaxScaler,PowerTransformer
import yaml
from yaml.loader import SafeLoader
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO
import itertools
from tqdm import tqdm
from auxiliary_functions import *

def disentangle_complex_parallel(name, group):
    complexMembers=len(group['proteinB'].unique())
    identifier_combinations = list(itertools.combinations(group['proteinB'].unique(), 2))
    new_rows = []
    for idA, idB in identifier_combinations:
         # treating experimentalMethod+complexID as separate experiments/assays
        new_rows.append([name[2]+"|"+name[0]+"|"+str(complexMembers), name[1], idA, idB])
    return pd.DataFrame(new_rows)

def compute_frequency(x):
    x_count = Counter(x)
    x_dict = dict([(p,val) for p,val in x_count.items()])
    return x_dict

def calculatePINscore(pin_row,gene_dict,pmid_dict):
    pmid = pin_row['metadata'].split(',')
    pmid_factor = [1/(np.log(pmid_dict[p]+1)) for p in pmid]
    degree_A = gene_dict[pin_row['proteinA_id']]
    degree_B = gene_dict[pin_row['proteinB_id']]
    avg_score = np.sum(pmid_factor) / len(pmid)
    interaction_score = avg_score / np.log(degree_A+degree_B)
    pin_row['score'] = interaction_score
    return pin_row

def getPINScoresForSpecies(linkDf,tax_id,evidenceConfig,proteomeConfig):
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    species_in_db = Species.objects.get(tax_id=tax_id)
    print('Extracting PIN for ' + tax_id)

    linkDf = linkDf.drop('taxA',axis=1).drop('taxB',axis=1)    

    print('# Collecting uniProt identifier or complex id')
    # Collecting uniProt identifier OR complex id (for complexes). Excluding rows without uniProt identifiers
    linkDf['proteinA']=linkDf['idA']+"|"+linkDf['idA_alt1']+"|"+linkDf['idA_alt2'] # Joining all potentail identifier locations for idA
    linkDf = linkDf[(linkDf.proteinA.str.contains("uniprotkb:") | linkDf.proteinA.str.contains("complex:"))] # Keeping rows having uniprot OR complex in idA
    linkDf['proteinA'] = linkDf['proteinA'].apply(lambda x: x.split("uniprotkb:")[1].split("|")[0] if 'complex:' not in x else x) # Extracting the uniprotID from the joined column (except for complexes)
    linkDf.proteinA = np.where(linkDf.idA.str.contains('complex:'), linkDf.idA, linkDf.proteinA) # Setting proteinA to complex ID for complexes
    linkDf['proteinB']=linkDf['idB']+"|"+linkDf['idB_alt1']+"|"+linkDf['idB_alt2'] # Joining all potential identifier locations for idB
    linkDf = linkDf[linkDf.proteinB.str.contains("uniprotkb:")] # Keeping wors having a uniprot id in idB
    linkDf['proteinB'] = linkDf['proteinB'].apply(lambda x: x.split("uniprotkb:")[1].split("|")[0]) # Extracting the uniprotID for proteinB
    linkDf = linkDf.drop(['idA','idA_alt1','idA_alt2','idB','idB_alt1','idB_alt2'],axis=1) # Drop cvolumns that are not needed anymore

    linkDf['experimentalMethod'] = linkDf['experimentalMethod'].apply(lambda x: x.split("(")[0].replace("psi-mi:","").replace("\"","")) # Clean up experimental method
    linkDf.experimentalMethod = np.where(linkDf.experimentalMethod == "-", "MI:0000", linkDf.experimentalMethod) # If experimental method is -, setting it to term for unspecified instead
    linkDf['pmid'] = linkDf['pmid'].apply(lambda x: x.split("|")[-1]) # cleaning up pubmed id
    linkDf = linkDf[linkDf.pmid != "-"] # Only using rows where a pubmed-id is specified

    # Extracting complexes into pairs
    complexDf = linkDf[linkDf.proteinA.str.contains("complex:")] # Extracting complex rows from the dataframe
    linkDf = linkDf[~linkDf.proteinA.str.contains("complex:")] # Removing complex rows from the pairwise dataframe
    groupedComplexes = complexDf.groupby(['proteinA','pmid','experimentalMethod']) # Group complexes by complexID-pmid-experimental method

    print('# Working on complexes')
    with tqdm_joblib(tqdm(desc='disentangle_complex', total=len(groupedComplexes))) as progress_bar:
        complex_df = pd.concat(parallel(delayed(disentangle_complex_parallel)(name,group) for name,group in groupedComplexes))
    complex_df.columns = ['experimentalMethod','pmid','proteinA','proteinB']
    linkDf = pd.concat([linkDf, complex_df]).reset_index(drop=True)

    # Getting mappings from the DB for the specified proteome
    print('# Mapping to protein_id')
    proteome = Proteome.objects.get(Q(version=proteomeConfig['genome']['version']),Q(species=species_in_db))
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    
    available_protein_df.columns = ['proteinA','protein_id']
    linkDf = pd.merge(linkDf, available_protein_df, on='proteinA', how='inner')
    linkDf.drop('proteinA',axis=1,inplace=True)
    linkDf.columns = ['experimentalMethod','pmid','proteinB','proteinA_id']
    available_protein_df.columns = ['proteinB','protein_id']
    linkDf = pd.merge(linkDf, available_protein_df, on='proteinB', how='inner')
    linkDf.drop('proteinB',axis=1,inplace=True)
    linkDf.columns = ['experimentalMethod','pmid','proteinA_id','proteinB_id']
    del available_protein_df
    
    # Put IDs in the correct order
    linkDf.proteinA_id, linkDf.proteinB_id = np.where(linkDf.proteinA_id < linkDf.proteinB_id, [linkDf.proteinB_id, linkDf.proteinA_id], [linkDf.proteinA_id, linkDf.proteinB_id])
    linkDf = linkDf.drop_duplicates() # Drop duplicates
    linkDf = linkDf[linkDf['proteinA_id']!=linkDf['proteinB_id']]
    
    ## Measure gene frequency
    ## How many distinct links? Remove duplicated lines on ['pmid','proteinA_id','proteinB_id'] and extract frequencies
    print('# Measuring gene and pmid frequencies')
    noduplicates_linkDf = linkDf.drop_duplicates(subset=['pmid', 'proteinA_id', 'proteinB_id'])
    genes = noduplicates_linkDf['proteinA_id'].to_list() + noduplicates_linkDf['proteinB_id'].to_list()
    gene_dict = compute_frequency(genes)
    pmids = noduplicates_linkDf['pmid'].to_list()
    pmid_dict = compute_frequency(pmids)

    linkDf.drop('experimentalMethod',axis=1,inplace=True)
    linkDf.columns = ['metadata','proteinA_id','proteinB_id']
    linkDf = linkDf.groupby(['proteinA_id','proteinB_id'])['metadata'].apply(lambda x: ','.join(x)).reset_index()
    
    # Scoring PIN
    print('# Scoring PIN')
    meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int', 'metadata':'object', 'score': 'float'}
    linkDDf = dd.from_pandas(linkDf, npartitions=288)  # Adjust the number of partitions as needed
    result_ddf = linkDDf.apply(lambda row: calculatePINscore(row,gene_dict,pmid_dict), axis=1, meta=meta_dict)
    with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
        linkDf = result_ddf.compute() 

    if len(linkDf)>0:
        scores = linkDf['score'].values.reshape(-1, 1)
        # Power Transformation https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.PowerTransformer.html#sklearn.preprocessing.PowerTransformer
        # Apply a power transform featurewise to make data more Gaussian-like.
        power_scaler = PowerTransformer(method='box-cox')
        power_scores = power_scaler.fit_transform(scores)
        # Min-Max scaling
        minmax_scaler = MinMaxScaler(feature_range=(0, 1))
        minmax_power_scores = minmax_scaler.fit_transform(power_scores)
        linkDf['score'] = np.round(minmax_power_scores,decimals=4)

        evidence_in_db = Evidence(type="PIN", species=species_in_db, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidence_in_db.save()

        linkDf = linkDf.sample(frac=1).reset_index(drop=True)
        linkDf.insert(0,'evidence',evidence_in_db.id)

        return linkDf


def getPINLinks(evidenceConfig, speciesToCollect,proteomeConfig):
    # Collecting PIN links from the iRefIndex file. Using complexes as well as pairwise interactions. Caclulates a score for each link
    # based on the number of papers, the number of assays, and the number of members in each assay supporting the link.
    iRefIndexFile = "data/tmp/"+evidenceConfig['url'].split("/")[-1]
    iRefIndex_df = pd.read_csv(iRefIndexFile, sep='\t',usecols=[0,1,2,3,6,8,9,10,38,39], skiprows=1, header=None, comment='#')
    print('Dataset ' + iRefIndexFile + ' open!')
    iRefIndex_df.columns = ['idA','idB','idA_alt1','idB_alt1','experimentalMethod','pmid','taxA','taxB','idA_alt2','idB_alt2']
    iRefIndex_df.taxA = np.where(iRefIndex_df.taxA == "-", iRefIndex_df.taxB, iRefIndex_df.taxA) # For complexes, taxa is -, then set it to the same as taxb
    iRefIndex_df = iRefIndex_df[iRefIndex_df['taxA'] == iRefIndex_df['taxB']].reset_index(drop=True) # Checking that the two interacting proteins belong to the same species

    ## TODO --> delete 
    # iRefIndex_df = iRefIndex_df[iRefIndex_df['taxA'].str.contains("taxid:"+tax_id+"\(", na=False)]
    # tax_id = '559292'
    for tax_id in speciesToCollect:
        subset_iRefIndex_df = iRefIndex_df[iRefIndex_df['taxA'].str.contains("taxid:"+tax_id+"\(", na=False)]
        links_pin_df = getPINScoresForSpecies(subset_iRefIndex_df,tax_id,evidenceConfig,proteomeConfig)
        if len(links_pin_df)==0: continue
        mem_csv = StringIO()
        links_pin_df.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        print("Writing %i PIN links to DB for %s" %(len(links_pin_df),tax_id))
        # Writing the csv to the evidenceLink table
        with closing(mem_csv) as csv_io:
            PIN.objects.from_csv(csv_io)

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
        with open('configFiles/exampleGenerate/trainingConfig.yml') as f:
            trainingConfig = yaml.load(f, Loader=SafeLoader)
    print(instanceConfig['instance']['species'])
    
    speciesToCollect = instanceConfig['instance']['species']
    evidenceConfig = evidenceConfig['PIN']
    getPINLinks(evidenceConfig, speciesToCollect, proteomeConfig)