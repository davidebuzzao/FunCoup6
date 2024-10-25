import os
from sklearn.preprocessing import MinMaxScaler
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

def pre_processGIN(evidenceConfig):
    raw_links_df=pd.DataFrame()
    for tax_id,species_name in zip(evidenceConfig['species'],evidenceConfig['species_name']):
        gin_file = 'data/tmp/BIOGRID-ORGANISM-%s-%s.mitab.txt' %(species_name, str(evidenceConfig['version'].split('_v')[1]))
        try:
            ## Read from GIN file only 2,3,6,14 columns
            all_gene_gin_df = pd.read_csv(gin_file, usecols=[2,3,6,14], sep='\t') # 7 for metadata
        except:
            print('Dataset ' + gin_file + ' not found!')
            raise SystemExit
        print('Dataset ' + gin_file + ' open!')
        all_gene_gin_df.columns = ['genea','geneb','interaction','score']
        ## Genes can be extracted as Gene_ID, Gene_Symbol or UniprotID.
        ## Different number of genes will be mapped accordingly.
        # all_gene_gin_df['genea'] = all_gene_gin_df['genea'].str.replace('entrez gene/locuslink:', '').astype(str)
        # all_gene_gin_df['geneb'] = all_gene_gin_df['geneb'].str.replace('entrez gene/locuslink:', '').astype(str)
        # all_gene_gin_df['genea'] = all_gene_gin_df['genea'].apply(lambda x: x.split('gene/locuslink:')[1].split('|')[0].upper())
        # all_gene_gin_df['geneb'] = all_gene_gin_df['geneb'].apply(lambda x: x.split('gene/locuslink:')[1].split('|')[0].upper())

        ## Select genetic interactions for which a score was precomputed. Use protein UniprotID (i.e. swiss-prot).
        all_gene_gin_df = all_gene_gin_df[(all_gene_gin_df['interaction'].str.contains('genetic')) \
                                          & (all_gene_gin_df['score']!='-') \
                                            & (all_gene_gin_df['genea'].str.contains('swiss-prot') \
                                                & (all_gene_gin_df['geneb'].str.contains('swiss-prot')))].reset_index(drop=True)
        all_gene_gin_df['genea'] = all_gene_gin_df['genea'].apply(lambda x: x.split('swiss-prot:')[1].split('|')[0].upper())
        all_gene_gin_df['geneb'] = all_gene_gin_df['geneb'].apply(lambda x: x.split('swiss-prot:')[1].split('|')[0].upper())

        all_gene_gin_df['score'] = all_gene_gin_df['score'].str.replace('score:', '').astype(float)
        all_gene_gin_df.drop('interaction',axis=1, inplace=True)
        all_gene_gin_df['tax_id'] = tax_id
        raw_links_df = pd.concat([raw_links_df,all_gene_gin_df], ignore_index=True)

    return raw_links_df

def compute_correlation(data):
    ## Compute correlations for pairs expressed in min 3 tissues,
    ## then extract upper matrix only and exclude diagonal.
    ## As the indices of data matrix are sorted in decreasing order, the stack()
    ## function will return a data frame with sorted proteins (i.e. A>B).
    r_df = data.corr(method="spearman",min_periods=5,numeric_only=True)\
        .where(np.triu(np.ones((data.shape[1],data.shape[1])), k=1).astype(bool))\
            .stack(dropna=True)\
                .reset_index()
    r_df.columns = ['proteinA_id','proteinB_id','score']
    return r_df

def measure_intersection(geneA,geneB):
    print(len(set(geneA)&set(geneB)))

def get_GINLinksForSpecies(linkDf, tax_id, evidenceConfig, proteomeConfig):
    print('Extracting GIN for ' + tax_id)

    ## Make sure that the identifier can be mapped in the database
    species_in_db = Species.objects.get(tax_id=tax_id)
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species=species_in_db))[0]
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    available_protein_df.columns = ['gene','protein_id']
    excluded_due_to_mapping = (set(linkDf['genea']) | set(linkDf['geneb'])) - set(available_protein_df['gene'])
    print("Excluded "+str(len(excluded_due_to_mapping))+" genes, as they couldnt be mapped. Extracting GIN pairs scores for "+str(len(set(linkDf['genea']) | set(linkDf['geneb'])))+" genes")

    linkDf.columns = ['gene','geneb','score']
    linkDf = pd.merge(linkDf,available_protein_df, on='gene', how='inner')
    linkDf.drop('gene',axis=1, inplace=True)
    linkDf.columns = ['gene','score','proteinA']
    linkDf = pd.merge(linkDf,available_protein_df, on='gene', how='inner')
    linkDf.drop('gene',axis=1, inplace=True)
    linkDf.columns = ['score','proteinA_id','proteinB_id']
    del available_protein_df
    # Sort identifiers in the correct tax_id
    linkDf.proteinA_id, linkDf.proteinB_id = np.where(linkDf.proteinA_id < linkDf.proteinB_id, [linkDf.proteinB_id, linkDf.proteinA_id], [linkDf.proteinA_id, linkDf.proteinB_id]) 
    linkDf = linkDf.groupby(['proteinA_id','proteinB_id']).mean(numeric_only=True).reset_index()

    # Create a symmetric matrix with both 'proteinA_id' and 'proteinB_id' as rows and columns
    unique_proteins = sorted(set(linkDf['proteinA_id'].to_list()+linkDf['proteinB_id'].to_list()), reverse=True)
    link_matrix = pd.DataFrame(index=unique_proteins, columns=unique_proteins)
    # Fill in the symmetric matrix with the 'score' values
    for index, row in linkDf.iterrows():
        protein_a = row['proteinA_id']
        protein_b = row['proteinB_id']
        score = row['score']
        link_matrix.at[protein_a, protein_b] = score
        link_matrix.at[protein_b, protein_a] = score
    
    ## Measure correlations
    link_matrix = link_matrix.astype(float)
    r_df = compute_correlation(link_matrix)

    scores = r_df['score'].values.reshape(-1, 1)
    min_score, max_score = str(np.round(np.min(scores), decimals=4)), str(np.round(np.max(scores), decimals=4))
    # MinMax scaling
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaled_scores = np.round(scaler.fit_transform(scores),decimals=4)
    r_df['score'] = scaled_scores

    evidence_in_db = Evidence(type='GIN', species=species_in_db, scoreRange='%s|%s' %(min_score,max_score), version=evidenceConfig['version'], scoringMethod = evidenceConfig['scoring_method'])
    evidence_in_db.save()
    
    r_df = r_df.sample(frac=1).reset_index(drop=True)
    r_df['evidence'] = evidence_in_db.id

    mem_csv = StringIO()
    r_df.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing "+str(len(r_df.index))+" GIN links to DB for "+tax_id)
    with closing(mem_csv) as csv_io:
        GIN.objects.from_csv(csv_io)

def get_GINLinks(evidenceConfig,proteomeConfig):
    all_gene_gin_df = pre_processGIN(evidenceConfig)
    all_gene_gin_df = all_gene_gin_df.drop_duplicates()
    all_gene_gin_df.columns = ['genea','geneb','score','tax_id']
    all_gene_gin_df[['tax_id','genea']].groupby('tax_id').count()
    cpus=np.min([os.cpu_count(),len(evidenceConfig['species'])])
    Parallel(n_jobs=cpus)(delayed(get_GINLinksForSpecies)(all_gene_gin_df[all_gene_gin_df['tax_id']==tax_id].reset_index(drop=True).drop('tax_id',axis=1), tax_id, evidenceConfig, proteomeConfig) for tax_id in set(all_gene_gin_df['tax_id'].to_list()))

if __name__ == '__main__':
    if os.path.exists('configFiles/exampleGenerateEvidenceSmall'):
        with open('configFiles/exampleGenerateEvidenceSmall/instanceConfig.yml') as f:
            instanceConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/goldStandardConfig.yml') as f:
            goldStandardConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/proteomeConfig.yml') as f:
            proteomeConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/evidenceConfig.yml') as f:
            evidenceConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/trainingConfig.yml') as f:
            trainingConfig = yaml.load(f, Loader=SafeLoader)
    print(instanceConfig['instance']['species'])
    evidenceConfig = evidenceConfig['GIN']
    instanceConfig = instanceConfig['instance']
    tax_id='9606'
    species_name='Homo sapiens'