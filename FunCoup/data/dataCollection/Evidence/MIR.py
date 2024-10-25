import os
from sklearn.preprocessing import MinMaxScaler,PowerTransformer
from scipy.stats import gmean
import yaml
from yaml.loader import SafeLoader
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
import math
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO

def calculateMIRscore(miRNAsA, miRNAsB):
    # Add 1 to the intersection and union to account for the size
    # of the union if the intersection is zero
    intersection = len(miRNAsA & miRNAsB) 
    union = len(miRNAsA | miRNAsB) 
    return intersection / union

def calculateMIRscore_InParallel(chunk,numChunks,geneMIRdict):

    # Calculating start and stop index for chunk based on the number of total calculations.
    numGenes=len(geneMIRdict)
    totalCalculations=(numGenes*(numGenes-1))/2
    calculationsPerChunk=round(totalCalculations/numChunks)
    startIndex=calculationsPerChunk*chunk
    stopIndex=(calculationsPerChunk*chunk)+calculationsPerChunk

    dataToAdd=[]
    pairIndex=0
    for index,geneA in enumerate(list(geneMIRdict.keys())[:-1]):
        for geneB in list(geneMIRdict.keys())[index+1:]:
            pairIndex+=1
            if pairIndex>startIndex and pairIndex<=stopIndex:
                if geneA!=geneB:
                    miRNAsA=set(geneMIRdict[geneA][0].split(','))
                    miRNAsB=set(geneMIRdict[geneB][0].split(','))
                    score = calculateMIRscore(miRNAsA,miRNAsB)
                    metadata = ','.join(miRNAsA & miRNAsB)
                    if len(metadata) == 0: metadata = None
                    newRecord=[geneA, geneB, score, metadata]
                    dataToAdd.append(newRecord)

    linkDf = pd.DataFrame.from_records(dataToAdd)
    linkDf.columns = ["proteinA","proteinB","score","metadata"]
    linkDf.proteinA, linkDf.proteinB = np.where(linkDf.proteinA < linkDf.proteinB, [linkDf.proteinB, linkDf.proteinA], [linkDf.proteinA, linkDf.proteinB])
    return linkDf

def calculateMIRscore_ForSpecies(mir_file,tax_id,evidenceConfig,proteomeConfig):
    try:
        ## Read from MIR xlsx file only 1,4,5 columns (it can take a while)
        gene_mir_df = pd.read_csv('data/tmp/' + mir_file,usecols=[0,2],sep='\t')
    except:
        print('Dataset ' + mir_file + ' not found!')
        return None

    print('Dataset ' + mir_file + ' open!')
    ## Remove duplicated rows, format the data frame
    gene_mir_df = gene_mir_df.drop_duplicates().reset_index(drop=True)
    gene_mir_df.columns = ['mirna','gene']
    gene_mir_df['gene'] = gene_mir_df['gene'].astype(str)

    print('Extracting MIR for ' + tax_id)
    ## Subset the MIR dataset
    ## Make sure that the identifier can be mapped in the database
    species = Species.objects.get(Q(tax_id=tax_id))
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species=species))[0]
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    available_protein_df.columns = ['gene','protein_id']
    excluded_due_to_mapping = set(gene_mir_df['gene']) - set(available_protein_df['gene'])
    gene_mir_df = pd.merge(gene_mir_df,available_protein_df, on='gene', how='inner')
    print("Excluded "+str(len(excluded_due_to_mapping))+" genes, as they couldnt be mapped. Calculating MIR pairs scores on "+str(len(set(gene_mir_df['gene'])))+" genes for "+str(tax_id))
    del available_protein_df
    gene_mir_df.drop('gene',axis=1,inplace=True)
    ## Group together tissue names for same protein id
    gene_mir_df[['protein_id','mirna']].groupby(['protein_id']).size().reset_index(name='Counts').sort_values(by='Counts', ascending=False).reset_index(drop=True)
    gene_mir_df = gene_mir_df.groupby('protein_id')['mirna'].apply(','.join).reset_index()

    if len(gene_mir_df.index)>1:
        geneMIRdict=gene_mir_df.set_index('protein_id').transpose().to_dict('list') # Turn into a dict for faster operations in score calculations
        cpus=os.cpu_count()
        if len(geneMIRdict)<=cpus:
            cpus=1

        parallel = Parallel(n_jobs=cpus)
        # Score gene pairs in parallell
        linkDf = pd.concat(parallel(delayed(calculateMIRscore_InParallel)(chunk,cpus,geneMIRdict) for chunk in range(0, cpus)))
        ## Drop all zero values
        linkDf = linkDf[linkDf['score']>0].reset_index(drop=True)

        scores = linkDf['score'].values.reshape(-1, 1)
        # Power Transformation https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.PowerTransformer.html#sklearn.preprocessing.PowerTransformer
        # Apply a power transform featurewise to make data more Gaussian-like.
        power_scaler = PowerTransformer(method='box-cox')
        power_scores = power_scaler.fit_transform(scores)
        # Min-Max scaling
        minmax_scaler = MinMaxScaler(feature_range=(0, 1))
        minmax_power_scores = minmax_scaler.fit_transform(power_scores)
        linkDf['score'] = np.round(minmax_power_scores,decimals=4)

        species_in_db = Species.objects.get(tax_id = tax_id)
        evidence_in_db = Evidence(type="MIR", species=species_in_db, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidence_in_db.save()

        linkDf = linkDf.sample(frac=1).reset_index(drop=True)
        linkDf.insert(0,'evidence',evidence_in_db.id)

        mem_csv = StringIO()
        linkDf.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        # Writing the csv to the evidenceLink table
        print("Writing "+str(len(linkDf.index))+" MIR links for "+str(tax_id))
        with closing(mem_csv) as csv_io:
            MIR.objects.from_csv(csv_io)

def get_MIRLinks(evidenceConfig,proteomeConfig):
    ## Distribute CPUs depending on number of genes
    for mir_file,tax_id in zip(evidenceConfig['files'],evidenceConfig['species']):
        calculateMIRscore_ForSpecies(mir_file,tax_id,evidenceConfig,proteomeConfig)

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
    print(instanceConfig['instance']['species'])
    evidenceConfig = evidenceConfig['MIR']
    tax_id='9606'
    mir_file = 'data/tmp/human_predictions_S_C_aug2010.txt'