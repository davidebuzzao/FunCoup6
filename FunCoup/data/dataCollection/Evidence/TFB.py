import os
from sklearn.preprocessing import MinMaxScaler,PowerTransformer
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
from auxiliary_functions import *

def calculateTFBScore(genea, geneb):
    # Add 1 to the intersection and union to account for the size
    # of the union if the intersection is zero
    # tfb_score = ((len(genea & geneb))**2) / (len(genea | geneb))
    tfb_score = (len(genea & geneb)) / (len(genea | geneb))
    return tfb_score

def calculateTFBScore_InParallel(chunk,numChunks,geneTFBdict):
    # Calculating start and stop index for chunk based on the number of total calculations.
    numGenes=len(geneTFBdict)
    totalCalculations=(numGenes*(numGenes-1))/2
    calculationsPerChunk=round(totalCalculations/numChunks)
    startIndex=calculationsPerChunk*chunk
    stopIndex=(calculationsPerChunk*chunk)+calculationsPerChunk

    print("Calculating TFB scores, chunk : "+str(chunk)+" start: "+str(startIndex)+" stop: "+str(stopIndex))
    dataToAdd=[]
    pairIndex=0
    for index,geneA in enumerate(list(geneTFBdict.keys())[:-1]):
        for geneB in list(geneTFBdict.keys())[index+1:]:
            pairIndex+=1
            if pairIndex>startIndex and pairIndex<=stopIndex:
                if geneA!=geneB:
                    genea=set(geneTFBdict[geneA][0].split(','))
                    geneb=set(geneTFBdict[geneB][0].split(','))
                    score = calculateTFBScore(genea,geneb)
                    metadata = ','.join(genea & geneb)
                    if len(metadata) == 0: metadata = None
                    newRecord=[geneA, geneB, score, metadata]
                    dataToAdd.append(newRecord)

    linkDf = pd.DataFrame.from_records(dataToAdd)
    linkDf.columns = ["proteinA","proteinB","score","metadata"]
    linkDf.proteinA, linkDf.proteinB = np.where(linkDf.proteinA < linkDf.proteinB, [linkDf.proteinB, linkDf.proteinA], [linkDf.proteinA, linkDf.proteinB])
    return linkDf

def calculateTFBScore_ForSpecies(url,evidenceConfig,proteomeConfig):
    tfb_file="data/tmp/"+url.split("/")[-1].split(".gz")[0]
    if not os.path.exists("data/tmp/"+url.split("/")[-1].split(".gz")[0]): # Might already have been downloaded.
        if '.gz' in url: downloadAndUnzipFile(url,format='gz')
        else: downloadFile(url, "data/tmp/"+url.split("/")[-1])
    try:
        ## Read from TFB file only 0,1,8 columns
        gene_tfb_df = pd.read_csv(tfb_file, usecols=[0,1,8], sep='\t')
    except:
        print('Dataset ' + tfb_file + ' not found!')
        raise SystemExit

    print('Dataset ' + tfb_file + ' open!')
    ## Remove duplicated rows, format the data frame
    gene_tfb_df = gene_tfb_df.drop_duplicates()
    gene_tfb_df.columns = ['tf','gene', 'species_name']
    species_name = list(gene_tfb_df['species_name'])[0]
    gene_tfb_df.drop('species_name',axis=1, inplace=True)
    gene_tfb_df['tf'] = gene_tfb_df['tf'].astype(str)
    gene_tfb_df['gene'] = gene_tfb_df['gene'].astype(str)

    species_in_db = Species.objects.get(species_name__contains=species_name)
    tax_id = species_in_db.tax_id
    print('Extracting TFB for ' + tax_id)

    ## Make sure that the identifier can be mapped in the database
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=tax_id))[0]
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    available_protein_df.columns = ['gene','protein_id']
    excluded_due_to_mapping = set(gene_tfb_df['gene']) - set(available_protein_df['gene'])
    gene_tfb_df = pd.merge(gene_tfb_df,available_protein_df, on='gene', how='inner')
    print("Excluded "+str(len(excluded_due_to_mapping))+" genes, as they couldnt be mapped. Calculating TFB pairs scores on "+str(len(set(gene_tfb_df['gene'])))+" genes")
    del available_protein_df

    ## Group together tf for same protein id
    gene_tfb_df = gene_tfb_df.groupby('protein_id')['tf'].apply(','.join).reset_index()

    if len(gene_tfb_df)>0:
        gene_tfb_df = gene_tfb_df.set_index('protein_id').transpose()
        geneTFBdict=gene_tfb_df.to_dict('list') # Turn into a dict for faster operations in score calculations

        cpus=os.cpu_count()
        if len(geneTFBdict)<=cpus:
            cpus=1
        parallel = Parallel(n_jobs=cpus)

        # Save gene pairs in parallell
        linkDf = pd.concat(parallel(delayed(calculateTFBScore_InParallel)(chunk,cpus,geneTFBdict) for chunk in range(0, cpus)))
        
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

        evidence_in_db = Evidence(type="TFB", species=species_in_db, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidence_in_db.save()

        linkDf = linkDf.sample(frac=1).reset_index(drop=True)
        linkDf.insert(0,'evidence',evidence_in_db.id)

        # Writing the csv to the evidenceLink table
        print("Done with "+str(len(linkDf.index))+" TFB links for "+str(tax_id))
        return linkDf

def get_TFBLinks(evidenceConfig,proteomeConfig):
    for s,url in zip(evidenceConfig['species'],evidenceConfig['url']):
        linkDf = calculateTFBScore_ForSpecies(url,evidenceConfig,proteomeConfig)
        if len(linkDf)==0: continue
        mem_csv = StringIO()
        linkDf.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        print("Writing %i TFB links to DB for %s" %(len(linkDf),s))
        # Writing the csv to the evidenceLink table
        with closing(mem_csv) as csv_io:
            TFB.objects.from_csv(csv_io)
