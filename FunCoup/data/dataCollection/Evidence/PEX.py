import os
import time 
import pandas as pd
import numpy as np
import yaml
from yaml.loader import SafeLoader
from sklearn.preprocessing import MinMaxScaler,PowerTransformer
from contextlib import closing
from io import StringIO
import operator
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed

############################
############################
### Jaccard Index
def pre_processPaxDb(evidenceConfig,available_species):

    raw_links=[]
    for species in available_species:
        for filein in os.listdir("data/tmp/"+ evidenceConfig['url'].split("/")[-1].split(".zip")[0] + "/" + species):
            if("integrated" not in filein):
                organ=""
                filename=""
                geneScoreDict={}
                sortedDict={}
                top75=1
                ## Parse raw QMS file
                with open("data/tmp/"+evidenceConfig['url'].split("/")[-1].split(".zip")[0]+"/"+species+"/"+filein, 'r') as in_file:
                    for line in in_file:
                        line = line.rstrip().split()
                        if line[0].startswith("#organ"):
                            organ=line[1]
                        elif line[0].startswith("#filename"):
                            filename=line[1].split("-")[1].split(".")[0]
                        elif not line[0].startswith("#") and organ not in ["WHOLE_ORGANISM","CELL_LINE"]:
                            geneScoreDict[line[0].split(".")[1]]=float(line[1])
                ## Sort tissues by decreasing protein expression, for each protein pick tissues with top expression.
                if len(geneScoreDict)>0:
                    sortedDict = sorted(geneScoreDict.items(), key=operator.itemgetter(1), reverse=True)
                    top75=int(round(len(sortedDict)*0.25))
                    counter=0
                    for s in sortedDict:
                        counter+=1
                        if counter <=top75:
                            raw_links.append([organ,s[0],species])

    raw_links_df = pd.DataFrame(raw_links)
    return raw_links_df

def calculateJaccardIndex(genea, geneb):
    # Add 1 to the intersection and union to account for the size
    # of the union if the intersection is zero
    # qms_score = ((len(genea & geneb))**2) / (len(genea | geneb))
    qms_score = (len(genea & geneb)) / (len(genea | geneb))
    return qms_score

def calculateJaccardIndex_InParallel(chunk,numChunks,geneQMSdict):
    # Calculating start and stop index for chunk based on the number of total calculations.
    numGenes=len(geneQMSdict)
    totalCalculations=(numGenes*(numGenes-1))/2
    calculationsPerChunk=round(totalCalculations/numChunks)
    startIndex=calculationsPerChunk*chunk
    stopIndex=(calculationsPerChunk*chunk)+calculationsPerChunk

    print("Calculating QMS scores, chunk : "+str(chunk)+" start: "+str(startIndex)+" stop: "+str(stopIndex))
    dataToAdd=[]
    pairIndex=0
    for index,geneA in enumerate(list(geneQMSdict.keys())[:-1]):
        for geneB in list(geneQMSdict.keys())[index+1:]:
            pairIndex+=1
            if pairIndex>startIndex and pairIndex<=stopIndex:
                if geneA!=geneB:
                    genea=set(geneQMSdict[geneA][0].split(','))
                    geneb=set(geneQMSdict[geneB][0].split(','))
                    score = calculateJaccardIndex(genea,geneb)
                    metadata = ','.join(genea & geneb)
                    if len(metadata) == 0: metadata = None
                    newRecord=[geneA, geneB, score, metadata]
                    dataToAdd.append(newRecord)

    evidenceLink = pd.DataFrame.from_records(dataToAdd)
    evidenceLink.columns = ["proteinA","proteinB","score","metadata"]
    evidenceLink.proteinA, evidenceLink.proteinB = np.where(evidenceLink.proteinA < evidenceLink.proteinB, [evidenceLink.proteinB, evidenceLink.proteinA], [evidenceLink.proteinA, evidenceLink.proteinB])
    return evidenceLink

def calculateJaccardIndex_ForSpecies(CPUSforWithinSpeciesCalc,gene_qms_df,tax_id,proteomeConfig,evidenceConfig):
    print('Extracting QMS for ' + tax_id)
    species_in_db = Species.objects.get(tax_id = tax_id)

    if len(gene_qms_df.index)>1:
        ## Make sure that the identifier can be mapped in the database
        gene_qms_df['gene'] = gene_qms_df['gene'].str.replace(rf'{tax_id}\.', '', regex=True)
        proteome = Proteome.objects.get(Q(version=proteomeConfig['genome']['version']),Q(species=species_in_db))
        available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
        available_protein_df.columns = ['gene','protein_id']
        available_protein_df['gene'] = available_protein_df['gene'].str.replace(r'\.\d+', '', regex=True)
        excluded_due_to_mapping = set(gene_qms_df['gene']) - set(available_protein_df['gene'])
        gene_qms_df = pd.merge(gene_qms_df,available_protein_df, on='gene', how='inner')
        print("Excluded "+str(len(excluded_due_to_mapping))+" genes, as they couldnt be mapped. Calculating QMS pairs scores on "+str(len(set(gene_qms_df['gene'])))+" genes")
        del available_protein_df

        ## Group together tissue names for same protein id
        gene_qms_df = gene_qms_df.groupby('protein_id')['tissue'].apply(','.join).reset_index()
        #print(gene_qms_df)

        if len(gene_qms_df.index)>1:
            gene_qms_df = gene_qms_df.set_index('protein_id').transpose()
            geneQMSdict=gene_qms_df.to_dict('list') # Turn into a dict for faster operations in score calculations

            if len(geneQMSdict)<=CPUSforWithinSpeciesCalc:
                CPUSforWithinSpeciesCalc=1
            # Score gene pairs in parallell
            parallel = Parallel(n_jobs=CPUSforWithinSpeciesCalc)
            evidenceLink = pd.concat(parallel(delayed(calculateJaccardIndex_InParallel)(chunk,CPUSforWithinSpeciesCalc,geneQMSdict) for chunk in range(0, CPUSforWithinSpeciesCalc)))  
            
            ## Drop all zero values
            evidenceLink = evidenceLink[evidenceLink['score']>0].reset_index(drop=True)

            scores = evidenceLink['score'].values.reshape(-1, 1)
            # Power Transformation https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.PowerTransformer.html#sklearn.preprocessing.PowerTransformer
            # Apply a power transform featurewise to make data more Gaussian-like.
            power_scaler = PowerTransformer(method='box-cox')
            power_scores = power_scaler.fit_transform(scores)
            # Min-Max scaling
            minmax_scaler = MinMaxScaler(feature_range=(0, 1))
            minmax_power_scores = minmax_scaler.fit_transform(power_scores)
            evidenceLink['score'] = np.round(minmax_power_scores,decimals=4)

            species_in_db = Species.objects.get(tax_id = tax_id)
            evidence_in_db = Evidence(type="PEX", species=species_in_db, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
            evidence_in_db.save()

            evidenceLink = evidenceLink.sample(frac=1).reset_index(drop=True)
            evidenceLink.insert(0,'evidence',evidence_in_db.id)

            mem_csv = StringIO()
            evidenceLink.to_csv(mem_csv, index=False)
            mem_csv.seek(0)
            # Writing the csv to the evidenceLink table
            print("Writing %i PEX links w/ method=%s for %s" %(len(evidenceLink),evidenceConfig['scoring_method'],str(tax_id)))
            with closing(mem_csv) as csv_io:
                PEX.objects.from_csv(csv_io)

def get_JaccardIndexLinks(evidenceConfig,proteomeConfig,available_species):
    qms_file = "data/tmp/"+evidenceConfig['url'].split("/")[-1].split(".zip")[0]
    try:
        all_gene_qms_df = pre_processPaxDb(evidenceConfig,available_species)
    except:
        print('Dataset ' + qms_file + ' not found!')
        raise SystemExit

    print('Dataset ' + qms_file + ' open!')
    ## Keep only data for species under study
    if len(all_gene_qms_df.index)>1: # No need to continue if dataframe is empty. Happens for some species that only have whole-organism info
        all_gene_qms_df = all_gene_qms_df.drop_duplicates().reset_index(drop=True)
        all_gene_qms_df.columns = ['tissue','gene','tax_id']
        all_gene_qms_df['gene'] = all_gene_qms_df['gene'].astype(str)
        all_gene_qms_df['tissue'] = all_gene_qms_df['tissue'].str.capitalize()
        
        # Group the data frame by "tax_id" and count the unique "tissue" values
        grouped = all_gene_qms_df.groupby("tax_id")["tissue"].nunique()
        # Filter the original data frame based on the selected tax_id
        all_gene_qms_df = all_gene_qms_df[all_gene_qms_df["tax_id"].isin(grouped[grouped > 3].index)].reset_index(drop=True)
        gene_counts = all_gene_qms_df[['gene','tax_id']].drop_duplicates().groupby(['tax_id']).size().reset_index(name='Counts').sort_values(by='Counts', ascending=False).reset_index(drop=True)

        ## Distribute CPUs depending on number of genes
        cpus=os.cpu_count()
        total_counts = gene_counts['Counts'].sum()
        # Calculate the proportion for each species
        gene_counts['Proportion'] = gene_counts['Counts'] / total_counts
        # Calculate the minimum cores (at least 1) for each species
        gene_counts['Cores'] = 1
        remaining_cpus = cpus - gene_counts['Cores'].sum()
        gene_counts['Cores'] = gene_counts['Counts'].apply(lambda x: max(1, round(remaining_cpus * x / total_counts)))
        # Calculate the remaining cores
        remaining_cpus = cpus - gene_counts['Cores'].sum()
        # Distribute the remaining cores to the species with the highest proportions
        top_species = gene_counts.nlargest(2, 'Proportion')
        gene_counts.loc[top_species.index, 'Cores'] += remaining_cpus//2
        gene_counts.drop(['Proportion','Counts'], axis=1, inplace=True)
        gene_counts_dict = gene_counts.set_index('tax_id').transpose().to_dict('list')

        Parallel(n_jobs=1)(delayed(calculateJaccardIndex_ForSpecies)(gene_counts_dict[tax_id][0],all_gene_qms_df[all_gene_qms_df['tax_id']==tax_id].reset_index(drop=True).drop('tax_id',axis=1),tax_id,proteomeConfig,evidenceConfig) for tax_id in evidenceConfig['species'])

############################
############################
### Spearman Correlation
def extract_PaxDbOrgan(filein):
    with open(filein, 'r') as in_file:
        for line in in_file:
            line = line.rstrip().split()
            if line[0].startswith("#organ"):
                return line[1]

def compute_correlation(data):
    ## Transpose the dataset
    data = data.T
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

def calculateCorrelation_ForSpecies(tax_id,evidenceConfig,proteomeConfig):
    print('Extracting PEX for ' + tax_id)
    species_in_db = Species.objects.get(tax_id = tax_id)

    indir = "data/tmp/%s/%s/" %(evidenceConfig['url'].split("/")[-1].split(".zip")[0], tax_id)
    file_list = os.listdir(indir)
    file_dictionary = dict([(i,[file_list[i],extract_PaxDbOrgan(indir + file_list[i])]) for i in range(len(file_list))])
    
    evidenceLink = pd.DataFrame()
    for col_name,file_info in file_dictionary.items():
        filein = file_info[0]
        if("integrated" not in filein):
            organ = file_info[1]
            if organ not in ["WHOLE_ORGANISM","CELL_LINE"]:
                tmp_evidenceLink = pd.read_csv(indir + filein, sep='\t', comment='#', usecols=[0,1],names=['gene',col_name], dtype={'gene':str,col_name:float})
                if len(tmp_evidenceLink)==0: continue
                if evidenceLink.empty: evidenceLink = tmp_evidenceLink.copy()
                else: evidenceLink = evidenceLink.merge(tmp_evidenceLink, on=['gene'], how='outer')

    if len(evidenceLink)>1:
        ## Make sure that the identifier can be mapped in the database
        evidenceLink['gene'] = evidenceLink['gene'].str.replace(rf'{tax_id}\.', '', regex=True)
        proteome = Proteome.objects.get(Q(version=proteomeConfig['genome']['version']),Q(species=species_in_db))
        available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
        available_protein_df.columns = ['gene','protein_id']
        available_protein_df['gene'] = available_protein_df['gene'].str.replace(r'\.\d+', '', regex=True)
        evidenceLink = pd.merge(evidenceLink,available_protein_df, on='gene', how='inner')
        del available_protein_df
        
        evidenceLink.drop('gene',inplace=True,axis=1)
        evidenceLink = evidenceLink.groupby('protein_id').mean(numeric_only=True)
        evidenceLink = evidenceLink.sort_index(ascending=False).rename_axis(None).rename_axis(None, axis=0)
        
        startTime = time.time()
        evidenceLink = compute_correlation(evidenceLink)
        endTime = time.time()
        print('Time: %.2f' %(endTime-startTime))

        if len(evidenceLink)>0:
            scores = evidenceLink['score'].values.reshape(-1, 1)
            min_score, max_score = str(np.round(np.min(scores), decimals=4)), str(np.round(np.max(scores), decimals=4))
            # Min-Max scaling
            minmax_scaler = MinMaxScaler(feature_range=(0, 1))
            minmax_scores = minmax_scaler.fit_transform(scores)
            evidenceLink['score'] = np.round(minmax_scores,decimals=4)
            # filein = 'PEX(min5)_%s.gz' %tax_id
            # evidenceLink.to_csv(filein, compression='gzip', sep='\t', index=False)
            evidence_in_db = Evidence(type="PEX", species=species_in_db, scoreRange='%s|%s' %(min_score,max_score), version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
            evidence_in_db.save()

            evidenceLink = evidenceLink.sample(frac=1).reset_index(drop=True)
            evidenceLink.insert(0,'evidence',evidence_in_db.id)
            
            mem_csv = StringIO()
            evidenceLink.to_csv(mem_csv, index=False)
            mem_csv.seek(0)
            # Writing the csv to the evidenceLink table
            print("Writing %i PEX links w/ method=%s for %s" %(len(evidenceLink),evidenceConfig['scoring_method'][0],str(tax_id)))
            with closing(mem_csv) as csv_io:
                PEX.objects.from_csv(csv_io)

def getCorrelationLinks(evidenceConfig,proteomeConfig,available_species):
    cpus=np.min([os.cpu_count(),len(available_species)])
    Parallel(n_jobs=cpus)(delayed(calculateCorrelation_ForSpecies)(tax_id,evidenceConfig,proteomeConfig) for tax_id in available_species)

############################
############################
### Main function
def getPEXLinks(evidenceConfig,proteomeConfig):
    ## Work only on data for species under study
    available_species = set(os.listdir("data/tmp/"+evidenceConfig['url'].split("/")[-1].split(".zip")[0]))
    available_species = available_species & set(evidenceConfig['species'])

    if 'jaccard' in evidenceConfig['scoring_method'].lower():
        ### Jaccard Index
        print('# Computing PEX with Jaccard Index')
        get_JaccardIndexLinks(evidenceConfig,proteomeConfig,available_species)

    elif 'correlation' in evidenceConfig['scoring_method'].lower():
        ### Spearman Correlation
        print('# Computing PEX with Spearman Correlation')
        getCorrelationLinks(evidenceConfig,proteomeConfig,available_species)


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
    evidenceConfig = evidenceConfig['PEX']
    tax_id='9606'
