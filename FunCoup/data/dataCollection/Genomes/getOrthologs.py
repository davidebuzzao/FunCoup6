import requests
import os
import pandas as pd
import yaml
from yaml.loader import SafeLoader
from joblib import Parallel, delayed
from contextlib import closing
from io import StringIO
from tqdm import tqdm
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from django.conf import settings
from auxiliary_functions import *

def getSpeciesOrder(speciesA, speciesB): # Always gets species in the same order
    if int(speciesA)>int(speciesB):
        return speciesA,speciesB
    else:
        return speciesB,speciesA

def getProteinIdFromDB(p):
    proteinIDInDB=[]
    try: # First, checking if protein is a primary identifier
        proteinInDB=Protein.objects.get(uniprot_id = p)
        proteinIDInDB.append(proteinInDB.id)
    except Protein.DoesNotExist: # If not, checking for a mapping
        mappedID = IdMapping.objects.filter(Q(mapped_id=p))
        if mappedID.exists() :
            for mapping in mappedID:
                if mapping.protein.id not in proteinIDInDB:
                    proteinIDInDB.append(mapping.protein.id)
    return proteinIDInDB

def getPairsFromGroup(spAGenes,spBGenes,spA,spB):
    orthologsToWrite=[]
    for gA in spAGenes:
        for gB in spBGenes:
            orth=gA+"|"+gB+"|"+spA+"|"+spB
            orthologsToWrite.append(orth)
    return orthologsToWrite

def readSqlTable(content,spA,spB):
    currentGroup=0
    spAGenes=[]
    spBGenes=[]
    orthologsToWrite=[]
    proteinsUsed=[]
    for line in content:
        if len(line)>1:
            groupId,bitscore,species,inparalogScore,proteinId,*seedscore=line.strip().split("\t")
            proteinsUsed.append(proteinId)
            if groupId!=currentGroup: # For each new group encountered, writing all pairs of orthologs to list
                orthologsToWrite.extend(getPairsFromGroup(spAGenes,spBGenes,spA,spB))
                currentGroup=groupId
                spAGenes=[]
                spBGenes=[]
            if species.split(".")[0]==spA:
                spAGenes.append(proteinId.strip())
            else:
                spBGenes.append(proteinId.strip())
    orthologsToWrite.extend(getPairsFromGroup(spAGenes,spBGenes,spA,spB))
    return orthologsToWrite,proteinsUsed

def getOrthologsFromFile(sp):
    spA,spB=sp.split("|")
    # Checking if file is available as species AB or BA
    if os.path.isfile("data/tmp/resultsInParanoid/SQLtable."+spA+".fa-"+spB+".fa"):
        sqltableFile=open("data/tmp/resultsInParanoid/SQLtable."+spA+".fa-"+spB+".fa","r")
    elif os.path.isfile("data/tmp/resultsInParanoid/SQLtable."+spB+".fa-"+spA+".fa"):
        sqltableFile=open("data/tmp/resultsInParanoid/SQLtable."+spB+".fa-"+spA+".fa","r")
    else:
        print("Cant find file: "+"data/tmp/resultsInParanoid/SQLtable."+spA+".fa-"+spB+".fa")
        exit()
    orthologsToWrite,proteinsUsed = readSqlTable(sqltableFile.readlines(),spA,spB)
    return orthologsToWrite,proteinsUsed

def getOrthologsFromUrl(sp,url):
    spA,spB=sp.split("|")
    orthologsRequest = requests.get(url.replace("TAXIDA",spA).replace("TAXIDB",spB))
    if orthologsRequest.status_code==200: # If inparanoidb has orthologs for this protein, returns a list of all ortholog groups
        orthologsToWrite,proteinsUsed = readSqlTable(orthologsRequest.content.decode("utf-8").split("\n"),spA,spB)
    elif orthologsRequest.status_code==500:
        orthologsRequest = requests.get(url.replace("TAXIDA",spB).replace("TAXIDB",spA))
        if orthologsRequest.status_code==200: # If inparanoidb has orthologs for this protein, returns a list of all ortholog groups
            orthologsToWrite, proteinsUsed = readSqlTable(orthologsRequest.content.decode("utf-8").split("\n"),spA,spB)
        else:
            print("Cant download file from: "+url.replace("TAXIDA",spB).replace("TAXIDB",spA)+" error code: "+str(orthologsRequest.status_code))
    else: # if error code is 500, it just means that the protein is not in inparanoiDB
        print("Cant download file from: "+url.replace("TAXIDA",spA).replace("TAXIDB",spB)+" error code: "+str(orthologsRequest.status_code))
        exit()
    return orthologsToWrite,proteinsUsed

def createOrthologEntry_parallel(uniprot_dict,chunk_orthologPairs):
    # start,end = 0,10; chunk_orthologPairs = orthologPairs[start:end]
    # Writing all collected orthologs, using lookup map to get protein ids
    orthologsToWrite=[]
    for o in chunk_orthologPairs:
        pA,pB,spA,spB=o.split("|")
        if pA in uniprot_dict:
            for idA in uniprot_dict[pA]:
                if pB in uniprot_dict:
                    for idB in uniprot_dict[pB]:
                        orth = [idA, idB, spA, spB]
                        orthologsToWrite.append(orth)
    return pd.DataFrame(orthologsToWrite,columns=['proteinA_id','proteinB_id','speciesA','speciesB'])


def updateOrthologs(species,proteomeConfig):
    # Getting orthologs for all species included in the isntance config, unless orthologs for the species
    # and version already exist in the DB. Has the opportunity to fetch orthologs either from the
    # InParanoidb website OR from a collection of all SQLtable files for the species (in this case, the
    # sqltable files must be available under data/tmp/resultsInParanoid). creates pairs of orthologs
    # From the ortholog groups, so that iformatioon can be transferred to both seed orthologs and
    # inparalogs. Saves ortholg relations in he DB where taxIDA is always > taxIDB

    species_to_update = ["9615","9031","7955"]

    proteomeVersion=proteomeConfig['genome']['version']
    orthologVersion=proteomeConfig['ortholog']['version']
    # Checking which species pairs already has orthologs in the DB
    speciesPairsToGet=set()
    for sA in species_to_update:
        for sB in species:
            if sA!=sB:
                speciesA,speciesB=getSpeciesOrder(sA,sB) # Always puts species in the correct order
                speciesPairsToGet.add(speciesA+"|"+speciesB)

    if len(speciesPairsToGet)>0:
        print("Getting Orthologs")
        ## STEP 1 --> Fast
        orthologPairs=[] # ["proteinA|proteinB|speciesA|speciesB",...]
        proteinsUsed=[] # ["proteinA","proteinB"]
        for sp in speciesPairsToGet:
            if proteomeConfig['ortholog']['orthologCollectionMethod']=="url": # If collecting from url
                orthologsToWrite,proteinsUsedInFile = getOrthologsFromUrl(sp,proteomeConfig['ortholog']['orthologsUrl'])
            elif proteomeConfig['ortholog']['orthologCollectionMethod']=="file": # If collecting from file
                orthologsToWrite,proteinsUsedInFile = getOrthologsFromFile(sp)
            else:
                print("ERROR: an orthologCollectionMethod must be defined in the evidenceConfig.yml file. Choose either url or file")
                exit()
            orthologPairs.extend(orthologsToWrite)
            proteinsUsed.extend(proteinsUsedInFile)
        
        ## Primary identifiers
        species_df = pd.DataFrame.from_records(Species.objects.filter(Q(tax_id__in=species)).values('id','tax_id'))
        species_df.columns = ['species_id','tax_id']
        proteome = Proteome.objects.filter(Q(version=proteomeVersion))
        proteome_df = pd.DataFrame.from_records(Proteome.objects.filter(Q(version=proteomeVersion)).values('id','species_id'))
        proteome_df.columns = ['proteome_id','species_id']
        proteome_df = pd.merge(species_df,proteome_df, on='species_id', how='inner')
        proteome_df = proteome_df.drop('species_id', axis=1)
        uniprot_df = pd.DataFrame.from_records(Protein.objects.filter(Q(proteome__in=proteome)).values('id','uniprot_id','proteome_id'))
        uniprot_df = pd.merge(uniprot_df,proteome_df, on='proteome_id', how='inner')
        uniprot_df = uniprot_df.drop('proteome_id', axis=1)
        uniprot_df.columns = ['protein_id','uniprot_id','tax_id']

        ## Secondary mappings
        uniprot_alt_df = pd.DataFrame.from_records(IdMapping.objects.filter().values('mapped_id','protein_id'))
        uniprot_alt_df = pd.merge(uniprot_alt_df,uniprot_df, on='protein_id', how='inner')
        uniprot_alt_df = uniprot_alt_df.drop('uniprot_id', axis=1)
        uniprot_alt_df.columns = ['uniprot_id','protein_id','tax_id']

        uniprot_df = pd.concat([uniprot_df,uniprot_alt_df],axis=0).drop_duplicates().reset_index(drop=True)
        
        print("Getting proteinID map for proteins: "+str(len(set(proteinsUsed))))
        uniprot_df = uniprot_df[uniprot_df['uniprot_id'].isin(set(proteinsUsed))].reset_index(drop=True)
        uniprot_df[uniprot_df['tax_id']=='9823'].reset_index(drop=True)
        
        # Group by 'protein_id' and aggregate 'uniprot_id' into lists
        uniprot_dict = uniprot_df.groupby('uniprot_id')['protein_id'].apply(list).to_dict()
        # uniprot_df[['uniprot_id','protein_id']].groupby('uniprot_id').count().sort_values('protein_id')

        # Calculate the chunk size
        num_cores=os.cpu_count()
        len_orthologPairs = len(orthologPairs)
        chunk_size = len_orthologPairs // num_cores
        # Initialize a list to store start and end indices for each chunk
        chunk_indices = []
        # Create chunks and calculate start and end indices
        for i in range(num_cores):
            start_index = i * chunk_size
            end_index = (i + 1) * chunk_size if i < num_cores - 1 else len_orthologPairs
            chunk_indices.append((start_index, end_index))
        
        with tqdm_joblib(tqdm(desc='Create Ortholog table', total=len(chunk_indices))) as progress_bar:
            orthologs_list = Parallel(n_jobs=num_cores)(delayed(createOrthologEntry_parallel)(uniprot_dict,orthologPairs[start:end]) for start,end in chunk_indices)
        orthologs_df = pd.concat(orthologs_list,ignore_index=True).drop_duplicates().reset_index(drop=True)
        orthologs_df['version'] = orthologVersion
        # orthologs_df[(orthologs_df['speciesA']=='9606') & (orthologs_df['speciesB']=='9031')][['proteinA_id']].drop_duplicates().reset_index(drop=True)#pig
        # orthologs_df[(orthologs_df['speciesA']=='9615') & (orthologs_df['speciesB']=='9606')][['proteinA_id']].drop_duplicates().reset_index(drop=True)#dog
        # orthologs_df[(orthologs_df['speciesA']=='9606') & (orthologs_df['speciesB']=='7955')][['proteinA_id']].drop_duplicates().reset_index(drop=True)#d.rerio
        

        ## delete old mappings
        old_orthologs_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA__in=species_to_update) | Q(proteinB_id__in=species_to_update)).values('id'))
        old_orthologs = ProteinOrtholog.objects.filter(Q(speciesA__in=species_to_update) | Q(proteinB_id__in=species_to_update))
        old_orthologs.count()
        if len(old_orthologs)>0:
            old_orthologs.delete()

        # Creating an in-memory csv for the data
        mem_csv = StringIO()
        orthologs_df.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        # Writing the csv to the evidenceLink table
        print("Writing %i links to DB" %len(orthologs_df))
        ## TODO --> add to table the csv_io module
        with closing(mem_csv) as csv_io:
            ProteinOrtholog.objects.from_csv(csv_io)


def getOrthologs(species,proteomeConfig):
    # Getting orthologs for all species included in the isntance config, unless orthologs for the species
    # and version already exist in the DB. Has the opportunity to fetch orthologs either from the
    # InParanoidb website OR from a collection of all SQLtable files for the species (in this case, the
    # sqltable files must be available under data/tmp/resultsInParanoid). creates pairs of orthologs
    # From the ortholog groups, so that iformatioon can be transferred to both seed orthologs and
    # inparalogs. Saves ortholg relations in he DB where taxIDA is always > taxIDB

    proteomeVersion=proteomeConfig['genome']['version']
    orthologVersion=proteomeConfig['ortholog']['version']
    # Checking which species pairs already has orthologs in the DB
    speciesPairsToGet=set()
    for sA in species:
        for sB in species:
            if sA!=sB:
                speciesA,speciesB=getSpeciesOrder(sA,sB) # Always puts species in the correct order
                speciesPairsToGet.add(speciesA+"|"+speciesB)
                ## TODO --> it will check for the existence of ≥1 ortholog, source of error!!!
                # orthologVersionInDB = ProteinOrtholog.objects.filter(Q(version=orthologVersion) &  Q(speciesA=speciesA) &  Q(speciesB=speciesB))
                # if not orthologVersionInDB.exists():
                #     speciesPairsToGet.add(speciesA+"|"+speciesB)
                # else:
                #     print("Orthologs already exist for "+str(sA)+" - "+str(sB))
    if len(speciesPairsToGet)>0:
        print("Getting Orthologs")
        ## STEP 1 --> Fast
        orthologPairs=[] # ["proteinA|proteinB|speciesA|speciesB",...]
        proteinsUsed=[] # ["proteinA","proteinB"]
        for sp in speciesPairsToGet:
            if proteomeConfig['ortholog']['orthologCollectionMethod']=="url": # If collecting from url
                orthologsToWrite,proteinsUsedInFile = getOrthologsFromUrl(sp,proteomeConfig['ortholog']['orthologsUrl'])
            elif proteomeConfig['ortholog']['orthologCollectionMethod']=="file": # If collecting from file
                orthologsToWrite,proteinsUsedInFile = getOrthologsFromFile(sp)
            else:
                print("ERROR: an orthologCollectionMethod must be defined in the evidenceConfig.yml file. Choose either url or file")
                exit()
            orthologPairs.extend(orthologsToWrite)
            proteinsUsed.extend(proteinsUsedInFile)
        
        ## Primary identifiers
        species_df = pd.DataFrame.from_records(Species.objects.filter(Q(tax_id__in=species)).values('id','tax_id'))
        species_df.columns = ['species_id','tax_id']
        proteome = Proteome.objects.filter(Q(version=proteomeVersion))
        proteome_df = pd.DataFrame.from_records(Proteome.objects.filter(Q(version=proteomeVersion)).values('id','species_id'))
        proteome_df.columns = ['proteome_id','species_id']
        proteome_df = pd.merge(species_df,proteome_df, on='species_id', how='inner')
        proteome_df = proteome_df.drop('species_id', axis=1)
        uniprot_df = pd.DataFrame.from_records(Protein.objects.filter(Q(proteome__in=proteome)).values('id','uniprot_id','proteome_id'))
        uniprot_df = pd.merge(uniprot_df,proteome_df, on='proteome_id', how='inner')
        uniprot_df = uniprot_df.drop('proteome_id', axis=1)
        uniprot_df.columns = ['protein_id','uniprot_id','tax_id']

        ## Secondary mappings
        uniprot_alt_df = pd.DataFrame.from_records(IdMapping.objects.filter().values('mapped_id','protein_id'))
        uniprot_alt_df = pd.merge(uniprot_alt_df,uniprot_df, on='protein_id', how='inner')
        uniprot_alt_df = uniprot_alt_df.drop('uniprot_id', axis=1)
        uniprot_alt_df.columns = ['uniprot_id','protein_id','tax_id']

        uniprot_df = pd.concat([uniprot_df,uniprot_alt_df],axis=0).drop_duplicates().reset_index(drop=True)
        
        print("Getting proteinID map for proteins: "+str(len(set(proteinsUsed))))
        uniprot_df = uniprot_df[uniprot_df['uniprot_id'].isin(set(proteinsUsed))].reset_index(drop=True)

        dog_notUsed = list(set([i for i in proteinsUsed if i not in uniprot_df['uniprot_id']]))

        uniprot_df[uniprot_df['tax_id']=='9823'].reset_index(drop=True)
        
        # Group by 'protein_id' and aggregate 'uniprot_id' into lists
        uniprot_dict = uniprot_df.groupby('uniprot_id')['protein_id'].apply(list).to_dict()
        # uniprot_df[['uniprot_id','protein_id']].groupby('uniprot_id').count().sort_values('protein_id')

        # Calculate the chunk size
        num_cores=os.cpu_count()
        len_orthologPairs = len(orthologPairs)
        chunk_size = len_orthologPairs // num_cores
        # Initialize a list to store start and end indices for each chunk
        chunk_indices = []
        # Create chunks and calculate start and end indices
        for i in range(num_cores):
            start_index = i * chunk_size
            end_index = (i + 1) * chunk_size if i < num_cores - 1 else len_orthologPairs
            chunk_indices.append((start_index, end_index))
        
        with tqdm_joblib(tqdm(desc='Create Ortholog table', total=len(chunk_indices))) as progress_bar:
            orthologs_list = Parallel(n_jobs=num_cores)(delayed(createOrthologEntry_parallel)(uniprot_dict,orthologPairs[start:end]) for start,end in chunk_indices)
        orthologs_df = pd.concat(orthologs_list,ignore_index=True).drop_duplicates().reset_index(drop=True)
        orthologs_df['version'] = orthologVersion
        orthologs_df[(orthologs_df['speciesA']=='9823') & (orthologs_df['speciesB']=='9606')][['proteinA_id']].drop_duplicates().reset_index(drop=True)#pig
        orthologs_df[(orthologs_df['speciesA']=='9615') & (orthologs_df['speciesB']=='9606')][['proteinA_id']].drop_duplicates().reset_index(drop=True)#dog
        orthologs_df[(orthologs_df['speciesA']=='9606') & (orthologs_df['speciesB']=='7955')][['proteinA_id']].drop_duplicates().reset_index(drop=True)#d.rerio
        
        # Creating an in-memory csv for the data
        mem_csv = StringIO()
        orthologs_df.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        # Writing the csv to the evidenceLink table
        print("Writing %i links to DB" %len(orthologs_df))
        ## TODO --> add to table the csv_io module
        with closing(mem_csv) as csv_io:
            ProteinOrtholog.objects.from_csv(csv_io)

        ## STEP 2 --> Slow!!!
        # Collecting protein ids for the unique set of proteins to be used.
        # proteinIdMap={}
        # for p in set(proteinsUsed):
        #     ids=getProteinIdFromDB(p)
        #     proteinIdMap[p]=ids
        
        # Writing all collected orthologs, using lookup map to get protein ids
        # orthologsToWrite=[]
        # orthologsWritten=0
        # for o in orthologPairs:
        #     pA,pB,spA,spB=o.split("|")
        #     for idA in uniprot_dict[pA]:
        #         for idB in uniprot_dict[pB]:
        #             orth = ProteinOrtholog(proteinA_id=idA, proteinB_id=idB,speciesA=spA, speciesB=spB, version=orthologVersion)
        #             orthologsToWrite.append(orth)
        #             orthologsWritten+=1
        #             if orthologsWritten % 100000==0: # Writing 100 000 orthologs at the time
        #                 ProteinOrtholog.objects.bulk_create(orthologsToWrite)
        #                 orthologsToWrite=[]
        # ProteinOrtholog.objects.bulk_create(orthologsToWrite)
        # print("Wrote "+str(orthologsWritten)+" orthologs to DB")
    

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
    getOrthologs(instanceConfig['instance']['species'],proteomeConfig)