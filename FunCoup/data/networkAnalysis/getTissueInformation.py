import requests
import zipfile
import csv
import os
import requests
import os
import gzip
import shutil
import numpy as np
import pandas as pd
from contextlib import closing
from io import StringIO
import yaml
from yaml.loader import SafeLoader

import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q

def downloadFile(urlToDownload, filename):
    url = urlToDownload
    makeRequest=True
    numberOfAttempts=0
    # Sometimes the request times out, here it retries to send request 10 times
    while makeRequest==True and numberOfAttempts<=10:
        numberOfAttempts+=1
        try:
            myfile = requests.get(url)
            if myfile.status_code==200:
                open(filename, "wb").write(myfile.content)
            elif myfile.status_code!=404 and myfile.status_code!=403: # Do not print if file is not available for this species
                print("Error when downloading file: "+url+" status code: "+ str(myfile.status_code))
            makeRequest=False
        except:
            print("Request for: "+url+" timed out. Trying again, attempts: "+str(numberOfAttempts))
            makeRequest=True

def gzipFile(filename):
    with gzip.open(filename+".gz", 'rb') as f_in:
        with open(filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def unzipFile(filename):
    print("unzipping file: "+filename+".zip")
    with zipfile.ZipFile(filename+".zip", "r") as zip_ref:
       zip_ref.extractall('.')

def extractTissueInformation(tissue_file, taxID, available_protein_df, tissueIDCounter):

    ## Read csv file and extract tissue info
    if taxID=='9606':
        ## HPA
        tissue_df = pd.read_csv(tissue_file, sep='\t', usecols=[0,2,4,5])
        tissue_df = tissue_df[(tissue_df['Reliability']!="Uncertain") & (tissue_df['Level']!="Not detected")].reset_index(drop=True)
        tissue_df['Gene'] = tissue_df['Gene'].str.strip()
        tissue_df['Tissue'] = tissue_df['Tissue']\
            .str.capitalize()\
            .str.strip()
        tissue_df = tissue_df.drop(['Reliability', 'Level'],axis=1)
    else:
        ## BGEE
        tissue_df = pd.read_csv(tissue_file, sep='\t', usecols=[0,3,4,5]) 
        tissue_df = tissue_df[(tissue_df['Expression']=='present') & (tissue_df['Call quality']=='gold quality')].reset_index(drop=True)
        tissue_df['Gene ID'] = tissue_df['Gene ID'].str.strip()
        tissue_df['Anatomical entity name'] = tissue_df['Anatomical entity name']\
            .str.replace(r'["\']', '',regex=True)\
            .str.split("(").str[0]\
            .str.capitalize()\
            .str.strip()
        tissue_df = tissue_df.drop(['Expression', 'Call quality'],axis=1)

    tissue_df.columns = ['gene','tissue']    
    tissue_size_cutoff = 100
    tissueCounts = tissue_df['tissue'].value_counts()
    filtered_tissue_df = tissue_df[tissue_df['tissue'].isin(tissueCounts[tissueCounts >= tissue_size_cutoff].index)].reset_index(drop=True)
    filtered_tissue_df = pd.merge(filtered_tissue_df,available_protein_df, on='gene', how='inner')
    
    filtered_tissue_df = filtered_tissue_df.drop('gene',axis=1)
    filtered_tissue_df.columns = ['tissue','protein_id']
    filtered_tissue_df = filtered_tissue_df.sort_values('tissue')
    filtered_tissue_df = filtered_tissue_df.drop_duplicates().reset_index(drop=True)

    tissue_groups = filtered_tissue_df.groupby('tissue')
    # Generate unique IDs for each tissue group
    tissue_ids = {tissue: i + tissueIDCounter for i, (tissue, _) in enumerate(tissue_groups)}

    # Assign tissue IDs to each row based on tissue group
    filtered_tissue_df['tissue_id'] = filtered_tissue_df['tissue'].map(tissue_ids)
    filtered_tissue_df['species'] = taxID

    print("Number of tissues: %i " %filtered_tissue_df['tissue'].nunique())
    print("Unique genes: %i" %filtered_tissue_df['protein_id'].nunique())

    new_order = ['protein_id', 'species', 'tissue', 'tissue_id']
    new_tissueIDCounter = np.max(filtered_tissue_df['tissue_id'])+1
    return filtered_tissue_df[new_order], new_tissueIDCounter

def getTissueInformation():
    print("Initiating retrieval of all tissue-gene relations")
    if os.path.exists('configFiles/websiteConfig.yml'):
        with open('configFiles/websiteConfig.yml') as f:
            webConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your websiteConfig.yml does not exist")
        exit()
    instanceName=webConfig['instanceName']
    if os.path.exists('configFiles/'+instanceName):
        if os.path.exists('configFiles/'+instanceName+'/proteomeConfig.yml'):
            with open('configFiles/'+instanceName+'/proteomeConfig.yml') as f:
                proteomeConfig = yaml.load(f, Loader=SafeLoader)
        else:
            print("Your config file does not exist")
            exit()

    urlBgee="https://www.bgee.org/ftp/current/download/calls/expr_calls/SPECIES_expr_simple.tsv.gz" #v15.2
    bgeeSpecies=["Danio rerio","Mus musculus", "Caenorhabditis elegans","Rattus norvegicus","Bos taurus", "Sus scrofa","Canis lupus familiaris","Gallus gallus","Drosophila melanogaster"]
    #bgeeSpecies=["Mus musculus", "Caenorhabditis elegans","Rattus norvegicus","Bos taurus", "Sus scrofa","Canis lupus familiaris","Gallus gallus","Drosophila melanogaster"]
    #bgeeSpecies=["Danio rerio"]
    urlHPA="https://www.proteinatlas.org/download/normal_tissue.tsv.zip" #v23

    speciesToNCBI={
    "Danio rerio":"7955",
    "Homo sapiens":"9606",
    "Mus musculus":"10090",
    "Caenorhabditis elegans":"6239",
    "Rattus norvegicus":"10116",
    "Bos taurus":"9913",
    "Sus scrofa":"9823",
    "Canis lupus familiaris":"9615",
    "Gallus gallus":"9031",
    "Drosophila melanogaster":"7227"
    }

    tissueIDCounter=1
    for species_name,taxID in speciesToNCBI.items():
        print("Getting tissues for %s" %species_name)

        ## Extract known proteins for species under study
        proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=taxID))[0]
        available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
        available_protein_df.columns = ['gene','protein_id']
        print("Got mappings")

        if taxID == '9606':
            url = urlHPA
            filename = "data/tmp/"+urlHPA.split("/")[-1] #.split(".zip")[0]
        else:
            url=urlBgee.replace("SPECIES", species_name.replace(" ", "_"))
            filename="data/tmp/"+"SPECIES_expr_simple.tsv.gz".replace("SPECIES", species_name.replace(" ", "_"))

        if not os.path.exists(filename):
            downloadFile(url, filename)

        print("Downloaded file")

        tissuesThisSpecies, tissueIDCounter = extractTissueInformation(filename, taxID, available_protein_df, tissueIDCounter)

        print("Final dataframe to write:")
        print(tissuesThisSpecies)
        # Creating an in-memory csv for the data
        mem_csv = StringIO()
        tissuesThisSpecies.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        # Writing the csv to the Tissue table
        print("Writing %i tissue-gene associations to DB" %(len(tissuesThisSpecies)))
        with closing(mem_csv) as csv_io:
            Tissue.objects.from_csv(csv_io)
