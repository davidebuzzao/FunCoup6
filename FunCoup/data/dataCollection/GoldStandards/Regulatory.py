import requests
import zipfile
import os
import yaml
import swifter
from yaml.loader import SafeLoader
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q

# Following direction map
# 1:no direction 2: A->B; 3:A<-B; 5:A<->B

def downloadFile(urlToDownload, filename):
    # Downloading file from url
    url = urlToDownload
    myfile = requests.get(url)
    if myfile.status_code==200:
        open(filename, "wb").write(myfile.content)
    else:
        print("Cant download file from: "+url)

def downloadFileNoSSLCertificateVerification(urlToDownload, filename):
    # Downloading file from url
    url = urlToDownload
    myfile = requests.get(url,verify=False)
    if myfile.status_code==200:
        open(filename , "wb").write(myfile.content)
    else:
        print("Cant download file from: "+url)

def unzipFile(path,filename):
    # Unzipping file to path
    print("unzipping file: "+filename)
    with zipfile.ZipFile(path+filename, "r") as zip_ref:
       zip_ref.extractall(path=path)

# Function to transform rows based on conditions
def transform_row(row):
    if row['proteinA'] > row['proteinB']:
        return row
    else:
        new_row = row.copy()
        new_row['proteinA'], new_row['proteinB'] = row['proteinB'], row['proteinA']
        if new_row['direction'] == 2:
            new_row['direction'] = 3
        return new_row

def parseRegNetwork(filename):
    directionMap={'-->':2,'--|':2,'-/-':2,'-p':0} # -p missing direction
    linksDf = pd.read_csv(filename, sep=' ', usecols=[0,2,4], comment='#')
    linksDf.columns = ['proteinA','proteinB','direction']
    linksDf['direction'] = linksDf['direction'].apply(lambda x: directionMap[x]) 
    return linksDf[linksDf['direction'] == 2]

def parseRegulonDB(filename):
    directionMap={'Strong':2,'Confirmed':2,'Weak':2,'ND':0} # Weak is about 1/3 of interactions
    linksDf = pd.read_csv(filename, sep='\t', usecols=[2,4,6], comment='#')
    linksDf.columns = ['proteinA','proteinB','direction']
    linksDf['direction'] = linksDf['direction'].apply(lambda x: directionMap[x]) 
    return linksDf[linksDf['direction'] == 2]

def parseTRRUST(filename):
    directionMap={'Activation':2,'Repression':2,'Unknown':2}
    linksDf = pd.read_csv(filename, sep='\t', usecols=[0,1,2], header=None)
    linksDf.columns = ['proteinA','proteinB','direction']
    linksDf['direction'] = linksDf['direction'].apply(lambda x: directionMap[x]) 
    return linksDf[linksDf['direction'] == 2]

def parseYeastract(filename):
    directionMap={'Positive':2,'Negative':2,'Unknown':2}
    linksDf = pd.read_csv(filename, sep='\t', usecols=[1,3,5]) #212467; 203509 if usecols=[0,2,5]
    linksDf.columns = ['proteinA','proteinB','direction']
    linksDf['direction'] = linksDf['direction'].fillna('Unknown')
    linksDf['direction'] = linksDf['direction'].apply(lambda x: directionMap[x]) 
    return linksDf[linksDf['direction'] == 2]

def getRegulatoryPairs(sp,goldStandardConfig,proteomeConfig):
    if sp in ['9606','10090']:
        print('Extract RegNetwork')
        regnetwork_url = goldStandardConfig['RegNetwork']['url']
        regnetwork_folder = regnetwork_url.split('/')[-1]
        if not os.path.exists('data/tmp/'+regnetwork_folder):
            downloadFileNoSSLCertificateVerification(regnetwork_url,'data/tmp/'+regnetwork_folder)
        unzipFile('data/tmp/',regnetwork_folder)
        regnetwork_filename = goldStandardConfig['RegNetwork']['data'][0]
        if sp=='10090': regnetwork_filename = goldStandardConfig['RegNetwork']['data'][1]
        regnetwork_version = goldStandardConfig['RegNetwork']['version']
        linksDf_regnetwork = parseRegNetwork("data/tmp/"+regnetwork_filename)

        print('Extract TRRUST')
        trrust_filename = goldStandardConfig['TRRUST']['data'][0]
        trrust_url = goldStandardConfig['TRRUST']['url'][0]
        if sp=='10090': 
            trrust_filename = goldStandardConfig['TRRUST']['data'][1]
            trrust_url = goldStandardConfig['TRRUST']['url'][1]
        if not os.path.exists('data/tmp/'+trrust_filename):
            downloadFile(trrust_url,"data/tmp/"+trrust_filename)
        trrust_version = goldStandardConfig['TRRUST']['version']
        linksDf_trrust = parseTRRUST("data/tmp/"+trrust_filename)

        linksDf = pd.concat([linksDf_regnetwork,linksDf_trrust], ignore_index=True)
        version = ','.join([regnetwork_version,trrust_version])
    
    elif sp in ['83333']:
        print('Extract RegulonDB')
        # TF-gene links
        regulondb_url = goldStandardConfig['RegulonDB']['url'][0]
        tf_gene = goldStandardConfig['RegulonDB']['data'][0] 
        if not os.path.exists('data/tmp/'+tf_gene):
            downloadFile(regulondb_url,"data/tmp/"+tf_gene)
        tf_geneDf = parseRegulonDB("data/tmp/"+tf_gene)
        # TF-TF links 
        regulondb_url = goldStandardConfig['RegulonDB']['url'][1]
        tf_tf = goldStandardConfig['RegulonDB']['data'][1]
        if not os.path.exists('data/tmp/'+tf_tf):
            downloadFile(regulondb_url,"data/tmp/"+tf_tf)
        tf_tfDf = parseRegulonDB("data/tmp/"+tf_tf)

        linksDf = pd.concat([tf_geneDf,tf_tfDf], ignore_index=True)
        version = goldStandardConfig['RegulonDB']['version']

    elif sp in ['559292']:
        print('Extract Yeastract')
        yeastract_filename = goldStandardConfig['Yeastract']['data']
        yeastract_url = goldStandardConfig['Yeastract']['url']
        if not os.path.exists('data/tmp/'+yeastract_filename):
            downloadFile(yeastract_url,"data/tmp/"+yeastract_filename)
        
        linksDf = parseYeastract("data/tmp/"+yeastract_filename)
        version = goldStandardConfig['Yeastract']['version']

    linksDf['proteinA'] = linksDf['proteinA'].str.lower()
    linksDf['proteinB'] = linksDf['proteinB'].str.lower()

    # Collecting mappings from DB
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=sp))[0]
    availableProteinsDf = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id'))
    availableProteinsDf.columns = ['proteinA','idA']
    availableProteinsDf['proteinA'] = availableProteinsDf['proteinA'].str.lower()

    # mapping the protein identifiers to the database ID
    linksDf = pd.merge(linksDf, availableProteinsDf, how='inner', on='proteinA',)
    availableProteinsDf.columns = ['proteinB','idB']
    linksDf = pd.merge(linksDf, availableProteinsDf, how='inner', on='proteinB')
    
    linksDf.drop(['proteinA','proteinB'],axis=1,inplace=True)
    linksDf.columns = ['direction','proteinA','proteinB']
    # Apply the transformation function to each row
    linksDf = linksDf.apply(transform_row, axis=1)\
                .drop_duplicates(subset=['proteinA', 'proteinB', 'direction'])\
                .reset_index(drop=True)
    linksDf.columns = ['direction','proteinA','proteinB']
    # ## Remove duplicated, add option 3 (bidirectionality) and remove ambiguities (A,B,1) and (A,B,2) with (A,B,3)    
    # linksDf = transformedlinksDf.sort_values('direction', ascending=False)\
    #     .drop_duplicates(subset=['proteinA', 'proteinB', 'direction'])\
    #     .groupby(['proteinA','proteinB'])['direction'].sum().reset_index()

    # Creating a new gold standard to relate the new links to
    speciesInDB = Species.objects.get(tax_id=sp)
    goldStandardInDB = GoldStandard(version=version, type='Regulatory', species=speciesInDB)
    goldStandardInDB.save()

    linksDf['goldStandard'] = goldStandardInDB.id
    # Creating an in-memory csv for the data
    mem_csv = StringIO()
    linksDf.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing "+str(len(linksDf.index))+" Regulatory links to DB for "+str(sp))
    with closing(mem_csv) as csv_io:
        GoldStandardLink.objects.from_csv(csv_io)

    return goldStandardInDB