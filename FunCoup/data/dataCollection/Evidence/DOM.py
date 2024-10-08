import os
# import matplotlib.pyplot as plt
import yaml
from yaml.loader import SafeLoader
from collections import Counter
from sklearn.preprocessing import MinMaxScaler,PowerTransformer
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
from auxiliary_functions import *
from tqdm import tqdm

#### DOM
def getPfam():
    # Checks if pfam exists in tmp, otherwise downloads and prepares the files
    pfamScanExist=False
    for file in os.listdir("data/tmp/"):
        if file=="pfam":
            pfamScanExist=True
    if pfamScanExist==False:
        print("Cant find pfam, downloading")
        os.makedirs("data/tmp/pfam")
        os.system("wget -P  data/tmp/pfam ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz")
        os.system("wget -P data/tmp/pfam ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.dat.gz")
        os.system("wget -P data/tmp/pfam ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/active_site.dat.gz")
        os.system("wget -P data/tmp/pfam ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz")
        os.system("cd data/tmp/pfam && tar xvfz PfamScan.tar.gz")
        os.system("cd data/tmp/pfam && gunzip -d Pfam-A.hmm.gz && gunzip -d Pfam-A.hmm.dat.gz && gunzip -d active_site.dat.gz && hmmpress Pfam-A.hmm")

def getProteomeConfigAndEbiIdentifiers(speciesToCollect,proteomeConfig):
    # Collecting the information necessary to download the proteome files for the species
    # the proteome files are necessary  for pfam to be able to extract domains for each protein
    if os.path.isfile("data/tmp/README")==False:
        downloadFile(proteomeConfig['genome']['urlUniProtReadme'], "README")
    ebiIdentifiers={}
    with open("data/tmp/README", 'r') as in_file:
        for line in in_file:
            if line.startswith("UP0"):
                l=line.split("\t")
                if len(l)> 4 and (l[1] in speciesToCollect):
                    ebiIdentifiers[l[1]]=[l[0], l[3]]
    return proteomeConfig,ebiIdentifiers

def runPfamForProteome(s, proteomeConfig,ebiIdentifiers):
    # If a pfam file does not already exist, runs pfam scan on the proteome file.
    # This is VERY SLOW, so if you can, please prepare these files in advance
    if os.path.isfile("data/tmp/"+s+"_pfam.tsv")==False:
        downloadFile(proteomeConfig['genome']['urlUniProtDownload']+ebiIdentifiers[s][1].capitalize()+"/"+ebiIdentifiers[s][0]+"/"+ebiIdentifiers[s][0]+"_"+s+".fasta.gz", s+".fasta.gz")
        unzipFile(s+".fasta")
        os.remove("data/tmp/"+s+".fasta.gz")
        print("Running pfam scan on "+s+".fasta")
        os.system("perl -I data/tmp/pfam/PfamScan/ data/tmp/pfam/PfamScan/pfam_scan.pl -fasta data/tmp/"+s+".fasta -dir data/tmp/pfam -outfile data/tmp/"+s+"_pfam.tsv -cpu 1")
        os.remove("data/tmp/"+s+".fasta")

def compute_frequency(x):
    x_count = Counter(x)
    x_dict = dict([(p,val) for p,val in x_count.items()])
    return x_dict

def getDOMScoresForSpecies(s, domIntDf,evidenceConfig, proteomeConfig):
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)

    # Collects all proteins sharing a domain from the pfam scan file for this species
    # Collecting the pfam domain predictions in a dataframe
    pfamDf = pd.read_csv("data/tmp/"+s+"_pfam.tsv", delim_whitespace=True,usecols=[0,5], header=None, comment='#')
    pfamDf.columns = ['proteinIdentifier','domain']
    pfamDf['proteinIdentifier'] = pfamDf['proteinIdentifier'].apply(lambda x: x.split('|')[1]) # Getting the protein id
    pfamDf['domain'] = pfamDf['domain'].apply(lambda x: x.split('.')[0]) # Excluding the . from the domain id
                  
    # Getting mappings from the DB for the specified proteome
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=s))[0]
    availableProteinsDf = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id'))
    availableProteinsDf.columns = ['proteinIdentifier','id']
                               
    # mapping the protein identifiers to the database ID
    mapping = pd.DataFrame(list(pfamDf['proteinIdentifier']))
    mapping.columns = ['proteinIdentifier']
    preMap = set(mapping['proteinIdentifier'])
    mappedPfamDf = pd.merge(pfamDf,availableProteinsDf, how='inner',on='proteinIdentifier')
    postMap = set(mappedPfamDf['proteinIdentifier'])
    excluded = preMap - postMap
    print("Excluded "+str(len(excluded))+" genes, as they couldnt be mapped. Getting DOM scores for "+str(len(postMap))+" genes for species: "+s)
    del availableProteinsDf
    pfamDf_count = mappedPfamDf.groupby(['id']).size().reset_index(name='count').sort_values('count',ascending=False)

    # Getting domain interactions per protein pair
    mappedPfamDf.columns = ['proteinA','domainA','proteinA_id']
    linkDf = pd.merge(domIntDf,mappedPfamDf, on='domainA', how='inner')
    mappedPfamDf.columns = ['proteinB','domainB','proteinB_id']
    linkDf = pd.merge(linkDf,mappedPfamDf, on='domainB', how='inner')

    linkDf = linkDf.drop(['proteinA','proteinB'],axis=1) # Drop columns that are no longer needed
    linkDf = linkDf[linkDf.proteinA_id != linkDf.proteinB_id] # Exclude protein pairs where Px and Py is the same protein
    linkDf.proteinA_id, linkDf.proteinB_id = np.where(linkDf.proteinA_id < linkDf.proteinB_id, [linkDf.proteinB_id, linkDf.proteinA_id], [linkDf.proteinA_id, linkDf.proteinB_id])
    linkDf = linkDf.sort_values("score", ascending=False)
    linkDf = linkDf.drop_duplicates(subset=['proteinA_id', 'proteinB_id','domainA','domainB'], keep="first").reset_index(drop=True)
    
    pfam_dict = pfamDf_count.set_index('id')['count'].to_dict()
    pfam_domains = domIntDf['domainA'].to_list() + domIntDf['domainB'].to_list()
    pfam_domains_dict = compute_frequency(pfam_domains)

    linkDf['NumDomProteinA'] = linkDf['proteinA_id'].map(pfam_dict)
    linkDf['NumDomProteinB'] = linkDf['proteinB_id'].map(pfam_dict)
    linkDf['fDomainA'] = linkDf['domainA'].map(pfam_domains_dict)
    linkDf['fDomainB'] = linkDf['domainB'].map(pfam_domains_dict)
    linkDf['fweight'] = 1/np.log(1+linkDf['fDomainA']+linkDf['fDomainB'])
    linkDf['wscore'] = linkDf['score']*linkDf['fweight']
    linkDf['combos'] = linkDf['NumDomProteinA']*linkDf['NumDomProteinB'] 
    linkDf['metadata'] = linkDf['domainA']+'-'+linkDf['domainB']

    grouped_df = linkDf.sort_values("wscore", ascending=False).groupby(['proteinA_id', 'proteinB_id']).agg({'wscore': 'mean', 'combos': 'first','metadata': lambda x: ','.join(x)}).reset_index()
    # grouped_df = linkDf.sort_values("wscore", ascending=False).groupby(['proteinA_id', 'proteinB_id']).agg({'wscore': 'max', 'combos': 'first', 'metadata': 'first'}).reset_index()
    # Calculate the count of rows for each group in the original DataFrame
    grouped_df['count'] = linkDf.groupby(['proteinA_id', 'proteinB_id']).size().reset_index(name='count')['count']
    # Calculate the score based on the formula: avg(score) * (count / combos)
    grouped_df['wscore'] = grouped_df['wscore'] * (grouped_df['count'] / grouped_df['combos'])
    linkDf = grouped_df.drop(['combos','count'], axis=1)
    linkDf.columns = ['proteinA','proteinB','score','metadata']
    
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

        # Preparing the dataframe to be added to the database
        speciesInDB = Species.objects.get(tax_id = s)
        evidenceInDB = Evidence(type="DOM", species=speciesInDB, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidenceInDB.save()
        
        linkDf = linkDf.sample(frac=1).reset_index(drop=True)
        linkDf.insert(0,'evidence',evidenceInDB.id)

        return linkDf

def getDOMLinks(evidenceConfig, speciesToCollect, proteomeConfig):
    # Collects domain interactions and scores from UniDomInt. The uniDom scores are used
    # as is presented in the UniDomInt file. Then iterates over all
    # Proteins within a species having domains from the interaction, and adds the score to a protein pair for
    # its highest scoring domain interaction. To be able to get the domains included in a
    # protein, pfam scan is run on each proteome file. As PfamScan is VERY SLOW, this can preferrably be done
    # in advance. If you already have pfamScan files for all you species, name them taxId_pfam.tsv and put them
    # in the data/tmp/ dir. If files are already there, pfam will be skipped.

    # Get pfam domains for all species
    getPfam()
    proteomeConfig,ebiIdentifiers = getProteomeConfigAndEbiIdentifiers(speciesToCollect,proteomeConfig)
    Parallel(n_jobs=os.cpu_count())(delayed(runPfamForProteome)(s,proteomeConfig,ebiIdentifiers) for s in speciesToCollect)

    # Collecting the domain interactions, and their scores from the uniDomInt file
    uniDomIntFilepath="data/tmp/"+evidenceConfig['url'].split("/")[-1]
    domIntDf = pd.read_csv(uniDomIntFilepath, sep='\t',usecols=[0,1,23], skiprows=10, header=None, comment='/')
    domIntDf.columns = ['domainA','domainB','score']
    domIntDf['score'] = domIntDf['score'].str.replace(",",".").astype(float)

    domIntDf = domIntDf[domIntDf.score > 0] # Excluding scores that are 0

    # Getting scores and writing links to the database for all proteins having the domains of an interaction
    for s in speciesToCollect:
        linkDf = getDOMScoresForSpecies(s, domIntDf,evidenceConfig, proteomeConfig)
        if len(linkDf)==0: continue
        # Creating an in-memory csv for the data
        mem_csv = StringIO()
        linkDf.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        # Writing the csv to the evidenceLink table
        print("Writing "+str(len(linkDf.index))+" DOM links to DB for "+s)
        with closing(mem_csv) as csv_io:
            DOM.objects.from_csv(csv_io)

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
    speciesToCollect = instanceConfig['instance']['species']
    evidenceConfig = evidenceConfig['DOM']
    getDOMLinks(evidenceConfig['DOM'], ['9606'], proteomeConfig)