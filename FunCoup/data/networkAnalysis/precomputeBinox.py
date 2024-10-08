import os
import pandas as pd
import numpy as np
import time
import networkx as nx
from joblib import Parallel, delayed
from io import StringIO
from tqdm import tqdm
from auxiliary_functions import *
### NETWORK
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from django.conf import settings
import yaml
from yaml.loader import SafeLoader
import subprocess

def getBinoxRandomNetwork(networkFile,networksPath,cutoff_ppv):

    if networkFile.endswith("_full.gz") and "TRANSFERRED" not in networkFile:

        spName = networkFile.replace("FC6.0_","").replace("_full.gz","").split(".")[1]
        species_object = Species.objects.filter(species_name__icontains=spName)
        taxid=species_object[0].tax_id
        print(taxid)

        randFile="website/static/website/binox/randNet_"+taxid+".csv"
        network_df = pd.read_csv(networksPath+"/"+networkFile, compression='gzip', sep='\t', usecols=[2,3,5])
        network_df.columns = ['#node1', 'node2', 'score']
        network_df = network_df[network_df['score']>=cutoff_ppv] # Filter network before running binox

        if os.path.exists(randFile)==False:
            tmpFile="website/static/website/binox/binoxNetwork_"+taxid+".csv"
            network_df.to_csv(tmpFile, index=False, sep ='\t')
            print("Made a BinoX-adjusted version of the network file, now running BinoX for network randomization...")

            binox="website/static/website/binox/BinoX"
            ## We filter out low confidence links above, here we write 0.5 but won't have any effect
            random_call = [binox, '-c', str(0.5), '-n', tmpFile, '-r', randFile, '-i', str(2000), '-s', str(taxid)]
            response = subprocess.run(random_call, stdout=subprocess.PIPE)
            print(response.stdout.decode('utf-8'))
            if "ERROR" in response.stdout.decode('utf-8'):
                exit()
            os.remove(tmpFile)
            print("Created randomized network file at: "+randFile)
        else:
            print("Randomized network for "+taxid+" already exist under website/static/website/binox/")

        pathwayFile="website/static/website/binox/pathways_"+taxid
        if os.path.exists(pathwayFile)==False:
            pathwaysDf = pd.DataFrame.from_records(Pathway.objects.filter(Q(species=taxid)).values('protein_id','pathway_id'))
            pathwaysDf.to_csv(pathwayFile, index=False, sep ='\t')
        else:
            print("Pathway file for "+taxid+" already exist under website/static/website/binox/")


def precomputeBinox():
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)

    print("Initiating precomputation of BINOX crosstalk")
    if os.path.exists('configFiles/websiteConfig.yml'):
        with open('configFiles/websiteConfig.yml') as f:
            webConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your websiteConfig.yml does not exist")
        exit()

    cutoff_ppv = webConfig['binoxCutoff']
    networksPath="website/static/website/networks/FunCoup6.0/"
    networkFile_list=os.listdir(networksPath)
    with tqdm_joblib(tqdm(desc='Transfer network', total=len(networkFile_list))) as progress_bar:
        parallel(delayed(getBinoxRandomNetwork)(networkFile,networksPath,cutoff_ppv)for networkFile in networkFile_list)


           