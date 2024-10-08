import requests
import os
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from Bio.KEGG.KGML import KGML_parser
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO
import yaml
from yaml.loader import SafeLoader


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

def updateKEGGPathwaysToUse(filename):
    # Getting new KEGG pathways from the old manually picked set of pathways. Adding
    # newly added pathways to the set, and removing pathways that has been deprecated.
    allIds={}
    with open(filename, 'r') as in_file:
        for line in in_file:
            if line.startswith("C"):
                allIds[line.split()[1]]=" ".join(line.split()[2:])

    return allIds

keggSpeciesNameMap={'9615':'cfa'}

def getKeggPathwayInformation():
    print("Initiating retrieval of all kegg pathway-gene relations")
    if os.path.exists('configFiles/websiteConfig.yml'):
        with open('configFiles/websiteConfig.yml') as f:
            webConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your websiteConfig.yml does not exist")
        exit()
    instanceName=webConfig['instanceName']
    if os.path.exists('configFiles/'+instanceName):
        if os.path.exists('configFiles/'+instanceName+'/instanceConfig.yml'):
            with open('configFiles/'+instanceName+'/instanceConfig.yml') as f:
                instanceConfig = yaml.load(f, Loader=SafeLoader)
        else:
            print("Your config file does not exist")
            exit()
        if os.path.exists('configFiles/'+instanceName+'/goldStandardConfig.yml'):
            with open('configFiles/'+instanceName+'/goldStandardConfig.yml') as f:
                goldStandardConfig = yaml.load(f, Loader=SafeLoader)
        else:
            print("Your config file does not exist")
            exit()
        if os.path.exists('configFiles/'+instanceName+'/proteomeConfig.yml'):
            with open('configFiles/'+instanceName+'/proteomeConfig.yml') as f:
                proteomeConfig = yaml.load(f, Loader=SafeLoader)
        else:
            print("Your config file does not exist")
            exit()

    filename="data/tmp/allKeggPathways" #v108.1
    downloadFile(goldStandardConfig['Metabolic']['urlPathways'], filename)
    allPathways=updateKEGGPathwaysToUse(filename)

    speciesToGet = instanceConfig['instance']['species']
    speciesToGet=['9606','10090','3702','224308','9913','10116','559292','9031','7227','7955','7719','9615','6239','36329','284812','39947','9823','44689','243232','273057','83332']
    #speciesToGet=['83333']
    pathwaysDf=pd.DataFrame(columns=('protein', 'species', 'pathway_id','pathway_name','pathway_db'))
    for sp in speciesToGet:
        allPathwaysForSp=0
        allGenesForSp=0
        speciesInDB= Species.objects.get(tax_id = sp)
        # get mappings
        proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=sp))[0]
        mappings = IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')
        mappingsDict={}
        for m in mappings:
            mappingsDict[m['mapped_id']]=m['protein_id']
    
        # Construct short kegg species name from the scientific name for the species in the DB
        if sp not in keggSpeciesNameMap:
            shortSpName=(speciesInDB.species_name.split()[0][0]+speciesInDB.species_name.split()[1][:2]).lower()
        else: # Special case for dog, as it can not be composed form the scientific name
            shortSpName=keggSpeciesNameMap[sp]
        # Iterating over each pathway in the set
        for pathway in allPathways:
            # Downloading pathway kgml file if it does not already exist
            if os.path.exists("data/tmp/"+shortSpName+pathway)==False:
                downloadFile(goldStandardConfig['Metabolic']['urlKGML'].replace("PATHWAYID",shortSpName+pathway), "data/tmp/"+shortSpName+pathway)
            # Reading KGML file if the pathway exist for this species
            genesInPathway=0
            if os.path.exists("data/tmp/"+shortSpName+pathway):
                allPathwaysForSp+=1
                pw = KGML_parser.read(open("data/tmp/"+shortSpName+pathway, 'r')) # Using kgml parser to read xml file
                # Iterating over all "gene" entries in the file to be able to use later when reading relations.
                for gene in pw.genes:
                    genes=gene.name.split()
                    for g in genes:
                        if g in mappingsDict:
                            allGenesForSp+=1
                            genesInPathway+=1
                            protId=mappingsDict[g]
                            pathwaysDf = pd.concat([pd.DataFrame([[protId,sp,pathway,allPathways[pathway],"KEGG"]], columns=pathwaysDf.columns), pathwaysDf], ignore_index=True)
                        elif g.split(":")[1] in mappingsDict:
                            allGenesForSp+=1
                            genesInPathway+=1
                            protId=mappingsDict[g.split(":")[1]]
                            pathwaysDf = pd.concat([pd.DataFrame([[protId,sp,pathway,allPathways[pathway],"KEGG"]], columns=pathwaysDf.columns), pathwaysDf], ignore_index=True)
            print("Got "+str(genesInPathway)+" genes for pathway "+pathway+" in species "+sp)
        print("Got "+str(allPathwaysForSp)+" pathways including a total of "+str(allGenesForSp)+" genes for species "+str(sp))
        print(pathwaysDf)

    pathwaysDf = pathwaysDf.drop_duplicates()
    print("Final dataframe to write:")
    print(pathwaysDf)
    # Creating an in-memory csv for the data
    mem_csv = StringIO()
    pathwaysDf.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the pathway table
    print("Writing "+str(len(pathwaysDf.index))+" pathway-gene associations to DB ")
    with closing(mem_csv) as csv_io:
        Pathway.objects.from_csv(csv_io)
