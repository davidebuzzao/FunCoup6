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


def getOrder(idA,idB):
    if idA.id>idB.id:
        return idA,idB
    else:
        return idB,idA

def getPair(protIdA,protIdB):
    if protIdA>protIdB:
        return protIdA+"||"+protIdB
    else:
        return protIdB+"||"+protIdA


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

def removeTempFile(filename):
    if os.path.isfile(filename):
        os.remove(filename)


def updateKEGGPathwaysToUse(filename, existingPathways, type):
    # Getting new KEGG pathways from the old manually picked set of pathways. Adding
    # newly added pathways to the set, and removing pathways that has been deprecated.
    ids={}
    allIds={}
    with open(filename, 'r') as in_file:
        if (type=="Metabolic"):
            pathwayType=""
            for line in in_file:
                if(line.startswith("A")):
                    pathwayType = line.split(">")[1].split("<")[0]
                if line.startswith("C"):
                    allIds[line.split()[1]]=" ".join(line.split()[2:])
                    # Only use pathways listed as Metabolic, that are not global maps(011) or chemical structure map (010)
                    if(pathwayType=="Metabolism" and not line.split()[1].startswith("011") and not line.split()[1].startswith("010")):
                        ids[line.split()[1]]=" ".join(line.split()[2:])
        elif(type == "Signaling"):
            # Find new pathways related to the signaling keywords
            signalingKeywords=["signal","transduction","potenti","junction","sensing","adhesio","interaction","contraction"]
            pathwayType=""
            for line in in_file:
                if(line.startswith("A")):
                    pathwayType = line.split(">")[1].split("<")[0]
                if line.startswith("C"):
                    allIds[line.split()[1]]=" ".join(line.split()[2:])
                    if any(ext in " ".join(line.split()[2:]) for ext in signalingKeywords):
                        # Exclude an pathways related to human diseases
                        if pathwayType !="Human Diseases":
                            ids[line.split()[1]]=" ".join(line.split()[2:])

    currentPathways = existingPathways.split(",")
    removed=[]
    added=[]
    for i in ids:
        if i not in currentPathways:
            added.append(i)
    # Add all pathways used in the previous version (that still exists)
    for c in currentPathways:
        if c not in ids:
            if c in allIds:
                ids[c]=allIds[c]
            else:
                removed.append(c)
    # Print the removed/added pathways
    print("\n###############################")
    print("Removed pathways: "+",".join(removed))
    print("Added pathways: "+",".join(added))
    print("###############################\n")
    return list(ids.keys())

keggSpeciesNameMap={'9615':'cfa'}

def getPathwayPairs(pathwaysToGet, goldStandardConfig, sp, pathwayType,proteomeConfig):
    # Download kgml files for a species and extract the relation links from it.
    # Then mapping the identifiers and writing the links to the DB
    pairs= [] # proteinA||proteinB: number of occurrences
    allGenesInPathway=[] # list of list of all genes in a pathway
    pathways=0
    speciesInDB= Species.objects.get(tax_id = sp)
    # Construct short kegg species name from the scientific name for the species in the DB
    if sp not in keggSpeciesNameMap:
        shortSpName=(speciesInDB.species_name.split()[0][0]+speciesInDB.species_name.split()[1][:2]).lower()
    else: # Special case for dog, as it can not be composed form the scientific name
        shortSpName=keggSpeciesNameMap[sp]
    # Iterating over each pathway in the set
    for pathway in pathwaysToGet:
        # Downloading pathway kgml file if it does not already exist
        if os.path.exists("data/tmp/"+shortSpName+pathway)==False:
            downloadFile(goldStandardConfig['urlKGML'].replace("PATHWAYID",shortSpName+pathway), "data/tmp/"+shortSpName+pathway)
        # Reading KGML file if the pathway exist for this species
        if os.path.exists("data/tmp/"+shortSpName+pathway):
            pathways+=1
            genesInPathway=set()
            pw = KGML_parser.read(open("data/tmp/"+shortSpName+pathway, 'r')) # Using kgml parser to read xml file
            # Iterating over all "gene" entries in the file to be able to use later when reading relations.
            genesMap={} # id in kgml file : list of gene name(s) as kegg ids
            for gene in pw.genes:
                genes=gene.name.split()
                genesMap[gene.id]=genes
                for g in genes:
                    genesInPathway.add(g)
            # Iterating over all relations in the pahway.
            for relation in pw.relations:
                if relation.type=="ECrel" or relation.type=="PPrel" or relation.type=="GErel": # Relation type must be one of these
                    genes1=[]
                    genes2=[]
                    # Each relation has two entries that can either be a gene or a group of genes
                    if relation.entry1.type=="group":
                        for component in relation.entry1.components:
                            genes1.extend(genesMap[component.id])
                    else:
                        if relation.entry1.id in genesMap:
                            genes1=genesMap[relation.entry1.id]
                    if relation.entry2.type=="group":
                        for component in relation.entry2.components:
                            genes2.extend(genesMap[component.id])
                    else:
                        if relation.entry2.id in genesMap:
                            genes2=genesMap[relation.entry2.id]
                    # Iterating over the genes of entry1 and entry2 and combining all possible pairs from them
                    if genes1!=[] and genes2!=[]:
                        for protIdA in genes1:
                            for protIdB in genes2:
                                if protIdA!=protIdB: # Excluding self-links
                                    pair = getPair(protIdA,protIdB)
                                    pairs.append(pair)
            allGenesInPathway.append(genesInPathway)
            # removeTempFile("data/tmp/"+shortSpName+pathway)
    # In case only using relations does not yield over 1000 links, all combinations of all genes in the pathway are used instead OR goldStandardConfig file says to not use the direct links
    if len(pairs)<1000 or goldStandardConfig['directLinks']==False:
        pairs=[]
        for genesInPathway in allGenesInPathway:
            for index,g1 in enumerate(list(genesInPathway)[:-1]):
                for g2 in list(genesInPathway)[index+1:]:
                    pair = getPair(g1,g2)
                    pairs.append(pair)

    # Creating a new gold standard to relate the new links to
    goldStandardInDB = GoldStandard(version=goldStandardConfig['version'], type=pathwayType, species=speciesInDB )
    goldStandardInDB.save()

    # Converting pairs to dataframe
    linksDf = pd.DataFrame(pairs, columns=['pairString'])
    linksDf['identifierA'] = linksDf['pairString'].apply(lambda x: x.split('||')[0])
    linksDf['identifierB'] = linksDf['pairString'].apply(lambda x: x.split('||')[1])
    linksDf = linksDf.drop('pairString',axis=1)

    # Getting all identifiers as originals AND with the prefix cut (hsa:123 and 123)
    identifiersDfA = linksDf[['identifierA']].copy()
    identifiersDfA['originalId']=linksDf['identifierA']
    identifiersDfB = pd.DataFrame(linksDf['identifierA'].apply(lambda x: x.split(':')[1]), columns=['identifierA'])
    identifiersDfB['originalId']=linksDf['identifierA']
    identifiersDfC = linksDf[['identifierB']].copy()
    identifiersDfC['originalId']=linksDf['identifierB']
    identifiersDfC.columns = ['identifierA','originalId']
    identifiersDfD = pd.DataFrame(linksDf['identifierB'].apply(lambda x: x.split(':')[1]))
    identifiersDfD['originalId']=linksDf['identifierB']
    identifiersDfD.columns = ['identifierA','originalId']
    allIdentifiers = pd.concat([identifiersDfA, identifiersDfB, identifiersDfC, identifiersDfD])
    allIdentifiers = allIdentifiers.drop_duplicates()

    # Collecting mappings from DB
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=sp))[0]
    availableProteinsDf = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id'))
    availableProteinsDf.columns = ['identifierA','id']
    # mapping the protein identifiers to the database ID
    allIdentifiers = pd.merge(allIdentifiers,availableProteinsDf, how='inner',on='identifierA')
    del availableProteinsDf
    allIdentifiers = allIdentifiers.drop('identifierA',axis=1)
    allIdentifiers.columns = ['identifierA','id']
    linksDf = pd.merge(linksDf,allIdentifiers, how='inner',on='identifierA')
    allIdentifiers.columns = ['identifierB','id']
    linksDf = pd.merge(linksDf,allIdentifiers, how='inner',on='identifierB')
    linksDf = linksDf.drop('identifierA',axis=1).drop('identifierB',axis=1)
    linksDf.columns = ['proteinA','proteinB']
    linksDf.proteinA, linksDf.proteinB = np.where(linksDf.proteinA < linksDf.proteinB, [linksDf.proteinB, linksDf.proteinA], [linksDf.proteinA, linksDf.proteinB]) # Put IDs in the correct order
    linksDf = linksDf.drop_duplicates()
    linksDf = linksDf[linksDf.proteinA != linksDf.proteinB]

    # getting proteinLinks
    #proteinLinksDf = pd.DataFrame.from_records(ProteinLink.objects.filter(Q(proteinA__proteome=proteome)).values('id', 'proteinA', 'proteinB'))
    #proteinLinksDf.columns = ['proteinLink','idA','idB']
    #linksDf = pd.merge(linksDf,proteinLinksDf, how='inner',on=['idA','idB'])
    linksDf['goldStandard']=goldStandardInDB.id
    linksDf.insert(0,'direction','0')
    # linksDf['direction']=0
    #linksDf = linksDf.drop('idA',axis=1).drop('idB',axis=1)

    # Creating an in-memory csv for the data
    mem_csv = StringIO()
    linksDf.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing "+str(len(linksDf.index))+" "+pathwayType+" links to DB for "+str(sp))
    with closing(mem_csv) as csv_io:
        GoldStandardLink.objects.from_csv(csv_io)
