import csv
import os
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

def getOrder(idA,idB):
    if idA.id>idB.id:
        return idA,idB
    else:
        return idB,idA

def getPair(protIdA,protIdB):
    if protIdA>protIdB: # making sure the pair ends up in the same order every time
        return protIdA+"||"+protIdB
    else:
        return protIdB+"||"+protIdA

def writeComplexesToDB(sp,pairs, version,proteomeConfig):
    # for sp in allPairsPerSpecies:
    # Making new gold standard for each species
    speciesInDB = Species.objects.get(tax_id = sp)
    goldStandardInDB = GoldStandard(version=version, type="Complex", species=speciesInDB )
    goldStandardInDB.save()
    # Converting pairs to dataframe
    linksDf = pd.DataFrame(pairs, columns=['pairString'])
    linksDf['identifierA'] = linksDf['pairString'].apply(lambda x: x.split('||')[0])
    linksDf['identifierB'] = linksDf['pairString'].apply(lambda x: x.split('||')[1])
    linksDf = linksDf.drop('pairString',axis=1)
    linksDf = linksDf.drop_duplicates()

    # Collecting mappings from DB
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=sp))[0]
    availableProteinsDf = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id'))
    availableProteinsDf.columns = ['identifierA','id']
    # mapping the protein identifiers to the database ID
    linksDf = pd.merge(linksDf,availableProteinsDf, how='inner',on='identifierA')
    availableProteinsDf.columns = ['identifierB','id']
    linksDf = pd.merge(linksDf,availableProteinsDf, how='inner',on='identifierB')
    del availableProteinsDf
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

    #linksDf = linksDf.drop('idA',axis=1).drop('idB',axis=1)

    # Creating an in-memory csv for the data
    mem_csv = StringIO()
    linksDf.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing "+str(len(linksDf.index))+" Complex links to DB for "+str(sp))
    with closing(mem_csv) as csv_io:
        GoldStandardLink.objects.from_csv(csv_io)

    return goldStandardInDB


def extractComplexFomiRefIndex(iref_file, speciesOfInterest,allPairsPerSpecies):
    # Extracting complexes from iRefIndex file. Since all species are in the same file,
    # a dictionary is created where each species to fetch gets its own pairs dict
    print("Getting iRefIndex complexes from file for all species")
    # Generating dictionaries for all species
    complexDicts=[]
    for s in speciesOfInterest:
        complexDict={}
        complexDicts.append(complexDict)
    # Iterating over iRefIndex file
    with open(iref_file, 'r') as in_file:
        for line in csv.reader((x.replace('\0', '') for x in in_file), delimiter='\t'):
            if "#uidA"!=line[0]: # Skipping header line
                taxID=""
                if('complex' in line[0]): # only checking for complexes
                    if ('uniprot' in line[1]):
                        if "taxid:" in line[10] :
                            taxID=line[10].split(":")[1].split("(")[0]
                            if taxID in speciesOfInterest: # checking for taxIds included in our species list
                                name = line[0].split(':')[1]
                                gene = line[1].split(':')[1].split("-")[0]# This maps protein ids XXX-1 to id XXX
                                dictComplex = complexDicts[speciesOfInterest.index(taxID)]
                                if name not in dictComplex:
                                    dictComplex[name] = [gene]
                                else:
                                    if gene not in dictComplex[name]:
                                        dictComplex[name].append(gene)
    # Iterating over all complexes and generating pairs for all combinations of members
    indexInSpeciesList=0
    print("# TaxID,iRefIndex complexes")
    for dictComplex in complexDicts:
        allPairs=[]
        print(speciesOfInterest[indexInSpeciesList]+","+str(len(dictComplex)))
        for i in dictComplex:
            if len(dictComplex[i])<100: # Only including complexes with less than 100 members
                for index,protIdA in enumerate(dictComplex[i][:-1]):
                    for protIdB in dictComplex[i][index+1:]:
                        pair=getPair(protIdA,protIdB)
                        allPairs.append(pair)
        allPairsPerSpecies[speciesOfInterest[indexInSpeciesList]]=allPairs
        indexInSpeciesList+=1
    return allPairsPerSpecies


def extractComplexesFromCorum(corum_file, speciesOfInterest, allPairsPerSpecies,corumSpeciesMap):
    # Extracting complexes from the corum file
    print("Getting Corum complexes from file for all species")
    # Generating dictionaries for all species to be fetched
    collect=False
    complexDicts={}
    for s in speciesOfInterest:
        if s in corumSpeciesMap.values():
            collect=True
            complexDict={}
            complexDicts[s]=complexDict
    if collect:
        # Iterating over the lines in the file
        with open(corum_file, 'r') as in_file:
            for line in csv.reader( in_file, delimiter='\t'):
                if line[2] in corumSpeciesMap and corumSpeciesMap[line[2]] in complexDicts: # checking that the species is among the species to fetch
                    species = corumSpeciesMap[line[2]]
                    complexName=line[1]
                    uniProtIds = line[5].split(';')
                    dictComplex = complexDicts[species]
                    if complexName not in dictComplex:
                        dictComplex[complexName] = uniProtIds
                    else: # if complex name alredy exist, keeps adding to the list of members
                        for id in uniProtIds:
                            dictComplex[complexName].append(id)

    # Iterating over all collected complexes and making pairs for all combinations of members
    print("# TaxID,Corum complexes")
    for sp in complexDicts:
        speciesDict=complexDicts[sp]
        allPairs=[]
        print(sp+","+str(len(speciesDict)))
        for i in speciesDict:
            for index,protIdA in enumerate(speciesDict[i][:-1]):
                for protIdB in speciesDict[i][index+1:]:
                    pair=getPair(protIdA,protIdB)
                    allPairs.append(pair)
        if sp not in allPairsPerSpecies:
            allPairsPerSpecies[sp]=allPairs
        else: # if this species already have complex links (from irefindex)
            allPairsPerSpecies[sp].extend(allPairs)
    return allPairsPerSpecies


def extractComplexesFromComplexPortal(speciesOfInterest, allPairsPerSpecies):
    # Extracting complexes from complex portal
    print("Getting CpmplexPortal complexes ")
    complexDicts={}
    for s in speciesOfInterest: # Reading the complexPortal file for each species
        complexDict={}
        with open("data/tmp/"+s+".tsv", 'r') as in_file:
            for line in csv.reader( in_file, delimiter='\t'):
                if not "#Complex" in line[0] and len(line)>0: # Skipping header line and empty lines
                    complexName = line[0]
                    uniProtIds = line[4].split('|')
                    uniProtIDsClean=[]
                    for id in uniProtIds: # cleaning up protein ids to be able to map them later
                        if (not "CHEBI" in id) and ("(" in id):
                            if "-PRO_" in id:
                                uniProtIDsClean.append(id.split("(")[0].split("-")[0])
                            else:
                                uniProtIDsClean.append(id.split("(")[0])
                    if len(uniProtIDsClean) > 0:
                        if(complexName not in complexDict):
                            complexDict[complexName] = uniProtIDsClean
                        else: # if complex name alredy exist, keeps adding to the list of members
                            for id in uniProtIDsClean:
                                complexDict[complexName].append(id)
        complexDicts[s]=complexDict

    # Iterating over all collected complexes and making pairs for all combinations of members
    print("# TaxID,complexPortal complexes")
    for sp in complexDicts:
        speciesDict=complexDicts[sp]
        allPairs=[]
        print(sp+","+str(len(speciesDict)))
        for i in speciesDict:
            for index,protIdA in enumerate(speciesDict[i][:-1]):
                for protIdB in speciesDict[i][index+1:]:
                    pair=getPair(protIdA,protIdB)
                    allPairs.append(pair)
        if sp not in allPairsPerSpecies:
            allPairsPerSpecies[sp]=allPairs
        else: # if this species already have complex links (from irefindex)
            allPairsPerSpecies[sp].extend(allPairs) # only adding new pairs, that does not already exist
    return allPairsPerSpecies
