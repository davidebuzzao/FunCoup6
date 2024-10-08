import os
import csv
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO
from joblib import Parallel, delayed


def getPair(protIdA,protIdB):
    if protIdA.startswith("CELE_"): # Correct for identifiers that are in the DB
        protIdA=protIdA.replace("CELE_","")
    if protIdB.startswith("CELE_"): # Correct for identifiers that are in the DB
        protIdB=protIdB.replace("CELE_","")
    if protIdA>protIdB: # Making sure the proteins are added in the same order every time
        return protIdA+"||"+protIdB
    else:
        return protIdB+"||"+protIdA

def writeOperonToDB(goldStandardConfig, sp, pairs, proteomeConfig):
    # Mapping and Writing operon pairs to the database.
    # First, making a new gold standard where all new links are attached
    speciesInDB = Species.objects.get(tax_id = sp)
    goldStandardInDB = GoldStandard(version=goldStandardConfig['version'], type="Operon", species=speciesInDB )
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
    linksDf.id_x, linksDf.id_y = np.where(linksDf.id_x < linksDf.id_y, [linksDf.id_y, linksDf.id_x], [linksDf.id_x, linksDf.id_y]) # Put IDs in the correct order
    linksDf = linksDf.drop_duplicates()
    linksDf = linksDf[linksDf.id_x != linksDf.id_y]
    linksDf.columns = ['proteinA','proteinB']

    # getting proteinLinks and prepping table for writing to the DB
    #proteinLinksDf = pd.DataFrame.from_records(ProteinLink.objects.filter(Q(proteinA__proteome=proteome)).values('id', 'proteinA', 'proteinB'))
    #proteinLinksDf.columns = ['proteinLink','id_x','id_y']
    #linksDf = pd.merge(linksDf,proteinLinksDf, how='inner',on=['id_x','id_y'])
    linksDf['goldStandard']=goldStandardInDB.id
    linksDf.insert(0,'direction','0')
    #linksDf = linksDf.drop('id_x',axis=1).drop('id_y',axis=1)

    # Creating an in-memory csv for the data
    mem_csv = StringIO()
    linksDf.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing "+str(len(linksDf.index))+" Operon links to DB for "+str(sp))
    with closing(mem_csv) as csv_io:
        GoldStandardLink.objects.from_csv(csv_io)

    return goldStandardInDB

def getOperonsFromFile(filename, goldStandardList, goldStandardConfig, speciesOfInterest, operonDbMap,proteomeConfig):
    # Getting operons from operonDB file
    allPairs={}
    for s in speciesOfInterest: # Defining a dict of operons for each species
        pairs=[]
        allPairs[s]=pairs
    # Iterating over lines in file
    with open(filename, 'r') as in_file:
        for line in csv.reader(in_file, delimiter='\t'):
            if line[1] in operonDbMap and operonDbMap[line[1]] in speciesOfInterest:
                pairs=allPairs[operonDbMap[line[1]]]
                prots=line[3].split(",")
                # Creating all possible pairs from operon members
                for index,protIdA in enumerate(prots[:-1]):
                    for protIdB in prots[index+1:]:
                        pair=getPair(protIdA,protIdB)
                        pairs.append(pair)
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    results = parallel(delayed(writeOperonToDB)(goldStandardConfig, sp, allPairs[sp],proteomeConfig) for sp in allPairs)

    goldStandardList.extend(results)
    return goldStandardList
