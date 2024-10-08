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


# Approved interaction types includes all subclasses of molecular association/association/physical association/direct interaction
approvedInteractionTypes=["MI:2232","MI:0914","MI:0915","MI:0407","MI:0195","MI:0408","MI:0556","MI:1126","MI:1127","MI:0414","MI:0192","MI:2356","MI:0557",
"MI:0193", "MI:1143","MI:1148","MI:1139","MI:0194","MI:1355","MI:0212","MI:0198","MI:0200","MI:0201","MI:0202","MI:0910","MI:0572","MI:0902","MI:0571","MI:0570",
"MI:1310","MI:0197","MI:2280","MI:0985","MI:1140","MI:0199","MI:0558","MI:0871","MI:0569","MI:0203","MI:0568","MI:0204","MI:1027","MI:0207","MI:0559","MI:2252",
"MI:0210","MI:1250","MI:1251","MI:1237","MI:0211","MI:0206","MI:0209","MI:0214","MI:0216","MI:0213","MI:0567","MI:0986","MI:0701","MI:0987","MI:0881","MI:0882",
"MI:0883","MI:0945","MI:1146","MI:0971","MI:0217","MI:0844","MI:2263","MI:1327","MI:0566","MI:2272","MI:2273","MI:0220","MI:1230"]


def getPair(idA,idB):
    if idA>idB: # Making sure the proteins are added in the same order every time
        return idA+"||"+idB
    else:
        return idB+"||"+idA

def writePPIsToDB(version, s, pairs, goldStandardsForInstance,proteomeConfig):
    # Creating a new gold standard
    speciesInDB = Species.objects.get(tax_id = s)
    goldStandardInDB = GoldStandard(version=version, type="PPI", species=speciesInDB )
    goldStandardInDB.save()

     # Converting pairs to dataframe
    linksDf = pd.DataFrame(list(pairs.keys()), columns=['pairString'])
    linksDf['identifierA'] = linksDf['pairString'].apply(lambda x: x.split('||')[0])
    linksDf['identifierB'] = linksDf['pairString'].apply(lambda x: x.split('||')[1])
    linksDf['numPublications'] = linksDf['pairString'].apply(lambda x: pairs[x])
    linksDf = linksDf.drop('pairString',axis=1)
    linksDf = linksDf.drop_duplicates()

    # Collecting mappings from DB
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=s))[0]
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
    linksDf = linksDf.rename(columns={"id_x": "proteinA", "id_y": "proteinB"})

    # getting proteinLinks and prepping table for writing to the DB
    #proteinLinksDf = pd.DataFrame.from_records(ProteinLink.objects.filter(Q(proteinA__proteome=proteome)).values('id', 'proteinA', 'proteinB'))
    #proteinLinksDf.columns = ['proteinLink','id_x','id_y']
    #linksDf = pd.merge(linksDf,proteinLinksDf, how='inner',on=['id_x','id_y'])
    linksDf['goldStandard']=goldStandardInDB.id
    linksDf.insert(0,'direction','0')

    #linksDf = linksDf.drop('id_x',axis=1).drop('id_y',axis=1)

    # Excluding links that has less than 2 publications AND does not overlap with existing gold standard links
    existingGoldStandardLinks= pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard__in=goldStandardsForInstance)).values('proteinA','proteinB'))
    linksToAddA = linksDf[linksDf.numPublications > 1].reset_index(drop=True) # Adding links that have more than one publication
    linksToAddB = pd.merge(linksDf,existingGoldStandardLinks, how='inner',on=['proteinA','proteinB']) # OR overlaps with existing gold standard links
    linksToAdd = pd.concat([linksToAddA, linksToAddB]).drop('numPublications',axis=1).drop_duplicates().reset_index(drop=True)

    # Creating an in-memory csv for the data
    mem_csv = StringIO()
    linksToAdd.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing "+str(len(linksToAdd.index))+" PPI links to DB for "+str(s))
    with closing(mem_csv) as csv_io:
        GoldStandardLink.objects.from_csv(csv_io)

    return goldStandardInDB

def extractFomAllFile(iref_file, speciesOfInterest,version,goldStandardsForInstance,proteomeConfig):
    # Extracting links from irefindex file. Filters on:
        # Requires edgetype X (pairwise interaction)
        # TaxIDA and TaxIDB must be the same species
        # Interaction type must be a physical interaction (see approvedInteractionTypes list)
        # Must have uniprotID as identifiers
        # Must have at least 2 different publications supporting the interaction OR exist among the already added gold standard links
        # must have lpr<100
    ppiDicts={}
    for s in speciesOfInterest:
        pairs={} # protA||protB : numPublications
        ppiDicts[s]=pairs
    print("Start reading from iRefIndex File")
    # Reading all lines in the irefindex file (this file includes all species available in irefindex)
    with open(iref_file, 'r') as in_file:
        for line in csv.reader((x.replace('\0', '') for x in in_file), delimiter='\t'):
            if "#uidA"!=line[0]: # Skipping header line
                edgetype=line[52]
                if edgetype=="X": # Only getting pairwise interactions, not complexes
                    taxa=line[9].split(":")[1].split("(")[0]
                    taxb=line[10].split(":")[1].split("(")[0]
                    # CAN NOT limit on host taxID, as too many interactions are excluded.
                    if taxa in speciesOfInterest and taxa==taxb: # Checking that the two interacting proteins belong to the same species
                        interactionType=line[11]
                        if "MI:" in interactionType:
                            interactionType=interactionType.split("\"")[1]
                        # Extracting uniprot identifiers, can be in multiple different columns
                        if "uniprotkb" in line[38]:
                            idA=line[38].split("uniprotkb:")[1].split("|")[0].split("-")[0]
                        elif "uniprotkb" in line[2]:
                            idA=line[2].split("uniprotkb:")[1].split("|")[0].split("-")[0]
                        elif "uniprotkb" in line[0]:
                            idA=line[0].split("uniprotkb:")[1].split("-")[0]
                        else:
                            idA=""
                        if "uniprotkb" in line[39]:
                            idB=line[39].split("uniprotkb:")[1].split("|")[0].split("-")[0]
                        elif "uniprotkb" in line[3]:
                            idB=line[3].split("uniprotkb:")[1].split("|")[0].split("-")[0]
                        elif "uniprotkb" in line[1]:
                            idB=line[1].split("uniprotkb:")[1].split("-")[0]
                        else:
                            idB=""
                        # Adding pairs that have two uniprot identifiers. This excludes some irefindex ids, however a vast majority of these seems depracated
                        # For each added pair, adding the max numPublications found for it in the iRefIndex file.
                        if idA!="" and idB!="":
                            lpr=1
                            numPublications=1
                            if line[14]!="-":
                                numPublications=int(line[14].split("np:")[1].split("|")[0])
                                lpr=int(line[14].split("lpr:")[1].split("|")[0])
                            pairs = ppiDicts[taxa]
                            pair=getPair(idA,idB)
                            if interactionType in approvedInteractionTypes and lpr <= 100:  # Making sure pair is direct association with lowest number of distinct interactions in an experiment to be lower than 100 -> more specific
                                if pair in pairs:
                                    if numPublications>pairs[pair]:
                                        pairs[pair]=numPublications
                                else:
                                    pairs[pair]=numPublications

    # Iterating over all collected pairs for all species
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    results = parallel(delayed(writePPIsToDB)(version, s, ppiDicts[s], goldStandardsForInstance,proteomeConfig) for s in speciesOfInterest)
    goldStandardsForInstance.extend(results)
    return goldStandardsForInstance
