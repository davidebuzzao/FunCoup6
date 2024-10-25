import os
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed

import random
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO
import math
import itertools

def checkForNegGSInDB(instanceSpecies, instanceGoldStandards, goldStandardConfig,goldStandardsForInstance):
    speciesAndNegGSToCollect={}
    combosToCollect=0
    for s in instanceSpecies:
        for g in instanceGoldStandards:
            negGSInDB = GoldStandard.objects.filter(Q(version=goldStandardConfig[g]['version']) &  Q(type=g+"_neg") &  Q(species__tax_id=s))
            if negGSInDB.exists(): # If the evidence already exist in the DB, using existing one for the instance
                goldStandardsForInstance.append(negGSInDB[0])
                print("Negative gold standard Exists: "+g+" v."+goldStandardConfig[g]['version']+" species: "+s)
            else: # Otherwise adding the species to fetch links
                posGSInDB = GoldStandard.objects.filter(Q(version=goldStandardConfig[g]['version']) &  Q(type=g) &  Q(species__tax_id=s))
                if posGSInDB.exists():
                    if s not in speciesAndNegGSToCollect:
                        speciesAndNegGSToCollect[s]=[]
                    speciesAndNegGSToCollect[s].append(posGSInDB[0])
                    combosToCollect+=1
    return speciesAndNegGSToCollect,goldStandardsForInstance,combosToCollect


def random_sampling(protein_id,seed_combinations,sizeList=int(1e6),directed=False):
    if directed:
        pass
        ## TODO --> proteins have to be tagged with prior knowledge (i.e. TF in TRRUST/Regnetwork)
        # tf_id_range = [i for i in range(protein_id[0],protein_id[1]+1)]
        # target_id_range = [i for i in range(protein_id[0],protein_id[1]+1)]
        # random.seed(seed_combinations[0])
        # proteinA = random.choices(tf_id_range, k=sizeList)
        # random.seed(seed_combinations[1])
        # proteinB = random.choices(target_id_range, k=sizeList)
        protein_df = pd.DataFrame({'proteinA':proteinA,'proteinB':proteinB})\
            .drop_duplicates().reset_index(drop=True)
    else:
        ## Sample (with replacement) proteinA and proteinB from the protein id space (i.e. range of protein id)
        protein_id_range = [i for i in range(protein_id[0],protein_id[1]+1)]
        random.seed(seed_combinations[0])
        proteinA = random.choices(protein_id_range, k=sizeList)
        random.seed(seed_combinations[1])
        proteinB = random.choices(protein_id_range, k=sizeList)
        ## Create a data frame with sorted proteins (i.e. A>B), removing possible duplicates due to sampling with replacement.
        protein_df = pd.DataFrame(\
            pd.DataFrame({'proteinA':proteinA,'proteinB':proteinB})\
                .apply(lambda r: sorted(r,reverse=True), axis = 1)\
                    .to_list(),columns=["proteinA","proteinB"])\
                    .drop_duplicates().reset_index(drop=True)

    return protein_df


def getNegativeGoldStandardsForThread(thread, combosPerThread, speciesAndNegGSToCollect,goldStandardsForInstance,seeds,proteomeConfig,instanceConfig,goldStandardConfig):
    goldStandardsForInstance=[]
    start=thread*combosPerThread
    stop=(thread*combosPerThread)+combosPerThread
    comboCount=0
    count = -1
    for speciesA in speciesAndNegGSToCollect:
        count +=1
        speciesA_name = Species.objects.filter(Q(tax_id=speciesA))[0].species_name
        proteome_A= Proteome.objects.filter(Q(species__tax_id=speciesA) & Q(version=proteomeConfig['genome']['version']))[0]
        ## Extract proteins_id from proteoms and make a range out of it
        first_protein_id, *_, last_protein_id = [p.id for p in list(Protein.objects.filter(Q(proteome=proteome_A)))]
        ## For each species sample 1M negative gold standard links (~1.3M is the largest positive gold standard).
        speciesA_negative_goldStandardLink = random_sampling((first_protein_id,last_protein_id),seeds[count],sizeList=int(1e6))
        ## For each positive gold standard, retrieve the species negative gold standard that contains no respective positive examples
        for goldStandard in speciesAndNegGSToCollect[speciesA]:
            comboCount+=1
            if comboCount>start and comboCount<=stop:
                # NOTE: This is temporary, to test if we can use Metabolic_direct instead of Metabolic_all (same for Signaling)
                ## If so, keep the tmp_goldstandard var, else remove it.
                tmp_goldstandard = goldStandard
                # if goldStandard.type=="Metabolic" and goldStandard.version == '2022-02': tmp_goldstandard = GoldStandard.objects.filter(Q(version='2022-02_ALL') &  Q(type=goldStandard.type) &  Q(species__tax_id=speciesA))[0]
                # if goldStandard.type=="Signaling" and goldStandard.version == '2022-02': tmp_goldstandard = GoldStandard.objects.filter(Q(version='2022-02_ALL') &  Q(type=goldStandard.type) &  Q(species__tax_id=speciesA))[0]

                ## Make a temporary copy of the negative gold standard
                negative_goldStandardLink = speciesA_negative_goldStandardLink.copy()
                #print('#################################################\nSpeciesA: %s\nProteins: %i\nGoldStandard: %s\nLinks-: %i\n' \
                #      %(speciesA_name,last_protein_id-first_protein_id,goldStandard.type,len(negative_goldStandardLink)))

                ## Filter out the positive links of the respective species, and also the species ortholog positive links.
                for speciesB in instanceConfig['instance']['species']:
                    spBGoldStandard = GoldStandard.objects.filter(Q(version=goldStandardConfig[goldStandard.type]['version']) &  Q(type=goldStandard.type) &  Q(species__tax_id=speciesB))
                    if not spBGoldStandard.exists(): continue
                    spBGoldStandard=spBGoldStandard[0]
                    tmp_spBgoldstandard = spBGoldStandard
                    # if spBGoldStandard.type=="Metabolic" and spBGoldStandard.version == '2022-02': tmp_spBgoldstandard = GoldStandard.objects.filter(Q(version='2022-02_ALL') &  Q(type=spBGoldStandard.type) &  Q(species__tax_id=speciesB))[0]
                    # if spBGoldStandard.type=="Signaling" and spBGoldStandard.version == '2022-02': tmp_spBgoldstandard = GoldStandard.objects.filter(Q(version='2022-02_ALL') &  Q(type=spBGoldStandard.type) &  Q(species__tax_id=speciesB))[0]

                    speciesB_name = Species.objects.filter(Q(tax_id=speciesB))[0].species_name

                    if speciesB == speciesA:
                        # Extract positive gold standard links from species A
                        positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=tmp_goldstandard)).values('proteinA','proteinB'))
                        if (len(positive_goldStandardLink)==0):
                            continue
                        else:
                            positive_goldStandardLink.columns = ['proteinA','proteinB']
                    else:
                        # Extract positive gold standard links from species B
                        tmp_positive_goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=tmp_spBgoldstandard)).values('proteinA','proteinB'))
                        if (len(tmp_positive_goldStandardLink)==0):
                            continue
                        else:
                            tmp_positive_goldStandardLink.columns = ['proteina','proteinb']
                            # Extract orthologs between speciesA and speciesB
                            if int(speciesA)>int(speciesB):
                                ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesA) & Q(speciesB=speciesB)).values('proteinA','proteinB'))
                            else:
                                ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesB) & Q(speciesB=speciesA)).values('proteinB','proteinA'))
                            ortholog_protein_df.columns = ['proteinA','proteina']

                            ## Extract ortholog links by:
                            ## - first filter out links with non present homolog
                            ## - second merge and convert to homologs
                            tmp_positive_goldStandardLink = tmp_positive_goldStandardLink[tmp_positive_goldStandardLink['proteina'].isin(ortholog_protein_df['proteina']) & \
                                                                                        tmp_positive_goldStandardLink['proteinb'].isin(ortholog_protein_df['proteina'])]

                            tmp_merged = pd.merge(tmp_positive_goldStandardLink, ortholog_protein_df, on="proteina", how='inner')
                            ortholog_protein_df.columns = ['proteinB','proteinb']
                            tmp_merged = pd.merge(tmp_merged, ortholog_protein_df, on="proteinb", how='inner')

                            positive_goldStandardLink = pd.DataFrame(\
                                tmp_merged.drop(['proteina','proteinb'], axis=1)\
                                    .apply(lambda r: sorted(r,reverse=True), axis = 1)\
                                        .to_list(),columns=["proteinA","proteinB"])\
                                            .drop_duplicates()\
                                                .reset_index(drop=True)

                    if (len(positive_goldStandardLink)==0):
                            #print('SpeciesB: %s\nGoldStandard: %s\nLinks+: %i\nLinks-: %i\n' \
                            #%(speciesB_name,goldStandard.type,0,len(negative_goldStandardLink)))
                            continue
                    else:
                        negative_goldStandardLink = \
                            pd.merge(negative_goldStandardLink,positive_goldStandardLink, indicator=True, how='outer')\
                                .query('_merge=="left_only"')\
                                .drop('_merge', axis=1)

                        #print('SpeciesB: %s\nGoldStandard: %s\nLinks+: %i\nLinks-: %i\n' \
                            #%(speciesB_name,goldStandard.type,len(positive_goldStandardLink),len(negative_goldStandardLink)))

                speciesInDB = Species.objects.get(tax_id = speciesA)
                negGoldStandardInDB = GoldStandard(version=goldStandard.version, type=goldStandard.type+"_neg", species=speciesInDB)
                negGoldStandardInDB.save()

                # getting proteinLinks and prepping table for writing to the DB
                #proteinLinksDf = pd.DataFrame.from_records(ProteinLink.objects.filter(Q(proteinA__proteome=proteome_A)).values('id', 'proteinA', 'proteinB'))
                #proteinLinksDf.columns = ['proteinLink','proteinA','proteinB']
                #linksDf = pd.merge(negative_goldStandardLink,proteinLinksDf, how='inner',on=['proteinA','proteinB'])
                negative_goldStandardLink['goldStandard']=negGoldStandardInDB.id
                negative_goldStandardLink['direction'] = 0
                #linksDf = linksDf.drop('proteinA',axis=1).drop('proteinB',axis=1)

                # Creating an in-memory csv for the data
                mem_csv = StringIO()
                negative_goldStandardLink.to_csv(mem_csv, index=False)
                mem_csv.seek(0)
                # Writing the csv to the negatieGoldStandardLink table
                print("THREAD: "+str(thread)+" Writing "+str(len(negative_goldStandardLink.index))+" Negative "+goldStandard.type+" Gold Standard Links for: "+str(speciesA))
                doWrite=True
                numberOfAttempts=0
                while doWrite==True and numberOfAttempts<=10:
                    numberOfAttempts+=1
                    try:
                        with closing(mem_csv) as csv_io:
                            NegativeGoldStandardLink.objects.from_csv(csv_io)
                        doWrite=False
                    except:
                        print("Unexpected error when attempting to write NGS to DB. Trying again, attempts: "+str(numberOfAttepmts))
                        doWrite=True
                goldStandardsForInstance.append(negGoldStandardInDB)
    return goldStandardsForInstance


def generateNegativeGoldStandard(instanceConfig,goldStandardConfig,proteomeConfig):
    # We can generate one big pool of negative examples and ajust that to positive gold standards.
    # To ensure reproducibility (as long as gold standard and proteomes do not change)
    # we generate as many tuples of random seeds as the number of species.
    # NOTE: Enqueue new species in the config file to maintain the results reproducible for old species.
    print("\nGenerating Negative gold standards\n")
    goldStandardsForInstance=[]
    speciesAndNegGSToCollect,goldStandardsForInstance,combosToCollect = checkForNegGSInDB(instanceConfig['instance']['species'], instanceConfig['instance']['goldStandard'], goldStandardConfig,goldStandardsForInstance)

    random.seed(12345)
    seedsA = random.sample(range(int(1e6)),k=len(instanceConfig['instance']['species']))
    random.seed(54321)
    seedsB = random.sample(range(int(1e6)),k=len(instanceConfig['instance']['species']))
    seeds = []
    for i,j in zip(seedsA,seedsB):
        seeds.append([i,j])

    CPUS=os.cpu_count()
    combosPerThread=math.ceil(combosToCollect/CPUS)
    goldStandardsForInstancePerThread = Parallel(n_jobs=CPUS)(delayed(getNegativeGoldStandardsForThread)(thread, combosPerThread, speciesAndNegGSToCollect,goldStandardsForInstance,seeds,proteomeConfig,instanceConfig,goldStandardConfig) for thread in range(0, CPUS))
    goldStandardsForInstance = list(itertools.chain.from_iterable(goldStandardsForInstancePerThread))

    return goldStandardsForInstance
