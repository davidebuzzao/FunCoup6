import os
from sklearn.preprocessing import MinMaxScaler,PowerTransformer,StandardScaler
import yaml
from yaml.loader import SafeLoader
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
import math
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
import requests
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import copy
from contextlib import closing
from io import StringIO

def getOrder(idA,idB):
    if idA>idB:
        return idA,idB
    else:
        return idB,idA

def readDistanceMatrix(distanceMatrix):
    inparanoidbSpecies=[] # list of all available inparanoidb species
    distanceMatrixFile=open(distanceMatrix, "r")
    lineCounter=0
    matrix=[] # lower triangular of the distance matrix as a list of lists
    for line in distanceMatrixFile.readlines():
        if len(inparanoidbSpecies)==0: # The first line includes all inparanoidb species
            for taxid in line.split("\t"):
                if len(taxid)>2: # excluding empty or short strings
                    inparanoidbSpecies.append(taxid)
        else:
            rowForMatrix=[]
            for index,dist in enumerate(line.strip().split("\t")):
                if index<=lineCounter:
                    rowForMatrix.append(float(dist))
            lineCounter+=1
            matrix.append(rowForMatrix)
    distanceMatrixFile.close()
    return inparanoidbSpecies, matrix

def getProfilesFromUrl(evidenceConfig,inparanoidbSpecies,sp,proteomeConfig):
    profiles={} # proteinName:[list with as many positions as there are inparanoidb species. 1 will be set for a species that has an ortholog to the gene and 0 to a species that does not have an ortholog]
    # collect proteome evidenceConfig to be able to extract all proteins for the proteome used in this instance
    allProteinNames=list(Protein.objects.filter(Q(proteome__version=proteomeConfig['genome']['version']) & Q(proteome__species__tax_id=sp)).values_list('protein_id', flat=True)) # Collecting all proteins in the DB for the species
    print("Collecting profiles from URL for species: "+sp+" for "+str(len(allProteinNames))+" proteins")
    for protein in allProteinNames: # Iterating over all proteins for the species
        orthologsRequest = requests.get(evidenceConfig['orthologsUrl'].replace("PROTEIN",protein))
        if orthologsRequest.status_code==200: # If inparanoidb has orthologs for this protein, returns a list of all ortholog groups
            profiles[protein]=[0] * len(inparanoidbSpecies) # Initializing list of zeroes for all inparanoid species
            firstline=1
            for line in orthologsRequest.content.decode("utf-8").split("\n"): # Iterating over all lines in the file
                if len(line)>0:
                    if firstline==1: # Skipping header line
                        firstline=0
                    else:
                        foundSp=line.split(",")[2] # each taxid in the file has an ortholog to the protein in question
                        spIndex=inparanoidbSpecies.index(foundSp) # collecting the index of the found species
                        profiles[protein][spIndex]=1 # Setting found species to 1 in the proteins profile
        elif orthologsRequest.status_code!=500: # if error code is 500, it just means that the protein is not in inparanoiDB
            print("Cant download file from: "+evidenceConfig['orthologsUrl'].replace("PROTEIN",protein)+" error code: "+str(orthologsRequest.status_code))
            exit()
    return profiles

def getProfilesFromFile(inparanoidbSpecies,sp):
    print("Collecting profiles from File for species: "+sp)
    profiles={} # proteinName:[list with as many positions as there are inparanoidb species. 1 will be set for a species that has an ortholog to the gene and 0 to a species that does not have an ortholog]
    for speciesB in inparanoidbSpecies: # Iterating over all inparanoidb species
        if speciesB!=sp:
            # Finding the file containing orthologs for this species combination
            if os.path.isfile("data/tmp/resultsInParanoid/SQLtable."+sp+".fa-"+speciesB+".fa"):
                sqltableFile=open("data/tmp/resultsInParanoid/SQLtable."+sp+".fa-"+speciesB+".fa","r")
            elif os.path.isfile("data/tmp/resultsInParanoid/SQLtable."+speciesB+".fa-"+sp+".fa"):
                sqltableFile=open("data/tmp/resultsInParanoid/SQLtable."+speciesB+".fa-"+sp+".fa","r")
            else:
                print("Cant find file: "+"data/tmp/resultsInParanoid/SQLtable."+sp+".fa-"+speciesB+".fa")
                exit()
            spIndexThis=inparanoidbSpecies.index(sp) # Collecting profile index of the species to be collected
            spIndexB=inparanoidbSpecies.index(speciesB) # Collecting the profile index of the inparanoid species
            for line in sqltableFile.readlines():
                groupId,bitscore,species,inparalogScore,proteinId,*seedscore=line.strip().split("\t")
                if species.split(".")[0]==sp: # If the Species to be collected has an ortholog in the inParanoid species, setting its profile index to 1
                    if proteinId not in profiles:
                        profiles[proteinId]=[0] * len(inparanoidbSpecies) # Initializing list of zeroes for all inparanoid species
                        profiles[proteinId][spIndexThis]=1 # Setting found species to 1 in the proteins profile
                    profiles[proteinId][spIndexB]=1 # Setting found species to 1 in the proteins profile
            sqltableFile.close()
    return profiles

def calculatePhpScore(profileA,profileB,totalBranchLength,allLeafBranchLengths):
    # profileA,profileB = profileDict[pA],profileDict[pB]
    visitedNodesCorrelating=[] # List of all nodes visited by genes having 1 in both species positions
    visitedNodesDiffering=[] # List of all nodes visited by genes having 1 in either of the species positions
    existsInBoth,existingInTotal = 0, 0
    for spIndex,occurrenceSpA in enumerate(profileA): # Iterating over all species positions in the profile
        occurrenceSpB=profileB[spIndex]
        if occurrenceSpA==1 and occurrenceSpB==1:
            existsInBoth+=1
            existingInTotal+=1
            visitedNodesCorrelating.extend(allLeafBranchLengths[spIndex])
        elif (occurrenceSpA==1 and occurrenceSpB==0) or (occurrenceSpA==0 and occurrenceSpB==1):
            existingInTotal+=1
            visitedNodesDiffering.extend(allLeafBranchLengths[spIndex])

    branchLengthCorrelating=1 # Initiating a branch length for 1:1 (positive)
    branchLengthDiffering=1 # Initiating a branch length for 1:0 or 0:1 (negative)
    visitedNodesCorrelatingSet=set(visitedNodesCorrelating) # Only counting each internal branch once
    for n in visitedNodesCorrelatingSet:
        branchLengthCorrelating+=n.branch_length
    visitedNodesDifferingSet=set(visitedNodesDiffering) # Only counting each internal branch once
    for n in visitedNodesDifferingSet:
        if n not in visitedNodesCorrelatingSet: # Excluding internal nodes used in the positive branch length
            branchLengthDiffering+=n.branch_length

    scorePos=float(branchLengthCorrelating/totalBranchLength)
    scoreNeg=float(branchLengthDiffering/totalBranchLength)
    # combinedScore=scorePos-scoreNeg #old way
    combinedScore=np.log(scorePos/scoreNeg)
    return combinedScore, existsInBoth, existingInTotal

def calculatePhpScoresInParallel(chunk,numChunks,profileDict,totalBranchLength,allLeafBranchLengths):
    # Scoring pairs of proteins
    numGenes=len(list(profileDict.keys()))
    totalCalculations=(numGenes*(numGenes-1))/2
    calculationsPerChunk=round(totalCalculations/numChunks)
    startIndex=calculationsPerChunk*chunk
    stopIndex=(calculationsPerChunk*chunk)+calculationsPerChunk
    print("Thread "+str(chunk)+" start: "+str(startIndex)+" stop: "+str(stopIndex))

    dataToAdd, pairIndex, createdLinks = [], 0, 0
    for index,pA in enumerate(list(profileDict.keys())[:-1]): # Iterating over all combinations of genes (only including A-B not B-A)
        for pB in list(profileDict.keys())[index+1:]:
            pairIndex+=1
            if pairIndex>startIndex and pairIndex<=stopIndex:
                combinedScore, existsInBoth, existingInTotal = calculatePhpScore(profileDict[pA],profileDict[pB],totalBranchLength,allLeafBranchLengths)
                if len(str(existsInBoth))>0 and len(str(existingInTotal))>0:
                    metadata=str(existsInBoth)+","+str(existingInTotal)
                else:
                    metadata= None
                newRecord=[pA, pB, combinedScore, metadata]
                dataToAdd.append(newRecord)
                createdLinks+=1
    print("Done collecting PHP links. Thread "+str(chunk))
    linkDf = pd.DataFrame.from_records(dataToAdd)
    linkDf.columns = ["proteinA","proteinB","score","metadata"]

    linkDf.proteinA, linkDf.proteinB = np.where(linkDf.proteinA < linkDf.proteinB, [linkDf.proteinB, linkDf.proteinA], [linkDf.proteinA, linkDf.proteinB])
    
    return linkDf

def calculatePhpScoreForSpecies(sp,tree,inparanoidbSpecies,evidenceConfig,proteomeConfig):
    # Creating phylogenetic profiles
    if evidenceConfig['orthologCollectionMethod']=="url": # If collecting from url
        profiles = getProfilesFromUrl(evidenceConfig,inparanoidbSpecies,sp,proteomeConfig)
    elif evidenceConfig['orthologCollectionMethod']=="file": # If collecting from file
        profiles = getProfilesFromFile(inparanoidbSpecies,sp)
    else:
        print("ERROR: an orthologCollectionMethod must be defined in the evidenceConfig.yml file. Choose either url or file")
        exit()

    # Convert dict to dataframe with proteins as rows
    profileDf = pd.DataFrame.from_dict(profiles, orient='index').reset_index()
    profileDf.rename(columns={'index':'proteinIdentifier'}, inplace=True)

    # Getting mappings from the DB for the specified proteome
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=sp))[0]
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    available_protein_df.columns = ['proteinIdentifier','id']

    # mapping the protein identifiers to the database ID
    mapping = pd.DataFrame(list(profileDf['proteinIdentifier']))
    mapping.columns = ['proteinIdentifier']
    preMap = set(mapping['proteinIdentifier'])
    profileDf = pd.merge(profileDf,available_protein_df, how='inner',on='proteinIdentifier')
    postMap = set(profileDf['proteinIdentifier'])
    excluded = preMap - postMap
    print("Excluded "+str(len(excluded))+" genes, as they couldnt be mapped. Getting PHP scores for "+str(len(postMap))+" genes, for species: "+sp)
    profileDf = profileDf.drop(['proteinIdentifier'],axis=1).set_index('id') # Dropping the proteinIdentifier column, no longer needed, using peorein id from now on
    del available_protein_df

    # In case duplicate appeared in the mapping, dropping duplicate rows with most 1s in the profile
    profileDf['numPresence'] = profileDf.sum(axis=1)
    profileDf = profileDf.sort_values('numPresence', ascending=False)
    profileDf = profileDf.reset_index().drop_duplicates(subset='id', keep='first').set_index('id').drop(['numPresence'],axis=1)
    ## To balance the load on cpus later on, we shuffle rows
    profileDf = profileDf.sample(frac=1)
    profileDict = profileDf.T.to_dict('list')

    # Making a new tree to reroot with sp as root
    rootedSpTree = copy.deepcopy(tree)
    outgroup_clade = rootedSpTree.find_clades(name=sp)
    rootedSpTree.root_with_outgroup(outgroup_clade)
    totalBranchLength = rootedSpTree.total_branch_length()
    # Collecting list of branch lengths/paths from root for all species in advance, to avoid doing multiple get_path() when calculating score.
    allLeafBranchLengths={} 
    for i,s in enumerate(inparanoidbSpecies):
        allLeafBranchLengths[i] = rootedSpTree.get_path(s)

    print("Done with re-rooting")
    # Scoring pairs of proteins in parallel, and writing resulting links to the DB
    cpus=os.cpu_count()
    linkDfList = Parallel(n_jobs=cpus)(delayed(calculatePhpScoresInParallel)(chunk,cpus,profileDict,totalBranchLength,allLeafBranchLengths) for chunk in range(0, cpus))
    linkDf = pd.concat(linkDfList)
    if len(linkDf)>0:
        # linkDf.to_csv('PHP_logScore_%s.gz' %sp, sep="\t", compression='gzip', index=False)
        # MinMax scaling
        scores = linkDf['score'].values.reshape(-1, 1)
        minmax_scaler = MinMaxScaler(feature_range=(0, 1))
        minmax_scores = minmax_scaler.fit_transform(scores)
        linkDf['score'] = np.round(minmax_scores,decimals=4)

        # Creating new PHP evidence for the species
        speciesInDB = Species.objects.get(tax_id = sp)
        evidenceInDB = Evidence(type="PHP", species=speciesInDB, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidenceInDB.save()

        linkDf = linkDf.sample(frac=1).reset_index(drop=True)
        linkDf.insert(0,'evidence',evidenceInDB.id)

    return linkDf

def getPHPLinks(distanceMatrix, evidenceConfig, speciesToCollect,proteomeConfig):
    # Collects PHP information from inParanoiDB orthologs (vi URL or File). Creates a phylogenetic profile
    # as a dict with a list for each gene, where each list position represent an inparanoidb species. occurrence
    # of an ortholog in a species is represented by 1, while no orthologs in a species is represented by 0
    # A distance matrix containing distances between all inparanoidb species is collected from the inparanoid website.
    # A tree is built from the distance matrix, which for each species in the instance, is rooted at the species in question.
    # A score is then calculated from the phylogenetic profile for each combination of genes in the species as the total branch length
    # of the subtree of all species having an ortholog for both genes over the total branch length of the tree (positive score) minus
    # the total branch length of the subtree of all species where one of the genes has an ortholog over the total branch length of
    # the tree (negative score). All links are then written to the database, and the evidences they are attached to are returned.
    inparanoidbSpecies, matrix = readDistanceMatrix(distanceMatrix)
    dm = DistanceMatrix(inparanoidbSpecies, matrix) # loading the distance matrix
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm) # Creating a neighbor-joining tree from the distance matrix
    print("Done constructing tree")

    # to minimize the time spent on indexing the table PHP
    # compute species in increasing proteome size order
    ordered_speciesToCollect = []
    for s in speciesToCollect:
        proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=s))[0]
        num_protein =  Protein.objects.filter(Q(proteome=proteome)).count()
        ordered_speciesToCollect.append((s,num_protein))
    # Sort the list by proteome size
    ordered_speciesToCollect = sorted(ordered_speciesToCollect, key=lambda x: x[1])

    for s,num_protein in ordered_speciesToCollect:
        linkDf = calculatePhpScoreForSpecies(s,tree,inparanoidbSpecies,evidenceConfig,proteomeConfig)
        if len(linkDf)==0: continue
        mem_csv = StringIO()
        linkDf.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        print("Writing %i PHP links to DB for %s" %(len(linkDf),s))
        # Writing the csv to the evidenceLink table
        with closing(mem_csv) as csv_io:
            PHP.objects.from_csv(csv_io)


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
    print(instanceConfig['instance']['species'])
    
    speciesToCollect = instanceConfig['instance']['species']
    speciesToCollect = ['9606']
    evidenceConfig = evidenceConfig['PHP']
    distanceMatrix = "data/tmp/"+evidenceConfig['distanceMatrixFileName']
    getPHPLinks(distanceMatrix, evidenceConfig,speciesToCollect,proteomeConfig)
