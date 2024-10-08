import requests
import pandas as pd
import os
import random
import time
import yaml
from yaml.loader import SafeLoader
from io import BytesIO
from zipfile import ZipFile
from urllib.request import urlopen
from joblib import Parallel, delayed
from io import StringIO
from tqdm import tqdm
from auxiliary_functions import *
# from Bio import Phylo
# from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import gzip

#fcSpeciesOrderedByGenomeCoverage=['559292','224308','9606','10090','83333','10116','9823','9615','9913','3702','7227','6239','284812','9031','7955','36329','44689','7719','39947','243232','273057','83332']

# Getting the FunCoup species used in the instance that is used for the website, as specified in the website config 
# And getting all inparanoid species to transfer networks to
def getSpecies():
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
    FCspecies=instanceConfig['instance']['species']

    speciesList = requests.get("https://inparanoidb.sbc.su.se/download/specieslist") # Getting species list from the inparanoid website
    taxIDsInInParanoiDB=[]
    speciesNameDict={}
    for speciesLine in speciesList.content.decode("utf-8").split("\n"):
        if speciesLine.startswith("#")==False and len(speciesLine)>0:
            taxIDsInInParanoiDB.append(speciesLine.split(",")[0])
            speciesNameDict[speciesLine.split(",")[0]]=speciesLine.split(",")[2]

    return taxIDsInInParanoiDB,FCspecies,speciesNameDict

# Getting the closest FunCoup species for all InParanoid species from a distance matrix
def getClosestFromDistanceMatrix(distanceMatrix, fcSpecies):
    inparanoidbSpecies=[] # list of all available inparanoidb species
    distanceMatrixFile=open(distanceMatrix, "r")
    closestFcSpecies={} # For storing closest
    distancesFCSpecies={} # For storing distances to all FC species from all inparanoid species
    for lineCounter, line in enumerate(distanceMatrixFile.readlines()):
        if len(inparanoidbSpecies)==0: # The first line includes all inparanoidb species
            for taxid in line.strip().split("\t"):
                if len(taxid)>2: # excluding empty or short strings
                    inparanoidbSpecies.append(taxid)
        else:
            otherSp=inparanoidbSpecies[lineCounter-1] # otherSP == InparanoidSpecies
            # TODO --> Transfer for FcSpecies for testing
            # if otherSp not in fcSpecies: # only transferring for non-FC species
            thisDistances={}
            closestSp=""
            closestDist=10000000
            for index,dist in enumerate(line.strip().split("\t")):
                dist=float(dist)
                s=inparanoidbSpecies[index]
                if s in fcSpecies: # Only considerng distances to the FC species
                    if s == otherSp: continue
                    thisDistances[s]=dist
                    if dist<closestDist:
                        closestSp=s
                        closestDist=dist
            closestFcSpecies[otherSp]=closestSp
            distancesFCSpecies[otherSp]=thisDistances
    distanceMatrixFile.close()
    return closestFcSpecies,distancesFCSpecies

# Getting information on taxonomy from NCBI. Obtaining all parent nodes for each given node in the taxonomy
def getNcbiTaxonomy(): 
    taxdump = urlopen("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip")
    taxdumpzip = ZipFile(BytesIO(taxdump.read()))
    parentNodes={}
    nodeRanks={}
    for line in taxdumpzip.open('nodes.dmp').readlines():
        line=line.decode('utf-8')
        parentNodes[line.split("|")[0].strip()]=line.split("|")[1].strip()
        nodeRanks[line.split("|")[0].strip()]=line.split("|")[2].strip()
    return parentNodes, nodeRanks

# Getting the lineage for a given species, by collecting its parent nodes
def getNcbiLineage(taxID, parentNodes):
    lineage=[]
    while taxID not in lineage and taxID in parentNodes:
        lineage.append(taxID)
        taxID=parentNodes[taxID]
    return lineage

# Getting ortholog pairs 
def getPairsFromGroup(spAGenes,spBGenes):
    orthologsToWrite=[]
    for gA in spAGenes:
        for gB in spBGenes:
            orthologsToWrite.append(gA+","+gB)
    return orthologsToWrite

# Reading SQLtable ortholog file
def readSqlTable(content,spA,spB):
    currentGroup=0
    spAGenes=[]
    spBGenes=[]
    orthologsToWrite=[]
    for line in content:
        if len(line)>1:
            groupId,bitscore,species,inparalogScore,proteinId,*seedscore=line.strip().split("\t")
            if groupId!=currentGroup: # For each new group encountered, writing all pairs of orthologs to list
                orthologsToWrite.extend(getPairsFromGroup(spAGenes,spBGenes))
                currentGroup=groupId
                spAGenes=[]
                spBGenes=[]
            if species.split(".")[0]==spA:
                spAGenes.append(proteinId.strip())
            else:
                spBGenes.append(proteinId.strip())
    orthologsToWrite.extend(getPairsFromGroup(spAGenes,spBGenes))
    return orthologsToWrite

# Reading distance matrix FOR TREE CONSTRUCTION (only used if distances from branch lengths in trees is applied) Not currently used
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

# Getting closes secies based on branch lengths in a tree (not currently used)
def getClosestFromTree(tree, fcSpecies,inparanoidSpecies):
    closest={}
    for sp in inparanoidSpecies:
        if sp not in fcSpecies:
            closestSp=""
            closestDist=1000000000
            for fcsp in fcSpecies:
                branchLengthDist = tree.distance(sp,fcsp)
                if branchLengthDist<closestDist:
                    closestSp=fcsp
                    closestDist=branchLengthDist
            closest[sp]=closestSp
    return closest

def getOrthologFiles(species,fcSpecies,speciesNameDict):
    fcSpeciesName=speciesNameDict[fcSpecies].split()[0][0].upper()+"."+speciesNameDict[fcSpecies].split()[1]
    speciesName=speciesNameDict[species].split()[0]+"_"+speciesNameDict[species].split()[1]

    orthoFile = "website/static/website/orthologFiles/"+species+":"+speciesNameDict[species].replace("/","").replace(" ","_")
    if not os.path.isfile(orthoFile):
        print(orthoFile)
        # time.sleep(random.randint(10,30))
        sqltableFile = requests.get("https://inparanoidb.sbc.su.se/download/sqltable/"+species+"&"+fcSpecies+"&prot")
        try: orthologsToWrite = readSqlTable(sqltableFile.content.decode("utf-8").split("\n"),species,fcSpecies)
        except: return
        fout = open(orthoFile, "w")
        fout.write("#"+fcSpecies+"\n")
        for o in orthologsToWrite:
            fout.write(o+"\n")
        fout.close()
    
    # if not os.path.isfile(orthoFile):
    #     print("Error, could not find ortholog file for FunCoup species: "+fcSpeciesName)
    

def getTransferNetwork(species,fcSpecies,speciesNameDict):
    # fcSpeciesName=speciesNameDict[fcSpecies].split()[0][0].upper()+"."+speciesNameDict[fcSpecies].split()[1]
    fcSpeciesName_name_clean = speciesNameDict[fcSpecies].split(' (')[0]
    fcSpeciesName = '%s.%s' %(fcSpeciesName_name_clean.split()[0][0],fcSpeciesName_name_clean.split()[-1])
    if fcSpeciesName=='O.japonica': fcSpeciesName = 'O.sativa'
    
    speciesName=speciesNameDict[species].split()[0]+"_"+speciesNameDict[species].split()[1]
    # speciesName_name_clean = speciesNameDict[species].split(' (')[0]
    # speciesName = '%s.%s' %(speciesName_name_clean.split()[0][0],speciesName_name_clean.split()[-1])

    orthoFile = "website/static/website/orthologFiles/"+species+":"+speciesNameDict[species].replace("/","").replace(" ","_")
    if not os.path.isfile(orthoFile):
        # print("Error, could not find ortholog file for species: "+species)
        return
    
    print('Working on %s, from %s' %(species,fcSpecies))

    # Read network for friend species, and translate to new network file - for FULL networks
    transferred_network_full_f = "website/static/website/networks/FunCoup6.0/FC6.0_TRANSFERRED_"+speciesName+"("+species+")"+"_from_"+fcSpeciesName+"_full.gz"        
    transferred_network_compact_f = "website/static/website/networks/FunCoup6.0/FC6.0_TRANSFERRED_"+speciesName+"("+species+")"+"_from_"+fcSpeciesName+"_compact.gz"        

    # network_f = 'data/tmp/orthologNetwork/FC6.0_s%s_t%s_full.gz' %(fcSpecies,species)
    network_f = "website/static/website/networks/FunCoup6.0/FC6.0_"+fcSpeciesName+"_full.gz"
    if all([os.path.isfile(transferred_network_full_f),os.path.isfile(transferred_network_compact_f)]):
        print("Already inferred species: "+speciesName)
    elif not all([os.path.isfile(transferred_network_full_f),os.path.isfile(transferred_network_compact_f)]) and os.path.isfile(network_f):
        ortho_df = pd.read_csv(orthoFile, comment='#', sep=',', header=None)
        ortho_df.columns = ['pA','ProteinA']

        network_df = pd.read_csv(network_f, comment='#', sep='\t')
        network_df.columns = [col.split(':')[1] for col in network_df.columns]
        network_df = network_df.drop(['FunCoupAid', 'FunCoupBid'],axis=1)

        network_df = pd.merge(network_df, ortho_df, on='ProteinA', how='inner')
        ortho_df.columns = ['pB','ProteinB']
        network_df = pd.merge(network_df, ortho_df, on='ProteinB', how='inner')
        network_df = network_df.drop(['ProteinA','ProteinB'],axis=1)
        
        # Rename the columns
        network_df = network_df.rename(columns={'pA': 'ProteinA', 'pB': 'ProteinB'})

        # Move the renamed columns to the first and second positions
        network_df.insert(0, 'ProteinA', network_df.pop('ProteinA'))
        network_df.insert(1, 'ProteinB', network_df.pop('ProteinB'))

        # Keep only maxPPV, remove self-interactions
        network_df = network_df.sort_values(by=['PPV','FBS_max'],ascending=False).reset_index(drop=True)
        network_df = network_df.drop_duplicates(subset=['ProteinA', 'ProteinB'], keep='first')
        network_df = network_df[network_df['ProteinA']!=network_df['ProteinB']].reset_index(drop=True)

        ## Rename columns
        col_names = network_df.columns.to_list()
        network_df.columns = [str(i) + ':' + col_name for i,col_name in enumerate(col_names)]

        ## Output full network
        network_df.to_csv(transferred_network_full_f, sep="\t", compression='gzip', index=False)
        
        ## Output compact network
        direction_col = [dir_col for dir_col in network_df.columns if 'direction' in dir_col]
        network_df.iloc[:,0:(9+len(direction_col))].to_csv(transferred_network_compact_f, sep="\t", compression='gzip', index=False)
    else:
        print("Error, could not find network file for FunCoup species: "+fcSpeciesName)
    
    return

def transferNetworks():
    print("Initiating network transfer")
    inparanoidSpecies,FCspecies,speciesNameDict = getSpecies()

    ## For getting distances from branch lengths from tree (NOT CURRENTLY USED)
    # inparanoidbSpecies, matrix = readDistanceMatrix("../data/tmp/distanceMatrix_inparanoidb9")
    # dm = DistanceMatrix(inparanoidbSpecies, matrix) # loading the distance matrix
    # constructor = DistanceTreeConstructor()
    # tree = constructor.upgma(dm) # Creating a tree from the distance matrix. can be either upgma or nj
    # closestTree = getClosestFromTree(tree, FCspecies,inparanoidSpecies)

    # Getting closest species from distance matrix
    closestInParanoid,distancesFCSpecies = getClosestFromDistanceMatrix("data/tmp/distanceMatrix_InParanoid8Method", FCspecies)

    # Getting closest species from NCBI lineages
    parentNodes, nodeRanks = getNcbiTaxonomy()
    lineagePerFcSpecies={}
    for species in FCspecies:
        lineage = getNcbiLineage(species,parentNodes)
        for l in lineage:
            lineagePerFcSpecies[l] = lineagePerFcSpecies.get(l,[])
            lineagePerFcSpecies[l].append(species)
    
    closest={}
    for species in inparanoidSpecies:
        if species in FCspecies:
            distances = distancesFCSpecies[species]
            shortestDistSp = min(distances, key=lambda k: distances[k])
            closest[species]=shortestDistSp
        else:
            ## Try to find the closest FC species according to NCBI lineage
            ## If not found, it will be added later
            lineage = getNcbiLineage(species,parentNodes)
            # [i for i in lineage if i in lineagePerFcSpecies]
            found=False
            for l in lineage:
                if found==False:
                    if l in lineagePerFcSpecies:
                        found=True # The lineage of parent nodes is reported in sorted order, from closest to furthest, that's why for loop stops ASAP
                        if len(lineagePerFcSpecies[l])>1: # If there are ties (there often is!!!) the closest species among the tied ones are selected from the distance-matrix-closest species
                            distances = {key: distancesFCSpecies[species][key] for key in lineagePerFcSpecies[l]}
                            shortestDistSp = min(distances, key=lambda k: distances[k])
                            
                            ## Alternatives to ortholog-based metrics
                            # indexesInOrderedFCSp=[]
                            # distance = tree.distance(species,candidate) # Use this in case ties should be resolved with branch distances from tree
                            #indexesInOrderedFCSp.append(fcSpeciesOrderedByGenomeCoverage.index(candidate)) # Use this in case ties should be resolved by taking the most complete FC network
                            #networkSP=fcSpeciesOrderedByGenomeCoverage[min(indexesInOrderedFCSp)] # Use this in case ties should be resolved by taking the most complete FC network
                            networkSP=shortestDistSp
                        else:
                            networkSP=lineagePerFcSpecies[l][0]
                        closest[species]=networkSP

    # Adding closest species from distanceMatrix for the species that are not in the NCBI taxonomy
    for species in closestInParanoid:
        if species not in closest:
            print("Missing from ncbi, adding closest from distanceMatrix: "+species+" -> "+closestInParanoid[species])
            closest[species]=closestInParanoid[species]

    # How to find out what species our 22 FC6 are closest to
    speciesNameDict_df = pd.DataFrame.from_dict(speciesNameDict, orient='index', columns=['spA_name'])
    speciesNameDict_df['spA'] =  speciesNameDict_df.index
    speciesNameDict_df = speciesNameDict_df[['spA','spA_name']].reset_index(drop=True)

    # Write orthologs to file - This is used when showing networks for transferred species on the website
    transferNetwork_f = 'data/tmp/orthologNetwork/NCBI&InParanoidSpAB.tsv'
    if os.path.isfile(transferNetwork_f):
        species_df = pd.read_csv(transferNetwork_f, sep="\t", dtype=str)
    else:
        species_df = pd.DataFrame.from_dict(closest, orient='index', columns=['spB'])
        species_df['spA'] = species_df.index
        species_df = species_df[['spA','spB']].reset_index(drop=True)
        species_df = pd.merge(species_df, speciesNameDict_df, on='spA', how='left')
        speciesNameDict_df.columns = ['spB','spB_name']
        species_df = pd.merge(species_df, speciesNameDict_df, on='spB', how='left')

        species_df['FC_species'] = species_df['spA'].isin(instanceConfig['instance']['species']).astype(int)
        species_df = species_df.sort_values('FC_species',ascending=True).reset_index(drop=True)
        species_df.to_csv('data/tmp/orthologNetwork/NCBI&InParanoidSpAB.tsv', sep="\t", index=False)
        # for sp in closest: # uncomment to print all the closest species
        #     print(" ".join(speciesNameDict[sp].split()[:2])+" ("+sp+") -> " + " ".join(speciesNameDict[closest[sp]].split()[:2])+" ("+closest[sp]+")")

    # species_df = species_df[species_df['FC_species']=='1'] #uncomment for assessment
    species_df = species_df[species_df['FC_species']=='0'] #uncomment for production
    closest = dict(zip(species_df['spA'], species_df['spB']))
    # closest={'122586':'83333','1123071':'83333','2043170':'83333','760142':'83333','1620215':'83333','387092':'83333'} # Using this for testing, comment or remove to get networks for ALL species

    for species in closest:
        getOrthologFiles(species,closest[species],speciesNameDict)

    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    # Extract ortholog files from InParanoiDB9 --> this does not work
    # with tqdm_joblib(tqdm(desc='Get orthologs', total=len(closest))) as progress_bar:
    #     parallel(delayed(getOrthologFiles)(species,closest[species],speciesNameDict) for species in closest)

    # Create transferred network from network files - This is used for downloading the entire network for a transferred species
    with tqdm_joblib(tqdm(desc='Transfer network', total=len(closest))) as progress_bar:
        parallel(delayed(getTransferNetwork)(species,closest[species],speciesNameDict) for species in closest)

    print("Done with network transfer!")



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
