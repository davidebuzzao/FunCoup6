import requests
import random
import os
from sklearn.preprocessing import MinMaxScaler
import networkx
import obonet
import math
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

# Code adapted from https://github.com/cbrl-nuces/GOntoSim
# Documentation on https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html

def downloadFile(urlToDownload, filename):
    # Downloading file from url
    url = urlToDownload
    myfile = requests.get(url)
    if myfile.status_code==200:
        open(filename, "wb").write(myfile.content)
    else:
        print("Cant download file from: "+url)

is_a = 0.8
part_of = 0.6

def intersection(lst1, lst2): 
	'''Helper Function to find intersecting terms from the two input lists of (term, s_Value)'''
	da = {v:k for v, k in lst1}
	db = {v:k for v, k in lst2} 
	return [(da[k],db[k]) for k in da.keys() & db.keys()]

def semanticValue(go_id, allPathsToRoot):
    all_all_paths = allPathsToRoot[go_id]
    S_values = list()
    for index, path in enumerate(all_all_paths):
        S_values.append([])
        path.reverse()
        for idx, term in enumerate(path):
            if idx == 0:
                S_values[index].append((go_id, 1))
            if idx < len(path)-1:
                if term['relationship'] != {}:
                    if 'part_of' in term['relationship']:
                        if path[idx+1] in term['relationship']['part_of']:
                            S_values[index].append((path[idx+1]['id'], S_values[index][idx][1] * part_of))
                        else: 
                            S_values[index].append((path[idx+1]['id'],S_values[index][idx][1] * is_a))
                    else: 
                        S_values[index].append((path[idx+1]['id'], S_values[index][idx][1] * is_a))
                else: 
                    S_values[index].append((path[idx+1]['id'],S_values[index][idx][1] * is_a))
    unique_terms_s_values = []
    for path in S_values:
        for term in path:
            unique_terms_s_values.append(term)
    unique_terms_s_values = sorted(unique_terms_s_values, key=lambda x: x[0])
    _s_values = {}
    for y, x in unique_terms_s_values: 
        if y in _s_values: 
            _s_values[y].append((y,x)) 
        else: 
            _s_values[y] = [(y, x)]
    final_s_values = []
    for node in _s_values:
        final_s_values.append(max(_s_values[node]))
    return final_s_values

def similarityOfGoTerms(go_id1, go_id2, allPathsToRoot):
    sv_a = semanticValue(go_id1, allPathsToRoot)
    sv_b = semanticValue(go_id2, allPathsToRoot)
    intersecting_terms = intersection(sv_a,sv_b)
    numerator = sum([x for t in intersecting_terms for x in t])
    denominator = sum(x for y,x in sv_a) + sum(x for y,x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
    if denominator > 0:
        Similarity = (numerator/denominator)
    else: 
        Similarity = 0
    return Similarity

def getSemanticSimilarity(termsA,termsB,allPathsToRoot):
    Sim1 = []
    Sim2 = []
    for idx, goterm in enumerate(termsA):
        Sim1.append([])
        for goid in termsB:
            Sim1[idx].append((goterm, goid,(similarityOfGoTerms(goterm, goid, allPathsToRoot))))
    for idx, goterm in enumerate(termsB):
        Sim2.append([])
        for goid in termsA:
            Sim2[idx].append((goterm, goid,(similarityOfGoTerms(goterm, goid, allPathsToRoot))))
    sem1 = []
    sem2 = []
    for index, goterm in enumerate(Sim1):	
        sem1.append((max(Sim1[index], key=lambda x: x[2])))
    for index, goterm in enumerate(Sim2):	
        sem2.append((max(Sim2[index], key=lambda x: x[2])))

    # Combine methods: Best-Match Average strategy
    similarity = (sum(x[2] for x in sem1)+ sum(x[2] for x in sem2) )/(len(termsA) + len(termsB))
    return similarity

def calculateSCLScore(geneAIndexes,geneBIndexes,weightLookup,numTerms):
    weightedFreqA={"0":0,"1":0} # Initializing all nonexisting(0) and existing(1) with 0 for geneA
    weightedFreqB={"0":0,"1":0} # Initializing all nonexisting(0) and existing(1) with 0 for geneB
    existingNonexistingCombos=["00","01","10","11"]
    weightedFreqAB={"00":0,"01":0,"10":0,"11":0} # Initializing all combinations of existing/nonexisting for geneA-geneB with 0
    ABindexes = np.unique(np.concatenate((geneAIndexes,geneBIndexes),0))
    for termIndex in ABindexes:
        existingA="1" if termIndex in geneAIndexes else "0"
        existingB="1" if termIndex in geneBIndexes else "0"
        wa=weightLookup[str(termIndex)+"|"+existingA]
        wb=weightLookup[str(termIndex)+"|"+existingB]
        weightedFreqA[existingA]+=wa
        weightedFreqB[existingB]+=wb
        weightedFreqAB[str(existingA)+str(existingB)]+=(wa*wb)

    # Divide all frequencies by the Total number of terms
    weightedFreqA["1"]=weightedFreqA["1"]/numTerms
    weightedFreqA["0"]=weightedFreqA["0"]/numTerms
    weightedFreqB["1"]=weightedFreqB["1"]/numTerms
    weightedFreqB["0"]=weightedFreqB["0"]/numTerms
    for c in existingNonexistingCombos:
        weightedFreqAB[c]=weightedFreqAB[c]/numTerms

    # Calculate final score for all combinations of expression levels
    score=0
    for c in existingNonexistingCombos:
        wfA=weightedFreqA[c[0]]
        wfB=weightedFreqB[c[1]]
        wfAB=weightedFreqAB[c]
        if wfAB > 0:
            score+= wfAB*(math.log10(wfAB/(wfA*wfB)))/math.log10(2) # Unclear why division by log10(2)
    return score

def calculateSclScoresInParallel(chunk,numChunks,genes,geneProfileDict,weightLookup,numTerms,termMap,calculateSemanticSimilarity,allPathsToRoot):
    createdLinks=0
    # Calculating start and stop index for chunk based on the total number of calculations to be made
    numGenes=len(genes)
    totalCalculations=(numGenes*(numGenes-1))/2
    calculationsPerChunk=round(totalCalculations/numChunks)
    startIndex=calculationsPerChunk*chunk
    stopIndex=(calculationsPerChunk*chunk)+calculationsPerChunk
    dataToAdd=[]
    pairIndex=0
    for index,geneA in enumerate(genes[:-1]):
        for geneB in genes[index+1:]:
            pairIndex+=1
            if pairIndex>startIndex and pairIndex<=stopIndex:
                if calculateSemanticSimilarity==False:
                    score = calculateSCLScore(geneProfileDict[geneA],geneProfileDict[geneB],weightLookup,numTerms)
                else:
                    termsA=[termMap[id] for id in geneProfileDict[geneA]]
                    termsB=[termMap[id] for id in geneProfileDict[geneB]]
                    score = getSemanticSimilarity(termsA,termsB,allPathsToRoot)
                intersectingTermIndexes = set(geneProfileDict[geneA]).intersection(set(geneProfileDict[geneB]))
                if len(intersectingTermIndexes)>0:
                    intersectingTerms=[]
                    for t in intersectingTermIndexes:
                        intersectingTerms.append(termMap[t])
                    metadata = ','.join(intersectingTerms)
                else:
                    metadata = None
                #print("score: "+str(score)+" metadata: "+metadata)
                newRecord=[geneA, geneB, score, metadata]
                dataToAdd.append(newRecord)
                createdLinks+=1

    linkDf = pd.DataFrame.from_records(dataToAdd)
    if len(linkDf)>0:
        linkDf.columns = ["proteinA","proteinB","score","metadata"]
        linkDf.proteinA, linkDf.proteinB = np.where(linkDf.proteinA < linkDf.proteinB, [linkDf.proteinB, linkDf.proteinA], [linkDf.proteinA, linkDf.proteinB])
    return linkDf

def getSCLLinksForSpecies(s,goTerms_cellularComponents,evidenceConfig,proteomeConfig,allParentTermsForAnnotations,calculateSemanticSimilarity,allPathsToRoot):
    print("Working on "+s)
    # Getting file
    spIndex=evidenceConfig['species'].index(s)
    url = evidenceConfig['goAnnotationUrl'].replace("SPECIES",evidenceConfig['goFilePerSpecies'][spIndex]) # Preparing url for gaf-file-download
    fileForSpecies="data/tmp/"+url.split("/")[-1]
    if os.path.exists(fileForSpecies)==False: # Downloads the gaf file if it does not already exist
        print("Downloading file: "+fileForSpecies)
        downloadFile(url,fileForSpecies)

    # Getting mappings from the DB for the specified proteome
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=s))[0]
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    available_protein_df.columns = ['proteinIdentifier','id']

    rowsToUse = []
    for chunk in  pd.read_csv(fileForSpecies, sep='\t',usecols=[1,3,4,6,12], skiprows=1, header=None, comment='!', chunksize=100000):
        rowsToUse.append(chunk[chunk[12] == "taxon:"+s]) # Only use rows with the required taxon
    termDf = pd.concat(rowsToUse, axis= 0)
    del rowsToUse

    termDf.columns = ['proteinIdentifier','association','goTerm','evidenceCode','species']
    termDf = termDf.drop('species',axis=1) # Drop species column as it is no longer needed
    termDf = termDf[termDf['evidenceCode'] != "IEA"] # Exclude electronic annotation evidence code (as this is not experimental at all)
    termDf = termDf[~termDf.association.str.startswith("NOT")] # Dont include associations starting with NOT
    termDf = termDf[termDf['goTerm'].isin(goTerms_cellularComponents)] # Only include terms that are a cellular component

    if calculateSemanticSimilarity==False:
        #mapping specific terms to more general terms
        termDf = pd.merge(termDf,allParentTermsForAnnotations, how='inner',on='goTerm')
        termDf = termDf.drop('goTerm',axis=1)
        termDf.rename(columns={'parentTerm':'goTerm'}, inplace=True)

    # mapping the protein identifiers to the database ID
    mapping = pd.DataFrame(list(termDf['proteinIdentifier']))
    mapping.columns = ['proteinIdentifier']
    preMap = set(mapping['proteinIdentifier'])
    termDf = pd.merge(termDf,available_protein_df, how='inner',on='proteinIdentifier')
    postMap = set(termDf['proteinIdentifier'])
    excluded = preMap - postMap
    print("Excluded "+str(len(excluded))+" genes, as they couldnt be mapped. Getting SCL scores for "+str(len(postMap))+" genes, for species: "+s)
    del available_protein_df
    del preMap
    del postMap
    del mapping
    del excluded

    # Cleaning up the data
    termDf = termDf.drop(['proteinIdentifier'],axis=1) # Dropping the proteinIdentifier column, no longer needed, using peorein id from now on
    termDf = termDf.drop(['association'],axis=1) # Dropping the association column, no longer needed
    termDf = termDf.drop(['evidenceCode'],axis=1) # Dropping the evidenceCode column, no longer needed
    termDf = termDf.drop_duplicates() # Drop duplicates
    termDf = termDf.groupby(['goTerm','id']).size().unstack(fill_value=0).reset_index().set_index('goTerm') # make protein ID as rows and goTerm as column, fills 1 or 0 for appearing/not appearing

    # Creating a weight lookup dict, so that calculations are not done more often than necessary
    numGenes=len(termDf.columns)
    numTerms=len(termDf.index)
    # print("Number of genes: "+str(numGenes))
    # print("Number of terms: "+str(numTerms))
    numOnes=0
    for termIndex,t in enumerate(list(termDf.index)):
        proteinIdsWithTerm = termDf.columns[termDf.loc[t] > 0].values
        numAppearancesOfTerm=len(proteinIdsWithTerm)
        numOnes+=numAppearancesOfTerm
    geneProfileIndexes={}
    termMap={}
    weightLookup={} # Getting all weights for all go-terms, either existing (1) or non-existing (0). Note that this is done for before calculating individual scores to save time
    for termIndex,t in enumerate(list(termDf.index) ):
        termMap[termIndex]=t
        proteinIdsWithTerm = termDf.columns[termDf.loc[t] > 0].values
        for p in proteinIdsWithTerm:
            if p not in geneProfileIndexes:
                geneProfileIndexes[p]=[]
            geneProfileIndexes[p].append(termIndex)
        numAppearancesOfTerm=len(proteinIdsWithTerm)
        weightLookup[str(termIndex)+"|1"]= 1-(numAppearancesOfTerm/numOnes) # weight of go term existing, 1-(number of appearances of the term/all genees(all possible appearances))
        weightLookup[str(termIndex)+"|0"]=1-((numOnes-numAppearancesOfTerm)/numOnes) # weight of go-term not existing, 1-(number of genes not having this go term/all genes (all possible appearances))
    del termDf
    genes = list(geneProfileIndexes.keys())
    random.seed(12345)
    genes = random.sample(genes,len(genes))
    
    # Score gene pairs in parallell
    cpus=os.cpu_count()
    linkDfList = Parallel(n_jobs=cpus)(delayed(calculateSclScoresInParallel)(chunk,cpus,genes,geneProfileIndexes,weightLookup,numTerms,termMap,calculateSemanticSimilarity,allPathsToRoot) for chunk in range(0, cpus))
    linkDf = pd.concat(linkDfList)
    if len(linkDf)>0:
        # print("Score range before normalization: "+str(linkDf['score'].min())+"-"+str(linkDf['score'].max()))    
        # MinMax scaling
        scores = linkDf['score'].values.reshape(-1, 1)
        scaler = MinMaxScaler(feature_range=(0, 1))
        scaled_scores = np.round(scaler.fit_transform(scores),decimals=4)
        linkDf['score'] = scaled_scores
        
        # Create a new evidence to attach new links to
        speciesInDB = Species.objects.get(tax_id = s)
        evidenceInDB = Evidence(type="SCL", species=speciesInDB, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidenceInDB.save()

        linkDf = linkDf.sample(frac=1).reset_index(drop=True)
        linkDf.insert(0,'evidence',evidenceInDB.id)

    return linkDf

def getSCLLinks(evidenceConfig, speciesToCollect, proteomeConfig):
    # Fetching GO annotations for all species to get and calculates a weighted mutual information score
    # for each gene pair in all species. As the information is divided in multiple fles, each taxID is mapped to
    # a gaf-file in dict taxIDToGOFile. For species that does not have their own file, the massive uniprot_all file is used
    if 'semanticSimilarity' in evidenceConfig:
        calculateSemanticSimilarity=evidenceConfig['semanticSimilarity']
    else:
        calculateSemanticSimilarity=False
    # Get graph from all oboterms, to be able to know which terms are cellular locations
    goTerms_cellularComponents=[]
    allParentTermsForAnnotations=[]
    graph = obonet.read_obo(evidenceConfig['goOntologyFile'])
    name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}
    allPathsToRoot={}
    for g in graph:
        if graph.nodes[g]['namespace']=="cellular_component":
            allPathsToRoot[g]=[]
            goTerms_cellularComponents.append(g)
            paths = networkx.all_simple_paths(graph,source=g,target=name_to_id['cellular_component'])
            addedTermsForG=[]
            for path in paths:
                pathWithValues=[]
                for p in path:
                    isa=[]
                    partof=[]
                    if 'is_a' in graph.nodes[p]:
                        isa = graph.nodes[p]['is_a']
                    if 'part_of' in graph.nodes[p]:
                        partof = graph.nodes[p]['part_of']
                    gObj={'id':p,'relationship': {'is_a':isa,'part_of':partof}}
                    pathWithValues.append(gObj)
                    if p not in addedTermsForG:
                        addedTermsForG.append(p)
                        newRecord=[str(g).strip(), str(p).strip()]
                        allParentTermsForAnnotations.append(newRecord)
                allPathsToRoot[g].append(pathWithValues)
    allParentTermsForAnnotations = pd.DataFrame.from_records(allParentTermsForAnnotations)
    allParentTermsForAnnotations.columns = ["goTerm", "parentTerm"]

    # CPUSforWithinSpeciesCalc=CPUS
    # CPUSforBetweenSpeciesCalc=1
    # if CPUS >= (2*len(speciesToCollect)):
    #     CPUSforWithinSpeciesCalc= math.floor(CPUS/len(speciesToCollect))
    #     CPUSforBetweenSpeciesCalc=len(speciesToCollect)
    for s in speciesToCollect:
        linkDf = getSCLLinksForSpecies(s,goTerms_cellularComponents,evidenceConfig,proteomeConfig,allParentTermsForAnnotations,calculateSemanticSimilarity,allPathsToRoot)
        if len(linkDf)==0: continue
        mem_csv = StringIO()
        linkDf.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        print("Writing %i SCL links to DB for %s" %(len(linkDf),s))
        # Writing the csv to the evidenceLink table
        with closing(mem_csv) as csv_io:
            SCL.objects.from_csv(csv_io)