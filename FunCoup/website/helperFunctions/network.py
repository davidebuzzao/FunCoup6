from data.models import *
from django.db.models import Q
from .colors import *
from .scoreCalculations import *
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from ..externalAlgorithms.topas import TOPAS


def ambiguousOrMissingMappings(geneQuery, genomeSelect, missingMapping, multiMapping, instanceConfigurations):
    # Gett mapping information for ambiguous or missing mappings using a more liberal search in the database

    # First, set up some parameters
    allQueryGenes=geneQuery.split(",")
    genesMissingMapping=missingMapping.split(",")
    genesWithMultipleMappings=multiMapping.split(",")
    normalMappings=[]
    missingMappings=[]
    multiMappings=[]
    proteomeInDB_id=Proteome.objects.filter(Q(species__tax_id=genomeSelect) & Q(version=instanceConfigurations['proteomeConfig']['genome']['version']))[0].id
    
    # Get all the identifiers that had no mapping issues, "normalMappings"
    for n in allQueryGenes:
        if n not in genesMissingMapping and n not in genesWithMultipleMappings:
            m = IdMapping.objects.filter(Q(mapped_id=n) & Q(protein__proteome=proteomeInDB_id)).values_list("protein__uniprot_id","mapped_id","protein__description","type").distinct()
            m=m[0]
            normalMappings.append({'input':m[1], 'output':m[0], 'type':m[3], 'description':m[2]})
            
    # Get all identifiers that initially had missing mappings, try looking up in DB using icontains
    for g in genesMissingMapping:
        min_length=3
        if g!="''" and len(g)>=min_length:
            mappingsPerGene={} # Using dict to get unique mapping rows for each protein
            mostSimilar="" # Using to get the identifier that is most similar to the query, to put first in the list
            caseInsensitiveMatches = IdMapping.objects.filter(Q(mapped_id__icontains=g) & Q(protein__proteome=proteomeInDB_id)).values_list("protein__uniprot_id","mapped_id","protein__description","type").distinct()
            if len(caseInsensitiveMatches)>0:
                for c in caseInsensitiveMatches:
                    if c[1].upper()==g.upper():
                        mostSimilar=c[0]
                    if c[0] not in mappingsPerGene:
                        mappingsPerGene[c[0]]={'input':g,'database':[c[1]], 'output':c[0], 'description':c[2]}
                    else:
                        mappingsPerGene[c[0]]['database'].append(c[1])
                listOfMappingsPerGene=list(mappingsPerGene.values())
                if mostSimilar!="": # Put the most similar match at the top of the list
                    listOfMappingsPerGene.remove(mappingsPerGene[mostSimilar])
                    listOfMappingsPerGene.insert(0,mappingsPerGene[mostSimilar])
                multiMappings.append(listOfMappingsPerGene)
            else:
                missingMappings.append(g)

    # Get all identifiers that has multiple uniprotIDs mapped to the identifier
    for g in genesWithMultipleMappings:
        if g!="''":
            mappingsPerGene={} # Using dict to get unique mapping rows for each protein
            mostSimilar=""  # Using to get the identifier that is most similar to the query, to put first in the list
            multi = IdMapping.objects.filter(Q(mapped_id=g) & Q(protein__proteome=proteomeInDB_id)).values_list("protein__uniprot_id","mapped_id","protein__description","type").distinct()
            for m in multi:
                if m[1].upper()==g.upper():
                    mostSimilar=m[0]
                if m[0] not in mappingsPerGene:
                    mappingsPerGene[m[0]]={'input':g,'database':[m[1]], 'output':m[0], 'description':m[2]}
                else:
                    mappingsPerGene[m[0]]['database'].append(m[1])
            if len(mappingsPerGene)>0:
                listOfMappingsPerGene=list(mappingsPerGene.values())
                if mostSimilar!="":  # Put the most similar match at the top of the list
                    listOfMappingsPerGene.remove(mappingsPerGene[mostSimilar])
                    listOfMappingsPerGene.insert(0,mappingsPerGene[mostSimilar])
                multiMappings.append(listOfMappingsPerGene)

    return normalMappings, multiMappings, missingMappings


def setUpTransferred(originalGenome, orthologSpecies, transferredNetwork, queryGenes, FCtoQueryMap, originalGeneList, genomeFileName, genomeSelect, taxIdToName, paramsForNetwork):
    # Setting up parameters and collecting data needed for a transferred network
    orthologyTransferMessage=""
    if originalGenome in orthologSpecies:
        transferredNetwork=True
        queryGenes=[] # overwriting query genes to be used in network construction (original is kept and nodes are labeled accordingly)
        # Reading orthologs from file
        file = open("website/static/website/orthologFiles/"+originalGenome+":"+orthologSpecies[originalGenome], "r")
        for line in file.readlines():
            if line.startswith("#"): # Header contains the taxid from where orthologs were transferred
                genomeSelect=line.strip().replace("#","")
            else:
                if line.split(",")[1].strip() in FCtoQueryMap:
                    FCtoQueryMap[line.split(",")[1].strip()].append(line.split(",")[0].strip())
                else:
                    FCtoQueryMap[line.split(",")[1].strip()]=[line.split(",")[0].strip()]
                if line.split(",")[0].strip() in originalGeneList:
                    queryGenes.append(line.split(",")[1].strip())
        if paramsForNetwork['individualEvidenceOnly']=="on" and paramsForNetwork['comparativeSpecies']!=["''"]: # Update individual sp evidence with the host network genome instead
            paramsForNetwork['constrainEv']="species"
            paramsForNetwork['constrainSpSelected']=[genomeSelect]
        genomeFileName="TRANSFERRED_"+orthologSpecies[originalGenome].split("_")[0]+"_"+orthologSpecies[originalGenome].split("_")[1]+"("+originalGenome+")_from_"
        orthologyTransferMessage="This network is transferred from "+taxIdToName[genomeSelect]['name']
    return transferredNetwork, genomeFileName, queryGenes, genomeSelect, FCtoQueryMap,orthologyTransferMessage, paramsForNetwork


def getOrthologInfoForQuery(genomeFileName, paramsForNetwork, allViewSpecies, genomeSelect, comp_species, taxIdToName, listOfGenes, proteinIdsToQuery, primaryColors):
    # Obtaining information on orthologs for comparative interactomics
    
    genomeFileName = taxIdToName[comp_species]['forFile']
    if paramsForNetwork['individualEvidenceOnly']=="on": # For individual evidence, the species restriction is updated to the current species
        paramsForNetwork['constrainSpSelected']=[comp_species]
    # Appending the comparative species to the species display
    allViewSpecies.append({'color':primaryColors[comp_species],'speciesName':taxIdToName[comp_species]['name'], 'short':taxIdToName[comp_species]['short']})
    # Collecting orthologs for the query proteins, using these for query in the comparative network
    orthologs_dict={'query':{}}
    
    if int(genomeSelect) > int(comp_species):
        # Orthologs for ALL network genes are collected in order to draw ortholines in the network
        orthologs=ProteinOrtholog.objects.filter(Q(speciesA=genomeSelect) & Q(speciesB=comp_species) & Q(proteinA__in=listOfGenes))
        for ortho in orthologs.filter(Q(proteinA__in=proteinIdsToQuery)):
            orthologs_dict['query'][ortho.proteinB.id] = orthologs_dict['query'].get(ortho.proteinB.id,{})
            orthologs_dict['query'][ortho.proteinB.id] = {'gene': ortho.proteinB.uniprot_id, 'uniprot': ortho.proteinB.uniprot_id}
    else:
        orthologs = ProteinOrtholog.objects.filter(Q(speciesB=genomeSelect) & Q(speciesA=comp_species) & Q(proteinB__in=listOfGenes))
        for ortho in orthologs.filter(Q(proteinB__in=proteinIdsToQuery)):
            orthologs_dict['query'][ortho.proteinA.id] = orthologs_dict['query'].get(ortho.proteinA.id,{})
            orthologs_dict['query'][ortho.proteinA.id]={'gene': ortho.proteinA.uniprot_id, 'uniprot': ortho.proteinA.uniprot_id}

    return genomeFileName, paramsForNetwork, allViewSpecies, orthologs_dict, orthologs 

def getOrthologInfoForQueryORTHOLOGSonly(genomeFileName, paramsForNetwork, allViewSpecies, genomeSelect, comp_species, taxIdToName, listOfGenes, proteinIdsToQuery, primaryColors):
    # Obtaining information on orthologs for comparative interactomics
    genomeFileName = taxIdToName[comp_species]['forFile']
    if paramsForNetwork['individualEvidenceOnly']=="on": # For individual evidence, the species restriction is updated to the current species
        paramsForNetwork['constrainSpSelected']=[comp_species]
    # Appending the comparative species to the species display
    allViewSpecies.append({'color':primaryColors[comp_species],'speciesName':taxIdToName[comp_species]['name'], 'short':taxIdToName[comp_species]['short']})
    
    # Collecting orthologs for the query proteins, using these for query in the comparative network
    orthologs_dict={'query':{},'other':{}}
    if int(genomeSelect) > int(comp_species):
        # Orthologs for ALL network genes are collected in order to draw ortholines in the network
        orthologs=ProteinOrtholog.objects.filter(Q(speciesA=genomeSelect) & Q(speciesB=comp_species) & Q(proteinA__in=listOfGenes))
        for ortho in orthologs:
            if ortho.proteinA.id in proteinIdsToQuery:
                orthologs_dict['query'][ortho.proteinB.id] = orthologs_dict['query'].get(ortho.proteinB.id,{})
                orthologs_dict['query'][ortho.proteinB.id] = {'gene': ortho.proteinB.uniprot_id, 'uniprot': ortho.proteinB.uniprot_id}
            else:
                orthologs_dict['other'][ortho.proteinB.id] = orthologs_dict['other'].get(ortho.proteinB.id,{})
                orthologs_dict['other'][ortho.proteinB.id] = {'gene': ortho.proteinB.uniprot_id, 'uniprot': ortho.proteinB.uniprot_id}
    else:
        orthologs = ProteinOrtholog.objects.filter(Q(speciesB=genomeSelect) & Q(speciesA=comp_species) & Q(proteinB__in=listOfGenes))
        for ortho in orthologs:
            if ortho.proteinB.id in proteinIdsToQuery:
                orthologs_dict['query'][ortho.proteinA.id] = orthologs_dict['query'].get(ortho.proteinA.id,{})
                orthologs_dict['query'][ortho.proteinA.id]={'gene': ortho.proteinA.uniprot_id, 'uniprot': ortho.proteinA.uniprot_id}
            else:
                orthologs_dict['other'][ortho.proteinA.id] = orthologs_dict['other'].get(ortho.proteinA.id,{})
                orthologs_dict['other'][ortho.proteinA.id]={'gene': ortho.proteinA.uniprot_id, 'uniprot': ortho.proteinA.uniprot_id}
    return genomeFileName, paramsForNetwork, allViewSpecies, orthologs_dict, orthologs 

def getPpvFunctions(paramsForNetwork, genomeSelect, instanceConfigurations):
    # Getting ppv functions for species for all goldstandards in case of constrain by evidence or species
    # Used to recvalculate the PPV

    ppvFunctions=[]
    if any(sub in paramsForNetwork['constrainEv'] for sub in ["species", "evidence"]): # Preparing LLR information needed to recompute Confidence scores
        for g in instanceConfigurations['trainingConfig']['Network']['GoldStandard_order']:
            gs = GoldStandard.objects.filter(Q(type=g) & Q(species__tax_id=genomeSelect))
            if gs.exists():
                gsid = gs[0].id
            else:
                gsid=0
            logPPV = LogisticPPV.objects.filter(Q(goldStandard=gsid))
            if logPPV.exists():
                ppvFunctions.append(logPPV[0].function)
            else:
                ppvFunctions.append('NONE')
    return ppvFunctions


def runTopas(genomeFileName, paramsForNetwork, proteinIdsToQuery):
    # This is run in case algorithm TOPAS is selected. For this, the depth iterations are disregarded
    network_f = "website/static/website/networks/FunCoup6.0/FC6.0_"+genomeFileName+"_compact.gz"
    network_df = pd.read_csv(network_f, compression='gzip', sep='\t', usecols=[0,1,5])
    network_df.columns = ['gene1','gene2','ppv']
    network_df['ppv'] = network_df['ppv'].astype(float)
    network_df = network_df[network_df['ppv'] >= paramsForNetwork['confidenceThreshold']]
    network_df = network_df.drop('ppv', axis=1).reset_index(drop=True)
    candidate_genes_df = TOPAS(network_df, [(proteinIdsToQuery[i]['uniprot']) for i in proteinIdsToQuery ], expansion_steps=paramsForNetwork['expansion_steps'], cores=1)
    listOfGenes=[]
    if candidate_genes_df is not None and len(candidate_genes_df.index)>0:
        candidate_genes=list(set(candidate_genes_df['gene1'].tolist()+candidate_genes_df['gene2'].tolist()))
        for g in (candidate_genes):
            listOfGenes.extend(list(IdMapping.objects.filter(Q(mapped_id=g)).values_list("protein_id",flat=True).distinct()))
    listOfGenes=list(set(listOfGenes))
    return listOfGenes


def updateForConstrainedEvidence(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames, linksDf):
    # Updates the linksDf to constrained evidences, recalculates fbs for the evidence types selected
    indexesToConsider=[]
    for ev in paramsForNetwork['constrainEvTypesSelected']:
        indexesToConsider.append(instanceConfigurations['trainingConfig']['Network']['Evidence_order'].index(ev))
    linksDf=linksDf.apply(lambda row: getNewFbs(row,indexesToConsider,ppvFunctions), axis=1)
    linksDf[gsColNames] = linksDf['ppv'].str.split(',',expand=True).astype(float) # Split ppv column by goldstandard to make easier to query 
    return linksDf


def updateForConstrainedSpecies(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames,linksDf):
    # Updates the linksDf to constrained species, recalculates fbs for the species selected

    indexesToConsider=[]
    for s in paramsForNetwork['constrainSpSelected']:
        indexesToConsider.append(instanceConfigurations['trainingConfig']['Network']['Species_order'].index(s))
    linksDf=linksDf.apply(lambda row: getNewFbsSpecies(row,indexesToConsider,ppvFunctions), axis=1)
    linksDf[gsColNames] = linksDf['ppv'].str.split(',',expand=True).astype(float) # Split ppv column by goldstandard to make easier to query 
    return linksDf


def updateForGroup(paramsForNetwork, linksDf, listOfGenes):
    # Updates linkDf to filter for the "group" expansion algorithm
    # This queries for all query proteins, and filter the highest scoring neighbors for the group
    # If prioritize common neighbors is activated, the proteins in the subnetwork with the most neighbors are kept
    if paramsForNetwork['prioritizeCommonNeighbors']=="on" and len(listOfGenes)>1:
        proteinsA=linksDf[['proteinA']].rename({'proteinA':'proteinB'}, axis='columns')
        proteinsB=linksDf[['proteinB']]
        tmpProteins=pd.concat([proteinsA,proteinsB]).groupby('proteinB').size().sort_values(ascending=False).reset_index()
    else:
        tmpProteinsA = linksDf[['proteinA',  paramsForNetwork['goldstandardColumn'] ]].copy().rename({'proteinA':'proteinB'}, axis='columns')
        tmpProteinsB = linksDf[['proteinB', paramsForNetwork['goldstandardColumn'] ]].copy()
        tmpProteins=pd.concat([tmpProteinsA,tmpProteinsB]).groupby('proteinB')[paramsForNetwork['goldstandardColumn']].max().reset_index().nlargest( paramsForNetwork['maxNodesPerStep']+len(listOfGenes), paramsForNetwork['goldstandardColumn'])
    tmpProteinsList=tmpProteins['proteinB'].tolist()
    initialGeneListLength=len(listOfGenes)
    for p in tmpProteinsList:
        if p not in listOfGenes and len(listOfGenes)<paramsForNetwork['maxNodesPerStep']+initialGeneListLength:
            listOfGenes.append(p)
    return listOfGenes


def updateForLocal(linksDf, paramsForNetwork, listOfGenes):
    # Updates linkDf to filter for the "local" expansion algorithm
    # here each query id treated individually so that each of its highest scoring neighbors are kept
    protsToAdd=[]
    for q in listOfGenes:
        queryGeneDf = linksDf[(linksDf['proteinA'] == q) | (linksDf['proteinB'] == q)].nlargest(paramsForNetwork['maxNodesPerStep'], paramsForNetwork['goldstandardColumn'])
        protsToAdd.extend(list(set(queryGeneDf['proteinA'].tolist()+queryGeneDf['proteinB'].tolist())))
    listOfGenes=list(set(protsToAdd))
    return listOfGenes


def updateForMaxlink(i, linksDf, paramsForNetwork,genomeSelect, proteinIdsToQuery, maxLinkAdditionalInfo, listOfGenes ):
    # update the dataframe of collected nodes according to maxlink
    # and save the additional maxlink information

    if i!=1: # First round (as maxlink is set to do two depth iterations)
        queryGeneDf = linksDf.nlargest(paramsForNetwork['maxNodesPerStep'], paramsForNetwork['goldstandardColumn'])
        listOfGenes=list(set(queryGeneDf['proteinA'].tolist()+queryGeneDf['proteinB'].tolist())) # Picking up top candidate genes
    else: # second round, at this stage, all candidate neighbours should be here as well
        allSubnetworkGenesList = list(linksDf['proteinA'].tolist()+linksDf['proteinB'].tolist())
        network_size = Stats.objects.filter(Q(tax_id=genomeSelect)).values_list('num_nodes', flat=True)[0]
        candidateGenes=listOfGenes.copy()
        listOfGenes=list(proteinIdsToQuery.keys())
        pvals={}
        degrees=[]
        indegrees=[]
        for g in candidateGenes: # for each candidate
            if g not in proteinIdsToQuery:
                degree=allSubnetworkGenesList.count(g)
                degrees.append(degree)
                queryLinks = linksDf[((linksDf['proteinA'] == g) & (linksDf['proteinB'].isin(proteinIdsToQuery)) ) | ((linksDf['proteinB'] == g) & (linksDf['proteinA'].isin(proteinIdsToQuery)) )]
                in_degree=queryLinks.shape[0]
                indegrees.append(in_degree)
                hypergeometric_p = hypergeom.sf(in_degree - 1, network_size, len(proteinIdsToQuery), degree)
                pvals[g]=hypergeometric_p               
        rejected_bonferroni, p_adjusted_bonferroni, _, alpha_corrected_bonferroni = multipletests(list(pvals.values()), alpha=paramsForNetwork['hypergeomCutoff'], method='bonferroni', is_sorted=False, returnsorted=False)
        rejected_bh, p_adjusted_bh, _, alpha_corrected_bh = multipletests(list(pvals.values()), alpha=paramsForNetwork['hypergeomCutoff'], method='fdr_bh', is_sorted=False, returnsorted=False)
        for j,cand in enumerate(pvals.keys()):
            if rejected_bh[j]==True:
                maxLinkAdditionalInfo[cand]={'p':"{:.2e}".format(hypergeometric_p), 'indegree':indegrees[j], 'totalLinks': degrees[j],'fdr':"{:.2e}".format(p_adjusted_bh[j]),'fwer':"{:.2e}".format(p_adjusted_bonferroni[j])}
                maxLinkAdditionalInfo[cand]['fdr']="{:.2e}".format(p_adjusted_bh[j])
                maxLinkAdditionalInfo[cand]['fwer']="{:.2e}".format(p_adjusted_bonferroni[j])
                listOfGenes.append(cand)
    return listOfGenes, maxLinkAdditionalInfo


def getMappingsPathwaysAndTissuesForSubnetworkGenes(listOfGenes, transferredNetwork, proteinIdsToQuery, mappingTypes, allNetworkPathways, allNetworkTissues, isAPI=False):
    # Mappings are already in place for the query genes, but mappings for all subnetwork genes has to be retrieved in order to be able to show
    # different label types.
    mappingDict={}
    pathwayIdDict={}
    tissueDict={}   
    
    # Dont collect mappings for orthology transferred networks
    if transferredNetwork: 
        pass
    else:

        # Query to get mappings, excluding alternatives
        allMappings = IdMapping.objects.filter(
            Q(protein_id__in=listOfGenes) &
            ~Q(type__startswith="Alt_") &
            ~Q(type__startswith="alt")
        )

        # Step 1: Collect unique mapping types
        mappingTypes = set()
        for m in allMappings:
            mappingTypes.add(m.type)

        # Step 2: Initialize data structures
        for p in listOfGenes:  # Initiate with empty data structures for all proteins
            pathwayIdDict[p] = []
            if p in proteinIdsToQuery.keys():
                mappingDict[p] = {'queryId': proteinIdsToQuery[p]['gene']}
            else:
                mappingDict[p] = {'queryId': 'None'}
            tissueDict[p] = []

        # Step 3: Process mappings
        for m in allMappings:
            mapping_type = m.type
            prot = m.protein_id

            # Replace incompatible characters for visualization?
            # mappedID = m.mapped_id.replace(".", "_").replace(":", "_").replace("\\", "_").replace("-", "_").replace("(", "_").replace(")", "_")
            mappedID = m.mapped_id

            # Store mappings for this protein
            if prot not in mappingDict:
                mappingDict[prot] = {}

            # Store the UniProtID immediately if encountered
            if mapping_type == "UniProtID":
                mappingDict[prot]['UniProtID'] = mappedID

            # Assign the mapped ID to the corresponding type
            if mapping_type not in mappingDict[prot]:
                mappingDict[prot][mapping_type] = mappedID

            if mapping_type=="Gene_Symbol" and (mappingDict[prot]['queryId']=='None' or mappingDict[prot]['queryId'].startswith("#")):
                mappingDict[prot]['queryId']=mappedID

        # Step 4: Ensure all required mapping types are present with UniProtID as fallback
        for prot, mappings in mappingDict.items():
            # Ensure all mapping types are in the dictionary
            for mapping_type in mappingTypes:
                if mapping_type != 'UniProtID' and mapping_type not in mappings:
                    # Fallback to UniProtID if available
                    if 'UniProtID' in mappings:
                        mappings[mapping_type] = mappings['UniProtID']

            # Ensure the queryId falls back to UniProtID if needed
            if mappings.get('queryId') == 'None' or mappings.get('queryId', '').startswith("#"):
                if 'UniProtID' in mappings:
                    mappings['queryId'] = mappings['UniProtID']

        pathwayAnnotations = Pathway.objects.filter(Q(protein_id__in=listOfGenes))
        for p in pathwayAnnotations:
            pathwayIdDict[p.protein_id].append(p.pathway_id)
            # allNetworkPathways.add(p.pathway_name+" (pathwayID:"+p.pathway_id+")") # for list of pathways to filter for
            allNetworkPathways[p.pathway_id]=p.pathway_name+" (pathwayID:"+p.pathway_id+")"
        tissueAnnotations = Tissue.objects.filter(Q(protein_id__in=listOfGenes))
        for p in tissueAnnotations:
            tissue=p.tissue.replace("/"," ")
            if isAPI:
                tissueDict[p.protein_id].append(tissue)
            else:
                tissueDict[p.protein_id].append(p.tissue_id) # each node keeps info on its own tissues
            allNetworkTissues[p.tissue_id]=tissue+" (tissueID:"+str(p.tissue_id)+")"  # for list of tissues to filter for
    
    return mappingDict, pathwayIdDict, tissueDict, mappingTypes, allNetworkPathways, allNetworkTissues


def addNodesToNetwork(G, listOfGenes, proteinIdsToQuery, nodePartners, transferredNetwork, FCtoQueryMap, maxLinkAdditionalInfo, tissueDict, mappingDict, pathwayIdDict, primaryColor, taxIdToName, originalGenome, compNetwork):
    # Adding the collected list of genes as nodes into a network, tohether with all information it needs for the frontend
    
    nodeProteins= Protein.objects.filter(Q(id__in=listOfGenes))
    excludedNodes=[]
    for node in nodeProteins:
        border = getQueryBorder() if node.id in proteinIdsToQuery.keys() else "None"
        borderWidth = "5px" if node.id in proteinIdsToQuery.keys() else "0px"
        queryNode = "True" if node.id in proteinIdsToQuery.keys() and compNetwork==False else "False"
        nodePartners[str(node.id)]=[]
        if transferredNetwork==True: # For transferred networks, node labels need to be translated back to the original query
            if node.uniprot_id in FCtoQueryMap:
                for orthologMapping in FCtoQueryMap[node.uniprot_id]:
                    maxlinkp='-'
                    maxlinkScore='-'
                    maxlinkDegree='-'
                    maxlinkFdr='0'
                    maxlinkFwer='-'
                    if node.id in maxLinkAdditionalInfo: # Add correct mxlink info to node if it exists
                        maxlinkp=maxLinkAdditionalInfo[node.id]['p']
                        maxlinkScore=maxLinkAdditionalInfo[node.id]['indegree']
                        maxlinkDegree=maxLinkAdditionalInfo[node.id]['totalLinks']
                        maxlinkFdr=maxLinkAdditionalInfo[node.id]['fdr']
                        maxlinkFwer=maxLinkAdditionalInfo[node.id]['fwer']
                    ## All other nodes
                    G.add_node(str(node.id), uniprotID=orthologMapping, queryNode=queryNode, tissues=[], pathwayIds=[], description="", shortDescription="", nodeBorder=border, nodeBorderWidth=borderWidth, mappings={'UniProtID':orthologMapping}, species={'color':primaryColor[originalGenome],'speciesName':taxIdToName[originalGenome]['name'],'taxId':originalGenome}, maxlinkp=maxlinkp, maxlinkScore=maxlinkScore, maxlinkDegree=maxlinkDegree,maxlinkFdr=maxlinkFdr,maxlinkFwer=maxlinkFwer)
            else: # Nodes that dont get an ortholog mapping are excluded, need to keep this info to know which links not to add
                excludedNodes.append(node.id)
        else: # If not transferred add node with its proper name and with/without maxlink info
            maxlinkp='-'
            maxlinkScore='-'
            maxlinkDegree='-'
            maxlinkFdr='0'
            maxlinkFwer='-'
            if node.id in maxLinkAdditionalInfo: # Add correct mxlink info to node if it exists
                maxlinkp=maxLinkAdditionalInfo[node.id]['p']
                maxlinkScore=maxLinkAdditionalInfo[node.id]['indegree']
                maxlinkDegree=maxLinkAdditionalInfo[node.id]['totalLinks']
                maxlinkFdr=maxLinkAdditionalInfo[node.id]['fdr']
                maxlinkFwer=maxLinkAdditionalInfo[node.id]['fwer']
            
            if len(primaryColor)>1:
                ## Host-pathogen nodes
                node_tax_id=str(node.proteome.species.tax_id)
                G.add_node(str(node.id), uniprotID=mappingDict[node.id]["queryId"].replace("#",""), queryNode=queryNode,  tissues=sorted(tissueDict[node.id]), pathwayIds=sorted(pathwayIdDict[node.id]), description=node.description, shortDescription=node.description[:20]+" ...", nodeBorder=border, nodeBorderWidth=borderWidth, mappings=mappingDict[node.id], species={'color':primaryColor[node_tax_id],'speciesName':taxIdToName[node_tax_id]['name'],'taxId':node_tax_id}, maxlinkp=maxlinkp, maxlinkScore=maxlinkScore, maxlinkDegree=maxlinkDegree,maxlinkFdr=maxlinkFdr,maxlinkFwer=maxlinkFwer)
            else:
                ## All other nodes
                G.add_node(str(node.id), uniprotID=mappingDict[node.id]["queryId"].replace("#",""), queryNode=queryNode,  tissues=sorted(tissueDict[node.id]), pathwayIds=sorted(pathwayIdDict[node.id]), description=node.description, shortDescription=node.description[:20]+" ...", nodeBorder=border, nodeBorderWidth=borderWidth, mappings=mappingDict[node.id], species={'color':primaryColor[originalGenome],'speciesName':taxIdToName[originalGenome]['name'],'taxId':originalGenome}, maxlinkp=maxlinkp, maxlinkScore=maxlinkScore, maxlinkDegree=maxlinkDegree,maxlinkFdr=maxlinkFdr,maxlinkFwer=maxlinkFwer)

    return G, nodePartners, excludedNodes


def addLinksToNetwork(G, linksDf, excludedNodes, instanceConfigurations, nodePartners, paramsForNetwork, taxIdToName, hasDirections, isAPI=False):
    # Iterating over all collected links, collects necessary link info and appends the link to the graph object
    for index, row in linksDf.iterrows():
        if row['proteinA'] not in excludedNodes and row['proteinB'] not in excludedNodes: # Make sure none of the nodes have been excluded due to missing mappings
            scoresPerGoldStandard=[] # Keeping list of score objects for each gold standard
            globalDirection=0 # Initiating a direction for the entire link
            ppv_list=row['ppv'].split(",")
            isGoldstandard=row['isGoldStandard']
            isPPVgoldstandard=row['isPPVgoldstandard']
            for i,gs in enumerate(instanceConfigurations['trainingConfig']['Network']['GoldStandard_order']): # Iterate over all gorld standards in the correct order
                ppv=ppv_list[i]
                isGoldstandardThis=int(isGoldstandard[i])
                isPPVgoldstandardThis=int(isPPVgoldstandard[i])
                # Only including info where conf score is >0 or node is a goldstandard node AND there is no restriction on GS or this is the restricted goldstandard in question
                if ((ppv!="" and float(ppv)>0.0) or isGoldstandardThis>0) and (paramsForNetwork['restrictGoldStandard']=="all" or gs==paramsForNetwork['restrictGoldStandard']):
                    fbs=row['fbs_goldstandard'].split(",")[i]
                    speciesLlrs=[]
                    if "evidence" not in paramsForNetwork['constrainEv']: # Not showing species breakdown if restriction on evidence (as the numbers then get misleading)
                        for j,instanceSp in enumerate(instanceConfigurations['trainingConfig']['Network']['Species_order']): # Iterating over the Correct species order for the instance
                            speciesLlrs.append({'type':taxIdToName[instanceSp]['short'],'llr':row['llr_evidence'].split(";")[i].split("|")[1].split(",")[j],'color':getLLRColor(float(row['llr_evidence'].split(";")[i].split("|")[1].split(",")[j]))})
                    evidenceLlrs=[]
                    grgLLR=0
                    if "species" not in paramsForNetwork['constrainEv']:  # Not showing evidence breakdown if restriction on species (as the numbers then get misleading)
                        for j,instanceEv in enumerate(instanceConfigurations['trainingConfig']['Network']['Evidence_order']): # Iterating over the Correct evidence order 
                            if instanceEv=="GRG": # Saving to check for link direction
                                grgLLR=float(row['llr_evidence'].split(";")[i].split("|")[0].split(",")[j])
                            evidenceLlrs.append({'type':instanceEv,'llr':row['llr_evidence'].split(";")[i].split("|")[0].split(",")[j],'color':getLLRColor(float(row['llr_evidence'].split(";")[i].split("|")[0].split(",")[j]))}) # For each evidence, add info on llr, color, and type (for hover)
                    
                    # Setting known cloupling info
                    knownCoupling="None"
                    if isGoldstandardThis>0:
                        knownCoupling="This"
                    elif isGoldstandard!="000000":
                        knownCoupling="Other"

                    # Setting direction for link in Gold standard network
                    globalDirection,localDirection = getDirection(int(row['GRG_direction'][i]), int(row['isGoldStandard'][i]), grgLLR, float(paramsForNetwork['grgLLRThreshold']), float(ppv),float(paramsForNetwork['confidenceThreshold']), globalDirection)

                    scoresPerGoldStandard.append({'type': gs,'ppv':float(ppv),'fbs':fbs,'evidenceLlrs':evidenceLlrs,'evidenceSpecies':speciesLlrs,'direction':localDirection,'known':knownCoupling,'isPPVgoldstandard':isPPVgoldstandardThis})
            
                    # Test JUN - CREB5
                    # if row['proteinA']==342209 and row['proteinB']==337024:
                    #     print('GlobalDir: ' + str(globalDirection))
                    #     print('LocalDir: ' + str(localDirection) + '\n')
            if globalDirection>0:
                hasDirections = True # Setting the global hasDirection variable in case of any dir in the network (to determine if the show directed links checkbox should be visible)
            
            scoresPerGoldStandard.sort(key=lambda x: (float(x['ppv']), float(x['fbs'])), reverse=True) # Sorting all score objects by ppv, highest first
            # Showing on API only link from highest scoring goldstandard network
            if isAPI: 
                G.add_edge(str(row['proteinA']),str(row['proteinB']),scoresPerGoldStandard=scoresPerGoldStandard[0],value="",direction=globalDirection)
            else:
                G.add_edge(str(row['proteinA']),str(row['proteinB']),scoresPerGoldStandard=scoresPerGoldStandard,value="",direction=globalDirection)
            
            nodePartners[str(row['proteinA'])].append(str(row['proteinB'])) # Adding all node partners for the nodes
            nodePartners[str(row['proteinB'])].append(str(row['proteinA']))

    return G, nodePartners, hasDirections

def getNetworkForGenes(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, proteinIdsToQuery, instanceConfigurations, genomeSelect,originalGenome, FCtoQueryMap, primaryColor, taxIdToName, genomeFileName, transferredNetwork, paramsForNetwork, compNetwork=False, isAPI=False):
    # Obtaining a species network for the given query parameters (can be either for the source or comparative species)

    # Getting names for the gold standard columns in the db
    gsColNames=[gs + "_ppv" for gs in instanceConfigurations['trainingConfig']['Network']['GoldStandard_order']] 
    listOfGenes=list(proteinIdsToQuery.keys())
    
    # Collecting all PPV functions, to use in case evidence or species are restricted to recalculate PPVs
    ppvFunctions = getPpvFunctions(paramsForNetwork, genomeSelect, instanceConfigurations)
    
    maxLinkAdditionalInfo={} # Keeping to store the additional maxlink information, in case of maxlink query
    if paramsForNetwork['expansionAlgorithm']=="topas":
        listOfGenes = runTopas(genomeFileName, paramsForNetwork, proteinIdsToQuery)
    else:
        # Iterating over depths
        for i in range(0,paramsForNetwork['depth']): 
            # First, collecting new nodes for each search depth
            linksDf = pd.DataFrame.from_records( Network.objects.filter((Q(proteinA_id__in=listOfGenes) | Q(proteinB_id__in=listOfGenes))).values('proteinA','proteinB','ppv','max_ppv','isGoldStandard','llr_evidence'))
            
            if len(linksDf)>0:
                linksDf[gsColNames] = linksDf['ppv'].str.split(',',expand=True).astype(float) # Split ppv column by goldstandard to make easier to query 
                if "evidence" in paramsForNetwork['constrainEv']: 
                    linksDf = updateForConstrainedEvidence(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames, linksDf)
                if "species" in paramsForNetwork['constrainEv']:
                    linksDf = updateForConstrainedSpecies(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames, linksDf)
                linksDf = linksDf[(linksDf[paramsForNetwork['goldstandardColumn']] >= paramsForNetwork['confidenceThreshold'])]
                
                if paramsForNetwork['restrictPathway']!="''" and paramsForNetwork['restrictPathway']!="": # Restrict search on given pathway
                    pid=paramsForNetwork['restrictPathway']
                    pathwayGenes=list(Pathway.objects.filter(Q(pathway_id=pid) & Q(species=genomeSelect)).values_list("protein", flat=True))
                    linksDf = linksDf[(linksDf['proteinA'].isin(pathwayGenes)) & (linksDf['proteinB'].isin(pathwayGenes))]
                elif paramsForNetwork['restrictTissue']!="": # Restrict search on given tissues
                    tissueGenes=list(Tissue.objects.filter(Q(tissue=paramsForNetwork['restrictTissue']) & Q(species=genomeSelect)).values_list("protein", flat=True))
                    linksDf = linksDf[(linksDf['proteinA'].isin(tissueGenes)) & (linksDf['proteinB'].isin(tissueGenes))]
                
                # Sort the DataFrame by the 'max_ppv' column in descending order
                linksDf = linksDf.sort_values(by='max_ppv', ascending=False).reset_index(drop=True)

                if paramsForNetwork['expansionAlgorithm']=="group":
                    listOfGenes = updateForGroup(paramsForNetwork, linksDf,listOfGenes)
                elif paramsForNetwork['expansionAlgorithm']=="local":
                    listOfGenes = updateForLocal(linksDf, paramsForNetwork,listOfGenes)
                elif paramsForNetwork['expansionAlgorithm']=="maxlink":
                    listOfGenes, maxlinkAdditionalInfo = updateForMaxlink(i, linksDf, paramsForNetwork,genomeSelect, proteinIdsToQuery, maxLinkAdditionalInfo, listOfGenes)

    # Done getting all new nodes, now getting all intermediate links between the added query nodes    
    linksDf = pd.DataFrame.from_records(Network.objects.filter((Q(proteinA_id__in=listOfGenes) & Q(proteinB_id__in=listOfGenes))).values('proteinA','proteinB','max_ppv','ppv','fbs_goldstandard','llr_evidence','GRG_direction','isGoldStandard','isPPVgoldstandard'))
    if len(linksDf)>0:
        linksDf[gsColNames] = linksDf['ppv'].str.split(',',expand=True).astype(float) # Split ppv column by goldstandard to make easier to query
        if "evidence" in paramsForNetwork['constrainEv']: 
            linksDf = updateForConstrainedEvidence(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames, linksDf)
        if "species" in paramsForNetwork['constrainEv']: 
            linksDf = updateForConstrainedSpecies(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames, linksDf)
        linksDf = linksDf[(linksDf[paramsForNetwork['goldstandardColumn']] >= paramsForNetwork['confidenceThreshold'])]
        
        # Sort the DataFrame by the 'max_ppv' column in descending order
        linksDf = linksDf.sort_values(by='max_ppv', ascending=False).reset_index(drop=True)
        
        # Getting mappings for all subnetwork proteins
        mappingDict, pathwayIdDict, tissueDict, mappingTypes, allNetworkPathways, allNetworkTissues = getMappingsPathwaysAndTissuesForSubnetworkGenes(listOfGenes, transferredNetwork, proteinIdsToQuery, mappingTypes, allNetworkPathways, allNetworkTissues, isAPI)
        
        # Adding nodes
        G, nodePartners, excludedNodes = addNodesToNetwork(G, listOfGenes, proteinIdsToQuery, nodePartners, transferredNetwork, FCtoQueryMap, maxLinkAdditionalInfo, tissueDict, mappingDict, pathwayIdDict, primaryColor, taxIdToName, originalGenome, compNetwork)
       
        # Adding links
        G, nodePartners, hasDirections = addLinksToNetwork(G, linksDf, excludedNodes, instanceConfigurations, nodePartners, paramsForNetwork, taxIdToName, hasDirections, isAPI)
       
    return G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes


def getNetworkForGenesORTHOLOGSonly(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, orthologQueryAndOther, instanceConfigurations, genomeSelect,originalGenome, FCtoQueryMap, primaryColor, taxIdToName, genomeFileName, transferredNetwork, paramsForNetwork, compNetwork=False, isAPI=False):
    # Obtaining a species network for the given query parameters (can be either for the source or comparative species)

    # Getting names for the gold standard columns in the db
    gsColNames=[gs + "_ppv" for gs in instanceConfigurations['trainingConfig']['Network']['GoldStandard_order']] 
    listOfGenes=list(orthologQueryAndOther['query'].keys())+list(orthologQueryAndOther['other'].keys())

    # Collecting all PPV functions, to use in case evidence or species are restricted to recalculate PPVs
    ppvFunctions = getPpvFunctions(paramsForNetwork, genomeSelect, instanceConfigurations)
    
    # Done getting all new nodes, now getting all intermediate links between the added query nodes    
    linksDf = pd.DataFrame.from_records(Network.objects.filter((Q(proteinA_id__in=listOfGenes) & Q(proteinB_id__in=listOfGenes))).values('proteinA','proteinB','max_ppv','ppv','fbs_goldstandard','llr_evidence','GRG_direction','isGoldStandard','isPPVgoldstandard'))
    if len(linksDf)>0:
        # Sort the DataFrame by the 'max_ppv' column in descending order
        #linksDf = linksDf.sort_values(by='max_ppv', ascending=False).reset_index(drop=True)
        linksDf[gsColNames] = linksDf['ppv'].str.split(',',expand=True).astype(float) # Split ppv column by goldstandard to make easier to query
        if "evidence" in paramsForNetwork['constrainEv']: 
            linksDf = updateForConstrainedEvidence(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames, linksDf)
        if "species" in paramsForNetwork['constrainEv']: 
            linksDf = updateForConstrainedSpecies(paramsForNetwork, instanceConfigurations, ppvFunctions, gsColNames, linksDf)
        # If no constraints or known couplings, always filter for confidenceThreshold
        linksDf = linksDf[(linksDf[paramsForNetwork['goldstandardColumn']] >= paramsForNetwork['confidenceThreshold'])]
        # Sort the DataFrame by the 'max_ppv' column in descending order
        linksDf = linksDf.sort_values(by='max_ppv', ascending=False).reset_index(drop=True)

        proteinIdsToQuery = orthologQueryAndOther['query']

        # Getting mappings for all subnetwork proteins
        mappingDict, pathwayIdDict, tissueDict, mappingTypes, allNetworkPathways, allNetworkTissues = getMappingsPathwaysAndTissuesForSubnetworkGenes(listOfGenes, transferredNetwork, proteinIdsToQuery, mappingTypes, allNetworkPathways, allNetworkTissues, isAPI)
        
        # Adding nodes
        maxLinkAdditionalInfo={} # Not used, just initialized
        G, nodePartners, excludedNodes = addNodesToNetwork(G, listOfGenes, proteinIdsToQuery, nodePartners, transferredNetwork, FCtoQueryMap, maxLinkAdditionalInfo, tissueDict, mappingDict, pathwayIdDict, primaryColor, taxIdToName, originalGenome, compNetwork)
       
        # Adding links
        G, nodePartners, hasDirections = addLinksToNetwork(G, linksDf, excludedNodes, instanceConfigurations, nodePartners, paramsForNetwork, taxIdToName, hasDirections, isAPI)
       
    return G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes
