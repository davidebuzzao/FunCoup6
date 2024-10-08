from django.shortcuts import render
from django.http import HttpResponse,HttpResponseServerError,JsonResponse
from django.shortcuts import redirect
from data.models import *
from django.db.models import Q
import os
from django.db import connection
import yaml
from yaml.loader import SafeLoader
from .forms import *
import re
import mimetypes
import networkx as nx
from networkx.readwrite import json_graph
import json
from .helperFunctions.colors import *
from .helperFunctions.configInfo import *
from .helperFunctions.scoreCalculations import *
from .helperFunctions.query import *
from .helperFunctions.linkDetails import *
from .helperFunctions.enrichment import *
from .helperFunctions.network import *
from django.apps import apps

# Here, all website views are defined. Each function is called from a specific url (see urls.py)
# and the code in the function is run to collect necessary information. All information collected is packaged in a dict 
# and is passed to the frontend, together with rendering of an html-file.
# All functions called from here are located under either 'externalAlgorithms' or 'helperFunctions'


def index(request):
    # Collecting information and rendering the index page, and the initial search form
    query,instanceConfigurations,parameters,orthologSpecies=getParametersForSearch(request)
    if request.method == 'POST':
        querystring=buildSearchQueryString(request,query,orthologSpecies)
        if querystring != None: # If an error was found, None was returned, then the form is loaded again, including the error text(s)
            return redirect('/network/'+querystring)

    context={'query':query,'parameters':parameters}
    return render(request, "index.html", context)

def network(request,geneQuery,genomeSelect,confidenceThreshold,grgLLRThreshold,depth,nodesPerStepOptions,expansionAlgorithm,prioritizeNeighbors,comparativeGenomes,individualEvidenceOnly,orthologsOnly,categoryID,constrainEvidence,evidenceSpecies,evidenceDataTypes,restriction,restrictPathway,restrictTissue,showAdvanced, missingMapping, multiMapping):
    # Getting all search parameters as an input
    # Rendering the gene selection page in case of unclarities in the mappings
    # Otherwise the network view is rendered. Collecting all necessary information to obtain the network with the given search parameters.
    
    # First, getting all information to be shown in the "modified search"-form
    query,instanceConfigurations,parameters,orthologSpecies=getParametersForSearch(request)

    # Collecting necessary information and rendering the gene selection page in case any mapping was ambiguous or missing, otherwise getting the network
    if multiMapping!="''" or missingMapping!="''": 
        # Preparing a new querystring with cleared up multi/missing mappings
        querystring="/network/###&"+genomeSelect+"&"+confidenceThreshold+"&"+grgLLRThreshold+"&"+depth+"&"+nodesPerStepOptions+"&"+expansionAlgorithm+"&"+prioritizeNeighbors+"&"+comparativeGenomes+"&"+individualEvidenceOnly+"&"+orthologsOnly+"&"+categoryID+"&"+constrainEvidence+"&"+evidenceSpecies+"&"+evidenceDataTypes+"&"+restriction+"&"+restrictPathway+"&"+restrictTissue+"&"+showAdvanced+"&''&''/"
        # Getting info on normal mappings, missing mappings and mappings to multiple uniprotIDs to render the gene selection page
        if genomeSelect not in orthologSpecies:
            normalMappings, multiMappings, missingMappings = ambiguousOrMissingMappings(geneQuery, genomeSelect, missingMapping, multiMapping, instanceConfigurations)
        ## TODO: fix ambiguous mappings for orthology transferred species
        else:
            geneQuery,multiMapping,missingMapping = geneQuery.strip(),multiMapping.strip(),missingMapping.strip()
            allQueryGenes=geneQuery.strip().split(",")
            multiMappings = [m for m in multiMapping.split(",") if m.strip() != "''" and m.strip()]
            missingMappings = [m for m in missingMapping.split(",") if m.strip() != "''" and m.strip()]
            normalMappings=[{'input':g, 'output':g, 'type':'UniProtID', 'description':''} for g in allQueryGenes if g not in multiMappings+missingMappings]
            multiMappings=[{'input':g, 'output':g, 'type':'UniProtID', 'description':''} for g in multiMappings]

        context={'querystring':querystring,'normalMappings':normalMappings,'multiMappings':multiMappings,'missingMappings':missingMappings}
        return render(request, "geneSelection.html", context)


    # If post from the form for "Modify search", redirect with new search parameters
    if request.method == 'POST': 
        querystring=buildSearchQueryString(request,query,orthologSpecies)
        if querystring != None:
            return redirect('/network/'+querystring)
        else: # If errors occured, redirect to network page for modify search with errors
            context={'query':query,'parameters':parameters,'data':{}}
            return render(request, "network.html", context)
    
    # If not a post, all fields for the search form needs to be populated with the data that was used for the search. Also updates input parametrs where necessary
    query, paramsForNetwork = fillSearchForm(query, parameters, geneQuery, genomeSelect, showAdvanced, confidenceThreshold, grgLLRThreshold, depth, nodesPerStepOptions, expansionAlgorithm, prioritizeNeighbors, comparativeGenomes, individualEvidenceOnly, orthologsOnly, categoryID, constrainEvidence, evidenceSpecies, evidenceDataTypes, restriction, restrictPathway,restrictTissue )

    # Defining some common variables for the entire subnetwork
    allNetworkPathways={} # Including all pathways in the subnetwork, used for filter
    allNetworkTissues={} # Including all tissues in the subnetwork, used for filter
    mappingTypes=set() # Keeping track of distinct mapping types (to show labels)
    queryGenes=geneQuery.split(",")
    G=nx.Graph()
    hasDirections=False # Keeping track of if there are ANY directions in the subnetwork. If not the checkbox "Show directions" will be hidden
    originalGeneList=queryGenes.copy() # Keeping original for transferred networks
    originalGenome=genomeSelect
    FCtoQueryMap={} # used for transferred networks
    genomeFileName="" # Used for Topas 
    transferredNetwork=False
    nodePartners={} # Dict for all node partners to be set as attributes to all nodes (used for frontend viz)

    taxIdToName={} # Creating lookup to obtain species short names (like hsa, mmu etc.), names, and name as in file ( like H.sapiens)
    for sp in query.fields['genomeSelect'].choices: # Using the species information already collected for the query
        if sp['name'].startswith("*"):
            shortName=sp['name'][1]+sp['name'].split()[1][:2].upper()
            forFile=sp['name'][1].upper()+"."+sp['name'].split('(')[0].split()[-1].lower()
        else:
            if sp['taxid']=='2697049': ##TODO: adjust in case there are more viruses
                shortName='SARS-CoV-2'
                forFile='SARS-CoV-2'
            else:
                shortName=sp['name'][0]+sp['name'].split()[1][:2].upper()
                forFile=sp['name'][0].upper()+"."+sp['name'].split('(')[0].split()[-1].lower()
        taxIdToName[sp['taxid']]={'short':shortName, 'name':sp['name'], 'forFile':forFile}

    # Setting up for orthology transferred network    
    transferredNetwork, genomeFileName, queryGenes, genomeSelect, FCtoQueryMap, orthologyTransferMessage, paramsForNetwork = setUpTransferred(originalGenome, orthologSpecies, transferredNetwork, queryGenes, FCtoQueryMap, originalGeneList, genomeFileName, genomeSelect, taxIdToName, paramsForNetwork)
    genomeFileName+=taxIdToName[genomeSelect]['forFile'] # Finish prep of file name to use for topas
    # Getting mappings for the query proteins, saving proteinID as key, uniprotID as well as original query symbol
    proteinIdsToQuery={}
    for gene in queryGenes:
        m=IdMapping.objects.filter(Q(mapped_id=gene) & Q(protein__proteome__species__tax_id=genomeSelect))
        for mapping in m:
            if mapping.protein_id not in proteinIdsToQuery:
                proteinIdsToQuery[mapping.protein.id]={}
            proteinIdsToQuery[mapping.protein.id]={'gene': gene, 'uniprot': mapping.protein.uniprot_id}

    if genomeSelect in parameters['instancePathogens']:
        # Getting network for the query pathogen-host species according to parameters
        hostSpecies=parameters['instancePathogens'][genomeSelect]
        primaryColors={hostSpecies: getNodeColors()[0],genomeSelect:getNodeColors()[-1]}
        allViewSpecies=[
            {'color':getNodeColors()[-1],'speciesName':taxIdToName[originalGenome]['name'], 'short':taxIdToName[originalGenome]['short']},
            {'color':getNodeColors()[0],'speciesName':taxIdToName[hostSpecies]['name'], 'short':taxIdToName[hostSpecies]['short']}
            ]
    else:
        # Getting network for the query species according to parameters
        primaryColors={originalGenome: getNodeColors()[0]}
        allViewSpecies=[{'color':getNodeColors()[0],'speciesName':taxIdToName[originalGenome]['name'], 'short':taxIdToName[originalGenome]['short']}]

    G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes = getNetworkForGenes(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, proteinIdsToQuery, instanceConfigurations, genomeSelect,originalGenome,FCtoQueryMap,primaryColors,taxIdToName,genomeFileName,transferredNetwork, paramsForNetwork, compNetwork=False, isAPI=False)
        
    # If no links was found for the query genes, stop here and return an error
    if G.number_of_edges()==0: # If no links could be found, an error is passed to the search form
        parameters['error']="No links could be found for your query. Try lowering the confidence cutoff."
        # Collect parameters for Modify Search
        parameters = collectQueryParameters(parameters,paramsForNetwork,genomeSelect,taxIdToName,originalGenome)
        return render(request, "network.html", {'query':query,'parameters':parameters,'data':{}})

    # In case of comparative interactomics, iterate over all comparative species, and get networks for each with orthologs to query proteins as query      
    if paramsForNetwork['comparativeSpecies']!=["''"]:
        transferredNetwork=False
        for i,comp_species in enumerate(paramsForNetwork['comparativeSpecies']):
            primaryColors[comp_species]=primaryColors.get(comp_species,getNodeColors()[i+1])
            if orthologsOnly=='on':
                ## Extract orthologs for all genes (query and others), find links between them
                genomeFileName,paramsForNetwork,allViewSpecies,orthologQueryAndOther,orthologs = getOrthologInfoForQueryORTHOLOGSonly(genomeFileName, paramsForNetwork, allViewSpecies, genomeSelect, comp_species, taxIdToName, listOfGenes, proteinIdsToQuery, primaryColors)
                if len(orthologQueryAndOther['query'])>0:
                    G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes_orthoSp = getNetworkForGenesORTHOLOGSonly(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, orthologQueryAndOther, instanceConfigurations,comp_species,comp_species,FCtoQueryMap,primaryColors,taxIdToName,genomeFileName,transferredNetwork,paramsForNetwork,compNetwork=True,isAPI=False)
                    for o in orthologs: # Adding ortholog links (green dashed) to the network
                        if str(o.proteinA_id) in nodePartners and str(o.proteinB_id) in nodePartners: # Makes sure that the orthologs are in the network
                            G.add_edge(str(o.proteinA_id),str(o.proteinB_id),scoresPerGoldStandard="",value="orthoLink",direction=0)
                            nodePartners[str(o.proteinA_id)].append(str(o.proteinB_id))
                            nodePartners[str(o.proteinB_id)].append(str(o.proteinA_id))
            else:
                ## Query the comparative species network independently on the query species with the query genes only
                genomeFileName,paramsForNetwork,allViewSpecies,orthologsQuery,orthologs = getOrthologInfoForQuery(genomeFileName, paramsForNetwork, allViewSpecies, genomeSelect, comp_species, taxIdToName, listOfGenes, proteinIdsToQuery, primaryColors)
                if len(orthologsQuery['query'])>0:
                    tmp_proteinIdsToQuery = {}
                    for q in orthologsQuery['query']:
                        tmp_proteinIdsToQuery[q] = orthologsQuery['query'][q]
                    G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes_orthoSp = getNetworkForGenes(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, tmp_proteinIdsToQuery, instanceConfigurations,comp_species,comp_species,FCtoQueryMap,primaryColors,taxIdToName,genomeFileName,transferredNetwork,paramsForNetwork,compNetwork=True,isAPI=False)
                    for o in orthologs: # Adding ortholog links (green dashed) to the network
                        if str(o.proteinA_id) in nodePartners and str(o.proteinB_id) in nodePartners: # Makes sure that the orthologs are in the network
                            G.add_edge(str(o.proteinA_id),str(o.proteinB_id),scoresPerGoldStandard="",value="orthoLink",direction=0)
                            nodePartners[str(o.proteinA_id)].append(str(o.proteinB_id))
                            nodePartners[str(o.proteinB_id)].append(str(o.proteinA_id))

    # Finally preparing all collected data for display in the network view
    degrees = dict(G.degree()) # Getting degrees for all nodes, appending to the node objects in the graph
    nx.set_node_attributes(G, degrees, "degree") # Appending degree information to each node object
    nx.set_node_attributes(G, nodePartners, "linkedNodes") # Appending neighbor nodes to each node object

    graphData=json_graph.node_link_data(G) # Making a json graph from the graph object, to be able to pass it to the front end
    linkDict={} # Organizing links in a dict where keys are proteinA|protienB
    for l in graphData['links']:
        linkDict[l['source']+"|"+l['target']]=l

    nodeDict={} # Organizing nodes in a dict where keys are nodeIDs
    sortedNodes={} # Sorting nodes to show sorted list of nodes in interactor view
    for n in graphData['nodes']:
        nodeDict[n['id']]=n
        sortedNodes[n['id']]=float(n['maxlinkFdr'])
    sortedNodes = sorted(sortedNodes.items(), key=lambda x:x[1])
    sortedNodes = list(map(lambda x: x[0], sortedNodes))
    networkData={'nodes':nodeDict,'links':linkDict}

    # Prepare parameters to pass to the frontend
    parameters['numQueryGenes']=len(queryGenes)
    parameters['subnetworkGenes']=G.number_of_nodes()
    parameters['subnetworkLinks']=G.number_of_edges()
    parameters['species']=allViewSpecies
    parameters['numSpecies']=len(allViewSpecies)
    parameters['mappingTypes']=sorted(mappingTypes)
    parameters['hasDirections']=hasDirections
    parameters['evidenceTypes']=[evidence[1] for evidence in query.fields['evidenceDataTypes'].choices]
    parameters['instanceSpecies']=[(taxIdToName[i]) for i in instanceConfigurations['trainingConfig']['Network']['Species_order'] ] # Used for interactions matrix
    parameters['llrColorsLegend']=[getLLRColor(10),getLLRColor(3.3),getLLRColor(0),getLLRColor(-1),getLLRColor(-3)]
    parameters['error']=""
    parameters['maxlink']="True" if expansionAlgorithm=="maxlink" else "False"
    parameters['orthologyTransferMessage']=orthologyTransferMessage
    parameters['allNetworkPathways']=allNetworkPathways # To be used for filter
    parameters['allNetworkTissues']=allNetworkTissues # To be used for filter
    parameters['hasTissues']= 'true' if len(list(allNetworkTissues))>0 else 'false'
    parameters['hasPathways']= 'true' if len(list(allNetworkPathways))>0 else 'false'
    parameters['sortedNodes']=sortedNodes
    
    parameters = collectQueryParameters(parameters,paramsForNetwork,genomeSelect,taxIdToName,originalGenome)

    context={'query':query,'parameters':parameters,'data':networkData}
    return render(request, "network.html", context)

def linkDetails(request, idA, idB, goldstandard):
    # Obtaining detailed information for a link between idA-idB in a given gold standard
    # To show information when clicking the details button in the interactions view.
    # The fetching of the information is triggered by an Ajax call
    instanceConfigurations = getConfig()

    # First, make sure proteinIDs are in the correct order
    if int(idA)<int(idB): 
        tmpA=idA
        idA=idB
        idB=tmpA
    
    # Collecting some necessary underlying information, species, goldstandardId
    idA_speciesInDb=Protein.objects.get(id=idA).proteome.species
    idA_species=idA_speciesInDb.tax_id
    idB_speciesInDb=Protein.objects.get(id=idB).proteome.species
    idB_species=idB_speciesInDb.tax_id

    ## Host-pathogen interactions are scored using host goldstandards
    host_pathogen_flag=False
    pathogens = dict([(pathogenA,hostA) for hostA,pathogenA,pathogenA_name in instanceConfigurations['instanceConfig']['instance']['hostPathogen']])
    if any(id in pathogens for id in [idA_species,idB_species]):
        host_pathogen_flag=True
        if idA_species in pathogens:
            species=pathogens[idA_species] #host species
            query_species=idA_species
        else: 
            species=pathogens[idB_species] #host species
            query_species=idB_species
    else:
        species=idA_species
        query_species=species

    speciesInDb=Species.objects.get(tax_id=species)

    gsVersion=instanceConfigurations['goldStandardConfig'][goldstandard]['version']
    if gsVersion!="":
        goldStandardInDb=GoldStandard.objects.filter(Q(type=goldstandard) & Q(species=speciesInDb) & Q(version=instanceConfigurations['goldStandardConfig'][goldstandard]['version']))[0]
    else:
        goldStandardInDb=GoldStandard.objects.filter(Q(type=goldstandard) & Q(species=speciesInDb))[0]
    goldStandardIndex=instanceConfigurations['trainingConfig']['Network']['GoldStandard_order'].index(goldstandard)
    
    # Collecting information on the link from the network table
    linkInfo = Network.objects.filter(Q(proteinA_id=idA) & Q(proteinB_id=idB))

    evidenceDetails=[]
    if linkInfo.exists(): # It should exist, but just to make sure
        linkInfo=linkInfo[0]
        evidenceSupportingLink=linkInfo.llr_evidence.split(";")[goldStandardIndex].split("|")[0].split(",")
        speciesSupportingLink=linkInfo.llr_evidence.split(";")[goldStandardIndex].split("|")[1].split(",")
        orthologTransfer=instanceConfigurations['trainingConfig']['OrthologyTransfer']['Evidences']
        orthologMap={} # taxID : {'A': [pA, pB], 'B':[pC,pD]}
        
        if host_pathogen_flag:
            host_species_index = instanceConfigurations['trainingConfig']['Network']['Species_order'].index(species)
            llrForSpecies=speciesSupportingLink[host_species_index]
            if llrForSpecies!="" and llrForSpecies!="0.0": 
                orthologMap[query_species]={'A':[idA],'B':[idB]}  
        else:

            # Get all orthologs for the genes, in order to fetch information on the evidence from different species
            for i,sp in enumerate(instanceConfigurations['trainingConfig']['Network']['Species_order']):
                llrForSpecies=speciesSupportingLink[i]
                if llrForSpecies!="" and llrForSpecies!="0.0": 
                    if sp==species: # If query species
                        orthoA=[idA]
                        orthoB=[idB]
                    else:
                        if int(species)>int(sp): # To get correct order of orthologs in DB query
                            orthoA=ProteinOrtholog.objects.filter(Q(speciesA=species) & Q(speciesB=sp) & Q(proteinA=idA)).values_list("proteinB_id", flat=True)
                            orthoB=ProteinOrtholog.objects.filter(Q(speciesA=species) & Q(speciesB=sp) & Q(proteinA=idB)).values_list("proteinB_id", flat=True)
                        else:
                            orthoA=ProteinOrtholog.objects.filter(Q(speciesB=species) & Q(speciesA=sp) & Q(proteinB=idA)).values_list("proteinA_id", flat=True)
                            orthoB=ProteinOrtholog.objects.filter(Q(speciesB=species) & Q(speciesA=sp) & Q(proteinB=idB)).values_list("proteinA_id", flat=True)
                    orthologMap[sp]={'A':list(orthoA),'B':list(orthoB)}

        dataForChart={} # Evidence type distribution to be put in bar chart
        # Iterating over all evidence types
        for i,evidenceType in enumerate(instanceConfigurations['trainingConfig']['Network']['Evidence_order']):
            llrForEv=evidenceSupportingLink[i]
            if llrForEv!="" and llrForEv!="0.0":
                if float(llrForEv)>0:
                    dataForChart[evidenceType]=float(llrForEv)
                evidenceDetails = getDetailsForEvidence(evidenceType,evidenceDetails,species,query_species,goldStandardInDb,orthologMap,orthologTransfer,host_pathogen_flag)
        # Sorting by LLR
        evidenceDetails.sort(key=lambda x: x['llr'], reverse=True)

    # Getting uniprotIDs for the protein ids
    uniprotIdA=Protein.objects.get(id=idA).uniprot_id
    uniprotIdB=Protein.objects.get(id=idB).uniprot_id
    
    # Collecting data for viewing in details table 
    context={'idA':idA, 'uniprotIdA':uniprotIdA, 'idB':idB, "uniprotIdB": uniprotIdB,'goldstandard':goldstandard,'species':species,'evidenceDetails':evidenceDetails,'dataForChart':dataForChart}
    return HttpResponse(json.dumps(context), content_type="application/json")


def enrichment(request, queryGenes, species):
    # Obtaining enriched pathways for the subnetwork genes
    # To show information in the enrichment tab
    # The fetching of the information is triggered by an Ajax call, which is initiated after page load of the network view
    # Note that javascript to handle this is in network.html. from there information is filled into the enrichment.html template
    instanceConfigurations = getConfig()

    if os.path.exists('configFiles/websiteConfig.yml'):
        with open('configFiles/websiteConfig.yml') as f:
            webConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your websiteConfig.yml does not exist")
    
    cutoff_confidence=float(webConfig['anubixCutoff'])
    fdr_cutoff,fwer_cutoff=0.05,0.05

    spList=species.split("|") # Passed as a string separated by | to be able to handle multiple species
    queryGenes=queryGenes.split("|") # Passed as a string separated by | to be able to handle multiple species

    enrichmentsPerSpecies=[]
    for i,sp in enumerate(spList): # Iterating over all species to get enrichments for
        if not sp in instanceConfigurations['instanceConfig']['instance']['species']: continue

        enrichedPathways={}
        species_objects = Species.objects.get(Q(tax_id=sp))
        speciesA_name = species_objects.species_name
        speciesA_name_short=speciesA_name.lower()[0]+speciesA_name.lower().split()[1][:2]
        query_set=list(set(list(map(int, queryGenes[i].split(",")))))
        
        pathways = Pathway.objects.filter(Q(species=sp)).values('pathway_id','pathway_name').distinct()
        pathwayMap={} # maps pathway ID to pathway name
        for p in pathways:
            pathwayMap[p['pathway_id']]=p['pathway_name']

        # ANUBIX
        enrichedPathways = getAnubix(query_set, sp, cutoff_confidence, fwer_cutoff, enrichedPathways, pathwayMap)  
        # BINOX TODO: faster python implementation, now ~2.30min on human
        # enrichedPathways = getBinox(query_set, sp, cutoff_confidence, fwer_cutoff, enrichedPathways, pathwayMap)
        # EASE
        enrichedPathways = getEase(query_set, sp, fwer_cutoff, enrichedPathways, pathwayMap)
        
        # Extracting and sorting all enrichments based on min FWER values
        enrichedPathways=enrichedPathways.values()
        sortedEnrichments = sorted(enrichedPathways, key=lambda d: d['min_fwer']) 

        # Appending all information per species
        enrichmentsPerSpecies.append({
            'pathways':sortedEnrichments,
            'speciesName':speciesA_name,
            'speciesNameShort':speciesA_name_short,
            'taxID':sp
        })

    context={'enrichmentsPerSpecies':enrichmentsPerSpecies}
    return HttpResponse(json.dumps(context), content_type="application/json")


def statistics(request):
    instanceConfigurations = getConfig()
    instanceSpecies = Species.objects.filter(tax_id__in=instanceConfigurations['instanceConfig']['instance']['species'])
    
    ## Evidence data
    evidenceTypes = set(instanceConfigurations['instanceConfig']['instance']['evidence'])
    instanceData=[]
    for e in evidenceTypes:
        evidence = Evidence.objects.filter(Q(species__in=instanceSpecies) & Q(type=e) & Q(version=instanceConfigurations['evidenceConfig'][e]['version']))
        if evidence.exists():
            for ev in evidence:
                instanceData.append({'species':ev.species.species_name,'type':ev.type,'name':ev.version,'url':instanceConfigurations['evidenceConfig'][e]['sourceUrl']})
    
    ## Gold Standard data
    goldStandardTypes = set(instanceConfigurations['instanceConfig']['instance']['goldStandard'])
    instanceGoldStandards=[]
    for e in goldStandardTypes:
        gs = GoldStandard.objects.filter(Q(species__in=instanceSpecies) & Q(type=e) & Q(version=instanceConfigurations['goldStandardConfig'][e]['version']))
        if gs.exists():
            for g in gs:
                instanceGoldStandards.append({'species':g.species.species_name,'type':g.type,'name':g.version,'url':instanceConfigurations['goldStandardConfig'][e]['sourceUrl']})

    ## Species data
    # goldStandardTypes = set(instanceConfigurations['instanceConfig']['instance']['goldStandard'])
    # instanceGoldStandards=[]

    context={'instanceData':instanceData,'instanceGoldStandards':instanceGoldStandards}
    return render(request, "statistics.html", context)


def contact(request):
    context={}
    return render(request, "contact.html", context)


def api(request):
    context={}
    return render(request, "api.html", context)


def help(request):
    context={}
    return render(request, "help.html", context)


def archive(request):
    files=[]
    for dir in os.listdir('website/static/website/networks/'):
        if dir!="FunCoup6.0" and dir.startswith("FunCoup"):
            filesForVersion=[]
            for file in os.listdir('website/static/website/networks/'+dir):
                size = str(round((os.path.getsize('website/static/website/networks/'+dir+"/"+file))/1000000, 2))+" MB"
                filesForVersion.append({'name':file,'size':size,'url':"/download/archive&"+file})
            files.append({'version': dir,'filesList': filesForVersion})
    
    # Example in Python
    files.sort(key=lambda x: x['version'], reverse=True)
    context={'files':files}
    return render(request, "archive.html", context)


def downloads(request):
    fileListFull=[]
    fileListCompact=[]
    fileListFullTransferred=[]
    fileListCompactTransferred=[]
    path_to_files = 'website/static/website/networks/FunCoup6.0/'
    file_list = sorted(os.listdir(path_to_files))
    for file in file_list:
        if file.endswith("_full.gz"):
            if "TRANSFERRED" not in file:
                size = str(round((os.path.getsize(path_to_files+file))/1000000, 2))+" MB"
                fileListFull.append({'name':file,'size':size,'url':"/download/network&"+file})
            else:
                size = str(round((os.path.getsize(path_to_files+file))/1000000, 2))+" MB"
                spName=file.split("TRANSFERRED_")[1].split("_from")[0].replace("_"," ").replace("("," (")
                fileListFullTransferred.append({'name':spName,'size':size,'url':"/download/network&"+file})
        elif file.endswith("_compact.gz"):
            if "TRANSFERRED" not in file:
                size = str(round((os.path.getsize(path_to_files+file))/1000000, 2))+" MB"
                fileListCompact.append({'name':file,'size':size,'url':"/download/network&"+file})
            else:
                size = str(round((os.path.getsize(path_to_files+file))/1000000, 2))+" MB"
                spName=file.split("TRANSFERRED_")[1].split("_from")[0].replace("_"," ").replace("("," (")
                fileListCompactTransferred.append({'name':spName,'size':size,'url':"/download/network&"+file})

    context={'filesFull':fileListFull,'filesCompact':fileListCompact, 'fileListCompactTransferred':fileListCompactTransferred, 'fileListFullTransferred':fileListFullTransferred}
    return render(request, "downloads.html", context)


def download(request, type, filename):
    if type=="network":
        filepath = 'website/static/website/networks/FunCoup6.0/'+filename
    elif type=="archive":
        for dir in os.listdir('website/static/website/networks/'):
            for file in os.listdir('website/static/website/networks/'+dir):
                if file==filename:
                    filepath='website/static/website/networks/'+dir+"/"+filename
    filename = filename
    path = open(filepath, 'rb')
    mime_type, _ = mimetypes.guess_type(filepath)
    response = HttpResponse(path, content_type=mime_type)
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    return response

###################################################
######## REST API
###################################################

## TSV to work with PathwAX
def apiGeneTSV(request,taxid,query,identifiersource=""):
    try:
        queryList=re.split('\W+', query)
        geneMatrix="#Keywords\tinternalGeneID\tgeneID\tidentifierType\n"
        if identifiersource=="":
            mappings = IdMapping.objects.filter(Q(mapped_id__in=queryList) & Q(protein__proteome__species__tax_id=taxid))
        else:
            mappings = IdMapping.objects.filter(Q(mapped_id__in=queryList) & Q(protein__proteome__species__tax_id=taxid) & Q(type=identifiersource))   
        
        ## To remove duplicates in IdMapping table
        existing_mapping = []
        for m in mappings:
            geneString = (m.mapped_id,str(m.protein_id),m.protein.uniprot_id,m.type)
            if geneString not in existing_mapping:
                existing_mapping.append(geneString)
                geneMatrix+=geneString[0]+"\t"+geneString[1]+"\t"+geneString[2]+"\t"+geneString[3]+"\n"

        response = HttpResponse(geneMatrix, content_type="text/plain")
        response['Content-Disposition'] = "filename='species.tsv'"
        return response
    except Exception as e:
        # Log the error for debugging purposes
        print(f"An error occurred: {e}")
        # Return an appropriate error response
        return HttpResponseServerError("An error occurred while processing the request.")

def apiSpeciesTSV(request):
    try:
        speciesId_df = pd.DataFrame.from_records(Stats.objects.filter(Q(num_links__gt=0)).values('tax_id','species_name','common_name','origin'))
        speciesId_df.columns = ['tax_id','name','common_name','origin']
        speciesDictionary = speciesId_df.set_index('tax_id').T.to_dict()
        
        speciesMatrix = "#Scientific_name\tCommon_name\tNCBI_taxonomy\tOrigin\n"
        for tax_id in speciesDictionary:
            speciesMatrix += f"{speciesDictionary[tax_id]['name']}\t{speciesDictionary[tax_id]['common_name']}\t{tax_id}\t{speciesDictionary[tax_id]['origin']}\n"

        response = HttpResponse(speciesMatrix, content_type="text/plain")
        response['Content-Disposition'] = "filename='species.tsv'"
        return response
    except Exception as e:
        # Log the error for debugging purposes
        print(f"An error occurred: {e}")
        # Return an appropriate error response
        return HttpResponseServerError("An error occurred while processing the request.")

###################################################
## JSON to work with Cytoscape app
def apiSpeciesJSON(request):
    try:
        speciesId_df = pd.DataFrame.from_records(Stats.objects.filter(Q(num_links__gt=0)).values('tax_id','species_name','common_name','origin'))
        speciesId_df.columns = ['tax_id','name','common_name','origin']
        speciesDictionary = speciesId_df.set_index('tax_id').T.to_dict()

        # Return JSON response
        return JsonResponse(speciesDictionary, safe=False)

    except Exception as e:
        # Log the error for debugging purposes
        print(f"An error occurred: {e}")
        # Return an appropriate error response
        return HttpResponseServerError("An error occurred while processing the request.")

def apiGeneJSON(request,geneQuery,genomeSelect):
    try:
        query,instanceConfigurations,parameters,orthologSpecies=getParametersForSearch(request, isAPI=True)
        geneQuery,missingMapping,multiMapping=buildSearchQueryStringAPI(geneQuery,genomeSelect,orthologSpecies)

        # Getting info on normal mappings, missing mappings and mappings to multiple uniprotIDs to render the gene selection page
        geneQuery = ','.join(geneQuery)
        missingMapping = ','.join(missingMapping)
        multiMapping = ','.join(multiMapping)
        if genomeSelect not in orthologSpecies:
            normalMappings, multiMappings, missingMappings = ambiguousOrMissingMappings(geneQuery,genomeSelect,missingMapping,multiMapping,instanceConfigurations)
        else: 
            geneQuery,multiMapping,missingMapping = geneQuery.strip(),multiMapping.strip(),missingMapping.strip()
            allQueryGenes=geneQuery.strip().split(",")
            multiMappings = [m for m in multiMapping.split(",") if m.strip() != "''" and m.strip()]
            missingMappings = [m for m in missingMapping.split(",") if m.strip() != "''" and m.strip()]
            normalMappings=[{'input':g, 'output':g, 'type':'UniProtID', 'description':''} for g in allQueryGenes if g not in multiMappings+missingMappings]
            multiMappings=[{'input':g, 'output':g, 'type':'UniProtID', 'description':''} for g in multiMappings]

        geneQueryMappings = {'successfulMappings':normalMappings,'ambiguousMappings':multiMappings,'unsuccessfulMappings':missingMappings}

        # Return JSON response
        return JsonResponse(geneQueryMappings, safe=False)
        
    except Exception as e:
        # Log the error for debugging purposes
        print(f"An error occurred: {e}")
        # Return an appropriate error response
        return HttpResponseServerError("An error occurred while processing the request.")


def apiNetworkJSON(request,geneQuery,genomeSelect,confidenceThreshold,grgLLRThreshold,depth,nodesPerStepOptions,expansionAlgorithm,prioritizeNeighbors,comparativeGenomes,individualEvidenceOnly,orthologsOnly):
    try:
        # Getting all search parameters as an input
        # Using first hit in the gene selection page in case of unclarities in the mappings
        # Collecting all necessary information to obtain the network with the given search parameters.
        # Outputting the network view as tsv.
        
        ## TEST
        # https://funcoup6.scilifelab.se/api/json/network/MYC&9606&0.95&1&15&group&on&off&''&off

        # First, getting all information to be shown in the "modified search"-form
        query,instanceConfigurations,parameters,orthologSpecies=getParametersForSearch(request)
    
        # If not a post, all fields for the search form needs to be populated with the data that was used for the search. Also updates input parametrs where necessary
        categoryID = "all"
        constrainEvidence= "all" 
        evidenceSpecies= ""
        evidenceDataTypes= ""
        restriction= "all"
        restrictPathway= ""
        restrictTissue = ""
        showAdvanced= "false"

        query, paramsForNetwork = fillSearchForm(query, parameters, geneQuery, genomeSelect, showAdvanced, confidenceThreshold, grgLLRThreshold, depth, nodesPerStepOptions, expansionAlgorithm, prioritizeNeighbors, comparativeGenomes, individualEvidenceOnly, orthologsOnly, categoryID, constrainEvidence, evidenceSpecies, evidenceDataTypes, restriction, restrictPathway,restrictTissue )

        # Defining some common variables for the entire subnetwork
        allNetworkPathways={} # Including all pathways in the subnetwork, used for filter
        allNetworkTissues={} # Including all tissues in the subnetwork, used for filter
        mappingTypes=set() # Keeping track of distinct mapping types (to show labels)
        queryGenes=geneQuery.split(",")
        G=nx.Graph()
        hasDirections=False # Keeping track of if there are ANY directions in the subnetwork. If not the checkbox "Show directions" will be hidden
        originalGeneList=queryGenes.copy() # Keeping original for transferred networks
        originalGenome=genomeSelect
        FCtoQueryMap={} # used for transferred networks
        genomeFileName="" # Used for Topas 
        transferredNetwork=False
        nodePartners={} # Dict for all node partners to be set as attributes to all nodes (used for frontend viz)

        taxIdToName={} # Creating lookup to obtain species short names (like hsa, mmu etc.), names, and name as in file ( like H.sapiens)
        # query_genomeSelect_choices = [{'taxid': '3702', 'name': 'Arabidopsis thaliana (taxid:3702)'}, {'taxid': '224308', 'name': 'Bacillus subtilis (strain 168) (taxid:224308)'}, {'taxid': '9913', 'name': 'Bos taurus (taxid:9913)'}, {'taxid': '6239', 'name': 'Caenorhabditis elegans (taxid:6239)'}, {'taxid': '9615', 'name': 'Canis lupus familiaris (taxid:9615)'}, {'taxid': '7719', 'name': 'Ciona intestinalis (taxid:7719)'}, {'taxid': '7227', 'name': 'Drosophila melanogaster (taxid:7227)'}, {'taxid': '44689', 'name': 'Dictyostelium discoideum (taxid:44689)'}, {'taxid': '7955', 'name': 'Danio rerio (taxid:7955)'}, {'taxid': '83333', 'name': 'Escherichia coli (strain K12) (taxid:83333)'}, {'taxid': '9031', 'name': 'Gallus gallus (taxid:9031)'}, {'taxid': '9606', 'name': 'Homo sapiens (taxid:9606)'}, {'taxid': '243232', 'name': 'Methanocaldococcus jannaschii (strain ATCC 43067 / DSM 2661 / JAL-1 / JCM 10045 / NBRC 100440) (taxid:243232)'}, {'taxid': '10090', 'name': 'Mus musculus (taxid:10090)'}, {'taxid': '83332', 'name': 'Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv) (taxid:83332)'}, {'taxid': '39947', 'name': 'Oryza sativa subsp. japonica (taxid:39947)'}, {'taxid': '36329', 'name': 'Plasmodium falciparum (taxid:36329)'}, {'taxid': '10116', 'name': 'Rattus norvegicus (taxid:10116)'}, {'taxid': '559292', 'name': 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (taxid:559292)'}, {'taxid': '284812', 'name': 'Schizosaccharomyces pombe (strain 972 / ATCC 24843) (taxid:284812)'}, {'taxid': '9823', 'name': 'Sus scrofa (taxid:9823)'}, {'taxid': '273057', 'name': 'Saccharolobus solfataricus (strain ATCC 35092 / DSM 1617 / JCM 11322 / P2) (taxid:273057)'}, {'taxid': 'SARSCOV2', 'name': 'SARS-CoV-2 - Homo sapiens (taxid:0)'}, {'taxid': '387092', 'name': '*Nitratiruptor sp. (strain SB155-2) (taxid:387092)'}, {'taxid': '760142', 'name': '*Hippea maritima (strain ATCC 700847  DSM 10411  MH2) (taxid:760142)'}, {'taxid': '2043170', 'name': '*Paremcibacter congregatus (taxid:2043170)'}, {'taxid': '1123071', 'name': '*Rubritalea squalenifaciens DSM 18772 (taxid:1123071)'}, {'taxid': '122586', 'name': '*Neisseria meningitidis serogroup B (strain MC58) (taxid:122586)'}, {'taxid': '1620215', 'name': '*Sulfuricaulis limicola (taxid:1620215)'}]

        taxIdToName={} # Creating lookup to obtain species short names (like hsa, mmu etc.), names, and name as in file ( like H.sapiens)
        for sp in query.fields['genomeSelect'].choices: # Using the species information already collected for the query
            if sp['name'].startswith("*"):
                shortName=sp['name'][1]+sp['name'].split()[1][:2].upper()
                forFile=sp['name'][1].upper()+"."+sp['name'].split()[1].lower()
            else:
                shortName=sp['name'][0]+sp['name'].split()[1][:2].upper()
                forFile=sp['name'][0].upper()+"."+sp['name'].split()[1].lower()
            taxIdToName[sp['taxid']]={'short':shortName, 'name':sp['name'],'forFile':forFile}

        # Setting up for orthology transferred network    
        transferredNetwork, genomeFileName, queryGenes, genomeSelect, FCtoQueryMap, orthologyTransferMessage, paramsForNetwork = setUpTransferred(originalGenome, orthologSpecies, transferredNetwork, queryGenes, FCtoQueryMap, originalGeneList, genomeFileName, genomeSelect, taxIdToName, paramsForNetwork)
        genomeFileName+=taxIdToName[genomeSelect]['forFile'] # Finish prep of file name to use for topas
        
        # Getting mappings for the query proteins, saving proteinID as key, uniprotID as well as original query symbol
        proteinIdsToQuery={}
        for gene in queryGenes:
            m=IdMapping.objects.filter(Q(mapped_id=gene) & Q(protein__proteome__species__tax_id=genomeSelect))
            for mapping in m:
                if mapping.protein_id not in proteinIdsToQuery:
                    proteinIdsToQuery[mapping.protein.id]={}
                proteinIdsToQuery[mapping.protein.id]={'gene': gene, 'uniprot': mapping.protein.uniprot_id}

        if genomeSelect in parameters['instancePathogens']:
            # Getting network for the query pathogen-host species according to parameters
            hostSpecies=parameters['instancePathogens'][genomeSelect]
            primaryColors={hostSpecies: getNodeColors()[0],genomeSelect:getNodeColors()[-1]}
            allViewSpecies=[
                {'color':getNodeColors()[0],'speciesName':taxIdToName[hostSpecies]['name'], 'short':taxIdToName[hostSpecies]['short']},
                {'color':getNodeColors()[-1],'speciesName':taxIdToName[originalGenome]['name'], 'short':taxIdToName[originalGenome]['short']}
                ]
        else:
            # Getting network for the query species according to parameters
            primaryColors={originalGenome: getNodeColors()[0]}
            allViewSpecies=[{'color':getNodeColors()[0],'speciesName':taxIdToName[originalGenome]['name'], 'short':taxIdToName[originalGenome]['short']}]
        
        G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes = getNetworkForGenes(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, proteinIdsToQuery, instanceConfigurations, genomeSelect,originalGenome,FCtoQueryMap,primaryColors,taxIdToName,genomeFileName,transferredNetwork, paramsForNetwork, compNetwork=False, isAPI=True)
        
        # If no links was found for the query genes, stop here and return an error
        if G.number_of_edges()==0: # If no links could be found, an error is passed to the search form
            # Log the error for debugging purposes
            print(f"An error occurred")
            # Return an appropriate error response
            return HttpResponseServerError("No links could be found for your query. Try lowering the confidence cutoff.")

        # In case of comparative interactomics, iterate over all comparative species, and get networks for each with orthologs to query proteins as query      
        if paramsForNetwork['comparativeSpecies']!=["''"]:
            transferredNetwork=False
            for i,comp_species in enumerate(paramsForNetwork['comparativeSpecies']):
                primaryColors[comp_species]=primaryColors.get(comp_species,getNodeColors()[i+1])
                
                if orthologsOnly=='on':
                    ## Extract orthologs for all genes (query and others), find links between them
                    genomeFileName,paramsForNetwork,allViewSpecies,orthologQueryAndOther,orthologs = getOrthologInfoForQueryORTHOLOGSonly(genomeFileName, paramsForNetwork, allViewSpecies, genomeSelect, comp_species, taxIdToName, listOfGenes, proteinIdsToQuery, primaryColors)
                    if len(orthologQueryAndOther['query'])>0:
                        G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes_orthoSp = getNetworkForGenesORTHOLOGSonly(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, orthologQueryAndOther, instanceConfigurations,comp_species,comp_species,FCtoQueryMap,primaryColors,taxIdToName,genomeFileName,transferredNetwork,paramsForNetwork,compNetwork=True,isAPI=True)
                        for o in orthologs: # Adding ortholog links (green dashed) to the network
                            if str(o.proteinA_id) in nodePartners and str(o.proteinB_id) in nodePartners: # Makes sure that the orthologs are in the network
                                G.add_edge(str(o.proteinA_id),str(o.proteinB_id),scoresPerGoldStandard="",value="orthoLink",direction=0)
                                nodePartners[str(o.proteinA_id)].append(str(o.proteinB_id))
                                nodePartners[str(o.proteinB_id)].append(str(o.proteinA_id))
                else:
                    ## Query the comparative species network independently on the query species with the query genes only
                    genomeFileName,paramsForNetwork,allViewSpecies,orthologsQuery,orthologs = getOrthologInfoForQuery(genomeFileName, paramsForNetwork, allViewSpecies, genomeSelect, comp_species, taxIdToName, listOfGenes, proteinIdsToQuery, primaryColors)        
                    if len(orthologsQuery['query'])>0:
                        tmp_proteinIdsToQuery = {}
                        for q in orthologsQuery['query']:
                            tmp_proteinIdsToQuery[q] = orthologsQuery['query'][q]

                        G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, listOfGenes_orthoSp = getNetworkForGenes(G, nodePartners, mappingTypes, hasDirections, allNetworkPathways, allNetworkTissues, tmp_proteinIdsToQuery, instanceConfigurations,comp_species,comp_species,FCtoQueryMap,primaryColors,taxIdToName,genomeFileName,transferredNetwork,paramsForNetwork,compNetwork=True,isAPI=True)
                        for o in orthologs: # Adding ortholog links (green dashed) to the network
                            if str(o.proteinA_id) in nodePartners and str(o.proteinB_id) in nodePartners: # Makes sure that the orthologs are in the network
                                G.add_edge(str(o.proteinA_id),str(o.proteinB_id),scoresPerGoldStandard="",value="orthoLink",direction=0)
                                nodePartners[str(o.proteinA_id)].append(str(o.proteinB_id))
                                nodePartners[str(o.proteinB_id)].append(str(o.proteinA_id))

                            
        # Finally preparing all collected data for display in the network view
        degrees = dict(G.degree()) # Getting degrees for all nodes, appending to the node objects in the graph
        nx.set_node_attributes(G, degrees, "degree") # Appending degree information to each node object
        nx.set_node_attributes(G, nodePartners, "linkedNodes") # Appending neighbor nodes to each node object

        graphData=json_graph.node_link_data(G) # Making a json graph from the graph object, to be able to pass it to the front end
        linkDict={} # Organizing links in a dict where keys are proteinA|protienB
        for l in graphData['links']:
            linkDict[l['source']+"|"+l['target']]=l
        nodeDict={} # Organizing nodes in a dict where keys are nodeIDs
        sortedNodes={} # Sorting nodes to show sorted list of nodes in interactor view
        for n in graphData['nodes']:
            nodeDict[n['id']]=n
            sortedNodes[n['id']]=float(n['maxlinkFdr'])
        sortedNodes = sorted(sortedNodes.items(), key=lambda x:x[1])
        sortedNodes = list(map(lambda x: x[0], sortedNodes))
        networkData={'nodes':nodeDict,'links':linkDict}
        
        # Return JSON response
        return JsonResponse(networkData, safe=False)

    except Exception as e:
        # Log the error for debugging purposes
        print(f"An error occurred: {e}")
        # Return an appropriate error response
        return HttpResponseServerError("An error occurred while processing the request.")


def linkDetailsJSON(request, idA, idB, goldstandard):
    # Obtaining detailed information for a link between idA-idB in a given gold standard
    # To show information when clicking the details button in the interactions view.
    # The fetching of the information is triggered by an Ajax call
    instanceConfigurations = getConfig()

    # First, make sure proteinIDs are in the correct order
    if int(idA)<int(idB): 
        tmpA=idA
        idA=idB
        idB=tmpA
    
    # Collecting some necessary underlying information, species, goldstandardId
    idA_speciesInDb=Protein.objects.get(id=idA).proteome.species
    idA_species=idA_speciesInDb.tax_id
    idB_speciesInDb=Protein.objects.get(id=idB).proteome.species
    idB_species=idB_speciesInDb.tax_id

    ## Host-pathogen interactions are scored using host goldstandards
    host_pathogen_flag=False
    pathogens = dict([(pathogenA,hostA) for hostA,pathogenA,pathogenA_name in instanceConfigurations['instanceConfig']['instance']['hostPathogen']])
    if any(id in pathogens for id in [idA_species,idB_species]):
        host_pathogen_flag=True
        if idA_species in pathogens:
            species=pathogens[idA_species] #host species
            query_species=idA_species
        else: 
            species=pathogens[idB_species] #host species
            query_species=idB_species
    else:
        species=idA_species
        query_species=species

    speciesInDb=Species.objects.get(tax_id=species)

    gsVersion=instanceConfigurations['goldStandardConfig'][goldstandard]['version']
    if gsVersion!="":
        goldStandardInDb=GoldStandard.objects.filter(Q(type=goldstandard) & Q(species=speciesInDb) & Q(version=instanceConfigurations['goldStandardConfig'][goldstandard]['version']))[0]
    else:
        goldStandardInDb=GoldStandard.objects.filter(Q(type=goldstandard) & Q(species=speciesInDb))[0]
    goldStandardIndex=instanceConfigurations['trainingConfig']['Network']['GoldStandard_order'].index(goldstandard)
    
    # Collecting information on the link from the network table
    linkInfo = Network.objects.filter(Q(proteinA_id=idA) & Q(proteinB_id=idB))

    evidenceDetails=[]
    if linkInfo.exists(): # It should exist, but just to make sure
        linkInfo=linkInfo[0]
        evidenceSupportingLink=linkInfo.llr_evidence.split(";")[goldStandardIndex].split("|")[0].split(",")
        speciesSupportingLink=linkInfo.llr_evidence.split(";")[goldStandardIndex].split("|")[1].split(",")
        orthologTransfer=instanceConfigurations['trainingConfig']['OrthologyTransfer']['Evidences']
        orthologMap={} # taxID : {'A': [pA, pB], 'B':[pC,pD]}
        
        ## TODO --> To adjust if species in pathogens
        if host_pathogen_flag:
            # Get all orthologs for the genes, in order to fetch information on the evidence from different species
            host_species_index = instanceConfigurations['trainingConfig']['Network']['Species_order'].index(species)
            llrForSpecies=speciesSupportingLink[host_species_index]
            if llrForSpecies!="" and llrForSpecies!="0.0": 
                orthologMap[query_species]={'A':[idA],'B':[idB]}           
        else:
            # Get all orthologs for the genes, in order to fetch information on the evidence from different species
            for i,sp in enumerate(instanceConfigurations['trainingConfig']['Network']['Species_order']):
                llrForSpecies=speciesSupportingLink[i]
                if llrForSpecies!="" and llrForSpecies!="0.0": 
                    if sp==species: # If query species
                        orthoA=[idA]
                        orthoB=[idB]
                    else:
                        if int(species)>int(sp): # To get correct order of orthologs in DB query
                            orthoA=ProteinOrtholog.objects.filter(Q(speciesA=species) & Q(speciesB=sp) & Q(proteinA=idA)).values_list("proteinB_id", flat=True)
                            orthoB=ProteinOrtholog.objects.filter(Q(speciesA=species) & Q(speciesB=sp) & Q(proteinA=idB)).values_list("proteinB_id", flat=True)
                        else:
                            orthoA=ProteinOrtholog.objects.filter(Q(speciesB=species) & Q(speciesA=sp) & Q(proteinB=idA)).values_list("proteinA_id", flat=True)
                            orthoB=ProteinOrtholog.objects.filter(Q(speciesB=species) & Q(speciesA=sp) & Q(proteinB=idB)).values_list("proteinA_id", flat=True)
                    orthologMap[sp]={'A':list(orthoA),'B':list(orthoB)}

        dataForChart={} # Evidence type distribution to be put in pie chart
        # Iterating over all evidence types
        for i,evidenceType in enumerate(instanceConfigurations['trainingConfig']['Network']['Evidence_order']):
            llrForEv=evidenceSupportingLink[i]
            if llrForEv!="" and llrForEv!="0.0":
                if float(llrForEv)>0:
                    dataForChart[evidenceType]=float(llrForEv)
                evidenceDetails = getDetailsForEvidence(evidenceType,evidenceDetails,species,query_species,goldStandardInDb,orthologMap,orthologTransfer,host_pathogen_flag)
        # Sorting by LLR
        evidenceDetails.sort(key=lambda x: x['llr'], reverse=True)

    # Getting uniprotIDs for the protein ids
    uniprotIdA=Protein.objects.get(id=idA).uniprot_id
    uniprotIdB=Protein.objects.get(id=idB).uniprot_id
    
    # Collecting data for viewing in details table 
    context={'idA':idA, 'uniprotIdA':uniprotIdA,'idB':idB, "uniprotIdB": uniprotIdB,'goldstandard':goldstandard,'species':species,'evidenceDetails':evidenceDetails,'dataForChart':dataForChart}
    return JsonResponse(context, safe=False)
