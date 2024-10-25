import re
from data.models import *
from django.db.models import Q
from django.utils.safestring import mark_safe
from .configInfo import *
from ..forms import *

def getParametersForSearch(request,isAPI=False):
    # Fetching all data and parameters needed to fill the search form with information
    instanceConfigurations = getConfig()

    instanceSpeciesUnordered = Species.objects.filter(tax_id__in=instanceConfigurations['instanceConfig']['instance']['species'])
    exampleSpecies = instanceSpeciesUnordered.filter(tax_id=instanceConfigurations['webConfig']['exampleSpecies'])[0].species_name+" ("+instanceSpeciesUnordered.filter(tax_id=instanceConfigurations['webConfig']['exampleSpecies'])[0].tax_id+")"
    
    # Getting all instance species in order
    instanceSpecies=[]
    for orderSp in instanceConfigurations['trainingConfig']['Network']['Species_order']:
        for sp in instanceSpeciesUnordered:
            if orderSp==sp.tax_id:
                instanceSpecies.append(sp)

    # Add here species info from instanceConfig file species section
    allSpecies = [{'taxid':species.tax_id, 'name':species.species_name + " (" + species.tax_id+")"} for species in instanceSpecies]
    # Add here HOST PATHOGEN info from instanceConfig file hostPathogen section
    pathogens = dict([(pathogenA,hostA) for hostA,pathogenA,pathogenA_name in instanceConfigurations['instanceConfig']['instance']['hostPathogen']])
    for hostA,pathogenA,pathogenA_name in instanceConfigurations['instanceConfig']['instance']['hostPathogen']:
        allSpecies.append({'taxid':pathogenA,'host':hostA,'name':pathogenA_name+" ("+pathogenA+")"})

    # Appending all networks transferred from orthologs to the species list
    orthologSpecies={} # Keeping transferred taxID:'species name'
    for orthologSp in os.listdir("website/static/website/orthologFiles"):
        orthologSp_taxid = orthologSp.split(":")[0]
        orthologSp_name = orthologSp.split(":")[1]
        orthologSp_name_gs = orthologSp_name.replace("_"," ")
        allSpecies.append({'taxid':orthologSp_taxid,'name':"*"+orthologSp_name_gs +" ("+orthologSp_taxid+")"})
        orthologSpecies[orthologSp_taxid]=orthologSp_name
    
    # Collecting evidence types
    orderedEvidenceTypes=[]
    for orderEv in instanceConfigurations['trainingConfig']['Network']['Evidence_order']:
        orderedEvidenceTypes.append((orderEv, orderEv))

    # Collecting all pathways
    allPathways=Pathway.objects.filter(Q(pathway_db="KEGG")).values("pathway_name","pathway_id").distinct()
    allPathways=[{'pathwayId':p['pathway_id'], 'name':p['pathway_name']+" (pathwayID:"+p['pathway_id']+")"} for p in allPathways]

    # Collecting all tissues per species
    allTissues=Tissue.objects.values("tissue", "species").distinct()
    allTissuesDict={}
    for t in allTissues:
        if t['species'] not in allTissuesDict:
            allTissuesDict[t['species']]=[]
        allTissuesDict[t['species']].append(t['tissue'])

    # Setting up parameters to be used in search form
    parameters={'exampleQuery':instanceConfigurations['webConfig']['exampleQuery'],
                'examplesPerSpecies': instanceConfigurations['webConfig']['examplesPerSpecies'],
                'exampleQueryMaxLink':instanceConfigurations['webConfig']['exampleQueryMaxLink'],
                'exampleQueryTOPAS':instanceConfigurations['webConfig']['exampleQueryTOPAS'],
                'exampleQuerySarsCov2':instanceConfigurations['webConfig']['exampleQuerySarsCov2'],
                'exampleSpecies':exampleSpecies,
                'showAdvanced':'false',
                'instanceSpecies':instanceConfigurations['instanceConfig']['instance']['species'],
                'instanceEvidences':instanceConfigurations['instanceConfig']['instance']['evidence'],
                'instanceAllGoldStandards':instanceConfigurations['instanceConfig']['instance']['goldStandard'],
                'instancePathogens':pathogens,
                'allSpecies':allSpecies,
                'allPathways':allPathways,
                'allTissues':allTissuesDict}

    # Prefilling the search form (from forms.py named Query) with information (the information that is not filled here, is filled in the javascript for the search page)
    query = Query(request.POST or None)
    query.fields['genomeSelect'].choices = allSpecies
    
    instanceSpeciesTaxIDnames = [(species.tax_id, species.species_name.split()[0][0] + "." + species.species_name.split(' (')[0].split()[-1] + " (" + species.tax_id+")") for species in instanceSpecies]
    query.fields['comparativeGenomes'].choices = instanceSpeciesTaxIDnames
    query.fields['evidenceSpecies'].choices = instanceSpeciesTaxIDnames
    categoryIDs=[(goldStandard, goldStandard+" links") for goldStandard in sorted(parameters['instanceAllGoldStandards'])]
    query.fields['categoryID'].choices = categoryIDs
    query.fields['categoryID'].initial = ['all']
    query.fields['evidenceDataTypes'].choices = orderedEvidenceTypes
    query.fields['constrainEvidence'].initial = ['all']
    return query,instanceConfigurations,parameters,orthologSpecies

def buildSearchQueryString(request,query,orthologSpecies):
    # Building a query string to use in the url from the data entered in the search form
    #     
    # First, checking for obvious errors, and in case one is found, add an error to the form, and dont redirect to the network view
    # Genes not available
    genes = request.POST['geneQuery'].strip()
    if genes=="":
        query.add_error('geneQuery', mark_safe('Please, input a gene identifier'))
        return None
    
    # Species not available
    allSpeciesNames = [species['name'] for species in query.fields['genomeSelect'].choices]
    if not request.POST['genomeSelect'] in allSpeciesNames:
        query.add_error('genomeSelect', mark_safe('Sorry, no network available for this species'))
        return None
    
    geneList = [item for sublist in [c.strip().split() for c in genes.split(',') if not c.isspace()] for item in sublist] # Split list on space AND comma
    originalGeneList=geneList.copy()

    # Getting query genome as taxid
    originalGenome=re.search(r'\((\d+)\)$', request.POST['genomeSelect'])
    originalGenome = originalGenome.group(1)
    genome=originalGenome # we need to keep track of the queried genome as well as the orthology transferred genome (in case of transferred search)
    
    # Trying to map the query proteins
    # Initialize dictionaries and lists to track mappings
    mappingDict = {}
    multiMapped = []
    missingMapping = []
    normallyMapped = []
    updatedGeneList = []

    # If the query species is orthology transferred, first collect the gene mappings for the query
    if originalGenome in orthologSpecies:
        geneList=[]
        successfullOriginalGeneList=[]
        file = open("website/static/website/orthologFiles/"+originalGenome+":"+orthologSpecies[originalGenome], "r")
        for line in file.readlines():
            if line.startswith("#"):
                genome=line.strip().replace("#","")
            elif line.split(",")[0].strip() in originalGeneList:
                geneList.append(line.split(",")[1].strip())
                successfullOriginalGeneList.append(line.split(",")[0].strip())
        ## Check what query genes have not been successfully mapped
        for g in originalGeneList:
            if g not in successfullOriginalGeneList:
                missingMapping.append(g)

    # if no ortholog mappings could be made for a transferred query, add an error and dont continue
    if originalGenome!=genome and len(geneList)<1: 
        query.add_error('geneQuery', 'Sorry, we could not find your identifier in the database.')
        return None
    
    # Iterate through each gene and query the database
    for gene in geneList:
        # Query for case-insensitive exact match
        results = IdMapping.objects.filter(
            Q(mapped_id__iexact=gene) &
            Q(protein__proteome__species__tax_id=genome)
        ).values_list("protein__uniprot_id", "mapped_id").distinct()

        # Process results
        if results:
            for uniprot_id, mapped_id in results:
                if mapped_id not in mappingDict:
                    mappingDict[mapped_id] = []
                mappingDict[mapped_id].append(uniprot_id)  # Store UniProt ID
                normallyMapped.append(mapped_id)
                # Replace the original gene with the mapped_id
                if mapped_id not in updatedGeneList:
                    updatedGeneList.append(mapped_id)
        else:
            # Store original gene if not found in the database
            missingMapping.append(gene)

    # Identify multi-mapped genes
    for mapped_id in mappingDict:
        if len(mappingDict[mapped_id]) > 1:
            multiMapped.append(mapped_id)

    # In case of only missing mappings, checking each to see if they can be matched with a more liberal check
    if len(normallyMapped)==0 and len(multiMapped)==0 and len(missingMapping)>0: 
        min_length=3
        canGetMappings=False
        for missing in missingMapping:
            if len(missing)>=min_length:
                caseInsensitiveMatches = IdMapping.objects.filter(Q(mapped_id__icontains=missing) & Q(protein__proteome__species__tax_id=genome))
                if caseInsensitiveMatches.exists():
                    canGetMappings=True # This is settled in the gene selection page before showing a network
        if canGetMappings==False: # If nothing can be mapped, return an error    
            query.add_error('geneQuery', mark_safe(
                'Sorry, we could not find your identifier in the database.<br>'
                'See <a href="/help/#Search">Help:Search</a> for supported identifiers and query tips.'
            ))
            return None

    # Starting to collect all parameters from the posted form and collecting them in a query string (for the url to the network)
    if 'depth' in request.POST:
        depth=request.POST['depth']
    else:
        depth="1"
    if len(geneList)<2 and depth=="0": # Passing error if depth of 0 and only one gene
        query.add_error('geneQuery', 'You can not use an expansion depth of zero when searching for a single protein. Please add more proteins or increase the Expansion depth.')
        return None
    
    if 'nodesPerStepOptions' in request.POST:
        nodesPerStepOptions=request.POST['nodesPerStepOptions']
    else:
        nodesPerStepOptions="15"

    if len(updatedGeneList)>0 and originalGenome==genome:
        querystring=",".join(updatedGeneList)+"&"+originalGenome+"&"+request.POST['confidenceThreshold']+"&"+request.POST['grgLLRThreshold']+"&"+depth+"&"+nodesPerStepOptions+"&"+request.POST['expansionAlgorithm']
    else:
        querystring=",".join(originalGeneList)+"&"+originalGenome+"&"+request.POST['confidenceThreshold']+"&"+request.POST['grgLLRThreshold']+"&"+depth+"&"+nodesPerStepOptions+"&"+request.POST['expansionAlgorithm']

    if 'prioritizeNeighbors' in request.POST:
        querystring+="&"+request.POST['prioritizeNeighbors']
    else:
        querystring+="&off"
    if 'comparativeGenomes' in request.POST:
        querystring+="&"+",".join(request.POST.getlist('comparativeGenomes'))
    else:
        querystring+="&''"
    if 'individualEvidenceOnly' in request.POST:
        querystring+="&"+request.POST['individualEvidenceOnly']
    else:
        querystring+="&off"
    if 'orthologsOnly' in request.POST:
        querystring+="&"+request.POST['orthologsOnly']
    else:
        querystring+="&off"

    if 'constrainGoldStandardLinks' in request.POST:
        querystring+="&"+request.POST['constrainGoldStandardLinks']
    else:
        querystring+="&"+request.POST['constrainGoldStandard']

    querystring+="&"+request.POST['constrainEvidence'] #constrainGoldStandard
    
    if 'evidenceSpecies' in request.POST:
        querystring+="&"+",".join(request.POST.getlist('evidenceSpecies'))
    else:
        querystring+="&''"
    if 'evidenceDataTypes' in request.POST:
        querystring+="&"+",".join(request.POST.getlist('evidenceDataTypes'))
    else:
        querystring+="&''"
    querystring+="&"+request.POST['constrainGene']
    if 'restrictPathway' in request.POST and len(request.POST['restrictPathway'])>0:
        querystring+="&"+request.POST['restrictPathway'].split("pathwayID:")[1][:-1]
    else:
        querystring+="&''"
    if 'restrictTissue' in request.POST and len(request.POST['restrictTissue'])>0:
        querystring+="&"+request.POST['restrictTissue']
    else:
        querystring+="&''"
    querystring+="&"+request.POST['showAdvanced']

    if len(missingMapping)>0:
        querystring+="&"+",".join(missingMapping)
    else:
        querystring+="&''"
    if len(multiMapped)>0:
        querystring+="&"+",".join(multiMapped)
    else:
        querystring+="&''"
    querystring+="/"
    return querystring


def buildSearchQueryStringAPI(genes,genomeSelect,orthologSpecies):
    # Building a query string to use in the url from the data entered in the search form
    geneList = [item for sublist in [c.strip().split() for c in genes.split(',') if not c.isspace()] for item in sublist] # Split list on space AND comma
    originalGeneList=geneList.copy()

    # Getting query genome as taxid
    originalGenome=genomeSelect
    genome=originalGenome # we need to keep track of the queried genome as well as the orthology transferred genome (in case of transferred search)

    # Trying to map the query proteins
    # Initialize dictionaries and lists to track mappings
    mappingDict = {}
    multiMapped = []
    missingMapping = []
    normallyMapped = []
    updatedGeneList = []

    # If the query species is orthology transferred, first collect the gene mappings for the query
    if originalGenome in orthologSpecies: 
        geneList=[]
        successfullOriginalGeneList=[]
        file = open("website/static/website/orthologFiles/"+originalGenome+":"+orthologSpecies[originalGenome], "r")
        for line in file.readlines():
            if line.startswith("#"):
                genome=line.strip().replace("#","")
            elif line.split(",")[0].strip() in originalGeneList:
                geneList.append(line.split(",")[1].strip())
                successfullOriginalGeneList.append(line.split(",")[0].strip())
        ## Check what query genes have not been successfully mapped
        for g in originalGeneList:
            if g not in successfullOriginalGeneList:
                missingMapping.append(g)
    # if no ortholog mappings could be made for a transferred query, add an error and dont continue
    if originalGenome!=genome and len(geneList)<1: 
        return [],originalGeneList,[]

    # Iterate through each gene and query the database
    for gene in geneList:
        # Query for case-insensitive exact match
        results = IdMapping.objects.filter(
            Q(mapped_id__iexact=gene) &
            Q(protein__proteome__species__tax_id=genome)
        ).values_list("protein__uniprot_id", "mapped_id").distinct()

        # Process results
        if results:
            for uniprot_id, mapped_id in results:
                if mapped_id not in mappingDict:
                    mappingDict[mapped_id] = []
                mappingDict[mapped_id].append(uniprot_id)  # Store UniProt ID
                normallyMapped.append(mapped_id)
                # Replace the original gene with the mapped_id
                if mapped_id not in updatedGeneList:
                    updatedGeneList.append(mapped_id)
        else:
            # Store original gene if not found in the database
            missingMapping.append(gene)
    
    # Identify multi-mapped genes
    for mapped_id in mappingDict:
        if len(mappingDict[mapped_id]) > 1:
            multiMapped.append(mapped_id)

    # In case of only missing mappings, checking each to see if they can be matched with a more liberal check
    if len(normallyMapped)==0 and len(multiMapped)==0 and len(missingMapping)>0: 
        min_length=3
        canGetMappings=False
        for missing in missingMapping:
            if len(missing)>=min_length:
                caseInsensitiveMatches = IdMapping.objects.filter(Q(mapped_id__icontains=missing) & Q(protein__proteome__species__tax_id=genome))
                if caseInsensitiveMatches.exists():
                    canGetMappings=True # This is settled in the gene selection page before showing a network
        if canGetMappings==False: # If nothing can be mapped, return an error    
            return [],originalGeneList,[]

    if len(updatedGeneList)>0 and originalGenome==genome:
        return updatedGeneList, missingMapping, multiMapped
    else:
        return originalGeneList, missingMapping, multiMapped


def fillSearchForm(query, parameters, geneQuery, genomeSelect, showAdvanced, confidenceThreshold, grgLLRThreshold, depth, nodesPerStepOptions, expansionAlgorithm, prioritizeNeighbors, comparativeGenomes, individualEvidenceOnly, orthologsOnly, categoryID, constrainEvidence, evidenceSpecies, evidenceDataTypes, restriction, restrictPathway,restrictTissue ):
    # Filling the search form in modify search with all information that was used for the query
    query.fields['geneQuery'].initial = geneQuery

    for s in parameters['allSpecies']: # Getting the species name (passed here as taxid but name needs to be in form)
        if s['taxid']==genomeSelect:
            searchedSpeciesname=s['name']

    query.fields['genomeSelect'].initial = searchedSpeciesname
    query.fields['showAdvanced'].initial = showAdvanced
    query.fields['confidenceThreshold'].initial = confidenceThreshold
    query.fields['grgLLRThreshold'].initial = grgLLRThreshold
    query.fields['depth'].initial = depth
    query.fields['nodesPerStepOptions'].initial = nodesPerStepOptions
    query.fields['expansionAlgorithm'].initial = expansionAlgorithm
    query.fields['prioritizeNeighbors'].initial = prioritizeNeighbors if prioritizeNeighbors!="off" else False
    query.fields['comparativeGenomes'].initial = comparativeGenomes.split(",")
    query.fields['individualEvidenceOnly'].initial = individualEvidenceOnly if individualEvidenceOnly!="off" else False
    query.fields['orthologsOnly'].initial = orthologsOnly if orthologsOnly!="off" else False
    query.fields['categoryID'].initial = [categoryID]
    query.fields['constrainEvidence'].initial = [constrainEvidence]
    query.fields['evidenceSpecies'].initial = evidenceSpecies.split(",")
    query.fields['evidenceDataTypes'].initial = evidenceDataTypes.split(",")
    query.fields['restriction'].initial = [restriction]
    searchedPathway=""
    if restrictPathway!="":
        for s in parameters['allPathways']: # Getting all pathway names (passed here as pathwayIDs but full name needs to be in form)
            if s['pathwayId']==restrictPathway:
                searchedPathway=s['name']
    query.fields['restrictPathway'].initial = searchedPathway
    query.fields['restrictTissue'].initial = "" if len(restrictTissue)<3 else restrictTissue

    # Now update and organize necessary parameters to work in the network collection
    paramsForNetwork={'confidenceThreshold':float(query['confidenceThreshold'].value()),
                      'grgLLRThreshold':float(query['grgLLRThreshold'].value()),
                      'maxNodesPerStep':int(query['nodesPerStepOptions'].value()),
                      'depth':query['depth'].value(),
                      'expansionAlgorithm':query['expansionAlgorithm'].value(),
                      'hypergeomCutoff':0.05,
                      'prioritizeCommonNeighbors':query['prioritizeNeighbors'].value(),
                      'restrictGoldStandard':query['categoryID'].value()[0],
                      'goldstandardColumn':"max_ppv",
                      'constrainEv': query['constrainEvidence'].value()[0],
                      'constrainEvTypesSelected': query['evidenceDataTypes'].value(),
                      'constrainSpSelected':query['evidenceSpecies'].value(),
                      'restrictPathway':restrictPathway,
                      'restrictTissue': query['restrictTissue'].value(),
                      'comparativeSpecies':query['comparativeGenomes'].value(),
                      'individualEvidenceOnly':query['individualEvidenceOnly'].value(),
                      'orthologsOnly':query['orthologsOnly'].value()}
    if paramsForNetwork['expansionAlgorithm']=="maxlink":
        paramsForNetwork['hypergeomCutoff']=float(paramsForNetwork['depth'])
        paramsForNetwork['depth']=2
    if paramsForNetwork['expansionAlgorithm']=="topas":
        paramsForNetwork['expansion_steps']=float(paramsForNetwork['depth'])
    if paramsForNetwork['restrictGoldStandard']!="all":
        paramsForNetwork['restrictGoldStandard']=paramsForNetwork['restrictGoldStandard'].split("-")[0]
        paramsForNetwork['restrictGoldStandard']
        paramsForNetwork['goldstandardColumn']=paramsForNetwork['restrictGoldStandard']+"_ppv"
    if paramsForNetwork['individualEvidenceOnly']=="on" and paramsForNetwork['comparativeSpecies']!=["''"]:
        paramsForNetwork['constrainEv']="species"
        paramsForNetwork['constrainSpSelected']=[genomeSelect]
    if paramsForNetwork['constrainEv'] == 'on':
        constrainEv=''
        if paramsForNetwork['constrainEvTypesSelected']!=["''"]:
            constrainEv+="evidence"
        if paramsForNetwork['constrainSpSelected']!=["''"]:
            constrainEv+="species"
        paramsForNetwork['constrainEv']=constrainEv

    if "species" in paramsForNetwork['constrainEv'] and paramsForNetwork['comparativeSpecies']==["''"] and len(paramsForNetwork['constrainSpSelected'])==len(parameters['instanceSpecies']):
        ## If all species are left selected, no recomputation of network is done
        ## The URL will still contain the species names
        paramsForNetwork['constrainEv']='off'
    if "evidence" in paramsForNetwork['constrainEv'] and paramsForNetwork['comparativeSpecies']==["''"] and len(paramsForNetwork['constrainEvTypesSelected'])==len(parameters['instanceEvidences']):
        ## If all evidences are left selected, no recomputation of network is done
        ## The URL will still contain the evidence names
        paramsForNetwork['constrainEv']='off'

    paramsForNetwork['depth']=int(paramsForNetwork['depth'])
    
    return query, paramsForNetwork

def collectQueryParameters(parameters,paramsForNetwork,genomeSelect,taxIdToName,originalGenome):

    ## To display in network information
    if paramsForNetwork['expansionAlgorithm']=="group":
        expansionAlgorithmLabel = 'Genes as group'
    elif paramsForNetwork['expansionAlgorithm']=="local":
        expansionAlgorithmLabel = 'Each gene independently'
    elif paramsForNetwork['expansionAlgorithm']=="maxlink":
        expansionAlgorithmLabel = 'MaxLink'
    elif paramsForNetwork['expansionAlgorithm']=="topas":
        expansionAlgorithmLabel = 'TOPAS'
    else: 
        expansionAlgorithmLabel = ''
    parameters['expansionAlgorithmLabel'] = expansionAlgorithmLabel

    ## To configure Modify Search
    parameters['genomeSelect']=taxIdToName[originalGenome]['name']
    if genomeSelect in parameters['instancePathogens']:
        parameters['hostPathogenNetwork']=True
    else: 
        parameters['hostPathogenNetwork']=False

    parameters['expansionAlgorithm'] = paramsForNetwork['expansionAlgorithm']
    parameters['confidenceThreshold'] = paramsForNetwork['confidenceThreshold']
    parameters['grgLLRThreshold'] = paramsForNetwork['grgLLRThreshold']
    parameters['nodesPerStepOptions'] = paramsForNetwork['maxNodesPerStep']
    parameters['maxNodesPerStep'] = paramsForNetwork['maxNodesPerStep']
    parameters['depth'] = paramsForNetwork['depth']
    parameters['hypergeomCutoff'] = paramsForNetwork['hypergeomCutoff']
    parameters['prioritizeCommonNeighbors'] = paramsForNetwork['prioritizeCommonNeighbors']
    parameters['restrictGoldStandard'] = paramsForNetwork['restrictGoldStandard']
    parameters['goldstandardColumn'] = paramsForNetwork['goldstandardColumn']
    parameters['constrainEv'] = paramsForNetwork['constrainEv']
    parameters['constrainEvTypesSelected'] = paramsForNetwork['constrainEvTypesSelected']
    parameters['constrainSpSelected'] = paramsForNetwork['constrainSpSelected']
    parameters['restrictGoldStandard'] = paramsForNetwork['restrictGoldStandard']
    parameters['restrictPathway'] = paramsForNetwork['restrictPathway']
    parameters['restrictTissue'] = paramsForNetwork['restrictTissue']
    parameters['comparativeSpecies'] = paramsForNetwork['comparativeSpecies']
    parameters['individualEvidenceOnly'] = paramsForNetwork['individualEvidenceOnly']
    parameters['orthologsOnly'] = paramsForNetwork['orthologsOnly']

    parameters['restrictPathway'] = paramsForNetwork['restrictPathway']
    parameters['restrictTissue'] = paramsForNetwork['restrictTissue']
    parameters['comparativeSpecies'] = paramsForNetwork['comparativeSpecies']
    parameters['individualEvidenceOnly'] = paramsForNetwork['individualEvidenceOnly']

    if parameters['expansionAlgorithm']=="maxlink":
        parameters['depth']=int(parameters['depth'])

    if parameters['expansionAlgorithm']=="topas":
        parameters['expansion_steps']=float(parameters['depth'])

    return parameters