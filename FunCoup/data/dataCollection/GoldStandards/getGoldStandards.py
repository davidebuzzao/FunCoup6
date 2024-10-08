import requests
import zipfile
import os
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed
from .PPI import *
from .Complex import *
from .MetabolicAndSignaling import *
from .Operon import *
from .Regulatory import *

def downloadFile(urlToDownload, filename):
    # Downloading file from url
    url = urlToDownload
    myfile = requests.get(url)
    if myfile.status_code==200:
        open(filename, "wb").write(myfile.content)
    else:
        print("Cant download file from: "+url)

def unzipIRefIndexFile(filename,path):
    # Special unzip for iRefIndex, as the file is broken (v.22-08), where the compressed file is named "-"
    print("unzipping file: "+filename+".zip")
    with zipfile.ZipFile(path+filename+".zip", "r") as zip_ref:
       zip_ref.extract("-", path=path)
    os.rename(path+"-",path+filename)

def unzipFile(filename,path):
    # Unzipping file to path
    print("unzipping file: "+filename+".zip")
    with zipfile.ZipFile(path+filename+".zip", "r") as zip_ref:
       zip_ref.extractall( path=path)

def removeTempFile(filename):
    # Removing file named filename
    if os.path.isfile(filename):
        os.remove(filename)

def getComplex(species,goldStandardConfig,proteomeConfig):
    # Getting complexes from the iRefIndex file OR collecting existing gold standards from DB
    goldStandardList=[]
    speciesToGet=[]
    for s in species:
        goldStandardInDB = GoldStandard.objects.filter(Q(version=goldStandardConfig['version']) &  Q(type="Complex") &  Q(species__tax_id=s))
        if goldStandardInDB.exists(): # Using existing gild standard if available in the DB
            goldStandardList.append(goldStandardInDB[0])
            print("Gold standard Exists: Complex v."+goldStandardConfig['version']+" species: "+s)
        else: # otherwise adding species to be collected
            speciesToGet.append(str(s))

    if len(speciesToGet)>0:
        corumSpeciesMap={'Human':'9606','Mouse':'10090','Rat':'10116'} # Defining corum species map with species names as defined in the corum file
        allPairsPerSpecies={} # taxID:allPairs

        # Downloading and extracting irefindex file if it does not already exist
        #if os.path.exists("data/tmp/"+goldStandardConfig['urlIrefIndex'].split("/")[-1].split(".zip")[0])==False: # Might already have been downloaded.
        #    downloadFile(goldStandardConfig['urlIrefIndex'], "data/tmp/"+goldStandardConfig['urlIrefIndex'].split("/")[-1])
        #    unzipIRefIndexFile(goldStandardConfig['urlIrefIndex'].split("/")[-1].split(".zip")[0],"data/tmp/")
        # Collecting complexes from irefindex file
        allPairsPerSpecies = extractComplexFomiRefIndex("data/tmp/"+goldStandardConfig['urlIrefIndex'].split("/")[-1].split(".zip")[0], speciesToGet,allPairsPerSpecies)

        # Downloading and extracting corum file (for human, rat and mouse) if it does not already exist
        if os.path.exists("data/tmp/"+goldStandardConfig['urlCorum'].split("/")[-1].split(".zip")[0])==False:
            downloadFile(goldStandardConfig['urlCorum'], "data/tmp/"+goldStandardConfig['urlCorum'].split("/")[-1])
            unzipFile(goldStandardConfig['urlCorum'].split("/")[-1].split(".zip")[0],"data/tmp/")
        # Collecting complexes from corum
        allPairsPerSpecies = extractComplexesFromCorum("data/tmp/"+goldStandardConfig['urlCorum'].split("/")[-1].split(".zip")[0], speciesToGet,allPairsPerSpecies,corumSpeciesMap)

        # Downloading and extracting complexPortal file if it does not already exist
        speciesForComplexPortal=[]
        for s in speciesToGet:
            if s not in corumSpeciesMap.values(): # Not getting for species that have corum
                myfile = requests.get(goldStandardConfig['urlComplexPortal']+s+".tsv")
                if myfile.status_code==200: # If species exist in complexPortal (url returns a valid file)
                    speciesForComplexPortal.append(s)
                    open("data/tmp/"+s+".tsv", "wb").write(myfile.content)
        # Collecting complexes from complexPortal
        allPairsPerSpecies = extractComplexesFromComplexPortal(speciesForComplexPortal,allPairsPerSpecies)

        # Writing all complex pairs to DB
        #print("# TaxID,LinksWrittenToDB,UnMappableProteins")
        cpus = os.cpu_count()
        parallel = Parallel(n_jobs=cpus)
        results = parallel(delayed(writeComplexesToDB)(sp, allPairsPerSpecies[sp], goldStandardConfig['version'],proteomeConfig) for sp in allPairsPerSpecies)

        goldStandardList.extend(results)
        # Deleting downloaded files
        # Remove .zip file, keep unzipped file, since it may be neede again.
        #removeTempFile("data/tmp/"+goldStandardConfig['urlIrefIndex'].split("/")[-1])
        # removeTempFile("data/tmp/"+goldStandardConfig['urlIrefIndex'].split("/")[-1].split(".zip")[0])
        #removeTempFile("data/tmp/"+goldStandardConfig['urlCorum'].split("/")[-1])
        #removeTempFile("data/tmp/"+goldStandardConfig['urlCorum'].split("/")[-1].split(".zip")[0])
        #for s in speciesForComplexPortal:
        #    removeTempFile("data/tmp/"+s+".tsv")
    return goldStandardList

def getMetabolic(species,goldStandardConfig,proteomeConfig):
    # Getting metabolic gold standard from DB OR downloading and extracting links from KEGG
    goldStandardList=[]
    speciesToGet=[]
    for s in species:
        goldStandardInDB = GoldStandard.objects.filter(Q(version=goldStandardConfig['version']) &  Q(type="Metabolic") &  Q(species__tax_id=s))
        if goldStandardInDB.exists(): # If gold standard already exist in the database, using that for this instance
            goldStandardList.append(goldStandardInDB[0])
            print("Gold standard Exists: Metabolic v."+goldStandardConfig['version']+" species: "+s)
        else: # otherwise adding species in species to fetch links for
            speciesToGet.append(str(s))

    if len(speciesToGet)>0:
        # Getting the metabolic pathways to use
        filename="data/tmp/allKeggPathways"
        downloadFile(goldStandardConfig['urlPathways'], filename)
        pathwaysToGet=updateKEGGPathwaysToUse(filename, goldStandardConfig['existingPathways'], "Metabolic")
        # Collectig links for the metabolic pathways
        cpus = os.cpu_count()
        parallel = Parallel(n_jobs=cpus)
        results = parallel(delayed(getPathwayPairs)(pathwaysToGet, goldStandardConfig, sp, "Metabolic",proteomeConfig) for sp in speciesToGet)

        goldStandardList.extend(results)
        removeTempFile(filename)

    return goldStandardList

def getOperon(species,goldStandardConfig,proteomeConfig):
    # Getting operons to use from the DB OR collecting new operon links from operonDB
    goldStandardList=[]
    speciesToGet=[]
    operonDbMap={"511145":"83333","224308":"224308","6239":"6239"} #OperonDB only contains these species,Note that ecoli needs to be mapped to a different ecoli taxID (for a different strain)

    for s in species:
        if s in operonDbMap.values():
            goldStandardInDB = GoldStandard.objects.filter(Q(version=goldStandardConfig['version']) &  Q(type="Operon") &  Q(species__tax_id=s))
            if goldStandardInDB.exists(): # Using existing gold standard if it already exist in teh DB
                goldStandardList.append(goldStandardInDB[0])
                print("Gold standard Exists: Operon v."+goldStandardConfig['version']+" species: "+s)
            else: # Otherwise adding species to fetch operon links for
                speciesToGet.append(str(s))

    if len(speciesToGet)>0:
        # Downloading operon file from operonDB
        filename="data/tmp/knownOperons"
        downloadFile(goldStandardConfig['url'], filename)
        # Collecting operon links from file
        goldStandardList=getOperonsFromFile(filename, goldStandardList, goldStandardConfig, speciesToGet, operonDbMap,proteomeConfig)
        removeTempFile(filename)

    return goldStandardList

def getPPI(species,goldStandardConfig,goldStandardList,proteomeConfig):
    # Getting PPI links from the irefindex file OR picks up the existing gold standards from the DB
    speciesToGet=[]
    for s in species:
        goldStandardInDB = GoldStandard.objects.filter(Q(version=goldStandardConfig['version']) &  Q(type="PPI") &  Q(species__tax_id=s))
        if goldStandardInDB.exists(): # If the gold standard already exist in the DB, using existing one for the instance
            goldStandardList.append(goldStandardInDB[0])
            print("Gold standard Exists: PPI v."+goldStandardConfig['version']+" species: "+s)
        else: # otherwise adding species to fetch new links for
            speciesToGet.append(str(s))

    if len(speciesToGet)>0:
        # Getting irefindex file if it does not already exist
        #if os.path.exists("data/tmp/"+goldStandardConfig['url'].split("/")[-1].split(".zip")[0])==False: # Might already have been downloaded for complex.
        #    downloadFile(goldStandardConfig['url'], "data/tmp/"+goldStandardConfig['url'].split("/")[-1])
        #    unzipIRefIndexFile(goldStandardConfig['url'].split("/")[-1].split(".zip")[0],"data/tmp/")
        # Extracting links from the file
        goldStandardList=extractFomAllFile("data/tmp/"+goldStandardConfig['url'].split("/")[-1].split(".zip")[0], speciesToGet,goldStandardConfig['version'],goldStandardList,proteomeConfig)
        # Deleting downloaded files
        #if os.path.exists("data/tmp/"+goldStandardConfig['url'].split("/")[-1]):# Remove .zip file, keep unzipped file, since it may be neede again.
            #removeTempFile("data/tmp/"+goldStandardConfig['url'].split("/")[-1])
            # removeTempFile("data/tmp/"+goldStandardConfig['url'].split("/")[-1].split(".zip")[0])
    return goldStandardList

def getRegulatory(species,goldStandardConfig,proteomeConfig):
    goldStandardList = []
    species_with_data = ['9606','10090','83333','559292']
    species = list(set(species) & set(species_with_data))
    if len(species)==0: return

    speciesToGet=[]
    for s in species:
        if s in ['9606','10090']: version = ','.join([goldStandardConfig['RegNetwork']['version'],goldStandardConfig['TRRUST']['version']])
        elif s in ['83333']: version = goldStandardConfig['RegulonDB']['version']
        elif s in ['559292']: version = goldStandardConfig['Yeastract']['version']
        goldStandardInDB = GoldStandard.objects.filter(Q(version=version) & Q(type="Regulatory") & Q(species__tax_id=s))
        if goldStandardInDB.exists(): # If the gold standard already exist in the DB, using existing one for the instance
            goldStandardList.append(goldStandardInDB[0])
            print("Gold standard Exists: Regulatory v."+goldStandardConfig['version']+" species: "+s)
        else: # Otherwise adding the species to fetch links
            speciesToGet.append(str(s))

    if len(speciesToGet)>0:
        # Collectig links for the singaling pathways
        cpus = os.cpu_count()
        parallel = Parallel(n_jobs=cpus)
        results = parallel(delayed(getRegulatoryPairs)(sp,goldStandardConfig,proteomeConfig) for sp in speciesToGet)
        goldStandardList.extend(results)

    return goldStandardList

def getSignaling(species,goldStandardConfig,proteomeConfig):
    # Getting signaling gold standard from DB OR downloading and extracting links from KEGG
    goldStandardList=[]
    speciesToGet=[]
    for s in species:
        goldStandardInDB = GoldStandard.objects.filter(Q(version=goldStandardConfig['version']) &  Q(type="Signaling") &  Q(species__tax_id=s))
        if goldStandardInDB.exists(): # If the gold standard already exist in the DB, using existing one for the instance
            goldStandardList.append(goldStandardInDB[0])
            print("Gold standard Exists: Signaling v."+goldStandardConfig['version']+" species: "+s)
        else: # Otherwise adding the species to fetch links
            speciesToGet.append(str(s))

    if len(speciesToGet)>0:
        # Getting the signaling pathways to use
        filename="data/tmp/allKeggPathways"
        downloadFile(goldStandardConfig['urlPathways'], filename)
        pathwaysToGet=updateKEGGPathwaysToUse(filename, goldStandardConfig['existingPathways'], "Signaling")
        # Collectig links for the singaling pathways
        cpus = os.cpu_count()
        parallel = Parallel(n_jobs=cpus)
        results = parallel(delayed(getPathwayPairs)(pathwaysToGet, goldStandardConfig, sp, "Signaling",proteomeConfig) for sp in speciesToGet)

        goldStandardList.extend(results)
        removeTempFile(filename)

    return goldStandardList

def getGoldStandards(instanceConfig,goldStandardConfig,proteomeConfig):
    # Getting gold standards defined in the generate instance goldStandardConfig file.
    # Returning the gold standards to be used for this instance

    goldStandardsForInstance=[]
    usingPPI=False # PPI needs to be added last, as it is intersected with other gold standards!
    for goldStandardType in instanceConfig['instance']['goldStandard']:
        if goldStandardType=="PPI":
            usingPPI=True
        elif goldStandardType=="Complex":
            print("\nGetting Complex")
            goldStandardsForInstance.extend(getComplex(instanceConfig['instance']['species'],goldStandardConfig['Complex'],proteomeConfig))
        elif goldStandardType=="Metabolic":
            print("\nGetting Metabolic")
            goldStandardsForInstance.extend(getMetabolic(instanceConfig['instance']['species'],goldStandardConfig['Metabolic'],proteomeConfig))
        elif goldStandardType=="Operon":
            print("\nGetting Operon")
            goldStandardsForInstance.extend(getOperon(instanceConfig['instance']['species'],goldStandardConfig['Operon'],proteomeConfig))
        elif goldStandardType=="Regulatory":
            print("\nGetting Regulatory")
            goldStandardsForInstance.extend(getRegulatory(instanceConfig['instance']['species'],goldStandardConfig['Regulatory'],proteomeConfig))
        elif goldStandardType=="Signaling":
            print("\nGetting Signaling")
            goldStandardsForInstance.extend(getSignaling(instanceConfig['instance']['species'],goldStandardConfig['Signaling'],proteomeConfig))
        else:
            print("Did not recognize gold standard: "+goldStandardType+". \nPlease choose from PPI, Complex, Signaling, Metabolic or Operon")
            exit()

    # Getting PPI last, to be able to intersect it with other goldstandards
    if usingPPI:
        print("\nGetting PPI")
        goldStandardsForInstance = getPPI(instanceConfig['instance']['species'],goldStandardConfig['PPI'],goldStandardsForInstance,proteomeConfig)

    return goldStandardsForInstance
