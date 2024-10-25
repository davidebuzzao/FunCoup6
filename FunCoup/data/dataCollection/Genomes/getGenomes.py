import requests
import os
import gzip
import shutil
import yaml
from yaml.loader import SafeLoader
from contextlib import closing
from io import StringIO
from joblib import Parallel, delayed
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed

def downloadFile(urlToDownload, filename):
    # Download file from url. Write contents to filename.
    if not os.path.exists("data/tmp/"+filename):
        print("Downloading: "+urlToDownload)
        url = urlToDownload
        myfile = requests.get(url)
        open("data/tmp/"+filename, "wb").write(myfile.content)

def unzipFile(filename):
    # Unzipping file named filename.gz to filename
    print("Unzipping: "+ filename)
    with gzip.open("data/tmp/"+filename+".gz", 'rb') as f_in:
        with open("data/tmp/"+filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def readREADME(sp):
    # Getting species information from readme
    # Saving species to the database ( with taxID, species name and common name)
    # Returning dictionary of taxID -> ebiIdentifier|domain|speciesObjectFromDB
    ebiIdentifiers={}
    with open("data/tmp/"+"./README", 'r') as in_file:
        for line in in_file:
            if line.startswith("UP0"):
                l=line.split("\t")
                if len(l)> 4 and (l[1] in sp):
                    ebiIdentifiers[l[1]]=[l[0], l[3]]
                    if len(l)<8: # This is necessary for the broken README in v.2022_02, which is not always properly tab separated
                        spName=l[5].strip()
                    else:
                        spName=l[7].strip()
                    if "(" in spName:
                        if "strain" in spName:
                            if spName.count("(")>1:
                                commonName=spName.split(") ")[1].replace("(","").replace(")","").strip()
                                spName=spName.split(") ")[0]+")".strip()
                            else: # If there is no common name, species name is set
                                spName=spName
                                commonName=spName
                        else:
                            commonName=spName.split("(")[1].replace(")","").strip()
                            spName=spName.split("(")[0].strip()
                    else: # If there is no common name, species name is set
                        commonName=spName
                    taxID=l[1].strip()
                    # Adding species to DB if it does not already exist
                    try:
                        species = Species.objects.get(tax_id = taxID)
                    except Species.DoesNotExist:
                        print("Adding species to db: "+taxID+" species name: "+spName+" common name: "+commonName)
                        species = Species(tax_id=taxID, species_name=spName, common_name=commonName)
                        species.save()
                    ebiIdentifiers[l[1]].append(species)
    # Checking is all requested species are available at EBI. If not, they are listed, and require manual download.
    if (len(ebiIdentifiers)==len(sp)):
        print("all species can be fetched from UniProt")
    else:
        print("############### Species missing from EBI, requires manual handling: ##################")
        for s in sp:
            if s not in ebiIdentifiers:
                print(s)
    return ebiIdentifiers

def removeTempFiles(path):
    if os.path.exists(path):
        os.remove(path)


def getProteomesAndMappings(s, ebiIdentifiers, proteomeConfig,mappingsOfInterest ):
    # Creating proteome in DB
    print("\nAdding proteome to DB for "+s+" version: "+proteomeConfig['genome']['version'])
    proteomeInDB = Proteome.objects.filter(version=proteomeConfig['genome']['version'], species=ebiIdentifiers[s][2])
    # proteomeInDB = Proteome(version=proteomeConfig['genome']['version'], species=ebiIdentifiers[s][2])
    if not proteomeInDB.exists():
        proteomeInDB = Proteome(version=proteomeConfig['genome']['version'], species=ebiIdentifiers[s][2])
        proteomeInDB.save()
        # Downloading and extracting  cannonical fasta. To avoid getting redundant info/proteins not included in the reference proteome
        downloadFile(proteomeConfig['genome']['urlUniProtDownload']+ebiIdentifiers[s][1].capitalize()+"/"+ebiIdentifiers[s][0]+"/"+ebiIdentifiers[s][0]+"_"+s+".fasta.gz", s+".fasta.gz")
        proteinsToUse={}
        fastaFile=gzip.open("data/tmp/"+s+".fasta.gz","r")
        for line in fastaFile:
            if line.startswith(b">"):
                line=line.decode()
                protId=line.split("|")[1].strip()
                description=line.split("|")[2].split("OS=")[0].strip()
                proteinsToUse[protId]=description # Storing protein name and desription to put in DB later

        # Downloading idmapping file to collect mappings per protein
        downloadFile(proteomeConfig['genome']['urlUniProtDownload']+ebiIdentifiers[s][1].capitalize()+"/"+ebiIdentifiers[s][0]+"/"+ebiIdentifiers[s][0]+"_"+s+".idmapping.gz", s+".idmapping.gz")
        idmappings={}
        altUniprotIDsByEnsembl={} # Keeping all uniprotIDs linked to an ensemblID in a dict, to be able to extract alternative uniprotIDs
        idmappingFile=gzip.open("data/tmp/"+s+".idmapping.gz","r")
        for line in idmappingFile:
            uniprotID, type, mapId = line.decode().strip().split("\t")
            if "-" in uniprotID:
                uniprotID=uniprotID.split("-")[0]
            if type in mappingsOfInterest:
                if mappingsOfInterest[type]=="Ensembl":
                    mapId=mapId.split(".")[0]
                    if mapId in altUniprotIDsByEnsembl:
                        altUniprotIDsByEnsembl[mapId].append(uniprotID)
                    else:
                        altUniprotIDsByEnsembl[mapId]=[uniprotID]
                if mappingsOfInterest[type]=="Ensembl_Protein":
                    mapId=mapId.split(".")[0]
                if uniprotID in idmappings:
                    idmappings[uniprotID].append(mapId+"||"+type)
                else:
                    idmappings[uniprotID]=[mapId+"||"+type]

        # Creating all proteins and mappings for the species
        countAddedProteins=0
        countAddedIdMappings=0
        proteinsToSaveInDB=[]
        mappingsToSaveInDB=[]
        for p in proteinsToUse:
            # Creating protein to be put in DB
            prot = Protein(uniprot_id=p, description=proteinsToUse[p],proteome=proteomeInDB)
            proteinsToSaveInDB.append(prot)
            countAddedProteins+=1
            idMapping = IdMapping(protein=prot, type="UniProtID",mapped_id=p) # Adding mapping to itself for faster lookup
            mappingsToSaveInDB.append(idMapping)
            ensembl=""
            if p not in idmappings: print('WARNING %s not in idmappings!!!' %p)
            else:
                for m in idmappings[p]:
                    # Creating ID mapping to be put in DB
                    mappingID, mappingType=m.split("||")
                    idMapping = IdMapping(protein=prot, type=mappingsOfInterest[mappingType],mapped_id=mappingID)
                    mappingsToSaveInDB.append(idMapping)
                    countAddedIdMappings+=1
                    if mappingsOfInterest[mappingType]=="Ensembl":
                        ensembl=mappingID
                # If the same ensembl is mapped to different uniprotIDs,
                # they can be viewed as alternative uniprot IDs, and are stored with mapping type alt_uniprot_id
                if ensembl!="" and ensembl in altUniprotIDsByEnsembl:
                    for m in altUniprotIDsByEnsembl[ensembl]:
                        if m!=p and m not in proteinsToUse:
                            numEnsemblMappings=0 # Making sure not to add alternative ids having multiple ensembl mappings, as those can be to ambiguous
                            for altMapping in idmappings[m]:
                                if mappingsOfInterest[altMapping.split("||")[1]]=="Ensembl":
                                    numEnsemblMappings+=1
                            if numEnsemblMappings<2:
                                idMapping = IdMapping(protein=prot, type="alt_uniprot_id",mapped_id=m)
                                mappingsToSaveInDB.append(idMapping)
                                countAddedIdMappings+=1
                                for altMapping in idmappings[m]:
                                    mappingID, mappingType=altMapping.split("||")
                                    idMapping = IdMapping(protein=prot, type="alt_"+mappingsOfInterest[mappingType],mapped_id=mappingID)
                                    mappingsToSaveInDB.append(idMapping)
                                    countAddedIdMappings+=1

        print("Got "+str(countAddedProteins)+" proteins and "+str(countAddedIdMappings)+" mappings for: "+s)
        # Writing all proteins and mappings to DB
        Protein.objects.bulk_create(proteinsToSaveInDB)
        IdMapping.objects.bulk_create(mappingsToSaveInDB)
    
    return proteomeInDB

def updateMappings(s, ebiIdentifiers, proteomeConfig,mappingsOfInterest):
    existingProteome=Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']) & Q(species=ebiIdentifiers[s][2]))[0]
    proteinsForSpecies=Protein.objects.filter(Q(proteome=existingProteome))
    # Downloading idmapping file to collect mappings per protein
    downloadFile(proteomeConfig['genome']['urlUniProtDownload']+ebiIdentifiers[s][1].capitalize()+"/"+ebiIdentifiers[s][0]+"/"+ebiIdentifiers[s][0]+"_"+s+".idmapping.gz", s+".idmapping.gz")
    idmappings={}
    altUniprotIDsByEnsembl={} # Keeping all uniprotIDs linked to an ensemblID in a dict, to be able to extract alternative uniprotIDs
    idmappingFile=gzip.open("data/tmp/"+s+".idmapping.gz","r")
    for line in idmappingFile:
        uniprotID, type, mapId = line.decode().strip().split("\t")
        if "-" in uniprotID:
            uniprotID=uniprotID.split("-")[0]
        if type in mappingsOfInterest:
            if mappingsOfInterest[type]=="Ensembl":
                mapId=mapId.split(".")[0]
                if mapId in altUniprotIDsByEnsembl:
                    altUniprotIDsByEnsembl[mapId].append(uniprotID)
                else:
                    altUniprotIDsByEnsembl[mapId]=[uniprotID]
            if mappingsOfInterest[type]=="Ensembl_Protein":
                mapId=mapId.split(".")[0]
            if uniprotID in idmappings:
                idmappings[uniprotID].append(mapId+"||"+type)
            else:
                idmappings[uniprotID]=[mapId+"||"+type]
    print("Done reading mapping file")
    countAddedIdMappings=0
    mappingsToSaveInDB=[]
    for p in proteinsForSpecies:
        ensembl=""
        # First, checking for self-mappings
        existingMapping=IdMapping.objects.filter(Q(protein=p) & Q(type="UniProtID") & Q(mapped_id=p.uniprot_id))
        if not existingMapping.exists():
            idMapping = IdMapping(protein=p, type="UniProtID",mapped_id=p.uniprot_id)
            mappingsToSaveInDB.append(idMapping)
            countAddedIdMappings+=1
        # Then going over all mappings that SHOULD be there
        for m in idmappings[p.uniprot_id]:
            # Creating ID mapping to be put in DB
            mappingID, mappingType=m.split("||")
            existingMapping=IdMapping.objects.filter(Q(protein=p) & Q(type=mappingsOfInterest[mappingType]) & Q(mapped_id=mappingID))
            if not existingMapping.exists():
                idMapping = IdMapping(protein=p, type=mappingsOfInterest[mappingType],mapped_id=mappingID)
                mappingsToSaveInDB.append(idMapping)
                countAddedIdMappings+=1
                if mappingsOfInterest[mappingType]=="Ensembl":
                    ensembl=mappingID
        # If the same ensembl is mapped to different uniprotIDs,
        # they can be viewed as alternative uniprot IDs, and are stored with mapping type alt_uniprot_id
        if ensembl!="" and ensembl in altUniprotIDsByEnsembl:
            for m in altUniprotIDsByEnsembl[ensembl]:
                if m!=p and m not in proteinsToUse :
                    numEnsemblMappings=0 # Making sure not to add alternative ids having multiple ensembl mappings, as those can be to ambiguous
                    for altMapping in idmappings[m]:
                        if mappingsOfInterest[altMapping.split("||")[1]]=="Ensembl":
                            numEnsemblMappings+=1
                    if numEnsemblMappings<2:
                        idMapping = IdMapping(protein=p, type="alt_uniprot_id",mapped_id=m)
                        mappingsToSaveInDB.append(idMapping)
                        countAddedIdMappings+=1
                        for altMapping in idmappings[m]:
                            mappingID, mappingType=altMapping.split("||")
                            idMapping = IdMapping(protein=p, type="alt_"+mappingsOfInterest[mappingType],mapped_id=mappingID)
                            mappingsToSaveInDB.append(idMapping)
                            countAddedIdMappings+=1

    print("Got "+str(countAddedIdMappings)+" new mappings for: "+s)
    # Writing all proteins and mappings to DB
    IdMapping.objects.bulk_create(mappingsToSaveInDB)
    return []

def getGenome(species, proteomeConfig):
    # Gets proteomes for a list of species, and checks if they exist in the database.
    # For unrecognized proteomes (species and version, defined in proteomeConfig.yml), the proteome file is fetched,
    # proteins are stored in the database, using uniprotID as a primary identifier.
    # Mappings are extracted if their mapping type is in the mappingsOfInterest and stored in the database.
    # Returning a list of proteomes to be used in the instance.

    # Checking if the requested species for version already exist in the DB. If not, adding to list
    speciesToGet=[]
    proteomesForInstance=[]
    for s in species:
        proteomeInDB = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']) &  Q(species__tax_id=s))
        if proteomeInDB.exists():
            proteomesForInstance.append(proteomeInDB[0])
            print("Proteome Exists, species: "+s+" version: "+proteomeConfig['genome']['version'])
        else:
            speciesToGet.append(str(s))

    if len(speciesToGet)>0:
        # Downloading README file
        downloadFile(proteomeConfig['genome']['urlUniProtReadme'], "README")
        # Find translations for taxIDs, and store species information in DB
        ebiIdentifiers=readREADME(speciesToGet)

        # Collecting proteome files in parallell
        cpus = os.cpu_count()
        parallel = Parallel(n_jobs=cpus)
        results = parallel(delayed(getProteomesAndMappings)(s, ebiIdentifiers, proteomeConfig, proteomeConfig['genome']['mappingsOfInterest']) for s in ebiIdentifiers)
        proteomesForInstance.extend(results)

    elif proteomeConfig['genome']['addMappings']==True:
        print("Collecting new mappings")
        # Downloading README file
        downloadFile(proteomeConfig['genome']['urlUniProtReadme'], "README")
        # Find translations for taxIDs, and store species information in DB
        ebiIdentifiers=readREADME(species)
        # Collecting proteome files in parallell
        cpus = os.cpu_count()
        parallel = Parallel(n_jobs=cpus)
        parallel(delayed(updateMappings)(s, ebiIdentifiers, proteomeConfig, proteomeConfig['genome']['mappingsOfInterest'] ) for s in ebiIdentifiers)

    # Remove temporary files (readme, xml.gz and .xml)
    #removeTempFiles("data/tmp/README")
    # for s in speciesToGet:
    #     removeTempFiles("data/tmp/"+s+"fasta.gz")
    #     removeTempFiles("data/tmp/"+s+"idmapping.gz")

    return proteomesForInstance


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

    ## Tmp fix dog, chicking and zebrafish
    import pandas as pd
    idMappingToAdd = pd.read_csv('data/tmp/orthologNetwork/InParanoiDB9_mappingTOadd.tsv', sep='\t')
    # idMappingToAdd_list = idMappingToAdd.values.tolist()
    # IdMapping.objects.bulk_create(idMappingToAdd_list)
    mem_csv = StringIO()
    idMappingToAdd.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing %i idmapping to DB" %len(idMappingToAdd))
    ## TODO --> add to table the csv_io module
    with closing(mem_csv) as csv_io:
        IdMapping.objects.from_csv(csv_io)
