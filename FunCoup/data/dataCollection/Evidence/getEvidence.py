import zipfile
import os
from joblib import Parallel, delayed
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed

from auxiliary_functions import *
from .getHomologs4MEX import *

from .DOM import *
from .GIN import *
from .GRG import *
from .MEX import *
from .MIR import *
from .PEX import *
from .PHP import *
from .PIN import *
from .PIN_HostPathogen import *
from .SCL import *
from .TFB import *

def unzipIRefIndexFile(filename,path):
    # Special unzip for iRefIndex, as the file is broken (v.22-08), where the compressed file is named "-"
    print("unzipping file: "+filename+".zip")
    with zipfile.ZipFile(path+filename+".zip", "r") as zip_ref:
       zip_ref.extract("-", path=path)
    os.rename(path+"-",path+filename)

def checkForEvidenceInDB(instanceSpecies, evidenceConfig, evidenceType):
    speciesToCollect=[]
    for s in instanceSpecies:
        if s in evidenceConfig['species']:
            evidenceInDB = Evidence.objects.filter(Q(type=evidenceType) & Q(species__tax_id=s) & Q(scoringMethod=evidenceConfig['scoring_method']) & Q(version=evidenceConfig['version']))
            if evidenceInDB.exists(): # If the evidence already exist in the DB, using existing one for the instance
                print("Evidence Exists: "+evidenceType+" v."+str(evidenceConfig['version'])+" species: "+s)
            else: # Otherwise adding the species to fetch links
                speciesToCollect.append(s)
    return speciesToCollect

def checkForHostPathogenEvidenceInDB(instanceSpecies, evidenceConfig, evidenceType):
    pathogenToCollect=[]
    for hostPathogen in instanceSpecies:
        pathogen_id = hostPathogen[1]
        if pathogen_id in evidenceConfig['pathogens']:
            evidenceInDB = Evidence.objects.filter(Q(type=evidenceType) & Q(species__tax_id=pathogen_id) & Q(scoringMethod=evidenceConfig['scoring_method']) & Q(version=evidenceConfig['version']))
            if evidenceInDB.exists(): # If the evidence already exist in the DB, using existing one for the instance
                print("Evidence Exists: "+evidenceType+" v."+str(evidenceConfig['version'])+" species: "+pathogen_id)
            else: # Otherwise adding the species to fetch links
                pathogenToCollect.append(pathogen_id)
    return pathogenToCollect

def getDOM(species,evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "DOM")
    if len(speciesToCollect)>0:
        if not os.path.exists("data/tmp/"+evidenceConfig['url'].split("/")[-1]): # Might already have been downloaded.
            downloadFile(evidenceConfig['url'], "data/tmp/"+evidenceConfig['url'].split("/")[-1])
        getDOMLinks(evidenceConfig,speciesToCollect,proteomeConfig)
        #removeTempFile("data/tmp/"+evidenceConfig['url'].split("/")[-1])

def getGIN(species,evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "GIN")
    if len(speciesToCollect)>0:
        indexes = [evidenceConfig['species'].index(x) for x in speciesToCollect]
        evidenceConfig['species'] = [evidenceConfig['species'][i] for i in indexes]
        evidenceConfig['species_name'] = [evidenceConfig['species_name'][i] for i in indexes]
        get_GINLinks(evidenceConfig,proteomeConfig)
        # evidenceConfig['file_name'] = [evidenceConfig['file_name'][i] for i in indexes]
        # evidenceConfig['file_columns'] = [evidenceConfig['file_columns'][i] for i in indexes]
        # evidenceConfig['version'] = [evidenceConfig['version'][i] for i in indexes]

        # for i in range(len(evidenceConfig['species'])):
        #     tax_id = evidenceConfig['species'][i]
        #     get_GINLinks(evidenceConfig['file_name'][i],evidenceConfig['file_columns'][i],\
        #                  evidenceConfig['version'][i],tax_id,evidenceConfig,proteomeConfig)


def getGRG(species, evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "GRG")
    if len(speciesToCollect)>0:
        indexes = [evidenceConfig['species'].index(x) for x in speciesToCollect]
        evidenceConfig['species'] = [evidenceConfig['species'][i] for i in indexes]
        evidenceConfig['species_name'] = [evidenceConfig['species_name'][i] for i in indexes]
        evidenceConfig['species_chromosomes'] = [evidenceConfig['species_chromosomes'][i] for i in indexes]
        evidenceConfig['url_annotation'] = [evidenceConfig['url_annotation'][i] for i in indexes]
        evidenceConfig['genome_version'] = [evidenceConfig['genome_version'][i] for i in indexes]
        evidenceConfig['encode_pipeline'] = [evidenceConfig['encode_pipeline'][i] for i in indexes]

        get_GRGLinks(evidenceConfig,proteomeConfig)


def getMEX_inParallel(runParameters,evidenceConfig,proteomeConfig,eval_cutoff=0.001):
    gse_name_list,tax_id,database,exp = runParameters
    evidence = get_MEX(gse_name_list,tax_id,evidenceConfig,proteomeConfig,evidenceConfig['geo']['gpl'],database,exp,eval_cutoff)
    return evidence

def getMEX(species,evidenceConfig,proteomeConfig,eval_cutoff=0.001):
    ## Extract gene pairs with BitScore > cutoff
    getHomologs(species,proteomeConfig,eval_cutoff=eval_cutoff)

    ## The check if evidence exist is done inside MEX.py
    if len(species)>0:
        indexes = [evidenceConfig['species'].index(x) for x in species]
        evidenceConfig['species'] = [evidenceConfig['species'][i] for i in indexes]
        evidenceConfig['geo']['ma'] = [evidenceConfig['geo']['ma'][i] for i in indexes]
        evidenceConfig['geo']['rseq'] = [evidenceConfig['geo']['rseq'][i] for i in indexes]
        evidenceConfig['ebi'] = [evidenceConfig['ebi'][i] for i in indexes]

        runlist = []
        for i,sp in enumerate(evidenceConfig['species']):
            exp_count=0
            num_dataset=5
            for j,ebi in enumerate(evidenceConfig['ebi'][i]):
                if j<num_dataset:
                    exp_count+=1
                    runlist.append([ebi,sp,'EBI','rseq'])
            if exp_count<num_dataset:
                for rseq in evidenceConfig['geo']['rseq'][i]:
                    if exp_count < num_dataset:
                        exp_count+=1
                        runlist.append([rseq,sp,'GEO','rseq'])
            if exp_count<num_dataset:
                for ma in evidenceConfig['geo']['ma'][i]:
                    if exp_count < num_dataset:
                        exp_count+=1
                        runlist.append([ma,sp,'GEO','ma'])

        cpus=os.cpu_count()
        Parallel(n_jobs=cpus)(delayed(getMEX_inParallel)(run,evidenceConfig,proteomeConfig) for run in runlist)


def getMIR(species,evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "MIR")
    if len(speciesToCollect)>0:
        indexes = [evidenceConfig['species'].index(x) for x in speciesToCollect]
        evidenceConfig['species'] = [evidenceConfig['species'][i] for i in indexes]
        evidenceConfig['files'] = [evidenceConfig['files'][i] for i in indexes]
        get_MIRLinks(evidenceConfig,proteomeConfig)


def getPEX(species,evidenceConfig,proteomeConfig):    
    ### Jaccard Index
    tmp_evidenceConfig = evidenceConfig.copy()
    tmp_evidenceConfig['scoring_method'] = evidenceConfig['scoring_method'][0]
    speciesToCollect = checkForEvidenceInDB(species,tmp_evidenceConfig, "PEX")
    if len(speciesToCollect)>0:
        tmp_evidenceConfig['species'] = speciesToCollect
        if not os.path.exists("data/tmp/"+tmp_evidenceConfig['url'].split("/")[-1]): # Might already have been downloaded.
            downloadAndUnzipFile(tmp_evidenceConfig['url'])
        getPEXLinks(tmp_evidenceConfig,proteomeConfig)
    
    ### Spearman Correlation
    tmp_evidenceConfig = evidenceConfig.copy()
    tmp_evidenceConfig['scoring_method'] = evidenceConfig['scoring_method'][1]
    speciesToCollect = checkForEvidenceInDB(species,tmp_evidenceConfig, "PEX")
    if len(speciesToCollect)>0:
        tmp_evidenceConfig['species'] = speciesToCollect
        if not os.path.exists("data/tmp/"+tmp_evidenceConfig['url'].split("/")[-1]): # Might already have been downloaded.
            downloadAndUnzipFile(tmp_evidenceConfig['url'])
        getPEXLinks(tmp_evidenceConfig,proteomeConfig)

def getPIN(species,evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "PIN")
    if len(speciesToCollect)>0:
        if os.path.exists("data/tmp/"+evidenceConfig['url'].split("/")[-1])==False: # Might already have been downloaded.
            downloadFile(evidenceConfig['url'], "data/tmp/"+evidenceConfig['url'].split("/")[-1])
        getPINLinks(evidenceConfig, speciesToCollect,proteomeConfig)
        #removeTempFile("data/tmp/"+evidenceConfig['url'].split("/")[-1])
    pathogenToCollect = checkForHostPathogenEvidenceInDB(species,evidenceConfig, "PIN")
    if len(pathogenToCollect)>0:
        indexes = [evidenceConfig['pathogens'].index(x) for x in species]
        evidenceConfig['pathogens'] = [evidenceConfig['pathogens'][i] for i in indexes]
        evidenceConfig['urlPathogen'] = [evidenceConfig['urlPathogen'][i] for i in indexes]
        evidenceConfig['filePathogen'] = [evidenceConfig['filePathogen'][i] for i in indexes]
        
        ## TODO --> To setup for more host-pathogen interaction
        for i,sp in enumerate(evidenceConfig['pathogens']):
            if os.path.exists("data/tmp/"+evidenceConfig['urlPathogen'][i].split("/")[-1])==False: # Might already have been downloaded.
                downloadFile(evidenceConfig['urlPathogen'][i], "data/tmp/"+evidenceConfig['filePathogen'][i])
        getPINLinks4Pathogen(pathogenToCollect, evidenceConfig, proteomeConfig)

def getPHP(species,evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "PHP")
    if len(speciesToCollect)>0:
        if not os.path.exists("data/tmp/"+evidenceConfig['distanceMatrixFileName']): # Might already have been downloaded.
            downloadFile(evidenceConfig['distanceMatrixURL'], "data/tmp/"+evidenceConfig['distanceMatrixFileName'])
        getPHPLinks("data/tmp/"+evidenceConfig['distanceMatrixFileName'], evidenceConfig,speciesToCollect,proteomeConfig)
        #removeTempFile("data/tmp/"+evidenceConfig['distanceMatrixFileName'])

def getSCL(species,evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "SCL")
    if len(speciesToCollect)>0:
        getSCLLinks(evidenceConfig,speciesToCollect,proteomeConfig)


def getTFB(species, evidenceConfig,proteomeConfig):
    speciesToCollect = checkForEvidenceInDB(species,evidenceConfig, "TFB")
    if len(speciesToCollect)>0:
        indexes = [evidenceConfig['species'].index(x) for x in speciesToCollect]
        evidenceConfig['species'] = [evidenceConfig['species'][i] for i in indexes]
        evidenceConfig['url'] = [evidenceConfig['url'][i] for i in indexes]
        get_TFBLinks(evidenceConfig,proteomeConfig)

def getEvidence(instanceConfig,evidenceConfig,proteomeConfig):
    # Getting evidence defined in the generate instance config file.
    # Returning the evidence to be used for this instance
    for evidenceType in instanceConfig['instance']['evidence']:
        if evidenceType=="DOM":
            print("Getting DOM")
            getDOM(instanceConfig['instance']['species'],evidenceConfig['DOM'],proteomeConfig)
        elif evidenceType=="GIN":
            print("Getting GIN")
            getGIN(instanceConfig['instance']['species'],evidenceConfig['GIN'],proteomeConfig)
        elif evidenceType=="GRG":
            print("Getting GRG")
            getGRG(instanceConfig['instance']['species'],evidenceConfig['GRG'],proteomeConfig)
        elif evidenceType=="MEX":
            print("Getting MEX")
            getMEX(instanceConfig['instance']['species'],evidenceConfig['MEX'],proteomeConfig)
        elif evidenceType=="MIR":
            print("Getting MIR")
            getMIR(instanceConfig['instance']['species'],evidenceConfig['MIR'],proteomeConfig)
        elif evidenceType=="PEX":
            print("Getting PEX")
            getPEX(instanceConfig['instance']['species'],evidenceConfig['PEX'],proteomeConfig)
        elif evidenceType=="PIN":
            print("Getting PIN")
            getPIN(instanceConfig['instance']['species'],evidenceConfig['PIN'],proteomeConfig)
        elif evidenceType=="PHP":
            print("Getting PHP")
            getPHP(instanceConfig['instance']['species'],evidenceConfig['PHP'],proteomeConfig)
        elif evidenceType=="SCL":
            print("Getting SCL")
            getSCL(instanceConfig['instance']['species'],evidenceConfig['SCL'],proteomeConfig)
        elif evidenceType=="TFB":
            print("Getting TFB")
            getTFB(instanceConfig['instance']['species'],evidenceConfig['TFB'],proteomeConfig)
        else:
            print("Did not recognize evidence: "+evidenceType)
            exit()
