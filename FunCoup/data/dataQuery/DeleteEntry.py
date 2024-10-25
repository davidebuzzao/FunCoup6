import pandas as pd
import os
import yaml
from yaml.loader import SafeLoader
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from joblib import Parallel, delayed

def delete_Homologs4MEX(evidence,eval_cutoff=0.001,bs_cutoff=100):
    tax_id = evidence.species.tax_id
    # diamond_out = 'data/tmp/Diamond_bs%i/diamond%s-%s' %(bs_cutoff,tax_id,tax_id)
    formatted_eval_cutoff = '{:.0e}'.format(eval_cutoff)
    diamond_out = 'data/tmp/Diamond_eval%s/diamond%s-%s' %(formatted_eval_cutoff,tax_id,tax_id)
    if os.path.exists(diamond_out): 
        diamond_df = pd.read_csv(diamond_out, sep='\t', header=None)
        diamond_df.columns = ['proteinA','proteinB']
    
    ProteinLinkScore = apps.get_model(app_label='data', model_name=evidence.type)
    evidenceLink = pd.DataFrame.from_records(ProteinLinkScore.objects.filter(Q(evidence=evidence)).values('id','proteinA','proteinB'))
    evidenceLink.columns = ['id','proteinA','proteinB']

    ## Delete 
    evidenceLink_todelete = pd.merge(evidenceLink, diamond_df, on=['proteinA','proteinB'], how='inner')
    if len(evidenceLink_todelete)>0:
        id_list = evidenceLink_todelete['id'].to_list()
        print('Deleting #%i links from %s(%s)' %(len(id_list),evidence.type,evidence.version))
        queryset = ProteinLinkScore.objects.filter(Q(id__in=id_list))
        queryset.delete()
    else:
        print('No links to delete from %s(%s)' %(evidence.type,evidence.version))
        return

def delete_EvidenceWith0Links(evidence):
    ProteinLinkScore = apps.get_model(app_label='data', model_name=evidence.type)
    nlinks = ProteinLinkScore.objects.filter(Q(evidence=evidence)).count()
    if nlinks==0:
        print('Removing %s (%s) for species %s' %(evidence.type,evidence.version,evidence.species.tax_id))
        evidence.delete()
    else:
        print('%s (%s) for species %s: %i links' %(evidence.type,evidence.version,evidence.species.tax_id,nlinks))

def delete_NetworkLinks(tax_id): #["9031","7955", "9615"]
    # tax_id = '9031'
    proteome = Proteome.objects.get(Q(species__tax_id=tax_id))
    proteins = [p.id for p in Protein.objects.filter(Q(proteome=proteome))]
    proteins_list = list(proteins)
    network_links = pd.DataFrame.from_records(Network.objects.filter(Q(proteinA_id__in=proteins_list) & Q(proteinB_id__in=proteins_list)).values('id'))
    network_links = Network.objects.filter(Q(proteinA_id__in=proteins_list) & Q(proteinB_id__in=proteins_list))
    network_links.count()
    if len(network_links)>0:
        network_links.delete()

def getIndex(evidence):
    ProteinLinkScore = apps.get_model(app_label='data', model_name=evidence.type)
    p = ProteinLinkScore.objects.filter(Q(evidence=evidence))

def deleteLinksScoreStats(instanceConfig,goldStandardConfig,proteomeConfig,evidenceConfig,opt='0links'):
    species_id = [species.id for species in Species.objects.filter(Q(tax_id__in=instanceConfig['instance']['species']))]
    cpus = os.cpu_count()

    ### Stats on goldstandard
    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species_id__in=species_id)) \
                         if goldstandard_t.type in goldStandardConfig and ((goldstandard_t.version==goldStandardConfig[goldstandard_t.type]['version']) or \
                                                                           (goldstandard_t.version == goldStandardConfig['Regulatory']['RegNetwork']['version'] + ',' + goldStandardConfig['Regulatory']['TRRUST']['version'] or \
                                                                            goldstandard_t.version == goldStandardConfig['Regulatory']['RegulonDB']['version'] or \
                                                                            goldstandard_t.version == goldStandardConfig['Regulatory']['Yeastract']['version']))]
    
    ### Stats on evidences
    evidence_list = [evidence_t for evidence_t in Evidence.objects.filter(Q(species_id__in=species_id) & Q(type__in=instanceConfig['instance']['evidence'])) \
        if (evidence_t.type in evidenceConfig and \
            ((evidence_t.version==evidenceConfig[evidence_t.type]['version']) or \
            ((evidence_t.version in sum(evidenceConfig['MEX']['geo']['ma'],[]) or \
            (evidence_t.version in sum(evidenceConfig['MEX']['geo']['rseq'],[])) or \
            (evidence_t.version in sum(evidenceConfig['MEX']['ebi'],[]))))))]

    if opt=='0links':
        print('Delete Evidence entries with 0 links')
        Parallel(n_jobs=cpus)(delayed(delete_EvidenceWith0Links)(evidence) for evidence in evidence_list)
    elif opt=='MEXhomologs':
        print('Delete MEX entries between homologs')
        Parallel(n_jobs=cpus)(delayed(delete_Homologs4MEX)(evidence) for evidence in evidence_list if evidence.type=='MEX')
        ##TODO --> Have deleted some DOM,MIR,QMS,SCL links by mistake

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
    print(instanceConfig['instance']['species'])