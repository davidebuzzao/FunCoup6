import os
import yaml
from yaml.loader import SafeLoader
import pandas as pd
import matplotlib.pyplot as plt
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from joblib import Parallel, delayed

def count_numProteins(tax_id):
    proteome_obj = Proteome.objects.get(Q(species__tax_id=tax_id))
    proteins = Protein.objects.filter(Q(proteome=proteome_obj))
    sp_name = Species.objects.get(Q(tax_id=tax_id)).species_name
    speciesA_name_compact = '%s.%s' %(sp_name.split()[0][0],sp_name.split()[1])

    if proteins.exists():
        nproteins = proteins.count()
        print('%s\t%i'%(speciesA_name_compact,nproteins))
    else:
        print('%s '%(speciesA_name_compact))

def count_goldstandardlinks(goldStandard):
    g = GoldStandardLink.objects.filter(Q(goldStandard=goldStandard))
    # species_name = goldStandard.species.species_name
    species_tax_id = goldStandard.species.tax_id
    goldStandard_name = '%s (%s)' %(goldStandard.type,goldStandard.version)
    if g.exists():
        nlinks = g.count()
        print('%s\t%s\t%i'%(species_tax_id,goldStandard_name,nlinks))
    else:
        print('%s\t%s\t '%(species_tax_id,goldStandard_name))

def count_evidencelinks(evidence):
    ProteinLinkScore = apps.get_model(app_label='data', model_name=evidence.type)
    p = ProteinLinkScore.objects.filter(Q(evidence=evidence))
    evidencelink = ProteinLinkScore.objects.filter(Q(evidence=evidence)).order_by('score')
    species = evidence.species.species_name
    evidence_name = '%s (%s)' %(evidence.type,evidence.version)
    if evidencelink.exists():
        nlinks = evidencelink.count()
        print('%s\t%s\t%i'%(species,evidence_name,nlinks))
        # minscore = evidencelink.first().score
        # maxscore = evidencelink.last().score
        # print('%s\t%s\t%i\t%.4f\t%.4f'%(species,evidence_name,nlinks,minscore,maxscore))
    else:
        print('%s\t%s\t \t \t '%(species,evidence_name))

def getIndex(evidence):
    ProteinLinkScore = apps.get_model(app_label='data', model_name=evidence.type)
    p = ProteinLinkScore.objects.filter(Q(evidence=evidence))
    # print(p.count(), p[0])

def countSuccessfulTraining(speciesA):
    ## Extract Gold Standard and Evidence
    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id=speciesA)) \
                        if goldstandard_t.type in instanceConfig['instance']['goldStandard']]
    gs_df = pd.DataFrame.from_records(GoldStandard.objects.filter(Q(species__tax_id=speciesA)).values('id','type'))
    pLLR = pd.DataFrame.from_records(PolynomialLLR.objects.filter(Q(goldStandard__in=goldstandard_list)).values('goldStandard','error'))
    pLLR.columns = ['id','error']
    gs_df = pd.merge(gs_df,pLLR, on='id')
    gs_df = gs_df.drop('id',axis=1)
    gs_df['error'].fillna('Successful', inplace=True)
    gs_df = gs_df.groupby(by=['type','error']).size().reset_index(name='count')
    # Sort the DataFrame so that the "Successful" error category is plotted first
    gs_df_sorted = gs_df.copy()
    gs_df_sorted['error_order'] = gs_df_sorted['error'].apply(lambda x: 0 if x == 'Successful' else 1)
    gs_df_sorted = gs_df_sorted.sort_values(by=['error_order'])

    # Define data
    types = gs_df_sorted['type'].unique()
    errors = gs_df_sorted['error'].unique()
    num_types = len(types)
    num_errors = len(errors)
    counts = gs_df_sorted.pivot_table(index='type', columns='error', values='count', fill_value=0)

    # Define colors for different types
    type_colors = {'Complex': '#5f8aee',
                'Metabolic': '#fd9e02',
                'PPI': '#a25fac',
                'Regulatory': '#ff9191',
                'Signaling': '#fff56a'}

    # Plot
    fig, axes = plt.subplots(nrows=num_errors, ncols=1, figsize=(10, 12), sharex=True)

    for i, error in enumerate(errors):
        ax = axes[i]
        for j, type_name in enumerate(types):
            ax.bar(type_name, counts[error][j], color=type_colors[type_name], label=type_name)

        ax.set_ylabel('Count')
        if error=='Successful':ax.set_title(error)
        else: ax.set_title(f'Error: {error}')
        ax.legend()

    plt.tight_layout()
    plt.savefig('%s_GScounts' %(speciesA))
    plt.close()

def computeLinksScoreStats(instanceConfig,goldStandardConfig,proteomeConfig,evidenceConfig,opt='delete'):
    species_id = [species.id for species in Species.objects.filter(Q(tax_id__in=instanceConfig['instance']['species']))]
    cpus = os.cpu_count()

    ### Stats on goldstandard
    
    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species_id__in=species_id)) \
                         if goldstandard_t.type in goldStandardConfig and ((goldstandard_t.version==goldStandardConfig[goldstandard_t.type]['version']) or \
                                                                           (goldstandard_t.version == goldStandardConfig['Regulatory']['RegNetwork']['version'] + ',' + goldStandardConfig['Regulatory']['TRRUST']['version'] or \
                                                                            goldstandard_t.version == goldStandardConfig['Regulatory']['RegulonDB']['version'] or \
                                                                            goldstandard_t.version == goldStandardConfig['Regulatory']['Yeastract']['version']))]
    
    if opt=='count_goldstandard':
        print('Species\tGoldStandard\tNr.ofLinks')
        Parallel(n_jobs=cpus)(delayed(count_goldstandardlinks)(goldStandard) for goldStandard in goldstandard_list)


    ### Stats on evidences
    evidence_list = [evidence_t for evidence_t in Evidence.objects.filter(Q(species_id__in=species_id) & Q(type__in=instanceConfig['instance']['evidence'])) \
        if (evidence_t.type in evidenceConfig and \
            ((evidence_t.version==evidenceConfig[evidence_t.type]['version']) or \
            ((evidence_t.version in sum(evidenceConfig['MEX']['geo']['ma'],[]) or \
            (evidence_t.version in sum(evidenceConfig['MEX']['geo']['rseq'],[])) or \
            (evidence_t.version in sum(evidenceConfig['MEX']['ebi'],[]))))))]

    if opt=='count_evidence':
        # print('Species\tEvidence\tNr.ofLinks\tMINscore\tMAXscore')
        print('Species\tEvidence\tNr.ofLinks')
        Parallel(n_jobs=cpus//2)(delayed(count_evidencelinks)(evidence) for evidence in evidence_list)
    elif opt=='index':
        print('Species\tEvidence\tNr.ofLinks\tSTART\tEND')
        Parallel(n_jobs=cpus)(delayed(getIndex)(evidence) for evidence in evidence_list)

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
    speciesA='9606'
    for tax_id in instanceConfig['instance']['species']:
        count_numProteins(tax_id)