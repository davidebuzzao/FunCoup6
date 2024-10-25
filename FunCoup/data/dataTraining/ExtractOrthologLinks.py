import os
import pandas as pd
import swifter
import yaml
from yaml.loader import SafeLoader
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.apps import apps
from django.db.models import Q

def extractOrthologGoldStandardLinks(goldStandard,speciesA):
    speciesB = goldStandard.species.tax_id

    # Extract orthologs between speciesA and speciesB
    if int(speciesA)>int(speciesB):
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesA) & Q(speciesB=speciesB)).values('proteinA','proteinB'))
    else:
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesB) & Q(speciesB=speciesA)).values('proteinB','proteinA'))
    ortholog_protein_df.columns = ['proteinA','proteina']

    # links = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=goldStandard)).values('proteinA_id','proteinB_id'))
    ## TODO --> could throw error 'pq: could not resize shared memory segment. No space left on device'
    linksDf = pd.DataFrame.from_records(\
        GoldStandardLink.objects.filter(Q(goldStandard=goldStandard) & \
                                    Q(proteinA__in=ortholog_protein_df['proteina']) & \
                                        Q(proteinB__in=ortholog_protein_df['proteina']))\
                                            .values('proteinA_id','proteinB_id'))
    if (linksDf.shape[0]==0): return None
    linksDf.columns = ['proteina','proteinb']

    ## Extract ortholog links by merge and convert to homologs
    tmp_merged = pd.merge(linksDf, ortholog_protein_df, on='proteina', how='inner')
    ortholog_protein_df.columns = ['proteinB','proteinb']
    tmp_merged = pd.merge(tmp_merged, ortholog_protein_df, on='proteinb', how='inner')

    OrthologLinks = pd.DataFrame(\
                tmp_merged.drop(['proteina','proteinb'], axis=1)\
                    .apply(lambda r: sorted(r,reverse=True), axis = 1)\
                        .to_list(),columns=['proteinA_id','proteinB_id'])\
                            .drop_duplicates()\
                                .reset_index(drop=True)

    return OrthologLinks


def extractOrthologLinks4Training(evidence,speciesA):
    speciesB = evidence.species.tax_id

    # Extract orthologs between speciesA and speciesB
    if int(speciesA)>int(speciesB):
        ## ProteinA - ProteinB
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesA) & Q(speciesB=speciesB)).values('proteinA','proteinB'))
    else:
        ## ProteinB - ProteinA
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesB) & Q(speciesB=speciesA)).values('proteinB','proteinA'))
    ortholog_protein_df.columns = ['proteinA','proteina']
    # links = pd.DataFrame.from_records(ProteinLinkScore.objects.filter(Q(evidence=evidence)).values('proteinA_id','proteinB_id', 'score'))
    ## TODO --> could throw error 'pq: could not resize shared memory segment. No space left on device'
    ProteinLinkScore = apps.get_model(app_label='data', model_name=evidence.type)
    if not ProteinLinkScore.objects.filter(Q(evidence=evidence)).exists(): return None, 0, 0
    p = ProteinLinkScore.objects.filter(Q(evidence=evidence)).order_by('score')
    abs_min, abs_max = p.first().score, p.last().score
    linksDf = pd.DataFrame.from_records(\
        ProteinLinkScore.objects.filter(Q(evidence=evidence) & \
                                    Q(proteinA__in=ortholog_protein_df['proteina']) & \
                                        Q(proteinB__in=ortholog_protein_df['proteina']))\
                                            .values('proteinA_id','proteinB_id', 'score'))
    if (linksDf.shape[0]==0): return None, 0, 0
    linksDf.columns = ['proteina','proteinb','score']

    ## Extract ortholog links by merge and convert to homologs
    tmp_merged = pd.merge(linksDf, ortholog_protein_df, on='proteina', how='inner')
    ortholog_protein_df.columns = ['proteinB','proteinb']
    tmp_merged = pd.merge(tmp_merged, ortholog_protein_df, on='proteinb', how='inner')
    tmp_merged = tmp_merged.drop(['proteina','proteinb'], axis=1)\
                    .drop_duplicates()\
                            .reset_index(drop=True)
    
    orthologLinks = pd.concat([
        pd.DataFrame(\
            tmp_merged.drop('score', axis=1)\
                .apply(lambda r: sorted(r,reverse=True), axis = 1)\
                    .to_list(),columns=['proteinA_id','proteinB_id']),tmp_merged['score']],axis=1)\
                        .drop_duplicates()\
                            .reset_index(drop=True)

    return orthologLinks, abs_min, abs_max

def extractOrthologLinks4Redundancy(evidenceBLink,speciesA,speciesB):
    # Extract orthologs between speciesA and speciesB
    if int(speciesA)>int(speciesB):
        ## ProteinA - ProteinB
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesA) & Q(speciesB=speciesB)).values('proteinA','proteinB'))
    else:
        ## ProteinB - ProteinA
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesB) & Q(speciesB=speciesA)).values('proteinB','proteinA'))
    ortholog_protein_df.columns = ['proteinA','proteina']
    if len(evidenceBLink)==0: return None
    evidenceBLink.columns = ['proteina','proteinb','score']
    # print(evidenceBLink)
    ## Extract ortholog links by merge and convert to homologs
    tmp_merged = pd.merge(evidenceBLink, ortholog_protein_df, on='proteina', how='inner')
    ortholog_protein_df.columns = ['proteinB','proteinb']
    tmp_merged = pd.merge(tmp_merged, ortholog_protein_df, on='proteinb', how='inner')
    tmp_merged = tmp_merged.drop(['proteina','proteinb'], axis=1)\
                    .drop_duplicates()\
                            .reset_index(drop=True)
    
    orthologBLinks = pd.concat([
        pd.DataFrame(\
            tmp_merged.drop('score', axis=1)\
                .apply(lambda r: sorted(r,reverse=True), axis = 1)\
                    .to_list(),columns=['proteinA_id','proteinB_id']),tmp_merged['score']],axis=1)\
                        .drop_duplicates()\
                            .reset_index(drop=True)
    
    return orthologBLinks

def extractOrthologLinks4Network(geneA,geneB,speciesA,speciesB):
    # Extract orthologs between speciesA and speciesB
    if int(speciesA)>int(speciesB):
        ## ProteinA - ProteinB
        pA = ProteinOrtholog.objects.filter(Q(speciesA=speciesA) & Q(speciesB=speciesB) & Q(proteinA_id=geneA))
        pB = ProteinOrtholog.objects.filter(Q(speciesA=speciesA) & Q(speciesB=speciesB) & Q(proteinA_id=geneB))
        if not pA.exists() or not pB.exists(): return None
        else:
            ortho_geneA = pA[0].proteinB_id
            ortho_geneB = pB[0].proteinB_id
    else:
        ## ProteinB - ProteinA
        pA = ProteinOrtholog.objects.filter(Q(speciesA=speciesB) & Q(speciesB=speciesA) & Q(proteinB_id=geneA))
        pB = ProteinOrtholog.objects.filter(Q(speciesA=speciesB) & Q(speciesB=speciesA) & Q(proteinB_id=geneB))

        if not pA.exists() or not pB.exists(): return None
        else:
            ortho_geneA = pA[0].proteinA_id
            ortho_geneB = pB[0].proteinA_id

    if ortho_geneA>ortho_geneB: return ortho_geneA, ortho_geneB
    else: return ortho_geneB, ortho_geneA

def extract_all_orthologs(speciesA):
    return pd.DataFrame.from_records(\
        ProteinOrtholog.objects.filter(Q(speciesA=speciesA) \
                                       | Q(speciesB=speciesA)))\
                                        .values('proteinA','proteinB','speciesA','speciesB')

def transform_row(row):
    if row['proteinA_id'] > row['proteinB_id']:
        return row
    else:
        new_row = row.copy()
        new_row['proteinA_id'], new_row['proteinB_id'] = row['proteinB_id'], row['proteinA_id']
        swap_direction = {'2': '3', '3': '2'}
        new_row['direction'] = swap_direction.get(new_row['direction'])
        return new_row

def convert_Blinks(evidenceBLink,e_id,speciesA,speciesB,prefix_col):
    if len(evidenceBLink)==0: return None
    tmp_prefix_col = prefix_col.copy()
    tmp_prefix_col[0:2] = ['proteina','proteinb']
    evidenceBLink.columns = tmp_prefix_col + [e_id]

    # Extract orthologs between speciesA and speciesB
    if int(speciesA)>int(speciesB):
        ## ProteinA - ProteinB
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesA) & Q(speciesB=speciesB)).values('proteinA','proteinB'))
    else:
        ## ProteinB - ProteinA
        ortholog_protein_df = pd.DataFrame.from_records(ProteinOrtholog.objects.filter(Q(speciesA=speciesB) & Q(speciesB=speciesA)).values('proteinB','proteinA'))
    if len(ortholog_protein_df)==0: return None
    ortholog_protein_df.columns = ['proteinA','proteina']

    ## Extract ortholog links by merge and convert to homologs
    tmp_merged = pd.merge(evidenceBLink, ortholog_protein_df, on='proteina', how='inner')
    ortholog_protein_df.columns = ['proteinB','proteinb']
    tmp_merged = pd.merge(tmp_merged, ortholog_protein_df, on='proteinb', how='inner')
    tmp_merged = tmp_merged.drop(['proteina','proteinb'], axis=1).drop_duplicates().reset_index(drop=True)

    if len(tmp_merged)==0: return None
    if 'direction' in prefix_col: 
        tmp_merged.columns = ['direction',e_id,'proteinA_id','proteinB_id']
        ## if protein_id order inverts and evidence is directed, invert direction 2--> to <--3
        orthologBLinks = pd.concat([
            pd.DataFrame(\
                tmp_merged.drop([e_id], axis=1)\
                    .apply(lambda row: transform_row(row), axis = 1),columns=prefix_col),tmp_merged[[e_id]]],axis=1)\
                            .drop_duplicates()\
                                .reset_index(drop=True)
    else: 
        tmp_merged.columns = [e_id,'proteinA_id','proteinB_id']
        orthologBLinks = pd.concat([
            pd.DataFrame(\
                tmp_merged.drop([e_id], axis=1)\
                    .apply(lambda row: sorted(row,reverse=True), axis = 1).to_list(),columns=prefix_col),tmp_merged[[e_id]]],axis=1)\
                            .drop_duplicates()\
                                .reset_index(drop=True)

    return orthologBLinks


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