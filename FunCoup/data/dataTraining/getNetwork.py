import os
import time
import pandas as pd
pd.options.mode.chained_assignment = None
import dask
import dask.dataframe as dd
import numpy as np
import yaml
from yaml.loader import SafeLoader
from functools import reduce
from sklearn.preprocessing import MinMaxScaler

# from pandarallel import pandarallel
from joblib import Parallel, delayed
from io import StringIO
from tqdm import tqdm

import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from django.conf import settings
from contextlib import closing

from data.dataTraining.ExtractOrthologLinks import *
from data.dataTraining.NegativeGoldStandard import *
from data.dataTraining.PolynomialFittingLLR import *
from data.dataTraining.PPVlogisticFitting import *
from auxiliary_functions import *
import warnings

'''
For each evidence type (e.g. MEX, QMS, TFB, …)
|     Read all evidence datasets (e.g. [A, B, MEX_hsa1, MEX_hsa2, …, MEX_mmu1] ) into memory
|     For each species, make a full dataset:
|     |     For each evidence in memory
|     |     |     Extract orthologs if species_evidence != species_goldstandard

For each species
|     For each gold standard
|     |     Train positive/negative gold standard, store LLR Polynomial function (for quality check)
|     |     Compute LLR of all raw score column-wise
|     |     Perform weighted redundancy integration row-wise 
|     |     Save total LLR_evidence_goldstandard
|     Save:
        Max_pfc
        Pfc as #Complex,Metabolic,Operon,PPI,Regulatory,Signaling
        Fbs_goldstandard = 
        Llr_evidence  = Complex:DOM,MEX,MIR,GIN,PEX,PHP,PIN,QMS,SCL,TFB|hsa,mmu,rno,cfa,...,osa;Metabolic:DOM,MEX,...||hsa,mmu,...;Operon:...
        Direction = # 0 no existing, 1 gs undirected, 2 ->, 3 <- , 5 <->
|     Write to database
'''

def createDictIndexes(speciesA_gs_order,ev_order,sp_order):
    counter = 0
    index_dict = {}
    for gs in speciesA_gs_order:
        index_dict[gs] = index_dict.get(gs,{})
        for ev in ev_order:
            index_dict[gs][ev] = index_dict[gs].get(ev,counter)
            counter += 1
        for sp in sp_order:
            index_dict[gs][sp] = index_dict[gs].get(sp,counter)
            counter += 1
    return index_dict

def tmpDictIndexes(index_dict,subset_tuples,first_ev,last_ev):
    tmp_index_dict = {}
    for sub in subset_tuples:
        g,e_type,e_id,s = sub
        tmp_index_dict[g] = tmp_index_dict.get(g,{})
        ## This works if funcoup is trained on EV at a time, otherwise use ev name
        tmp_index_dict[g]['EV'] = tmp_index_dict[g].get('EV',index_dict[g][e_type])
        tmp_index_dict[g][e_id] = tmp_index_dict[g].get(e_id,index_dict[g][s])

        ## For speed up always include first and last evidence index
        tmp_index_dict[g][first_ev] = tmp_index_dict[g].get(first_ev,index_dict[g][first_ev])
        tmp_index_dict[g][last_ev] = tmp_index_dict[g].get(last_ev,index_dict[g][last_ev])

    return tmp_index_dict

###########################
## Reading input data
def readGoldStandard(g_t,gs_order,gs_directed):
    if gs_directed[gs_order.index(g_t.type)]:
        goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=g_t)).values('proteinA_id','proteinB_id','direction'))
        goldStandardLink['direction'] = goldStandardLink['direction'].astype(float)
        goldStandardLink.columns = ['proteinA_id','proteinB_id',g_t.type]
    else: 
        goldStandardLink = pd.DataFrame.from_records(GoldStandardLink.objects.filter(Q(goldStandard=g_t)).values('proteinA_id','proteinB_id'))
        goldStandardLink[g_t.type] = '1'
    goldStandardLink = goldStandardLink[goldStandardLink.notnull()].reset_index(drop=True) # Same as inner join? but might be necessary for index reasons
    if len(goldStandardLink)>0:
        goldStandardLink = goldStandardLink[goldStandardLink['proteinA_id'] != goldStandardLink['proteinB_id']].reset_index(drop=True)
        goldStandardLink['proteinA_id'] = goldStandardLink['proteinA_id'].astype(int)
        goldStandardLink['proteinB_id'] = goldStandardLink['proteinB_id'].astype(int)
        goldStandardLink = goldStandardLink.sort_values(['proteinA_id','proteinB_id'])
        return goldStandardLink
    else: return None

def orthologEvidence(df,e_id,speciesA,speciesB,prefix_col):
    df = convert_Blinks(df,e_id,speciesA,speciesB,prefix_col)
    if df is None: return None
    ## Many-to-one relationships can indipendently support final LLRs
    df = df.groupby(prefix_col, as_index=False).mean(numeric_only=True)
    # evidenceLink = evidenceLink.groupby(['proteinA_id','proteinB_id'], as_index=False).max(numeric_only=True)
    return df

def convertToRawScore(score,original_min,original_max):
    return score * (original_max - original_min) + original_min

def readEvidence(e_t,prefix_col,speciesA,ev_order,ev_directed):
    e_id,e_type,speciesB = e_t.id,e_t.type,e_t.species.tax_id
    ProteinLinkScore = apps.get_model(app_label='data', model_name=e_type)
    if ev_directed[ev_order.index(e_type)]:
        evidenceLink = pd.DataFrame.from_records(ProteinLinkScore.objects.filter(Q(evidence=e_t)).values('proteinA','proteinB','direction','score'))
        evidenceLink['direction'] = evidenceLink['direction'].astype(float)
    else:
        evidenceLink = pd.DataFrame.from_records(ProteinLinkScore.objects.filter(Q(evidence=e_t)).values('proteinA','proteinB','score'))
    evidenceLink.columns = prefix_col + [e_id]

    ### Convert to raw values and set cutoff for MEX
    if e_type=='MEX':
        # print(evidenceLink)
        scoreRange = e_t.scoreRange.split('|')
        original_min,original_max = float(scoreRange[0]),float(scoreRange[1])
        evidenceLink[e_id] = convertToRawScore(evidenceLink[e_id],original_min,original_max)
        evidenceLink = evidenceLink[evidenceLink[e_id]>0.5]
        # print(evidenceLink)
        if evidenceLink is not None and len(evidenceLink)>0:
            scaler = MinMaxScaler(feature_range=(0, 1))
            evidenceLink[e_id] = np.round(scaler.fit_transform(evidenceLink[e_id].values.reshape(-1, 1)), decimals=4)
        else:
            return evidenceLink

    if evidenceLink is not None and len(evidenceLink)>0:
        if speciesA!=speciesB:
            evidenceLink = orthologEvidence(evidenceLink,e_id,speciesA,speciesB,prefix_col)
        if evidenceLink is not None and len(evidenceLink)>0:
            evidenceLink['proteinA_id'] = evidenceLink['proteinA_id'].astype(int)
            evidenceLink['proteinB_id'] = evidenceLink['proteinB_id'].astype(int)
            evidenceLink = evidenceLink.sort_values(['proteinA_id','proteinB_id'])
    return evidenceLink

###########################
## LLR extraction
def extractLLR(tmp_speciesA_evidence_df,e_id,positive_goldStandardLink_score,goldstandard_t,instanceConfig,goldStandardConfig,trainingConfig,saveToDatabase=False):
    # Negative gold standard
    min_n_sample = trainingConfig['nSample'][0]
    tax_id = Evidence.objects.get(Q(id=e_id)).species.tax_id
    negative_goldStandardLink = sampleNegativeGoldstandard(tmp_speciesA_evidence_df,tax_id,trainingConfig)
    if negative_goldStandardLink is None and saveToDatabase: 
        write_pLLRToDatabase(Evidence.objects.get(Q(id=e_id)),goldstandard_t,['Nan'],motivation='|NGS|<%i' %min_n_sample)
        return #print('negative_goldStandardLink = None')
    negative_goldStandardLink = filterNegativeGoldStandard(negative_goldStandardLink,goldstandard_t,e_id,instanceConfig)
    if negative_goldStandardLink is None and saveToDatabase: 
        write_pLLRToDatabase(Evidence.objects.get(Q(id=e_id)),goldstandard_t,['Nan'],motivation='|NGS|<%i' %min_n_sample)
        return
    # Train positive/negative gold standard, store LLR Polynomial function (for quality check)
    # Run KDE + Polynomial fitting training
    polynomialLLR(e_id,positive_goldStandardLink_score,goldstandard_t,negative_goldStandardLink,trainingConfig,saveToDatabase=saveToDatabase)

def trainGoldStandard(positive_goldStandardLink,tmp_speciesA_evidence_df,e_id,goldstandard_t,ev_directed,gs_directed,ev_order,gs_order,e_type,g_type,instanceConfig,goldStandardConfig,trainingConfig,saveToDatabase=False):
    pLLR = PolynomialLLR.objects.filter(Q(evidence_id=e_id) & Q(goldStandard=goldstandard_t))
    if not pLLR.exists() or not saveToDatabase:
        if ev_directed[ev_order.index(e_type)] & gs_directed[gs_order.index(g_type)]:
            positive_goldStandardLink_score = pd.merge(positive_goldStandardLink,tmp_speciesA_evidence_df,on=['proteinA_id','proteinB_id','direction'],how='inner')[e_id].to_list()
        else:
            positive_goldStandardLink_score = pd.merge(positive_goldStandardLink,tmp_speciesA_evidence_df,on=['proteinA_id','proteinB_id'],how='inner')[e_id].to_list()
        
        min_n_sample = trainingConfig['nSample'][0]
        if (positive_goldStandardLink_score is None or len(positive_goldStandardLink_score)<min_n_sample) and saveToDatabase:
            write_pLLRToDatabase(Evidence.objects.get(Q(id=e_id)),goldstandard_t,['Nan'],motivation='|PGS|<%i' %min_n_sample)
            return
        
        extractLLR(tmp_speciesA_evidence_df,e_id,positive_goldStandardLink_score,goldstandard_t,instanceConfig,goldStandardConfig,trainingConfig,saveToDatabase=saveToDatabase)

def convertToLLR_parallel(scores, goldstandard_t, prefix_col):
    ev_id = scores.name
    if ev_id not in prefix_col:
        pLLR = PolynomialLLR.objects.get(Q(evidence_id=ev_id) & Q(goldStandard=goldstandard_t))
        coefficients = pLLR.function
        if 'Nan' not in coefficients:
            min_score, max_score = pLLR.scoreRange.split('|')
            ## TODO: if min_score>max_score: return None ## to skip decreasing LLR functions
            g_type = GoldStandard.objects.get(Q(id=pLLR.goldStandard_id)).type
            coefficients = [float(i) for i in coefficients.split(',')]
            result_col = np.round(np.polyval(coefficients, np.clip(scores, float(min_score), float(max_score))), 3)
            col_name = f'{ev_id}_{g_type}'  
            return pd.Series(result_col,name=col_name).astype("Sparse[float]")          
            # return pd.Series(result_col,name=col_name).astype(pd.SparseDtype("float", np.nan)) # return sparse array
            # return pd.Series(result_col,name=col_name)
    else: return None

def format_correlation_rows(row,gs_id):
    #correlation | evidenceA_id | evidenceB_id |         pvalue          |   size    | goldStandard_id
    A,B = sortPair2(row['A'],row['B'])
    gsA,gsB = A.split('_')[1],B.split('_')[1]
    if gsA==gsB:
        row['A'],row['B'] = A.split('_')[0],B.split('_')[0]
        row['C'] = gs_id[gsA]
        return row[['A','B','C','corr']]

def compute_correlation(data,gs_id,llr_collected,method_redundancy,size_cutoff):
    ## Compute correlations for pairs expressed with a min nr. samples,
    ## then extract upper matrix only and exclude diagonal.
    ## As the indices of data matrix are sorted in decreasing order, the stack()
    ## function will return a data frame with sorted proteins (i.e. A>B).
    for gs_llr in llr_collected:
        time_start = time.perf_counter()
        if method_redundancy == 'spearman': 
            r_df = data[gs_llr].corr(method="spearman",min_periods=size_cutoff,numeric_only=True)
        elif method_redundancy == 'pearson': 
            r_df = data[gs_llr].corr(method="pearson",min_periods=size_cutoff,numeric_only=True)
        
        shape_data = len(gs_llr)
        r_df = r_df.where(np.triu(np.ones((shape_data,shape_data)), k=1).astype(bool))\
                .stack(dropna=True)\
                    .reset_index()
        if len(r_df)>0:
            r_df.columns = ['A','B','corr']
            r_df = r_df.apply(lambda row: format_correlation_rows(row,gs_id),axis=1)
            r_df.columns = ['evidenceA_id','evidenceB_id','goldStandard_id','correlation']
            # r_df['correlation'] = r_df['correlation'].fillna(0)
            ## We only store positive correlations given max(0,corr) used later on
            r_df = r_df[r_df['correlation']>0]
        
        if len(r_df)>0:
            r_df = r_df.sort_values(by=['evidenceA_id', 'evidenceB_id']).reset_index(drop=True)
            print(r_df)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Done %s in %5.1f secs' %(gs_llr[0].split('_')[1],time_elapsed))

            mem_csv = StringIO()
            r_df.to_csv(mem_csv, index=False)
            mem_csv.seek(0)
            print("Writing %i Redundancy entries to DB" %len(r_df))
            # Writing the csv to the evidenceLink table
            with closing(mem_csv) as csv_io:
                Redundancy.objects.from_csv(csv_io)
        else:
            time_elapsed = (time.perf_counter() - time_start) 
            print('### Done %s in %5.1f secs' %(gs_llr[0].split('_')[1],time_elapsed))

###########################
## LLR integration
def weightedSequentialAdditiveRedundancy4Species(LLR,LLR_perEv_inGs,r_dict,tmp_gs_index_dict,alpha):
    '''
    The function computes the weighted integration of LLR for an input evidence type.
    Output for goldstandard Complex: DOM,GIN,GRG,MEX,MIR,PEX,PHP,PIN,SCL,TFB|hsa,mmu,rno,cfa,...,osa;
    '''
    num_llr = len(LLR_perEv_inGs)
    if num_llr==0: return LLR
    if num_llr>1:
        sortedByLLR=sorted(LLR_perEv_inGs, key=LLR_perEv_inGs.get, reverse=True) # Sorting by value, smallest values first, just getting sorted keys as list to iterate over
        weight,sp_llr = 1,[]
        ## Sum LLRs by starting from the 2nd largest one, and down-weight it with larger correlated LLRs. 
        ## The Sequential Additive mode operates such that 2,1 | 2,1*3,2 | 2,1*3,2*4,3 | …
        if not r_dict:
            ## If there is no correlation between evidences, multiply weight by alpha
            for i,evi_id in enumerate(sortedByLLR):
                wLLRe = float(LLR_perEv_inGs[evi_id])
                wLLRe *= weight
                LLR[tmp_gs_index_dict[evi_id]] += wLLRe
                sp_llr.append(wLLRe)
                weight *= alpha
        else:
            ## Otherwise, measure distance as alpha*(1-r). 
            ## r=max(0,correlation), and is already done in measureSimilarity
            for i,evi_id in enumerate(sortedByLLR):
                wLLRe = float(LLR_perEv_inGs[evi_id])
                evj_id = sortedByLLR[i-1]
                wLLRe *= weight
                LLR[tmp_gs_index_dict[evi_id]] += wLLRe
                sp_llr.append(wLLRe)
                if not evi_id in r_dict or not evj_id in r_dict[evi_id]: weight *= alpha
                else: weight *= alpha*(1-r_dict[evi_id][evj_id])
        LLR[tmp_gs_index_dict['EV']] = np.sum(sp_llr)
    else:
        evi_id = list(LLR_perEv_inGs.keys())[0]
        uwLLRe = float(LLR_perEv_inGs[evi_id])
        LLR[tmp_gs_index_dict[evi_id]] += uwLLRe
        LLR[tmp_gs_index_dict['EV']] = uwLLRe

    return LLR

def updateMultipleLLREvidenceAndSpeciesWithRedundancy(LLR_row,r_dict,r_flag,direction_col,tmp_index_dict,first_ev,last_ev,arraySize,alpha):
    ## For testing:
    # LLR_row = speciesA_df.dropna(thresh=4).reset_index(drop=True).iloc[0,]

    ## LLRs are ranked by their absolute value in increasing order
    LLR_row = LLR_row.dropna()
    if 'LLR' in LLR_row:
        tmpLLR_row = LLR_row[3+len(direction_col)::]
        LLR = LLR_row['LLR']
    else:
        tmpLLR_row = LLR_row[2+len(direction_col)::]
        LLR = np.zeros(arraySize)
    
    if tmpLLR_row.empty: return LLR_row
    LLR_row = LLR_row[['proteinA_id','proteinB_id']+direction_col].sparse.to_dense() #using SparseData
    # LLR_row = LLR_row[['proteinA_id','proteinB_id']+direction_col]

    ## Convert row into dictionary for faster operations
    tmpLLRDict = tmpLLR_row.to_dict()
    tmpLLRDict_byGS={}
    for t in tmpLLRDict:
        eid,gs=t.split("_")
        if gs not in tmpLLRDict_byGS:
            tmpLLRDict_byGS[gs]={}
        tmpLLRDict_byGS[gs][eid]=tmpLLRDict[t]

    ## Iterate over gold standards and run the Sequential Additive redundancy schema
    ## Update LLRs, separated by evidence and species
    gs_fbs = []
    for g_type in tmpLLRDict_byGS:
        evForGs=tmpLLRDict_byGS[g_type]
        ## The Sequential Additive schema 2,1 | 2,1*3,2 | 2,1*3,2*4,3 | …
        if r_flag or not g_type in r_dict: LLR = weightedSequentialAdditiveRedundancy4Species(LLR,evForGs,{},tmp_index_dict[g_type],alpha)
        else: LLR = weightedSequentialAdditiveRedundancy4Species(LLR,evForGs,r_dict[g_type],tmp_index_dict[g_type],alpha)
        ## This is to save up space later on
        gs_fbs.append(np.sum(LLR[tmp_index_dict[g_type][first_ev]:tmp_index_dict[g_type][last_ev]+1])>0)

    LLR_row['LLR'] = np.round(LLR,3)
    LLR_row['hasTrue'] = int(any(gs_fbs))
    return LLR_row

def sumLLR4Species(LLR,LLR_perEv_inGs,tmp_gs_index_dict):
    num_llr = len(LLR_perEv_inGs)
    if num_llr==0: return LLR
    sp_llr = []
    if num_llr>0:
        for evi_id in LLR_perEv_inGs:
            wLLRe = float(LLR_perEv_inGs[evi_id])
            sp_llr.append(wLLRe)
            LLR[tmp_gs_index_dict[evi_id]] += wLLRe
    LLR[tmp_gs_index_dict['EV']] = np.sum(sp_llr)

    return LLR

def updateMultipleLLREvidenceAndSpeciesWithOrthology(LLR_row,direction_col,tmp_index_dict,first_ev,last_ev,arraySize):
    ## For testing:
    # LLR_row = speciesA_df.dropna(thresh=4).reset_index(drop=True).iloc[0,]
    
    ## LLRs are ranked by their absolute value in increasing order
    LLR_row = LLR_row.dropna()
    if 'LLR' in LLR_row:
        tmpLLR_row = LLR_row[3+len(direction_col)::]
        LLR = LLR_row['LLR']
    else:
        tmpLLR_row = LLR_row[2+len(direction_col)::]
        LLR = np.zeros(arraySize)
    
    if tmpLLR_row.empty: return LLR_row
    LLR_row = LLR_row[['proteinA_id','proteinB_id']+direction_col].sparse.to_dense() #using SparseData
    # LLR_row = LLR_row[['proteinA_id','proteinB_id']+direction_col]

    ## Convert row into dictionary for faster operations
    tmpLLRDict = tmpLLR_row.to_dict()
    tmpLLRDict_byGS={}
    for t in tmpLLRDict:
        eid,gs=t.split("_")
        if gs not in tmpLLRDict_byGS:
            tmpLLRDict_byGS[gs]={}
        tmpLLRDict_byGS[gs][eid]=tmpLLRDict[t]

    ## Iterate over gold standards and run LLR sum
    ## Update LLRs, separated by evidence and species
    gs_fbs = []
    for g_type in tmpLLRDict_byGS:
        ## Simply keep track of species
        LLR = sumLLR4Species(LLR,tmpLLRDict_byGS[g_type],tmp_index_dict[g_type])
        gs_fbs.append(np.sum(LLR[tmp_index_dict[g_type][first_ev]:tmp_index_dict[g_type][last_ev]+1])>0)
    
    LLR_row['LLR'] = np.round(LLR,3)
    LLR_row['hasTrue'] = int(any(gs_fbs))
    return LLR_row

def updateMultipleLLREvidenceAndSpecies(LLR_row,direction_col,tmp_index_dict,first_ev,last_ev,llr_collected,arraySize):
    ## For testing:
    # LLR_row = speciesA_df.dropna(thresh=4).reset_index(drop=True).iloc[0,]
    ## LLRs are ranked by their absolute value in increasing order
    LLR_row = LLR_row.dropna()
    if 'LLR' in LLR_row:
        tmpLLR_row = LLR_row[3+len(direction_col)::]
        LLR = LLR_row['LLR']
    else:
        tmpLLR_row = LLR_row[2+len(direction_col)::]
        LLR = np.zeros(arraySize)
    
    if tmpLLR_row.empty: return LLR_row
    LLR_row = LLR_row[['proteinA_id','proteinB_id']+direction_col].sparse.to_dense() #using SparseData
    # LLR_row = LLR_row[['proteinA_id','proteinB_id']+direction_col]
    
    col_names = [item for sublist in llr_collected for item in sublist if item in tmpLLR_row]
    gs_fbs = []
    for col_name in col_names:
        e_i,g_type=col_name.split("_")
        score = float(tmpLLR_row[col_name])
        LLR[tmp_index_dict[g_type][e_i]] += score
        LLR[tmp_index_dict[g_type]['EV']] = score
        gs_fbs.append(np.sum(LLR[tmp_index_dict[g_type][first_ev]:tmp_index_dict[g_type][last_ev]+1])>0)

    LLR_row['LLR'] = np.round(LLR,3)
    LLR_row['hasTrue'] = int(any(gs_fbs))
    return LLR_row

def updateSingleLLREvidenceAndSpecies(LLR_row,sp_i,et_i,tmp_index_dict,first_ev,last_ev,direction_col,llr_collected,arraySize):
    # LLR_row = speciesA_df.dropna(thresh=10).reset_index(drop=True).iloc[0,]
    ## LLRs are ranked by their absolute value in increasing order
    LLR_row = LLR_row.dropna()
    if 'LLR' in LLR_row:
        tmpLLR_row = LLR_row[3+len(direction_col)::]
        LLR = LLR_row['LLR']
    else:
        tmpLLR_row = LLR_row[2+len(direction_col)::]
        LLR = np.zeros(arraySize)
    
    if tmpLLR_row.empty: return LLR_row
    LLR_row = LLR_row[['proteinA_id','proteinB_id']+direction_col].sparse.to_dense()
    
    col_name = llr_collected[0][0]
    score = np.round(float(tmpLLR_row[col_name]),3)
    LLR[et_i] = score
    LLR[sp_i] += score

    LLR_row['LLR'] = LLR

    gs_fbs = np.sum(LLR[tmp_index_dict[first_ev]:tmp_index_dict[last_ev]+1])>0
    LLR_row['hasTrue'] = int(gs_fbs)

    return LLR_row

def updateLLRParallel(speciesA_df,e_type,evidence_list_species,evidences_redundancy,evidences_orthology,ev_order,merged_llrs,llr_collected,r_dict,r_flag,direction_col,index_dict,arraySize,alpha):
    ## Subset index_dict
    #  merged_llrs = ['752_Complex','751_Metabolic','748_Complex','749_PPI','750_Complex'] 
    
    subset_tuples = []
    for t in merged_llrs:
        ev_id,gs = t.split('_')
        sp = evidence_list_species[int(ev_id)]
        subset_tuples.append((gs,e_type,ev_id,sp))
    first_ev,last_ev = ev_order[0],ev_order[-1]
    tmp_index_dict = tmpDictIndexes(index_dict,subset_tuples,first_ev,last_ev)
    
    # if len(direction_col)>0: meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int', 'direction':'float', 'LLR': 'object'}
    # else: meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int', 'LLR': 'object'}
    meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int'}
    for dir_col in direction_col:
        meta_dict[dir_col] = meta_dict.get(dir_col,'int')
        if dir_col in speciesA_df: speciesA_df[dir_col] = speciesA_df[dir_col].fillna(0)
    meta_dict['LLR'] = meta_dict.get('LLR','object')
    meta_dict['hasTrue'] = meta_dict.get('hasTrue','int8')
    print(speciesA_df)
    if len(merged_llrs)>1:
        if e_type in evidences_redundancy:
            speciesA_ddf = dd.from_pandas(speciesA_df, npartitions=288)  # Adjust the number of partitions as needed
            result_ddf = speciesA_ddf.apply(lambda row: updateMultipleLLREvidenceAndSpeciesWithRedundancy(row,r_dict,r_flag,direction_col,tmp_index_dict,first_ev,last_ev,arraySize,alpha), axis=1, meta=meta_dict)
            with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
                speciesA_df = result_ddf.compute() 
        else:
            if e_type in evidences_orthology:
                speciesA_ddf = dd.from_pandas(speciesA_df, npartitions=288)  # Adjust the number of partitions as needed
                result_ddf = speciesA_ddf.apply(lambda row: updateMultipleLLREvidenceAndSpeciesWithOrthology(row,direction_col,tmp_index_dict,first_ev,last_ev,arraySize), axis=1, meta=meta_dict)
                with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
                    speciesA_df = result_ddf.compute()
            else:
                speciesA_ddf = dd.from_pandas(speciesA_df, npartitions=288)  # Adjust the number of partitions as needed
                result_ddf = speciesA_ddf.apply(lambda row: updateMultipleLLREvidenceAndSpecies(row,direction_col,tmp_index_dict,first_ev,last_ev,llr_collected,arraySize), axis=1, meta=meta_dict)
                with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
                    speciesA_df = result_ddf.compute() 
    else:
        col_name = merged_llrs[0]
        ev_id,g_type=col_name.split("_")
        sp_i = tmp_index_dict[g_type][ev_id]
        et_i = tmp_index_dict[g_type]['EV']
        speciesA_ddf = dd.from_pandas(speciesA_df, npartitions=288)  # Adjust the number of partitions as needed
        result_ddf = speciesA_ddf.apply(lambda row: updateSingleLLREvidenceAndSpecies(row,sp_i,et_i,tmp_index_dict[g_type],first_ev,last_ev,direction_col,llr_collected,arraySize), axis=1, meta=meta_dict)
        with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
            speciesA_df = result_ddf.compute() 
    return speciesA_df

###########################
## FBS extraction
def FBSonly(LLR_row,index_dict,gs_order,ev_order,sp_order):
    """
    Extracts final bayesian score (FBS) and formats the output for a given LLR_row.

    Parameters:
    - LLR_row (pandas.Series): Row containing LLR information for a protein pair.
    - index_dict (dict): Dictionary mapping goldstandard IDs to their corresponding indices in the LLR array.
    - direction_col (list): List of columns representing the directions (e.g., ['GRG_direction',...]).
    - gs_order (list): Order of goldstandards.
    - ev_order (list): Order of evidence types.
    - sp_order (list): Order of species.

    Returns:
    - pandas.Series: Formatted output containing protein IDs, maximum PFC, PFC values, FBS for each goldstandard,
                     and LLR evidence information.

    Example:
    result = extractFBSandFormatOutput(LLR_row, index_dict, ['GRG_direction'],
                                       ['GS1', 'GS2', 'GS3'], ['EV1', 'EV2'], ['SP1', 'SP2'])
    """
    fbs,llr = [],''
    for gs in gs_order:
        if gs in index_dict:
            ev_llr = [LLR_row['LLR'][index_dict[gs][ev]] for ev in ev_order]
            sp_llr = [LLR_row['LLR'][index_dict[gs][sp]] for sp in sp_order]
            tmp_ev_llr = np.sum(ev_llr)
            fbs.append(tmp_ev_llr)
            LLR_row[gs] = np.round(tmp_ev_llr,3)
            llr += ','.join(map(str,np.round(ev_llr,3))) + '|' + ','.join(map(str,np.round(sp_llr,3))) + ';'
        else: 
            fbs.append(0)
            LLR_row[gs] = 0
            llr += ','.join(map(str,[0] * len(ev_order))) + '|' + ','.join(map(str,[0] * len(sp_order))) + ';'
                        
    fbs = np.round(fbs,3)
    LLR_row['max_fbs'] = np.max(fbs)
    return LLR_row[['proteinA_id','proteinB_id','max_fbs'] + gs_order]

def FBS_pfc(LLR_row,index_dict,direction_col,gs_order,ev_order,sp_order,pfc_constant=0.01):
    """
    Extracts final bayesian score (FBS) and formats the output for a given LLR_row.

    Parameters:
    - LLR_row (pandas.Series): Row containing LLR information for a protein pair.
    - index_dict (dict): Dictionary mapping goldstandard IDs to their corresponding indices in the LLR array.
    - direction_col (list): List of columns representing the directions (e.g., ['GRG_direction',...]).
    - gs_order (list): Order of goldstandards.
    - ev_order (list): Order of evidence types.
    - sp_order (list): Order of species.
    - pfc_constant (float, optional): Constant used in the calculation of probability of functional coupling (PFC).
                                     Defaults to 0.001.

    Returns:
    - pandas.Series: Formatted output containing protein IDs, maximum PFC, PFC values, FBS for each goldstandard,
                     and LLR evidence information.

    Example:
    result = extractFBSandFormatOutput(LLR_row, index_dict, ['GRG_direction'],
                                       ['GS1', 'GS2', 'GS3'], ['EV1', 'EV2'], ['SP1', 'SP2'])
    """
    fbs,llr = [],''
    for gs in gs_order:
        if gs in index_dict:
            ev_llr = [LLR_row['LLR'][index_dict[gs][ev]] for ev in ev_order]
            sp_llr = [LLR_row['LLR'][index_dict[gs][sp]] for sp in sp_order]
            tmp_ev_llr = np.sum(ev_llr)
            fbs.append(tmp_ev_llr)
            llr += ','.join(map(str,np.round(ev_llr,3))) + '|' + ','.join(map(str,np.round(sp_llr,3))) + ';'
        else: 
            fbs.append(0)
            llr += ','.join(map(str,[0] * len(ev_order))) + '|' + ','.join(map(str,[0] * len(sp_order))) + ';'
                        
    fbs = np.round(fbs,3)
    LLR_row['fbs_goldstandard'] = ','.join(map(str,fbs))
    LLR_row['pfc'] = np.round(1/(1+np.exp(-np.log(pfc_constant)-fbs)),3)
    LLR_row['max_pfc'] = np.max(LLR_row['pfc'])
    LLR_row['pfc'] = ','.join(map(str,LLR_row['pfc']))
    LLR_row['llr_evidence'] = llr[:-1]
    return LLR_row[['proteinA_id','proteinB_id','max_pfc','pfc','fbs_goldstandard','llr_evidence'] + direction_col]

def FBS(LLR_row,index_dict,direction_col,gs_order,ev_order,sp_order):
    """
    Extracts final bayesian score (FBS) and formats the output for a given LLR_row.

    Parameters:
    - LLR_row (pandas.Series): Row containing LLR information for a protein pair.
    - index_dict (dict): Dictionary mapping goldstandard IDs to their corresponding indices in the LLR array.
    - direction_col (list): List of columns representing the directions (e.g., ['GRG_direction',...]).
    - gs_order (list): Order of goldstandards.
    - ev_order (list): Order of evidence types.
    - sp_order (list): Order of species.

    Returns:
    - pandas.Series: Formatted output containing protein IDs, maximum PFC, PFC values, FBS for each goldstandard,
                     and LLR evidence information.

    Example:
    result = extractFBSandFormatOutput(LLR_row, index_dict, ['GRG_direction'],
                                       ['GS1', 'GS2', 'GS3'], ['EV1', 'EV2'], ['SP1', 'SP2'])
    """
    fbs,llr = [],''
    for gs in gs_order:
        if gs in index_dict:
            ev_llr = [LLR_row['LLR'][index_dict[gs][ev]] for ev in ev_order]
            sp_llr = [LLR_row['LLR'][index_dict[gs][sp]] for sp in sp_order]
            tmp_ev_llr = np.sum(ev_llr)
            fbs.append(tmp_ev_llr)
            LLR_row[gs] = np.round(tmp_ev_llr,3)
            llr += ','.join(map(str,np.round(ev_llr,3))) + '|' + ','.join(map(str,np.round(sp_llr,3))) + ';'
        else: 
            fbs.append(0)
            LLR_row[gs] = 0
            llr += ','.join(map(str,[0] * len(ev_order))) + '|' + ','.join(map(str,[0] * len(sp_order))) + ';'
                        
    fbs = np.round(fbs,3)
    LLR_row['fbs_goldstandard'] = ','.join(map(str,fbs))
    LLR_row['llr_evidence'] = llr[:-1]
    return LLR_row[['proteinA_id','proteinB_id','fbs_goldstandard','llr_evidence'] + gs_order + direction_col]

def computeFBS(speciesA_df,index_dict,direction_col,gs_order,ev_order,sp_order):
    ## The order of metadata works if GRG is not first evidence processed
    # meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int','max_fbs': 'float'}
    meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int', 'fbs_goldstandard': 'object', 'llr_evidence': 'object'}
    for gs_col in gs_order:
        meta_dict[gs_col] = meta_dict.get(gs_col,'float')
    for dir_col in direction_col:
        meta_dict[dir_col] = meta_dict.get(dir_col,'int')
    speciesA_df = speciesA_df.sample(frac=1).reset_index(drop=True)
    speciesA_ddf = dd.from_pandas(speciesA_df, npartitions=288)  # Adjust the number of partitions as needed
    # result_ddf = speciesA_ddf.apply(lambda row: FBSonly(row,index_dict,gs_order,ev_order,sp_order), axis=1, meta=meta_dict)
    result_ddf = speciesA_ddf.apply(lambda row: FBS(row,index_dict,direction_col,gs_order,ev_order,sp_order), axis=1, meta=meta_dict)
    with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
        speciesA_df = result_ddf.compute()

    return speciesA_df

def FBStoPPVgoldStandard(tmp_speciesA_evidence_df,g_type,goldstandard_t,instanceConfig,goldStandardConfig,trainingConfig):
    # Negative gold standard
    # tmp_speciesA_evidence_df=speciesA_df[prefix_col+[g_type]].dropna(how='any',axis=0)
    # trainingConfig = trainingConfig['Network']
    min_n_sample = trainingConfig['nSample'][0]
    speciesA = goldstandard_t.species.tax_id
    ## Only positive values
    tmp_speciesA_evidence_df = tmp_speciesA_evidence_df[tmp_speciesA_evidence_df[g_type]>0]
    # tmp_speciesA_evidence_df = tmp_speciesA_evidence_df[tmp_speciesA_evidence_df[g_type]>1]

    positive_goldStandardLink,negative_goldStandardLink = separateNegativePositiveGoldStandard(tmp_speciesA_evidence_df,goldstandard_t,g_type,instanceConfig)
    if positive_goldStandardLink is None or negative_goldStandardLink is None: 
        write_logPPVtoDatabase(goldstandard_t,['Nan'],motivation='|NGS|<%i' %min_n_sample)
        return
    lpos,lneg = len(positive_goldStandardLink),len(negative_goldStandardLink)
    # np.random.seed(int(speciesA))
    # negative_goldStandardLink = random.sample(negative_goldStandardLink,np.min([int(1e6),len(negative_goldStandardLink)]))
    # negative_goldStandardLink = random.sample(negative_goldStandardLink,np.min([lpos,lneg]))
    ## Randomly sample 6 times more NGS than PGS
    # negative_goldStandardLink = random.sample(negative_goldStandardLink,np.min([lpos,lneg]))
    # lneg = len(negative_goldStandardLink)
    print('Nr. of PGS:\t%i\nNr. of NGS:\t%i' %(lpos,lneg))
    logisticPPV(positive_goldStandardLink,negative_goldStandardLink,goldstandard_t,speciesA,g_type,trainingConfig,saveToDatabase=True)

def FBStoPPV(FBS,logPPV):
    coefficients = logPPV.function
    if 'Nan' not in coefficients:
        coefficients = [float(i) for i in coefficients.split(',')]
        print(coefficients)
        a, b, c = coefficients
        return np.round(np.clip(logistic(np.clip(FBS, 0, np.max(FBS)),a, b, c),0.5,1),3)
    return 0

def computePPV(speciesA_df,speciesA_goldstandard_list,prefix_col,speciesA,gs_order,instanceConfig,goldStandardConfig,trainingConfig):
    for goldstandard_t in speciesA_goldstandard_list:
        g_type = goldstandard_t.type
        time_start = time.perf_counter()
        ## Compute PPV
        logPPV = LogisticPPV.objects.filter(Q(goldStandard=goldstandard_t))
        if not logPPV.exists():
            time_start = time.perf_counter()
            FBStoPPVgoldStandard(speciesA_df[prefix_col+[g_type]].dropna(how='any',axis=0),g_type,goldstandard_t,instanceConfig,goldStandardConfig,trainingConfig)
            time_elapsed = (time.perf_counter() - time_start)
            print('### %s PPV extracted in %5.1f secs' %(g_type,time_elapsed))
            logPPV = LogisticPPV.objects.filter(Q(goldStandard=goldstandard_t))
        
        if 'Nan' not in logPPV[0].function:
            speciesA_df[g_type] = FBStoPPV(speciesA_df[g_type],logPPV[0])
        else:
            print('Error is: ' + str(logPPV[0].error))
            speciesA_df[g_type] = 0
        ## To test PFC, comment the previous code
        # speciesA_df[g_type] = np.round(1/(1+np.exp(-np.log(0.01)-speciesA_df[g_type])),3)
        time_elapsed = (time.perf_counter() - time_start)
        print('### %s PPV extracted in %5.1f secs' %(g_type,time_elapsed))
        
    # Extract max_ppv
    time_start = time.perf_counter()
    speciesA_df['max_ppv'] = speciesA_df[gs_order].max(axis=1)    
    speciesA_df['ppv'] = speciesA_df[gs_order].astype(str).agg(','.join, axis=1)
    speciesA_df = speciesA_df.drop(gs_order,axis=1)
    time_elapsed = (time.perf_counter() - time_start)
    print('### maxPPV extracted in %5.1f secs' %time_elapsed)
    print(speciesA_df)
    return speciesA_df

###########################
## Gold Standard Links
def replaceEmptyValue(df,gs_order,ev_order,sp_order):
    ## Otherwise do this once at the beginning of evidence calculation
    llr = ''.join([','.join(map(str,[0] * len(ev_order))) + '|' + ','.join(map(str,[0] * len(sp_order))) + ';'] * len(gs_order))      
    df['llr_evidence'] = df['llr_evidence'].fillna(llr[:-1])
    df['max_ppv'] = df['max_ppv'].fillna(0)
    df['ppv'] = df['ppv'].fillna(','.join(map(str,[0] * len(gs_order))))
    df['fbs_goldstandard'] = df['fbs_goldstandard'].fillna(','.join(map(str,[0] * len(gs_order))))
    df['isGoldStandard'] = df['isGoldStandard'].fillna('0'*len(gs_order))
    df['isPPVgoldstandard'] = df['isPPVgoldstandard'].fillna('0'*len(gs_order))
    return df

def gsPPV(nr_gold_standard,min_ppv):
    # This function takes in input GoldStandard links and 
    # assigns goldstandardPPV based on NrGoldStandard
    supplement = (1-min_ppv)/3
    if nr_gold_standard == 1: return min_ppv
    elif nr_gold_standard == 2: return np.round(min_ppv + supplement,3)
    elif nr_gold_standard == 3: return np.round(min_ppv + 2*supplement,3)
    else: return 1

def scoreGoldStandardLink(speciesA_goldstandard_row,min_ppv,direction_col):
    isGS = speciesA_goldstandard_row['isGoldStandard']
    nrGS,indexGS = 0,[]
    for i,val in enumerate(isGS):
        if val!='0': 
            nrGS += 1 
            indexGS.append(i)

    # Apply the function to create the 'inferred_ppv' column
    goldstandardPPV = gsPPV(nrGS,min_ppv)
    ppv = [float(i) for i in speciesA_goldstandard_row['ppv'].split(',')]
    gsPpv = ''
    for i,val in enumerate(ppv):
        if i in indexGS:
            if val<goldstandardPPV: 
                ppv[i] = goldstandardPPV
                gsPpv += '1'
            else: gsPpv += '0'
        else: gsPpv += '0'
        
    speciesA_goldstandard_row['max_ppv'] = np.max(ppv)
    speciesA_goldstandard_row['ppv'] = ','.join(map(str,ppv))
    speciesA_goldstandard_row['isPPVgoldstandard'] = gsPpv

    return speciesA_goldstandard_row[['proteinA_id','proteinB_id','isGoldStandard','isPPVgoldstandard','fbs_goldstandard','llr_evidence','max_ppv','ppv'] + direction_col]

def reformatGoldStandardLink(speciesA_goldstandard_df,speciesA_gs_missing,gs_order,ev_order,sp_order,min_ppv):
    # Fill 0 values to non existing goldstandard
    for i,val in enumerate(speciesA_gs_missing):
        if val == -1 and gs_order[i] not in speciesA_goldstandard_df.columns: 
            speciesA_goldstandard_df.insert(i+2, gs_order[i], 0)
    
    # Cast to int
    speciesA_goldstandard_df[gs_order] = speciesA_goldstandard_df[gs_order].fillna(0).astype(int)
    speciesA_goldstandard_df = speciesA_goldstandard_df.sort_values(by=['max_ppv'],ascending=False).reset_index(drop=True)
    speciesA_goldstandard_df = speciesA_goldstandard_df.drop_duplicates(subset=['proteinA_id', 'proteinB_id','Regulatory'], keep='first')

    direction_col = [dir_col for dir_col in speciesA_goldstandard_df.columns if 'direction' in dir_col]
    column_dict = {'Complex': 'first',
        'Operon': 'first',
        'PPI': 'first',
        'Signaling': 'first',
        'Metabolic': 'first',
        'Regulatory': 'sum',
        'fbs_goldstandard': 'first', 'llr_evidence': 'first',
        'max_ppv': 'first', 'ppv': 'first'
        }
    for dir_col in direction_col:
        column_dict[dir_col] = 'first'

    speciesA_goldstandard_df = \
        speciesA_goldstandard_df\
            .groupby(['proteinA_id', 'proteinB_id'])\
                .agg(column_dict).reset_index()

    # set(speciesA_goldstandard_df['Regulatory'].to_list())
    speciesA_goldstandard_df[gs_order] = speciesA_goldstandard_df[gs_order].astype(str).fillna('0') 

    # Create a new 'isGoldStandard' column by concatenating the values of the specified columns
    speciesA_goldstandard_df['isGoldStandard'] = speciesA_goldstandard_df[gs_order].astype(str).agg(''.join, axis=1)
    # Create a new 'isPPVgoldstandard' column by 
    speciesA_goldstandard_df['isPPVgoldstandard'] = '0'*len(gs_order)
    speciesA_goldstandard_df = replaceEmptyValue(speciesA_goldstandard_df,gs_order,ev_order,sp_order)
    for dir_col in direction_col:
        speciesA_goldstandard_df[dir_col] = speciesA_goldstandard_df[dir_col].fillna(0).astype(int)

    # Drop the individual columns
    speciesA_goldstandard_df = speciesA_goldstandard_df.drop(gs_order, axis=1)

    # Assign goldstandardPPV to GS links based on FBS distributions
    speciesA_goldstandard_df = speciesA_goldstandard_df[['proteinA_id', 'proteinB_id','isGoldStandard','isPPVgoldstandard','fbs_goldstandard','llr_evidence','max_ppv','ppv']+direction_col]
    meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int','isGoldStandard': 'object', 'isPPVgoldstandard': 'object', 'fbs_goldstandard': 'object', 'llr_evidence': 'object', 'max_ppv': 'float', 'ppv': 'object'}
    for dir_col in direction_col:
        meta_dict[dir_col] = meta_dict.get(dir_col,'int')
    speciesA_goldstandard_df = speciesA_goldstandard_df.sample(frac=1).reset_index(drop=True)
    speciesA_goldstandard_ddf = dd.from_pandas(speciesA_goldstandard_df, npartitions=288)  # Adjust the number of partitions as needed
    result_ddf = speciesA_goldstandard_ddf.apply(lambda row: scoreGoldStandardLink(row,min_ppv,direction_col), axis=1, meta=meta_dict)
    with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
        speciesA_goldstandard_df = result_ddf.compute()

    return speciesA_goldstandard_df

###########################
## Database and tsv output
def write_network_to_database(speciesA_df,speciesA):
    # Creating an in-memory csv for the data
    mem_csv = StringIO()
    speciesA_df.to_csv(mem_csv, index=False)
    mem_csv.seek(0)
    # Writing the csv to the evidenceLink table
    print("Writing %i links to DB for %s" %(len(speciesA_df),speciesA))
    with closing(mem_csv) as csv_io:
        Network.objects.from_csv(csv_io)

def expand_col_to_tsv(network_row,direction_col,direction_dir,gs_order,gs_columns,ev_columns,sp_columns):
    ## Extract max ppv
    max_ppv_list = [float(i) for i in network_row['ppv'].split(',')]
    max_fbs_list = [float(fbs) for i,fbs in enumerate(network_row['fbs_goldstandard'].split(','))]

    # Using list comprehension to find the index(es) of the maximum value(s)
    max_index_ppv_fbs = [(index,ppv,fbs) for index,ppv,fbs in zip(range(len(max_ppv_list)),max_ppv_list,max_fbs_list)]
    max_index_ppv_fbs = sorted(max_index_ppv_fbs, key=lambda x: (x[1], x[2]), reverse=True)
    max_ppv_fbs_index = max_index_ppv_fbs[0][0]
    max_ppv = max_index_ppv_fbs[0][1]
    max_fbs = max_index_ppv_fbs[0][2]

    ## In the TSV output we write combo of GS network with max PPV, e.g. CMP (Complex, Metabolic and PPI with max PPV)
    max_ppv_index = [index for index,ppv in enumerate(max_ppv_list) if ppv==max_ppv]
    max_gs_list = [gs_order[i][0].upper() for i in max_ppv_index]
    max_gs = ''.join(max_gs_list)

    # Split the values in the "llr_evidence" column by ';'
    # Expand the resulting list by '|' and take only one of the five lists (e.g., the first list)    
    split_llr_evidence = network_row['llr_evidence'].split(';')[max_ppv_fbs_index].split('|')[0].split(',')
    split_llr_species = network_row['llr_evidence'].split(';')[max_ppv_fbs_index].split('|')[1].split(',')
    split_gs = network_row['fbs_goldstandard'].split(',')

    # Create separate DataFrames for the lists
    data_dict = {'ProteinA': network_row['proteinA'],
                'ProteinB': network_row['proteinB'],
                'FunCoupAid': network_row['FCAid'],
                'FunCoupBid': network_row['FCBid'],
                'GoldStandard': max_gs,
                'PPV': network_row['max_ppv'],
                'FBS_max': max_fbs,
                'isGoldStandard': network_row['isGoldStandard'][max_ppv_fbs_index],
                'isPPVgoldstandard': network_row['isPPVgoldstandard'][max_ppv_fbs_index]}

    for dir_col in direction_col:
        data_dict[dir_col] = direction_dir[network_row[dir_col][max_ppv_fbs_index]]
    
    for col_name,col_val in zip(gs_columns,split_gs):
        data_dict[col_name] = col_val
    for col_name,col_val in zip(ev_columns,split_llr_evidence):
        data_dict[col_name] = col_val
    for col_name,col_val in zip(sp_columns,split_llr_species):
        data_dict[col_name] = col_val

    return pd.Series(data_dict)

def write_network_to_tsv(speciesA_df,speciesA,direction_col,gs_order,ev_order,sp_order):
    sp_name_order = []
    for sp in sp_order:
        sp_name = Species.objects.get(Q(tax_id=sp)).species_name
        sp_name_clean = sp_name.split(' (')[0]
        sp_name_compact = '%s%s' %(sp_name_clean.split()[0][0],sp_name_clean.split()[-1][0:2].upper())
        sp_name_order.append(sp_name_compact)

    # Prepare input
    direction_dir = {'0': '--', '1': '--', '2': '-->', '3': '<--'}
    meta_dict = {'ProteinA':'object','ProteinB':'object','FunCoupAid':'int','FunCoupBid':'int','GoldStandard':'object','PPV':'float','FBS_max':'float','isGoldStandard':'object','isPPVgoldstandard':'object'}
    for col_name in direction_col:
        meta_dict[col_name] = meta_dict.get(col_name,'object') 
    gs_columns = [f'FBS_{gs}' for gs in gs_order]
    ev_columns = [f'LLR_{ev}' for ev in ev_order]
    sp_columns = [f'LLR_{sp}' for sp in sp_name_order]
    for col_name in gs_columns + ev_columns + sp_columns:
        meta_dict[col_name] = meta_dict.get(col_name,'float') 

    speciesA_ddf = dd.from_pandas(speciesA_df, npartitions=288)  # Adjust the number of partitions as needed
    result_ddf = speciesA_ddf.apply(lambda row: expand_col_to_tsv(row,direction_col,direction_dir,gs_order,gs_columns,ev_columns,sp_columns), axis=1, meta=meta_dict)
    with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
        speciesA_df = result_ddf.compute()
        
    ## Sort by descending max_ppv
    # speciesA_df = speciesA_df.sort_values(by=['FBS_max'],ascending=False).reset_index(drop=True)
    speciesA_df = speciesA_df.sort_values(by=['PPV', 'FBS_max'],ascending=False).reset_index(drop=True)
    ## Rename columns
    col_names = speciesA_df.columns.to_list()
    speciesA_df.columns = [str(i) + ':' + col_name for i,col_name in enumerate(col_names)]

    ## Export full network
    f = 'website/static/website/networks/FunCoup6.0/FC6.0_%s_full.gz' %speciesA
    speciesA_df.to_csv(f, sep="\t", compression='gzip', index=False)
    ## Export compact network
    f = 'website/static/website/networks/FunCoup6.0/FC6.0_%s_compact.gz' %speciesA
    # f = 'website/static/website/networks/FunCoup6.0/FC6.0_%s_PFCmax_compact.gz' %speciesA
    speciesA_df.iloc[:,0:(9+len(direction_col))].to_csv(f, sep="\t", compression='gzip', index=False)

###########################
## Training framework
###########################
def unpackParams(instanceConfig,trainingConfig):
    ### KDE config
    visualize_LLR = trainingConfig['KDE']['VisualizeLLR']
    if visualize_LLR: 
        for tx in instanceConfig['instance']['species']:
            if not os.path.exists('data/tmp/LLR/%s' %tx): os.makedirs('data/tmp/LLR/%s' %tx)
    min_n_sample = trainingConfig['KDE']['nSample'][0]

    ## Redundancy config
    evidences_redundancy = trainingConfig['Redundancy']['Evidences']
    method_redundancy = trainingConfig['Redundancy']['Method']
    size_redundancy = trainingConfig['Redundancy']['MinSize']
    alpha = trainingConfig['Redundancy']['Alpha']

    ## Orthology config
    evidences_orthology = trainingConfig['OrthologyTransfer']['Evidences']

    ## Network config
    min_ppv = trainingConfig['Network']['MinPPV']
    visualize_PPV = trainingConfig['Network']['VisualizePPV']
    if visualize_PPV: 
        for tx in instanceConfig['instance']['species']:
            if not os.path.exists('data/tmp/PPV/%s' %tx): os.makedirs('data/tmp/PPV/%s' %tx)
    gs_order = trainingConfig['Network']['GoldStandard_order']
    ev_order = trainingConfig['Network']['Evidence_order']
    sp_order = trainingConfig['Network']['Species_order']
    gs_directed = trainingConfig['Network']['GoldStandard_directed']
    ev_directed = trainingConfig['Network']['Evidence_directed']

    params = [visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
            evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed]
    
    return params

def getGoldStandardInfo(speciesA,goldstandard_list,params,parallel):
    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params

    ## GoldStandard
    ## Subset species-specific gs info
    speciesA_goldstandard_list = [gt for gt in goldstandard_list if gt.species.tax_id==speciesA]
    speciesA_gs_order = [gs for gs in gs_order if gs in [gt.type for gt in speciesA_goldstandard_list]]
    speciesA_gs_missing = [0 if gs in speciesA_gs_order else -1 for gs in gs_order]
    speciesA_gs_directed = [gs_directed[i] for i in [gs_order.index(gs) for gs in speciesA_gs_order]]

    # For each gold standard
    with tqdm_joblib(tqdm(desc='Read GS', total=len(speciesA_goldstandard_list))) as progress_bar:
        goldstandard_df_list = parallel(delayed(readGoldStandard)(goldstandard_t,speciesA_gs_order,speciesA_gs_directed) for goldstandard_t in speciesA_goldstandard_list)
    speciesA_goldstandard_df = reduce(lambda df_left,df_right: pd.merge(df_left, df_right, 
                        on=['proteinA_id','proteinB_id'], how='outer'), 
                    [i for i in goldstandard_df_list if i is not None])
    speciesA_goldstandard_dict_count = dict([(g_type,speciesA_goldstandard_df[g_type].count()) for g_type in speciesA_goldstandard_df.columns if g_type in gs_order])
    print(speciesA_goldstandard_df)

    ## Update goldstandard info after reading
    speciesA_goldstandard_list = [gt for gt in speciesA_goldstandard_list if (gt.type in speciesA_goldstandard_df.columns and speciesA_goldstandard_dict_count[gt.type]>=min_n_sample)]
    if len(speciesA_goldstandard_list)==0: 
        print('### No network for %s' %speciesA)
        return None,None,None,None
    speciesA_gs_order = [gs for gs in gs_order if gs in [gt.type for gt in speciesA_goldstandard_list]]
    speciesA_gs_missing = [0 if gs in speciesA_gs_order else -1 for gs in gs_order]

    return [speciesA_goldstandard_df,speciesA_goldstandard_list,speciesA_gs_order,speciesA_gs_missing]

def getLLR(speciesA,goldStandardInfo,evidence_list,params,instanceConfig,goldStandardConfig,trainingConfig):
    '''
    For each evidence type (e.g. MEX, QMS, TFB, …)
    |     Read all evidence datasets (e.g. [A, B, MEX_hsa1, MEX_hsa2, …, MEX_mmu1] ) into memory
    |     For each species
    |     |     For each evidence in memory
    |     |     |     Extract orthologs if species_evidence != species_goldstandard
    |     |     For each gold standard
    |     |     |     Train positive/negative gold standard, store LLR Polynomial function (for quality check)
    |     |     |     Compute LLR of all raw score column-wise
    |     |     |     Perform weighted redundancy integration row-wise 
    |     |     |     Save total LLR_evidence_goldstandard
    |     |     Save FBS as sum(LLR_ev1,LLR_ev2,…,LLR_evN)
    |     |     Write to database
    ''' 
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)

    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params

    speciesA_goldstandard_df,speciesA_goldstandard_list,speciesA_gs_order,speciesA_gs_missing = goldStandardInfo

    evidence_list_species = dict([(et.id,et.species.tax_id) for et in evidence_list])

    ev_training_order = ['GRG','PIN','DOM','GIN','PEX','SCL','MIR','TFB','PHP','MEX']
    evidence_type_list = [et for et in ev_training_order if et in list({et.type for et in evidence_list})]

    speciesA_f1 = "data/tmp/network/%s_llr.pkl" % speciesA

    ## GoldStandard
    gs_id = {}
    arraySize = len(speciesA_gs_order) * (len(ev_order) + len(sp_order))
    index_dict = createDictIndexes(speciesA_gs_order,ev_order,sp_order)

    # Evidence
    if os.path.exists(speciesA_f1):
        print('### Reading pickled network for %s' %speciesA)
        speciesA_df = joblib.load(speciesA_f1)
        speciesA_df['proteinA_id'] = speciesA_df['proteinA_id'].astype(int)
        speciesA_df['proteinB_id'] = speciesA_df['proteinB_id'].astype(int)
        evidence_type_list = []
        print(speciesA_df)
    else:
        speciesA_df = pd.DataFrame()

    for e_type in evidence_type_list:
        prefix_col = ['proteinA_id','proteinB_id']
        if ev_directed[ev_order.index(e_type)]:
            prefix_col = ['proteinA_id','proteinB_id','direction']

        ## Check if orthology transfer is allowed for the evidence
        if e_type in evidences_orthology:
            tmp_evidence_t_list = [et for et in evidence_list if et.type == e_type]
        else:
            tmp_evidence_t_list = [et for et in evidence_list if (et.type == e_type and et.species.tax_id==speciesA)]
        
        ## If there is no evidence for speciesA
        if len(tmp_evidence_t_list)==0: continue
        
        ## Read and merge all evidence-specific datasets
        time_start = time.perf_counter()
        if len(tmp_evidence_t_list)>1:
            with tqdm_joblib(tqdm(desc='Read %s' %e_type, total=len(tmp_evidence_t_list))) as progress_bar:
                evidence_df_list = parallel(delayed(readEvidence)(tmp_evidence_t,prefix_col,speciesA,ev_order,ev_directed) for tmp_evidence_t in tmp_evidence_t_list)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Read %s in %5.1f secs' %(e_type,time_elapsed))
            
            ## Merge evidences with speciesA protein_id
            time_start = time.perf_counter()
            # speciesA_evidence_df = pd.DataFrame()
            speciesA_evidence_df = reduce(lambda df_left,df_right: pd.merge(df_left, df_right, 
                                on=prefix_col, how='outer'), 
                            [i for i in evidence_df_list if i is not None])
            ## Free up memory
            del evidence_df_list
            speciesA_evidence_df = speciesA_evidence_df.reset_index(drop=True)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Merged %s in %5.1f secs' %(e_type,time_elapsed))
        else:
            print('### Reading %s' %e_type)
            speciesA_evidence_df = readEvidence(tmp_evidence_t_list[0],prefix_col,speciesA,ev_order,ev_directed)
            speciesA_evidence_df = speciesA_evidence_df.reset_index(drop=True)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Read %s in %5.1f secs' %(e_type,time_elapsed))
        print(speciesA_evidence_df)
        # Set data type to SparseDtype to achieve memory efficiency?
        # speciesA_evidence_df.iloc[:, len(prefix_col)::] = speciesA_evidence_df.iloc[:, len(prefix_col)::].astype(pd.SparseDtype("float", np.nan))
        num_ev = speciesA_evidence_df.shape[1]
        if num_ev==0: continue
        else: 
            # Reorder the columns in the DataFrame
            suffix_col = [col for col in speciesA_evidence_df.columns if col not in prefix_col]
            speciesA_evidence_df = speciesA_evidence_df[prefix_col + suffix_col]

            # speciesA_evidence_df = speciesA_evidence_df[speciesA_evidence_df['proteinA_id'] != speciesA_evidence_df['proteinB_id']].reset_index(drop=True)
            speciesA_evidence_df['proteinA_id'] = speciesA_evidence_df['proteinA_id'].astype(int)
            speciesA_evidence_df['proteinB_id'] = speciesA_evidence_df['proteinB_id'].astype(int)

            if 'direction' in speciesA_evidence_df.columns: 
                speciesA_evidence_df['direction'] = speciesA_evidence_df['direction'].astype(float)

            ## GoldStandard,Ev training
            ev_collected = speciesA_evidence_df.columns[len(prefix_col):num_ev].to_list()
            llrs_converted = []

            for goldstandard_t in speciesA_goldstandard_list:
                g_type=goldstandard_t.type
                gs_id[g_type] = gs_id.get(g_type,goldstandard_t.id)
                # Subset positive gold standard links
                positive_goldStandardLink = speciesA_goldstandard_df[['proteinA_id','proteinB_id',goldstandard_t.type]].dropna().reset_index(drop=True)
                positive_goldStandardLink.columns.values[-1] = 'direction'
                positive_goldStandardLink['proteinA_id'] = positive_goldStandardLink['proteinA_id'].astype(int)
                positive_goldStandardLink['proteinB_id'] = positive_goldStandardLink['proteinB_id'].astype(int)
                
                ### Parallelize over evidences id
                trained_e_id = [e_id for e_id in ev_collected if PolynomialLLR.objects.filter(Q(evidence_id=e_id) & Q(goldStandard=goldstandard_t)).exists()]
                toTrain_e_id = [e_id for e_id in ev_collected if not e_id in trained_e_id]
                if len(toTrain_e_id)>0:
                    time_start = time.perf_counter()
                    with tqdm_joblib(tqdm(desc='Training %s' %g_type, total=len(toTrain_e_id))) as progress_bar:
                        _ = parallel(delayed(trainGoldStandard)(positive_goldStandardLink,speciesA_evidence_df[prefix_col+[e_id]].dropna(how='any',axis=0),e_id,goldstandard_t,ev_directed,gs_directed,ev_order,gs_order,e_type,g_type,instanceConfig,goldStandardConfig,trainingConfig['KDE'],saveToDatabase=True) \
                                            for e_id in toTrain_e_id)
                    time_elapsed = (time.perf_counter() - time_start)
                    print('### Training %s_%s in %5.1f secs' %(g_type,e_type,time_elapsed))
                
                ### Collect all ev,gs combos
                time_start = time.perf_counter()
                llr_converted = speciesA_evidence_df.apply(lambda column: convertToLLR_parallel(column, goldstandard_t, prefix_col), axis=0)
                llrs_converted.append(llr_converted)
                time_elapsed = (time.perf_counter() - time_start)
                print('### Converting raw score to LLR for %s_%s in %5.1f secs' %(g_type,e_type,time_elapsed))

            #################################################
            # Extract results and update DataFrame
            speciesA_evidence_df = speciesA_evidence_df.drop(ev_collected,axis=1)
            llr_collected,concatenated_llrs = [],[]
            for gs,_ in enumerate(llrs_converted):
                llr_gs_collected = []
                for ev_gs in llrs_converted[gs]:
                    if ev_gs is not None:
                        concatenated_llrs.append(ev_gs)
                        llr_gs_collected.append(ev_gs.name)
                if len(llr_gs_collected)>0: 
                    llr_collected.append(llr_gs_collected)
                        
            # Merge the sublists into one flat list
            merged_llrs = [item for sublist in llr_collected for item in sublist]
            if len(merged_llrs)>0: 
                concatenated_llrs = pd.concat(concatenated_llrs,axis=1)
                if e_type in evidences_redundancy:
                    print('### Extracting Redundancy factors for %s' %e_type)
                    time_start = time.perf_counter()
                    # Extract a set of integers by splitting on underscores
                    ev_list = set(int(i.split('_')[0]) for i in merged_llrs)
                    precomputed_r = Redundancy.objects.filter(Q(goldStandard__in=speciesA_goldstandard_list) & Q(evidenceA__in=ev_list) & Q(evidenceB__in=ev_list))
                    if not precomputed_r.exists():
                        if len(concatenated_llrs)>int(50e6):
                            random.seed(speciesA)
                            compute_correlation(concatenated_llrs.sample(frac=0.25),gs_id,llr_collected,method_redundancy,size_redundancy)
                        else:
                            compute_correlation(concatenated_llrs,gs_id,llr_collected,method_redundancy,size_redundancy)
                    time_elapsed = (time.perf_counter() - time_start)
                    print('### Done in %5.1f secs' %time_elapsed)

                speciesA_evidence_df = pd.concat([speciesA_evidence_df,concatenated_llrs],axis=1)
                # print(speciesA_evidence_df.dtypes)
                # If you want to keep the filtered rows in the original dataframe
                if 'direction' in prefix_col:
                    ev_direction_col = e_type + '_direction'
                    prefix_col = ['proteinA_id','proteinB_id',ev_direction_col]
                    speciesA_evidence_df.columns.values[2] = ev_direction_col
                
                speciesA_evidence_df = speciesA_evidence_df[prefix_col+merged_llrs]

                ## Drop links with no ev,gs successfull combo
                speciesA_evidence_df = speciesA_evidence_df.loc[speciesA_evidence_df.dropna(subset=merged_llrs, how='all').index]

                # Perform integration of LLR with weighted redundancy schema
                direction_col = [dir_col for dir_col in speciesA_evidence_df.columns if 'direction' in dir_col]
                r_flag,r_dict = False,{}
                if e_type in evidences_redundancy:
                    # Extract a set of integers by splitting on underscores
                    ev_list = set(int(i.split('_')[0]) for i in merged_llrs)
                    r_df = pd.DataFrame.from_records(\
                        Redundancy.objects.filter(Q(goldStandard__in=speciesA_goldstandard_list) & Q(evidenceA__in=ev_list) & Q(evidenceB__in=ev_list) & Q(correlation__gt=0))\
                            .values('goldStandard','evidenceA','evidenceB','correlation'))
                    r_df = r_df.dropna(subset=['correlation'])
                    if r_df.empty: r_flag = True
                    else: 
                        ## Convert dataframe into dictionary to speed up computations
                        for _,row in r_df.iterrows():
                            gold_standard = int(row['goldStandard'])
                            r_dict[gold_standard] = r_dict.get(gold_standard,{})
                            evidence_a = int(row['evidenceA'])
                            evidence_b = int(row['evidenceB'])
                            correlation = row['correlation']
                            r_dict[gold_standard][evidence_a] = r_dict[gold_standard].get(evidence_a,{})
                            r_dict[gold_standard][evidence_a][evidence_b] = correlation
                            r_dict[gold_standard][evidence_b] = r_dict[gold_standard].get(evidence_b,{})
                            r_dict[gold_standard][evidence_b][evidence_a] = correlation
                    print('### Computing weighted redundancy LLR integration')
                else: print('### Computing LLR integration')

                ## We shuffle the dataframe to homogenize amount of work
                speciesA_evidence_df = speciesA_evidence_df.sample(frac=1).reset_index(drop=True)
                time_start = time.perf_counter()
                speciesA_evidence_df = updateLLRParallel(speciesA_evidence_df,e_type,evidence_list_species,evidences_redundancy,evidences_orthology,ev_order,merged_llrs,llr_collected,r_dict,r_flag,direction_col,index_dict,arraySize,alpha)
                speciesA_evidence_df['proteinA_id'] = speciesA_evidence_df['proteinA_id'].astype(int)
                speciesA_evidence_df['proteinB_id'] = speciesA_evidence_df['proteinB_id'].astype(int)
                time_elapsed = (time.perf_counter() - time_start)
                print('### Integrated %s LLR in %5.1f secs' %(e_type,time_elapsed))
                
                if speciesA_df.empty:
                    speciesA_df = speciesA_evidence_df.copy()
                    speciesA_df['proteinA_id'] = speciesA_df['proteinA_id'].astype(int)
                    speciesA_df['proteinB_id'] = speciesA_df['proteinB_id'].astype(int)
                    del speciesA_evidence_df
                else:
                    time_start = time.perf_counter()
                    speciesA_df = pd.merge(speciesA_df,speciesA_evidence_df, on=['proteinA_id','proteinB_id'],how='outer')
                    time_elapsed = (time.perf_counter() - time_start)
                    print('### Merged in %5.1f secs' %time_elapsed)
                    del speciesA_evidence_df
                    time_start = time.perf_counter()
                    speciesA_df['LLR'] = speciesA_df[["LLR_x", "LLR_y"]].sum(axis=1)
                    speciesA_df = speciesA_df.drop(['LLR_x','LLR_y'],axis=1)
                    speciesA_df['hasTrue'] = speciesA_df[["hasTrue_x", "hasTrue_y"]].sum(axis=1)
                    speciesA_df = speciesA_df.drop(["hasTrue_x", "hasTrue_y"],axis=1)
                    time_elapsed = (time.perf_counter() - time_start)
                    print('### LLR summed up in %5.1f secs' %time_elapsed)
                
                ########################################
                ## Memory usage of pandas dataframe
                speciesA_df.info(verbose = False, memory_usage = 'deep')
                ## Dataframe sparsity: % of values that have not been “compressed”
                # if 'LLR' in speciesA_df.columns: print(speciesA_df.iloc[:,len(prefix_col)+1::].sparse.density)
                # else: print(speciesA_df.iloc[:,len(prefxix_col)::].sparse.density)
                ########################################
                print(speciesA_df)
            else:
                del speciesA_evidence_df 
                continue

            # Intermediate step, save to pkl. Change to save more intermediary steps.
            # if len(speciesA_df)>0 and e_type in ['TFB']: 
            #     print('Dumping speciesA_df into disk after %s' %e_type)
            #     joblib.dump(speciesA_df, "data/tmp/network/%s_llr_noMEX.pkl" % speciesA)
            if len(speciesA_df)>0 and e_type in ['MEX']: 
                if speciesA != '9606':
                    speciesA_df = speciesA_df[speciesA_df['hasTrue']>0].reset_index(drop=True)
                    speciesA_df = speciesA_df.drop(['hasTrue'],axis=1)
                    print(speciesA_df)
                print('Dumping speciesA_df into disk after %s' %e_type)
                joblib.dump(speciesA_df, "data/tmp/network/%s_llr.pkl" % speciesA)
    
    return speciesA_df

def getFBS(speciesA_df,speciesA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig):
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    
    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    speciesA_goldstandard_df,speciesA_goldstandard_list,speciesA_gs_order,speciesA_gs_missing = goldStandardInfo
    index_dict = createDictIndexes(speciesA_gs_order,ev_order,sp_order)

    #######################################
    #######################################
    prefix_col = ['proteinA_id','proteinB_id']
    direction_col = [dir_col for dir_col in speciesA_df.columns if 'direction' in dir_col]
    for dir_col in direction_col:
        prefix_col.append(dir_col)
        speciesA_df[dir_col] = speciesA_df[dir_col].fillna(0).astype(int)

    #######################################
    ## FBS extraction
    #######################################
    time_start = time.perf_counter()
    speciesA_df = computeFBS(speciesA_df,index_dict,direction_col,gs_order,ev_order,sp_order)
    time_elapsed = (time.perf_counter() - time_start)
    print('### FBS done in %5.1f secs' %time_elapsed)
    print(speciesA_df)
    # fbs_f = "data/tmp/network/%s_ONLYfbs.pkl" %speciesA 
    fbs_f = "data/tmp/network/%s_fbs.pkl" %speciesA 
    if not os.path.exists(fbs_f):
        joblib.dump(speciesA_df, fbs_f)

    return speciesA_df

def getPPV(speciesA_df,speciesA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig):
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    
    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    speciesA_goldstandard_df,speciesA_goldstandard_list,speciesA_gs_order,speciesA_gs_missing = goldStandardInfo
    index_dict = createDictIndexes(speciesA_gs_order,ev_order,sp_order)

    #######################################
    #######################################
    prefix_col = ['proteinA_id','proteinB_id']
    direction_col = [dir_col for dir_col in speciesA_df.columns if 'direction' in dir_col]
    for dir_col in direction_col:
        prefix_col.append(dir_col)
        speciesA_df[dir_col] = speciesA_df[dir_col].fillna(0).astype(int)

    #######################################
    ## ppv extraction
    #######################################
    time_start = time.perf_counter()
    speciesA_df = computePPV(speciesA_df,speciesA_goldstandard_list,prefix_col,speciesA,gs_order,instanceConfig,goldStandardConfig,trainingConfig['Network'])
    time_elapsed = (time.perf_counter() - time_start)
    print('### PPV done in %5.1f secs' %time_elapsed)
    print(speciesA_df)
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))  # Adjust figsize as needed
    fig.suptitle(speciesA, fontsize=14)
    # Histogram of ppv_max
    axes[0].hist(speciesA_df['max_ppv'], bins=100, color='#a7c957', edgecolor='black')
    axes[0].grid(True, linestyle='--', alpha=0.7)
    axes[0].set_yscale('log')
    axes[0].set_xlabel('PPV max')
    # axes[0].set_xlabel('PFC max')
    axes[0].set_ylabel('Nr. of links')
    axes[0].set_xlim(0.5, 1)
    axes[0].set_xticks(np.arange(0.5, 1.01, 0.1))

    cdf_values, bin_edges, _ = axes[1].hist(speciesA_df['max_ppv'], bins=1000, color='#a7c957', edgecolor='black', cumulative=True, density=False, alpha=0)
    # cdf_values, bin_edges, _ = axes[1].hist(np.clip(fbs_to_pfc, 0, np.max(fbs_to_pfc)), bins=50, color='skyblue', edgecolor='black', cumulative=True, density=False, alpha=0)
    complement_cdf_values = len(speciesA_df) - cdf_values
    # Plot the complement of the CDF
    ## To align the center of the bins with the points --> bin_edges[1:]
    axes[1].bar(bin_edges[1:], complement_cdf_values, width=bin_edges[1] - bin_edges[0], color='#e63946', edgecolor='#a7c957')
    axes[1].grid(True, linestyle='--', alpha=0.7)
    axes[1].set_yscale('log')
    axes[1].set_xlabel('PPV max')
    axes[1].set_ylabel('Nr. of links (cumulative)')
    axes[1].set_xlim(0.5, 1)
    axes[1].set_xticks(np.arange(0.5, 1.01, 0.1))

    # Adjust spacing between subplots
    plt.tight_layout()
    plt.savefig('data/tmp/PPV/%s/%s_PPVmax.png' %(speciesA,speciesA))
    plt.close()

    return speciesA_df

def getDatabaseOutput(speciesA_df,speciesA,goldStandardInfo,params):
    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    speciesA_goldstandard_df,speciesA_goldstandard_list,speciesA_gs_order,speciesA_gs_missing = goldStandardInfo
    speciesA_name = Species.objects.get(Q(tax_id=speciesA)).species_name
    speciesA_name_clean = speciesA_name.split(' (')[0]
    speciesA_name_compact = '%s.%s' %(speciesA_name_clean.split()[0][0],speciesA_name_clean.split()[-1])

    #######################################
    ## Gold standard links
    #######################################
    time_start = time.perf_counter()
    speciesA_goldstandard_df = speciesA_goldstandard_df.merge(speciesA_df,on=['proteinA_id','proteinB_id'],how='left')
    speciesA_df = speciesA_df[speciesA_df['max_ppv']>=min_ppv].reset_index(drop=True)
    # speciesA_duplicated_df = speciesA_df[speciesA_df.duplicated(subset=['proteinA_id', 'proteinB_id','Regulatory'], keep=False)]
    # speciesA_df[(speciesA_df['proteinA_id']==349000) & (speciesA_df['proteinB_id']==337527)][['proteinA_id','proteinB_id','max_ppv','GRG_direction']]

    speciesA_goldstandard_df = reformatGoldStandardLink(speciesA_goldstandard_df,speciesA_gs_missing,gs_order,ev_order,sp_order,min_ppv)

    ## Remove GS links from data, then put it back with inferred PPV
    speciesA_df = pd.merge(speciesA_df,speciesA_goldstandard_df[['proteinA_id','proteinB_id']], on=['proteinA_id','proteinB_id'], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
    speciesA_df['isGoldStandard'] = '0'*len(gs_order)
    speciesA_df = pd.concat([speciesA_df,speciesA_goldstandard_df],axis=0).reset_index(drop=True)
    time_elapsed = (time.perf_counter() - time_start)
    print('### Add GS links done in %5.1f secs' %time_elapsed)

    #######################################
    ## Storing the network
    #######################################

    #######################################
    ## Replace NaN values
    speciesA_df = replaceEmptyValue(speciesA_df,gs_order,ev_order,sp_order)
    direction_col = [dir_col for dir_col in speciesA_df.columns if 'direction' in dir_col]
    for dir_col in direction_col:
        speciesA_df[dir_col] = speciesA_df[dir_col].fillna(0).astype(int)

    ### TODO: GRG might support A-->B, B-->A, what to do?!
    ### Current approach: Max out duplicated links
    # speciesA_duplicated_df = speciesA_df[speciesA_df.duplicated(subset=['proteinA_id', 'proteinB_id'], keep=False)]
    # if len(speciesA_duplicated_df)>0:
    #     # Group by 'proteinA_id' and 'proteinB_id' and get the index of rows with maximum 'max_ppv'
    #     max_ppv_indices = speciesA_duplicated_df.groupby(['proteinA_id', 'proteinB_id'])['max_ppv'].idxmax()
    #     # Select the rows with maximum 'max_ppv' using the indices
    #     speciesA_duplicated_df = speciesA_duplicated_df.loc[max_ppv_indices]
    #     speciesA_df = pd.concat([speciesA_df,speciesA_duplicated_df],axis=0).reset_index(drop=True)

    ## Sort by max_ppv
    speciesA_df = speciesA_df.sort_values(by=['max_ppv'],ascending=False).reset_index(drop=True)
    speciesA_df = speciesA_df.drop_duplicates(subset=['proteinA_id', 'proteinB_id'], keep='first')

    #######################################
    ## Remove self-loops
    speciesA_df = speciesA_df[speciesA_df['proteinA_id']!=speciesA_df['proteinB_id']].reset_index(drop=True)

    ## Network model has a pre-set number of direction columns
    network_columns = [f.get_attname() for f in Network._meta.fields]
    network_direction_columns = [dir_col for dir_col in network_columns if 'direction' in dir_col]
    isG_id = speciesA_df.columns.to_list().index('isGoldStandard')
    for dir_col in network_direction_columns:
        if dir_col in direction_col:
            speciesA_df[dir_col] = speciesA_df[dir_col].astype(int).astype(str) * len(gs_order)
        else: 
            speciesA_df.insert(isG_id,dir_col, str(0) * len(gs_order))
            direction_col.append(dir_col)
            isG_id += 1

    direction_col = [dir_col for dir_col in speciesA_df.columns if 'direction' in dir_col]
    speciesA_df['proteinA_id'] = speciesA_df['proteinA_id'].astype(int)
    speciesA_df['proteinB_id'] = speciesA_df['proteinB_id'].astype(int)
    write_network_to_database(speciesA_df,speciesA_name_compact)

    return speciesA_df

def getTSVOutput(speciesA_df,speciesA,params):
    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    speciesA_name = Species.objects.get(Q(tax_id=speciesA)).species_name
    speciesA_name_clean = speciesA_name.split(' (')[0]
    speciesA_name_compact = '%s.%s' %(speciesA_name_clean.split()[0][0],speciesA_name_clean.split()[-1])

    direction_col = [dir_col for dir_col in speciesA_df.columns if 'direction' in dir_col]

    # Extract uniprot_id
    proteome = Proteome.objects.get(Q(species__tax_id=speciesA))
    uniprot_df = pd.DataFrame.from_records(Protein.objects.filter(Q(proteome=proteome)).values('id','uniprot_id'))
    uniprot_df['id'] = uniprot_df['id'].astype(int)

    new_order = ['proteinA_id','proteinB_id','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard']
    speciesA_df = speciesA_df[new_order]

    speciesA_df.columns = ['id','proteinB_id','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard']
    speciesA_df = pd.merge(speciesA_df, uniprot_df, on='id', how='left')

    speciesA_df.columns = ['FCAid','id','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard','proteinA']
    speciesA_df = pd.merge(speciesA_df, uniprot_df, on='id', how='left')

    ## TODO check if direction_col and 'isGoldStandard' are swapped here
    speciesA_df.columns = ['FCAid','FCBid','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard','proteinA','proteinB']
    print(speciesA_df)

    write_network_to_tsv(speciesA_df,speciesA_name_compact,direction_col,gs_order,ev_order,sp_order)

    return speciesA_df

###########################
## Main function =)
###########################
def trainNetwork(instanceConfig,goldStandardConfig,evidenceConfig,proteomeConfig,trainingConfig):
    '''
    '''
    # Filter out the specific warning related to DataFrame.min
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)

    # Unpack parameters
    params = unpackParams(instanceConfig,trainingConfig)
    
    ## Extract Gold Standard and Evidence
    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id__in=instanceConfig['instance']['species'])) \
                        if goldstandard_t.type in instanceConfig['instance']['goldStandard']]

    evidence_list = [evidence_t for evidence_t in Evidence.objects.filter(Q(type__in=instanceConfig['instance']['evidence']) \
                        & Q(species__tax_id__in=instanceConfig['instance']['species'])).order_by('species_id') \
                        if (evidence_t.type in evidenceConfig and \
                            (evidence_t.scoringMethod in str(evidenceConfig[evidence_t.type]['scoring_method'])))]
    
    # For each evidence type (e.g. MEX, QMS, TFB, …)
    ## Sort species by proteome size
    sorted_species = []
    for s in instanceConfig['instance']['species']:
        proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=s))[0]
        num_protein = Protein.objects.filter(Q(proteome=proteome)).count()
        sorted_species.append((s,num_protein))
    
    # Sort the list by proteome size
    sorted_species = sorted(sorted_species, key=lambda x: x[1])
    # sorted_species = [(speciesA,num_proteins) for speciesA,num_proteins in sorted_species if speciesA in []] 
    for speciesA,num_proteins in sorted_species:
        species_time_start = time.perf_counter()

        speciesA_name = Species.objects.get(Q(tax_id=speciesA)).species_name
        speciesA_name_clean = speciesA_name.split(' (')[0]
        speciesA_name_compact = '%s.%s' %(speciesA_name_clean.split()[0][0],speciesA_name_clean.split()[-1])
        print_frame('%s, %i proteins' %(speciesA_name_compact,num_proteins))

        # If network exists in TSV format, skip all computations
        network_f = 'website/static/website/networks/FunCoup6.0/FC6.0_%s_compact.gz' %speciesA_name_compact
        if os.path.exists(network_f): 
            print('### Network already exists for %s' %speciesA)
            continue

        # If network exists in the database, then only extract TSV format
        proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=speciesA))[0]
        proteins = Protein.objects.filter(Q(proteome=proteome))
        speciesA_network = Network.objects.filter(Q(proteinA_id__in=proteins))
        if speciesA_network.exists():
            col_names = tuple([col.name for col in Network._meta.get_fields() if col.name!='id'])
            speciesA_df = pd.DataFrame.from_records(speciesA_network.values(*col_names))
            speciesA_df = speciesA_df.rename(columns={'proteinA': 'proteinA_id', 'proteinB': 'proteinB_id'})

            # Output networks
            species_output_start = time.perf_counter()
            print('### STEP 4b: TSV Output')
            # goldStandardInfo[-1] = speciesA_gs_missing
            speciesA_df = getTSVOutput(speciesA_df,speciesA,params)
            species_output_elapsed = (time.perf_counter() - species_output_start)
            print('### TSV Output network for %s in %5.1f secs' %(speciesA,species_output_elapsed))
            continue

        ## OTHERWISE
        ## Gold standard info
        species_GS_start = time.perf_counter()
        goldStandardInfo = getGoldStandardInfo(speciesA,goldstandard_list,params,parallel)
        species_GS_elapsed = (time.perf_counter() - species_GS_start)
        print('### Extracted GS for %s in %5.1f secs' %(speciesA,species_GS_elapsed))

        fbs_f = "data/tmp/network/%s_fbs.pkl" %speciesA
        speciesA_df = []
        if os.path.exists(fbs_f):
            species_FBS_start = time.perf_counter()
            print('### STEP 2: FBS')
            speciesA_df = joblib.load(fbs_f)
            species_FBS_elapsed = (time.perf_counter() - species_FBS_start)
            print('### Extracted FBS for %s in %5.1f secs' %(speciesA,species_FBS_elapsed))
        else:
            ## LLR integrated
            species_LLR_start = time.perf_counter()
            print('### STEP 1: LLR')
            ### TODO --> remember to edit getLLR 
            speciesA_df = getLLR(speciesA,goldStandardInfo,evidence_list,params,instanceConfig,goldStandardConfig,trainingConfig)
            if 'hasTrue' in speciesA_df.columns:
                speciesA_df = speciesA_df[speciesA_df['hasTrue']>0].reset_index(drop=True)
                speciesA_df = speciesA_df.drop(['hasTrue'],axis=1)
            print(speciesA_df)
            species_LLR_elapsed = (time.perf_counter() - species_LLR_start)
            print('### Extracted LLR for %s in %5.1f secs' %(speciesA,species_LLR_elapsed))
            
            # Extract FBS
            species_FBS_start = time.perf_counter()
            print('### STEP 2: FBS')
            speciesA_df = getFBS(speciesA_df,speciesA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig)
            species_FBS_elapsed = (time.perf_counter() - species_FBS_start)
            print('### Extracted FBS for %s in %5.1f secs' %(speciesA,species_FBS_elapsed))

        # Extract PPV
        species_PPV_start = time.perf_counter()
        print('### STEP 3: PPV')
        speciesA_df = getPPV(speciesA_df,speciesA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig)
        species_PPV_elapsed = (time.perf_counter() - species_PPV_start)
        print('### Extracted ppv for %s in %5.1f secs' %(speciesA,species_PPV_elapsed))

        # Output networks
        species_output_start = time.perf_counter()
        print('### STEP 4a: Database Output')
        # goldStandardInfo[-1] = speciesA_gs_missing
        speciesA_df = getDatabaseOutput(speciesA_df,speciesA,goldStandardInfo,params)
        species_output_elapsed = (time.perf_counter() - species_output_start)
        print('### Database Output network for %s in %5.1f secs' %(speciesA,species_output_elapsed))

        # Output networks
        species_output_start = time.perf_counter()
        print('### STEP 4b: TSV Output')
        # goldStandardInfo[-1] = speciesA_gs_missing
        speciesA_df = getTSVOutput(speciesA_df,speciesA,params)
        species_output_elapsed = (time.perf_counter() - species_output_start)
        print('### TSV Output network for %s in %5.1f secs' %(speciesA,species_output_elapsed))

        species_time_elapsed = (time.perf_counter() - species_time_start)
        print('### Extracted network for %s in %5.1f secs' %(speciesA,species_time_elapsed))

        del speciesA_df

    # Reset the warning filter to its default behavior (optional)
    warnings.resetwarnings()

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

    # Filter out the specific warning related to DataFrame.min
    warnings.simplefilter(action='ignore', category=FutureWarning)
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus,prefer="threads")

    # Unpack parameters
    params = unpackParams(instanceConfig,trainingConfig)
    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    ## Extract Gold Standard and Evidence
    goldstandard_list = [goldstandard_t for goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id__in=instanceConfig['instance']['species'])) \
                        if goldstandard_t.type in instanceConfig['instance']['goldStandard']]

    evidence_list = [evidence_t for evidence_t in Evidence.objects.filter(Q(type__in=instanceConfig['instance']['evidence']) \
                        & Q(species__tax_id__in=instanceConfig['instance']['species'])).order_by('species_id') \
                        if (evidence_t.type in evidenceConfig and \
                            (evidence_t.scoringMethod in str(evidenceConfig[evidence_t.type]['scoring_method'])))]    

    speciesA = '559292'
    speciesA = '243232'
    speciesA = '9606'
    speciesA = '83333'

    speciesA_name = Species.objects.get(Q(tax_id=speciesA)).species_name
    speciesA_name_compact = '%s.%s' %(speciesA_name.split()[0][0],speciesA_name.split()[1])
    # speciesA_df = pd.read_csv('website/static/website/networks/FunCoup6.0/FC6.0_%s_full.gz' %speciesA_name_compact, usecols=[2,3,8,9,10,11,12,13], sep='\t', header=0)
    # speciesA_df.columns = ['proteinA_id','proteinB_id'] + gs_order
    print_frame(speciesA_name_compact)

    ## Gold standard info
    species_GS_start = time.perf_counter()
    goldStandardInfo = getGoldStandardInfo(speciesA,goldstandard_list,params,parallel)
    species_GS_elapsed = (time.perf_counter() - species_GS_start)
    # goldStandardInfo[0].to_csv('data/tmp/%s_goldStandard.gz' %speciesA, sep="\t", compression='gzip', index=False)
    print('### Extracted GS for %s in %5.1f secs' %(speciesA,species_GS_elapsed))