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
from contextlib import closing
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
from data.dataTraining.getNetwork import *
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
###########################
## LLR extraction
def convertToLLR_parallel(scores, e_type, hostA, host_goldstandard_t, prefix_col):
    pathogen_ev_id = scores.name
    host_ev = Evidence.objects.filter(Q(type=e_type) & Q(species__tax_id=hostA))
    if pathogen_ev_id not in prefix_col and host_ev.exists():
        host_ev_id = host_ev[0].id
        pLLR = PolynomialLLR.objects.get(Q(evidence_id=host_ev_id) & Q(goldStandard=host_goldstandard_t))
        coefficients = pLLR.function
        if 'Nan' not in coefficients:
            min_score, max_score = pLLR.scoreRange.split('|')
            ## TODO: if min_score>max_score: return None ## to skip decreasing LLR functions
            norm_scores = scores.values.reshape(-1, 1)
            # Min-Max scaling to the min,max of polyLLR in host species
            minmax_scaler = MinMaxScaler(feature_range=(float(min_score), float(max_score)))
            minmax_scores = minmax_scaler.fit_transform(norm_scores)
            scores = np.round(minmax_scores.flatten(),decimals=4)
            g_type = GoldStandard.objects.get(Q(id=pLLR.goldStandard_id)).type
            coefficients = [float(i) for i in coefficients.split(',')]
            result_col = np.round(np.polyval(coefficients, np.clip(scores, float(min_score), float(max_score))), 3)
            col_name = f'{host_ev_id}_{g_type}'  
            return pd.Series(result_col,name=col_name).astype("Sparse[float]")          
            # return pd.Series(result_col,name=col_name).astype(pd.SparseDtype("float", np.nan)) # return sparse array
            # return pd.Series(result_col,name=col_name)
    else: return None

###########################
## LLR integration
def updateLLRParallel(pathogenA_df,e_type,evidence_list_host,evidences_redundancy,evidences_orthology,ev_order,merged_llrs,llr_collected,r_dict,r_flag,direction_col,index_dict,arraySize,alpha):
    ## Subset index_dict
    #['726_Complex', '726_Signaling', '726_PPI', '726_Metabolic', '726_Regulatory']
    
    subset_tuples = []
    for t in merged_llrs:
        ev_id,gs = t.split('_')
        sp = evidence_list_host[int(ev_id)]
        subset_tuples.append((gs,e_type,ev_id,sp))
    first_ev,last_ev = ev_order[0],ev_order[-1]
    tmp_index_dict = tmpDictIndexes(index_dict,subset_tuples,first_ev,last_ev)
    
    # if len(direction_col)>0: meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int', 'direction':'float', 'LLR': 'object'}
    # else: meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int', 'LLR': 'object'}
    meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int'}
    for dir_col in direction_col:
        meta_dict[dir_col] = meta_dict.get(dir_col,'int')
        if dir_col in pathogenA_df: pathogenA_df[dir_col] = pathogenA_df[dir_col].fillna(0)
    meta_dict['LLR'] = meta_dict.get('LLR','object')
    meta_dict['hasTrue'] = meta_dict.get('hasTrue','int8')
    print(pathogenA_df)
    if len(merged_llrs)>1:
        if e_type in evidences_redundancy:
            speciesA_ddf = dd.from_pandas(pathogenA_df, npartitions=288)  # Adjust the number of partitions as needed
            result_ddf = speciesA_ddf.apply(lambda row: updateMultipleLLREvidenceAndSpeciesWithRedundancy(row,r_dict,r_flag,direction_col,tmp_index_dict,first_ev,last_ev,arraySize,alpha), axis=1, meta=meta_dict)
            with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
                pathogenA_df = result_ddf.compute() 
        else:
            if e_type in evidences_orthology:
                speciesA_ddf = dd.from_pandas(pathogenA_df, npartitions=288)  # Adjust the number of partitions as needed
                result_ddf = speciesA_ddf.apply(lambda row: updateMultipleLLREvidenceAndSpeciesWithOrthology(row,direction_col,tmp_index_dict,first_ev,last_ev,arraySize), axis=1, meta=meta_dict)
                with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
                    pathogenA_df = result_ddf.compute()
            else:
                speciesA_ddf = dd.from_pandas(pathogenA_df, npartitions=288)  # Adjust the number of partitions as needed
                result_ddf = speciesA_ddf.apply(lambda row: updateMultipleLLREvidenceAndSpecies(row,direction_col,tmp_index_dict,first_ev,last_ev,llr_collected,arraySize), axis=1, meta=meta_dict)
                with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
                    pathogenA_df = result_ddf.compute() 
    else:
        col_name = merged_llrs[0]
        ev_id,g_type=col_name.split("_")
        sp_i = tmp_index_dict[g_type][ev_id]
        et_i = tmp_index_dict[g_type]['EV']
        speciesA_ddf = dd.from_pandas(pathogenA_df, npartitions=288)  # Adjust the number of partitions as needed
        result_ddf = speciesA_ddf.apply(lambda row: updateSingleLLREvidenceAndSpecies(row,sp_i,et_i,tmp_index_dict[g_type],first_ev,last_ev,direction_col,llr_collected,arraySize), axis=1, meta=meta_dict)
        with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
            pathogenA_df = result_ddf.compute() 
    return pathogenA_df

###########################
## PPV extraction
def computePPV(pathogenA_df,host_goldstandard_list,prefix_col,pathogenA,gs_order,instanceConfig,goldStandardConfig,trainingConfig):
    for host_goldstandard_t in host_goldstandard_list:
        g_type = host_goldstandard_t.type
        time_start = time.perf_counter()
        ## Compute PPV
        logPPV = LogisticPPV.objects.filter(Q(goldStandard=host_goldstandard_t))
        if 'Nan' not in logPPV[0].function:
            pathogenA_df[g_type] = FBStoPPV(pathogenA_df[g_type],logPPV[0])
        else:
            print('Error is: ' + str(logPPV[0].error))
            pathogenA_df[g_type] = 0
        ## To test PFC, comment the previous code
        # pathogenA_df[g_type] = np.round(1/(1+np.exp(-np.log(0.01)-pathogenA_df[g_type])),3)
        time_elapsed = (time.perf_counter() - time_start)
        print('### %s PPV extracted in %5.1f secs' %(g_type,time_elapsed))
        
    # Extract max_ppv
    time_start = time.perf_counter()
    pathogenA_df['max_ppv'] = pathogenA_df[gs_order].max(axis=1)    
    pathogenA_df['ppv'] = pathogenA_df[gs_order].astype(str).agg(','.join, axis=1)
    pathogenA_df = pathogenA_df.drop(gs_order,axis=1)
    time_elapsed = (time.perf_counter() - time_start)
    print('### maxPPV extracted in %5.1f secs' %time_elapsed)
    print(pathogenA_df)
    return pathogenA_df

###########################
## Gold Standard Links
def replaceEmptyValue(df,gs_order,ev_order,sp_order):
    ## Otherwise do this once at the beginning of evidence calculation
    llr = ''.join([','.join(map(str,[0] * len(ev_order))) + '|' + ','.join(map(str,[0] * len(sp_order))) + ';'] * len(gs_order))      
    df['llr_evidence'] = df['llr_evidence'].fillna(llr[:-1])
    df['max_ppv'] = df['max_ppv'].fillna(0)
    df['ppv'] = df['ppv'].fillna(','.join(map(str,[0] * len(gs_order))))
    df['fbs_goldstandard'] = df['fbs_goldstandard'].fillna(','.join(map(str,[0] * len(gs_order))))
    df['isGoldStandard'] = '0'*len(gs_order)
    df['isPPVgoldstandard'] = '0'*len(gs_order)
    return df

###########################
## Training framework
###########################
def unpackParams(instanceConfig,trainingConfig):
    ### KDE config
    min_n_sample = trainingConfig['KDE']['nSample'][0]
    
    ## Redundancy config
    evidences_redundancy = trainingConfig['Redundancy']['Evidences']
    method_redundancy = trainingConfig['Redundancy']['Method']
    size_redundancy = trainingConfig['Redundancy']['MinSize']
    alpha = trainingConfig['Redundancy']['Alpha']

    ## Orthology config
    evidences_orthology = trainingConfig['OrthologyTransfer']['Evidences']
    
    ## Network config
    visualize_PPV = trainingConfig['Network']['VisualizePPV']
    if visualize_PPV: 
        for hostA_pathogenA in instanceConfig['instance']['hostPathogen']:
            pathogenA = hostA_pathogenA[1]
            if not os.path.exists('data/tmp/PPV/%s' %pathogenA): os.makedirs('data/tmp/PPV/%s' %pathogenA)
    min_ppv = trainingConfig['Network']['MinPPV']
    gs_order = trainingConfig['Network']['GoldStandard_order']
    ev_order = trainingConfig['Network']['Evidence_order']
    sp_order = trainingConfig['Network']['Species_order']
    gs_directed = trainingConfig['Network']['GoldStandard_directed']
    ev_directed = trainingConfig['Network']['Evidence_directed']

    params = [min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed]
    
    return params

def getGoldStandardInfo(pathogenA,goldstandard_list,params,parallel):
    min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed = params

    ## GoldStandard
    ## Subset species-specific gs info
    host_goldstandard_list = [gt for gt in goldstandard_list if gt.species.tax_id==pathogenA]
    host_gs_order = [gs for gs in gs_order if gs in [gt.type for gt in host_goldstandard_list]]
    host_gs_missing = [0 if gs in host_gs_order else -1 for gs in gs_order]
    speciesA_gs_directed = [gs_directed[i] for i in [gs_order.index(gs) for gs in host_gs_order]]

    # For each gold standard
    with tqdm_joblib(tqdm(desc='Read GS', total=len(host_goldstandard_list))) as progress_bar:
        goldstandard_df_list = parallel(delayed(readGoldStandard)(host_goldstandard_t,host_gs_order,speciesA_gs_directed) for host_goldstandard_t in host_goldstandard_list)
    host_goldstandard_df = reduce(lambda df_left,df_right: pd.merge(df_left, df_right, 
                        on=['proteinA_id','proteinB_id'], how='outer'), 
                    [i for i in goldstandard_df_list if i is not None])
    speciesA_goldstandard_dict_count = dict([(g_type,host_goldstandard_df[g_type].count()) for g_type in host_goldstandard_df.columns if g_type in gs_order])
    print(host_goldstandard_df)

    ## Update goldstandard info after reading
    host_goldstandard_list = [gt for gt in host_goldstandard_list if (gt.type in host_goldstandard_df.columns and speciesA_goldstandard_dict_count[gt.type]>=min_n_sample)]
    if len(host_goldstandard_list)==0: 
        print('### No network for %s' %pathogenA)
        return None,None,None,None
    host_gs_order = [gs for gs in gs_order if gs in [gt.type for gt in host_goldstandard_list]]
    host_gs_missing = [0 if gs in host_gs_order else -1 for gs in gs_order]

    return [host_goldstandard_list,host_gs_order,host_gs_missing]

def getLLR(pathogenA,hostA,goldStandardInfo,evidence_list,params,instanceConfig,goldStandardConfig,trainingConfig,parallel):
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

    min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    host_goldstandard_list,host_gs_order,host_gs_missing = goldStandardInfo

    # evidence_list_pathogen = dict([(et.id,et.species.tax_id) for et in evidence_list if et.species.tax_id==pathogenA])
    evidence_list_host = dict([(et.id,et.species.tax_id) for et in evidence_list if et.species.tax_id==hostA])
    evidence_list_pathogen = [et for et in evidence_list if et.species.tax_id==pathogenA]

    ev_training_order = ['PIN']
    evidence_type_list = [et for et in ev_training_order if et in list({et.type for et in evidence_list_pathogen})]

    speciesA_f1 = "data/tmp/network/%s_llr.pkl" % pathogenA

    ## GoldStandard
    gs_id = {}
    arraySize = len(host_gs_order) * (len(ev_order) + len(sp_order))
    index_dict = createDictIndexes(host_gs_order,ev_order,sp_order)

    # Evidence
    if os.path.exists(speciesA_f1):
        print('### Reading pickled network for %s' %pathogenA)
        pathogenA_df = joblib.load(speciesA_f1)
        pathogenA_df['proteinA_id'] = pathogenA_df['proteinA_id'].astype(int)
        pathogenA_df['proteinB_id'] = pathogenA_df['proteinB_id'].astype(int)
        evidence_type_list = []
        print(pathogenA_df)
    else:
        pathogenA_df = pd.DataFrame()

    for e_type in evidence_type_list:
        prefix_col = ['proteinA_id','proteinB_id']
        if ev_directed[ev_order.index(e_type)]:
            prefix_col = ['proteinA_id','proteinB_id','direction']

        ## Check if orthology transfer is allowed for the evidence
        # TODO --> activate if more host-pathogen interactomes are added
        # if e_type in evidences_orthology:
        #     tmp_evidence_t_list = [et for et in evidence_list_pathogen if et.type == e_type]
        # else:
            # tmp_evidence_t_list = [et for et in evidence_list_pathogen if (et.type == e_type and et.species.tax_id==pathogenA)]
        tmp_evidence_t_list = [et for et in evidence_list_pathogen if et.type == e_type]
        
        ## If there is no evidence for pathogenA
        if len(tmp_evidence_t_list)==0: continue
        
        ## Read and merge all evidence-specific datasets
        time_start = time.perf_counter()
        if len(tmp_evidence_t_list)>1:
            with tqdm_joblib(tqdm(desc='Read %s' %e_type, total=len(tmp_evidence_t_list))) as progress_bar:
                evidence_df_list = parallel(delayed(readEvidence)(tmp_evidence_t,prefix_col,pathogenA,ev_order,ev_directed) for tmp_evidence_t in tmp_evidence_t_list)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Read %s in %5.1f secs' %(e_type,time_elapsed))
            
            ## Merge evidences with pathogenA protein_id
            time_start = time.perf_counter()
            # pathogenA_evidence_df = pd.DataFrame()
            pathogenA_evidence_df = reduce(lambda df_left,df_right: pd.merge(df_left, df_right, 
                                on=prefix_col, how='outer'), 
                            [i for i in evidence_df_list if i is not None])
            ## Free up memory
            del evidence_df_list
            pathogenA_evidence_df = pathogenA_evidence_df.reset_index(drop=True)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Merged %s in %5.1f secs' %(e_type,time_elapsed))
        else:
            pathogenA_evidence_df = readEvidence(tmp_evidence_t_list[0],prefix_col,pathogenA,ev_order,ev_directed)
            pathogenA_evidence_df = pathogenA_evidence_df.reset_index(drop=True)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Read %s in %5.1f secs' %(e_type,time_elapsed))
        print(pathogenA_evidence_df)
        # pathogenA_evidence_df[(pathogenA_evidence_df['proteinA_id']==351481) & (pathogenA_evidence_df['proteinB_id']==334772)]
        # pathogenA_evidence_df[(pathogenA_evidence_df['proteinA_id']==351481)].sort_values(771)

        # Set data type to SparseDtype to achieve memory efficiency?
        # pathogenA_evidence_df.iloc[:, len(prefix_col)::] = pathogenA_evidence_df.iloc[:, len(prefix_col)::].astype(pd.SparseDtype("float", np.nan))
        num_ev = pathogenA_evidence_df.shape[1]
        if num_ev==0: continue
        else: 
            # Reorder the columns in the DataFrame
            suffix_col = [col for col in pathogenA_evidence_df.columns if col not in prefix_col]
            pathogenA_evidence_df = pathogenA_evidence_df[prefix_col + suffix_col]

            # pathogenA_evidence_df = pathogenA_evidence_df[pathogenA_evidence_df['proteinA_id'] != pathogenA_evidence_df['proteinB_id']].reset_index(drop=True)
            pathogenA_evidence_df['proteinA_id'] = pathogenA_evidence_df['proteinA_id'].astype(int)
            pathogenA_evidence_df['proteinB_id'] = pathogenA_evidence_df['proteinB_id'].astype(int)

            if 'direction' in pathogenA_evidence_df.columns: 
                pathogenA_evidence_df['direction'] = pathogenA_evidence_df['direction'].astype(float)

            ## GoldStandard,Ev training
            ev_collected = pathogenA_evidence_df.columns[len(prefix_col):num_ev].to_list()
            llrs_converted = []

            for host_goldstandard_t in host_goldstandard_list:
                g_type=host_goldstandard_t.type
                gs_id[g_type] = gs_id.get(g_type,host_goldstandard_t.id)
                
                ### Collect all ev,gs combos
                time_start = time.perf_counter()
                llr_converted = pathogenA_evidence_df.apply(lambda column: convertToLLR_parallel(column,e_type,hostA,host_goldstandard_t,prefix_col), axis=0)
                llrs_converted.append(llr_converted)
                time_elapsed = (time.perf_counter() - time_start)
                print('### Converting raw score to LLR for %s_%s in %5.1f secs' %(g_type,e_type,time_elapsed))

            #################################################
            # Extract results and update DataFrame
            pathogenA_evidence_df = pathogenA_evidence_df.drop(ev_collected,axis=1)
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
                    precomputed_r = Redundancy.objects.filter(Q(goldStandard__in=host_goldstandard_list) & Q(evidenceA__in=ev_list) & Q(evidenceB__in=ev_list))
                    if not precomputed_r.exists():
                        if len(concatenated_llrs)>int(50e6):
                            random.seed(pathogenA)
                            compute_correlation(concatenated_llrs.sample(frac=0.25),gs_id,llr_collected,method_redundancy,size_redundancy)
                        else:
                            compute_correlation(concatenated_llrs,gs_id,llr_collected,method_redundancy,size_redundancy)
                    time_elapsed = (time.perf_counter() - time_start)
                    print('### Done in %5.1f secs' %time_elapsed)

                pathogenA_evidence_df = pd.concat([pathogenA_evidence_df,concatenated_llrs],axis=1)
                # print(pathogenA_evidence_df.dtypes)
                # If you want to keep the filtered rows in the original dataframe
                if 'direction' in prefix_col:
                    ev_direction_col = e_type + '_direction'
                    prefix_col = ['proteinA_id','proteinB_id',ev_direction_col]
                    pathogenA_evidence_df.columns.values[2] = ev_direction_col
                
                pathogenA_evidence_df = pathogenA_evidence_df[prefix_col+merged_llrs]

                ## Drop links with no ev,gs successfull combo
                pathogenA_evidence_df = pathogenA_evidence_df.loc[pathogenA_evidence_df.dropna(subset=merged_llrs, how='all').index]

                # Perform integration of LLR with weighted redundancy schema
                direction_col = [dir_col for dir_col in pathogenA_evidence_df.columns if 'direction' in dir_col]
                r_flag,r_dict = False,{}
                if e_type in evidences_redundancy:
                    # Extract a set of integers by splitting on underscores
                    ev_list = set(int(i.split('_')[0]) for i in merged_llrs)
                    r_df = pd.DataFrame.from_records(\
                        Redundancy.objects.filter(Q(goldStandard__in=host_goldstandard_list) & Q(evidenceA__in=ev_list) & Q(evidenceB__in=ev_list) & Q(correlation__gt=0))\
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
                pathogenA_evidence_df = pathogenA_evidence_df.sample(frac=1).reset_index(drop=True)
                time_start = time.perf_counter()
                pathogenA_evidence_df = updateLLRParallel(pathogenA_evidence_df,e_type,evidence_list_host,evidences_redundancy,evidences_orthology,ev_order,merged_llrs,llr_collected,r_dict,r_flag,direction_col,index_dict,arraySize,alpha)
                pathogenA_evidence_df['proteinA_id'] = pathogenA_evidence_df['proteinA_id'].astype(int)
                pathogenA_evidence_df['proteinB_id'] = pathogenA_evidence_df['proteinB_id'].astype(int)
                time_elapsed = (time.perf_counter() - time_start)
                print('### Integrated %s LLR in %5.1f secs' %(e_type,time_elapsed))
                
                if pathogenA_df.empty:
                    pathogenA_df = pathogenA_evidence_df.copy()
                    pathogenA_df['proteinA_id'] = pathogenA_df['proteinA_id'].astype(int)
                    pathogenA_df['proteinB_id'] = pathogenA_df['proteinB_id'].astype(int)
                    del pathogenA_evidence_df
                else:
                    time_start = time.perf_counter()
                    pathogenA_df = pd.merge(pathogenA_df,pathogenA_evidence_df, on=['proteinA_id','proteinB_id'],how='outer')
                    time_elapsed = (time.perf_counter() - time_start)
                    print('### Merged in %5.1f secs' %time_elapsed)
                    del pathogenA_evidence_df
                    time_start = time.perf_counter()
                    pathogenA_df['LLR'] = pathogenA_df[["LLR_x", "LLR_y"]].sum(axis=1)
                    pathogenA_df = pathogenA_df.drop(['LLR_x','LLR_y'],axis=1)
                    pathogenA_df['hasTrue'] = pathogenA_df[["hasTrue_x", "hasTrue_y"]].sum(axis=1)
                    pathogenA_df = pathogenA_df.drop(["hasTrue_x", "hasTrue_y"],axis=1)
                    time_elapsed = (time.perf_counter() - time_start)
                    print('### LLR summed up in %5.1f secs' %time_elapsed)
                
                ########################################
                ## Memory usage of pandas dataframe
                pathogenA_df.info(verbose = False, memory_usage = 'deep')
                ## Dataframe sparsity: % of values that have not been “compressed”
                # if 'LLR' in pathogenA_df.columns: print(pathogenA_df.iloc[:,len(prefix_col)+1::].sparse.density)
                # else: print(pathogenA_df.iloc[:,len(prefxix_col)::].sparse.density)
                ########################################
                print(pathogenA_df)
            else:
                del pathogenA_evidence_df 
                continue

            # Intermediate step, save to pkl. Change to save more intermediary steps.
            if len(pathogenA_df)>0: 
                print('Dumping pathogenA_df into disk after %s' %e_type)
                joblib.dump(pathogenA_df, "data/tmp/network/%s_llr.pkl" % pathogenA)
    
    return pathogenA_df

def getFBS(pathogenA_df,pathogenA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig):
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    
    min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    host_goldstandard_list,host_gs_order,host_gs_missing = goldStandardInfo
    index_dict = createDictIndexes(host_gs_order,ev_order,sp_order)

    #######################################
    #######################################
    prefix_col = ['proteinA_id','proteinB_id']
    direction_col = [dir_col for dir_col in pathogenA_df.columns if 'direction' in dir_col]
    for dir_col in direction_col:
        prefix_col.append(dir_col)
        pathogenA_df[dir_col] = pathogenA_df[dir_col].fillna(0).astype(int)

    #######################################
    ## FBS extraction
    #######################################
    time_start = time.perf_counter()
    pathogenA_df = computeFBS(pathogenA_df,index_dict,direction_col,gs_order,ev_order,sp_order)
    time_elapsed = (time.perf_counter() - time_start)
    print('### FBS done in %5.1f secs' %time_elapsed)
    print(pathogenA_df)
    # fbs_f = "data/tmp/network/%s_ONLYfbs.pkl" %pathogenA 
    fbs_f = "data/tmp/network/%s_fbs.pkl" %pathogenA 
    if not os.path.exists(fbs_f):
        joblib.dump(pathogenA_df, fbs_f)

    return pathogenA_df

def getPPV(pathogenA_df,pathogenA,hostA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig):
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    
    min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    host_goldstandard_list,host_gs_order,host_gs_missing = goldStandardInfo
    index_dict = createDictIndexes(host_gs_order,ev_order,sp_order)

    #######################################
    #######################################
    prefix_col = ['proteinA_id','proteinB_id']
    direction_col = [dir_col for dir_col in pathogenA_df.columns if 'direction' in dir_col]
    for dir_col in direction_col:
        prefix_col.append(dir_col)
        pathogenA_df[dir_col] = pathogenA_df[dir_col].fillna(0).astype(int)

    #######################################
    ## ppv extraction
    #######################################
    time_start = time.perf_counter()
    pathogenA_df = computePPV(pathogenA_df,host_goldstandard_list,prefix_col,pathogenA,gs_order,instanceConfig,goldStandardConfig,trainingConfig['Network'])
    time_elapsed = (time.perf_counter() - time_start)
    print('### PPV done in %5.1f secs' %time_elapsed)
    print(pathogenA_df)
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))  # Adjust figsize as needed
    fig.suptitle(pathogenA, fontsize=14)
    # Histogram of ppv_max
    axes[0].hist(pathogenA_df['max_ppv'], bins=100, color='#a7c957', edgecolor='black')
    axes[0].grid(True, linestyle='--', alpha=0.7)
    axes[0].set_yscale('log')
    axes[0].set_xlabel('PPV max')
    # axes[0].set_xlabel('PFC max')
    axes[0].set_ylabel('Nr. of links')
    axes[0].set_xlim(0.5, 1)
    axes[0].set_xticks(np.arange(0.5, 1.01, 0.1))

    cdf_values, bin_edges, _ = axes[1].hist(pathogenA_df['max_ppv'], bins=1000, color='#a7c957', edgecolor='black', cumulative=True, density=False, alpha=0)
    # cdf_values, bin_edges, _ = axes[1].hist(np.clip(fbs_to_pfc, 0, np.max(fbs_to_pfc)), bins=50, color='skyblue', edgecolor='black', cumulative=True, density=False, alpha=0)
    complement_cdf_values = len(pathogenA_df) - cdf_values
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
    plt.savefig('data/tmp/PPV/%s/%s_PPVmax.png' %(pathogenA,pathogenA))
    plt.close()

    return pathogenA_df

def getDatabaseOutput(pathogenA_df,pathogenA,hostA,goldStandardInfo,params):
    min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    host_goldstandard_list,host_gs_order,host_gs_missing = goldStandardInfo
    speciesA_name = Species.objects.get(Q(tax_id=pathogenA)).species_name
    # speciesA_name_compact = '%s.%s' %(speciesA_name.split()[0][0],speciesA_name.split()[1])
    speciesA_name_compact = 'SARS-CoV-2'

    #######################################
    ## Storing the network
    #######################################

    #######################################
    ## Replace NaN values
    pathogenA_df = replaceEmptyValue(pathogenA_df,gs_order,ev_order,sp_order)
    direction_col = [dir_col for dir_col in pathogenA_df.columns if 'direction' in dir_col]
    for dir_col in direction_col:
        pathogenA_df[dir_col] = pathogenA_df[dir_col].fillna(0).astype(int)

    ### TODO: GRG might support A-->B, B-->A, what to do?!
    ### Current approach: Max out duplicated links
    # speciesA_duplicated_df = pathogenA_df[pathogenA_df.duplicated(subset=['proteinA_id', 'proteinB_id'], keep=False)]
    # if len(speciesA_duplicated_df)>0:
    #     # Group by 'proteinA_id' and 'proteinB_id' and get the index of rows with maximum 'max_ppv'
    #     max_ppv_indices = speciesA_duplicated_df.groupby(['proteinA_id', 'proteinB_id'])['max_ppv'].idxmax()
    #     # Select the rows with maximum 'max_ppv' using the indices
    #     speciesA_duplicated_df = speciesA_duplicated_df.loc[max_ppv_indices]
    #     pathogenA_df = pd.concat([pathogenA_df,speciesA_duplicated_df],axis=0).reset_index(drop=True)

    ## Sort by max_ppv
    pathogenA_df = pathogenA_df.sort_values(by=['max_ppv'],ascending=False).reset_index(drop=True)
    pathogenA_df = pathogenA_df.drop_duplicates(subset=['proteinA_id', 'proteinB_id'], keep='first')

    #######################################
    ## Remove self-loops
    pathogenA_df = pathogenA_df[pathogenA_df['proteinA_id']!=pathogenA_df['proteinB_id']].reset_index(drop=True)

    ## Network model has a pre-set number of direction columns
    network_columns = [f.get_attname() for f in Network._meta.fields]
    network_direction_columns = [dir_col for dir_col in network_columns if 'direction' in dir_col]
    isG_id = pathogenA_df.columns.to_list().index('isGoldStandard')
    for dir_col in network_direction_columns:
        if dir_col in direction_col:
            pathogenA_df[dir_col] = pathogenA_df[dir_col].astype(int).astype(str) * len(gs_order)
        else: 
            pathogenA_df.insert(isG_id,dir_col, str(0) * len(gs_order))
            direction_col.append(dir_col)
            isG_id += 1

    direction_col = [dir_col for dir_col in pathogenA_df.columns if 'direction' in dir_col]
    pathogenA_df['proteinA_id'] = pathogenA_df['proteinA_id'].astype(int)
    pathogenA_df['proteinB_id'] = pathogenA_df['proteinB_id'].astype(int)
    write_network_to_database(pathogenA_df,speciesA_name_compact)

    return pathogenA_df

def getTSVOutput(pathogenA_df,pathogenA,hostA,params):
    min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed = params

    # speciesA_name = Species.objects.get(Q(tax_id=pathogenA)).species_name
    # speciesA_name_compact = '%s.%s' %(speciesA_name.split()[0][0],speciesA_name.split()[1])
    speciesA_name_compact = 'SARS-CoV-2'
    direction_col = [dir_col for dir_col in pathogenA_df.columns if 'direction' in dir_col]
    
    # Extract uniprot_id
    proteome = Proteome.objects.filter(Q(species__tax_id__in=[pathogenA,hostA]))
    uniprot_df = pd.DataFrame.from_records(Protein.objects.filter(Q(proteome__in=proteome)).values('id','uniprot_id'))
    uniprot_df['id'] = uniprot_df['id'].astype(int)

    new_order = ['proteinA_id','proteinB_id','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard']
    pathogenA_df = pathogenA_df[new_order]

    pathogenA_df.columns = ['id','proteinB_id','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard']
    pathogenA_df = pd.merge(pathogenA_df, uniprot_df, on='id', how='left')

    pathogenA_df.columns = ['FCAid','id','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard','proteinA']
    pathogenA_df = pd.merge(pathogenA_df, uniprot_df, on='id', how='left')

    ## TODO check if direction_col and 'isGoldStandard' are swapped here
    pathogenA_df.columns = ['FCAid','FCBid','max_ppv','ppv','fbs_goldstandard','llr_evidence'] + direction_col + ['isGoldStandard','isPPVgoldstandard','proteinA','proteinB']

    print(pathogenA_df)
    write_network_to_tsv(pathogenA_df,speciesA_name_compact,direction_col,gs_order,ev_order,sp_order)

    return pathogenA_df

###########################
## Main function =)
###########################
def trainHostPathogenNetwork(instanceConfig,goldStandardConfig,evidenceConfig,proteomeConfig,trainingConfig):
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
    hostA = [id[0] for id in instanceConfig['instance']['hostPathogen']]
    pathogenA = [id[1] for id in instanceConfig['instance']['hostPathogen']]
    goldstandard_list = [host_goldstandard_t for host_goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id__in=hostA)) \
                        if host_goldstandard_t.type in instanceConfig['instance']['hostPathogenGoldStandard']]

    evidence_list = [evidence_t for evidence_t in Evidence.objects.filter(Q(type__in=instanceConfig['instance']['hostPathogenEvidence']) \
                        & Q(species__tax_id__in=pathogenA+hostA)).order_by('species_id') \
                        if (evidence_t.type in evidenceConfig and \
                            (evidence_t.scoringMethod in str(evidenceConfig[evidence_t.type]['scoring_method'])))]
    
    # For each evidence type (e.g. MEX, QMS, TFB, …)
    ## Sort species by proteome size
    sorted_species = []
    for hostA_pathogenA in instanceConfig['instance']['hostPathogen']:
        proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=hostA_pathogenA[1]))[0]
        num_protein =  Protein.objects.filter(Q(proteome=proteome)).count()
        sorted_species.append((hostA_pathogenA,num_protein))
    
    # Sort the list by proteome size
    sorted_species = sorted(sorted_species, key=lambda x: x[1])
    # sorted_species = [(pathogenA,num_proteins) for pathogenA,num_proteins in sorted_species if pathogenA in ['9606']]
    for hostA_pathogenA,num_proteins in sorted_species:
        species_time_start = time.perf_counter()
        hostA = hostA_pathogenA[0]
        pathogenA = hostA_pathogenA[1]
        host_name = Species.objects.get(Q(tax_id=hostA)).species_name
        host_name_compact = '%s.%s' %(host_name.split()[0][0],host_name.split()[1])
        pathogen_name = Species.objects.get(Q(tax_id=pathogenA)).species_name
        # pathogen_name_compact = '%s.%s' %(pathogen_name.split()[0][0],pathogen_name.split()[1])
        ## TODO --> change when adding more pathogens
        pathogen_name_compact = 'SARS-CoV-2'
        print_frame('%s-%s, %i proteins' %(host_name_compact,pathogen_name_compact,num_proteins))
        
        # If network exists in TSV format, skip all computations
        network_f = 'website/static/website/networks/FunCoup6.0/FC6.0_%s_compact.gz' %pathogen_name_compact
        if os.path.exists(network_f): 
            print('### Network already exists for %s' %pathogenA)
            continue

        # If network exists in the database, then only extract TSV format
        proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=pathogenA))[0]
        proteins = Protein.objects.filter(Q(proteome=proteome))
        pathogenA_network = Network.objects.filter(Q(proteinA_id__in=proteins))
        if pathogenA_network.exists():
            col_names = tuple([col.name for col in Network._meta.get_fields() if col.name!='id'])
            pathogenA_df = pd.DataFrame.from_records(pathogenA_network.values(*col_names))
            pathogenA_df = pathogenA_df.rename(columns={'proteinA': 'proteinA_id', 'proteinB': 'proteinB_id'})

            # Output networks
            species_output_start = time.perf_counter()
            print('### STEP 4b: TSV Output')
            # goldStandardInfo[-1] = speciesA_gs_missing
            pathogenA_df = getTSVOutput(pathogenA_df,pathogenA,hostA,params)
            species_output_elapsed = (time.perf_counter() - species_output_start)
            print('### TSV Output network for %s in %5.1f secs' %(pathogenA,species_output_elapsed))
            continue

        species_GS_start = time.perf_counter()
        goldStandardInfo = getGoldStandardInfo(hostA,goldstandard_list,params,parallel)
        species_GS_elapsed = (time.perf_counter() - species_GS_start)
        print('### Extracted GS for %s in %5.1f secs' %(pathogenA,species_GS_elapsed))

        fbs_f = "data/tmp/network/%s_fbs.pkl" %pathogenA
        pathogenA_df = []
        if os.path.exists(fbs_f):
            species_FBS_start = time.perf_counter()
            print('### 2: FBS')
            pathogenA_df = joblib.load(fbs_f)
            species_FBS_elapsed = (time.perf_counter() - species_FBS_start)
            print('### Extracted FBS for %s in %5.1f secs' %(pathogenA,species_FBS_elapsed))
        else:
            ## LLR integrated
            species_LLR_start = time.perf_counter()
            print('### 1: LLR')
            pathogenA_df = getLLR(pathogenA,hostA,goldStandardInfo,evidence_list,params,instanceConfig,goldStandardConfig,trainingConfig,parallel)
            if 'hasTrue' in pathogenA_df.columns:
                pathogenA_df = pathogenA_df[pathogenA_df['hasTrue']>0].reset_index(drop=True)
                pathogenA_df = pathogenA_df.drop(['hasTrue'],axis=1)
            print(pathogenA_df)
            species_LLR_elapsed = (time.perf_counter() - species_LLR_start)
            print('### Extracted LLR for %s in %5.1f secs' %(pathogenA,species_LLR_elapsed))
            
            # Extract FBS
            species_FBS_start = time.perf_counter()
            print('### 2: FBS')
            pathogenA_df = getFBS(pathogenA_df,pathogenA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig)
            species_FBS_elapsed = (time.perf_counter() - species_FBS_start)
            print('### Extracted FBS for %s in %5.1f secs' %(pathogenA,species_FBS_elapsed))

        # Extract PPV
        species_PPV_start = time.perf_counter()
        print('### 3: PPV')
        pathogenA_df = getPPV(pathogenA_df,pathogenA,hostA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig)
        species_PPV_elapsed = (time.perf_counter() - species_PPV_start)
        print('### Extracted ppv for %s in %5.1f secs' %(pathogenA,species_PPV_elapsed))

        # pathogenA_df[(pathogenA_df['proteinA_id']==351481) & (pathogenA_df['proteinB_id']==334772)]
        # pathogenA_df[(pathogenA_df['proteinA_id']==351481)].sort_values('max_ppv')

        # Output networks
        species_output_start = time.perf_counter()
        print('### STEP 4a: Database Output')
        # goldStandardInfo[-1] = speciesA_gs_missing
        pathogenA_df = getDatabaseOutput(pathogenA_df,pathogenA,hostA,goldStandardInfo,params)
        species_output_elapsed = (time.perf_counter() - species_output_start)
        print('### Database Output network for %s in %5.1f secs' %(pathogenA,species_output_elapsed))

        species_output_start = time.perf_counter()
        print('### STEP 4b: TSV Output')
        # goldStandardInfo[-1] = speciesA_gs_missing
        pathogenA_df = getTSVOutput(pathogenA_df,pathogenA,hostA,params)
        species_output_elapsed = (time.perf_counter() - species_output_start)
        print('### TSV Output network for %s in %5.1f secs' %(pathogenA,species_output_elapsed))

        species_time_elapsed = (time.perf_counter() - species_time_start)
        print('### Extracted network for %s in %5.1f secs' %(pathogenA,species_time_elapsed))

        del pathogenA_df

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
    min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,evidences_orthology,min_ppv,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    ## Extract Gold Standard and Evidence
    hostA = [id[0] for id in instanceConfig['instance']['hostPathogen']]
    pathogenA = [id[1] for id in instanceConfig['instance']['hostPathogen']]
    goldstandard_list = [host_goldstandard_t for host_goldstandard_t in GoldStandard.objects.filter(Q(species__tax_id__in=hostA)) \
                        if host_goldstandard_t.type in instanceConfig['instance']['hostPathogenGoldStandard']]

    evidence_list = [evidence_t for evidence_t in Evidence.objects.filter(Q(type__in=instanceConfig['instance']['hostPathogenEvidence']) \
                        & Q(species__tax_id__in=[pathogenA,hostA])).order_by('species_id') \
                        if (evidence_t.type in evidenceConfig and \
                            (evidence_t.scoringMethod in str(evidenceConfig[evidence_t.type]['scoring_method'])))]

    hostA = hostA[0]
    pathogenA = pathogenA[0]
    speciesA_name = Species.objects.get(Q(tax_id=pathogenA)).species_name
    # speciesA_name_compact = '%s.%s' %(speciesA_name.split()[0][0],speciesA_name.split()[1])
    speciesA_name_compact = 'SARS-CoV-2'
    # pathogenA_df.columns = ['proteinA_id','proteinB_id'] + gs_order
    print_frame(speciesA_name_compact)

    ## Gold standard info
    species_GS_start = time.perf_counter()
    goldStandardInfo = getGoldStandardInfo(hostA,goldstandard_list,params,parallel)
    species_GS_elapsed = (time.perf_counter() - species_GS_start)
    # goldStandardInfo[0].to_csv('data/tmp/%s_goldStandard.gz' %pathogenA, sep="\t", compression='gzip', index=False)
    print('### Extracted GS for %s in %5.1f secs' %(pathogenA,species_GS_elapsed))