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
## FBS extraction
def FBStoPPVgoldStandard(tmp_speciesA_evidence_df,g_type,goldstandard_t,instanceConfig,goldStandardConfig,trainingConfig):
    # Negative gold standard
    # tmp_speciesA_evidence_df=speciesA_df[prefix_col+[g_type]].dropna(how='any',axis=0)
    # trainingConfig = trainingConfig['Network']
    speciesA = goldstandard_t.species.tax_id
    ## Only positive values
    tmp_speciesA_evidence_df = tmp_speciesA_evidence_df[tmp_speciesA_evidence_df[g_type]>0]

    positive_goldStandardLink,negative_goldStandardLink = separateNegativePositiveGoldStandard(tmp_speciesA_evidence_df,goldstandard_t,g_type,instanceConfig)
    if positive_goldStandardLink is None or negative_goldStandardLink is None: return ['Nan']
    lpos,lneg = len(positive_goldStandardLink),len(negative_goldStandardLink)

    print('Nr. of PGS:\t%i\nNr. of NGS:\t%i' %(lpos,lneg))
    coefficients = logisticPPV(positive_goldStandardLink,negative_goldStandardLink,goldstandard_t,speciesA,g_type,trainingConfig,saveToDatabase=False)
    return coefficients

def FBStoPPV(FBS,coefficients):
    if 'Nan' not in coefficients:
        coefficients = [float(i) for i in coefficients.split(',')]
        print(coefficients)
        a, b, c = coefficients
        return np.round(np.clip(logistic(np.clip(FBS, 0, np.max(FBS)),a, b, c),0.5,1),3)
    return 0

def computePPV(speciesA_df,speciesA_goldstandard_list,prefix_col,speciesA,gs_order,instanceConfig,goldStandardConfig,trainingConfig):
    for goldstandard_t in speciesA_goldstandard_list:
        g_type = goldstandard_t.type
        ## Compute PPV
        time_start = time.perf_counter()
        coefficients = FBStoPPVgoldStandard(speciesA_df[prefix_col+[g_type]].dropna(how='any',axis=0),g_type,goldstandard_t,instanceConfig,goldStandardConfig,trainingConfig)
        
        if 'Nan' not in coefficients:
            speciesA_df[g_type] = FBStoPPV(speciesA_df[g_type],coefficients)
        else:
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
## tsv output
def write_network_to_tsv(speciesA_df,source_speciesA,target_speciesA,direction_col,gs_order,ev_order,sp_order):
    sp_name_order = []
    for sp in sp_order:
        sp_name = Species.objects.get(Q(tax_id=sp)).species_name
        sp_name_compact = '%s%s' %(sp_name.split()[0][0],sp_name.split()[1][0:2].upper())
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
    f = 'data/tmp/orthologNetwork/FC6.0_s%s_t%s_full.gz' %(source_speciesA,target_speciesA)
    speciesA_df.to_csv(f, sep="\t", compression='gzip', index=False)
    ## Export compact network
    f = 'data/tmp/orthologNetwork/FC6.0_s%s_t%s_compact.gz' %(source_speciesA,target_speciesA)
    speciesA_df.iloc[:,0:(9+len(direction_col))].to_csv(f, sep="\t", compression='gzip', index=False)

###########################
## Training framework
###########################
def getLLR(source_speciesA,target_speciesA,goldStandardInfo,evidence_list,params,instanceConfig,goldStandardConfig,trainingConfig):
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

    speciesA_f1 = "data/tmp/orthologNetwork/s%s_t%s_llr.pkl" %(source_speciesA,target_speciesA)

    ## GoldStandard
    gs_id = {}
    arraySize = len(speciesA_gs_order) * (len(ev_order) + len(sp_order))
    index_dict = createDictIndexes(speciesA_gs_order,ev_order,sp_order)

    # Evidence
    if os.path.exists(speciesA_f1):
        print('### Reading pickled network for s:%s-t:%s' %(source_speciesA,target_speciesA))
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

        ####################################################################
        #### HERE MAKE SURE THAT NO EVIDENCE FROM target_speciesA IS USED
        ####################################################################
        ## Check if orthology transfer is allowed for the evidence
        if e_type in evidences_orthology:
            tmp_evidence_t_list = [et for et in evidence_list if (et.type == e_type and et.species.tax_id!=target_speciesA)]
        else:
            tmp_evidence_t_list = [et for et in evidence_list if (et.type == e_type and et.species.tax_id==source_speciesA)]
        
        if any([et.species.tax_id==target_speciesA for et in tmp_evidence_t_list]):
            print('###!!!! ERROR !!!!###')
            raise SystemExit

        ## If there is no evidence for source_speciesA
        if len(tmp_evidence_t_list)==0: continue
        
        ## Read and merge all evidence-specific datasets
        time_start = time.perf_counter()
        if len(tmp_evidence_t_list)>1:
            with tqdm_joblib(tqdm(desc='Read %s' %e_type, total=len(tmp_evidence_t_list))) as progress_bar:
                evidence_df_list = parallel(delayed(readEvidence)(tmp_evidence_t,prefix_col,source_speciesA,ev_order,ev_directed) for tmp_evidence_t in tmp_evidence_t_list)
            time_elapsed = (time.perf_counter() - time_start)
            print('### Read %s in %5.1f secs' %(e_type,time_elapsed))
            
            ## Merge evidences with source_speciesA protein_id
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
            speciesA_evidence_df = readEvidence(tmp_evidence_t_list[0],prefix_col,source_speciesA,ev_order,ev_directed)
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
                        _ = parallel(delayed(trainGoldStandard)(positive_goldStandardLink,speciesA_evidence_df[prefix_col+[e_id]].dropna(how='any',axis=0),e_id,goldstandard_t,ev_directed,gs_directed,ev_order,gs_order,e_type,g_type,instanceConfig,goldStandardConfig,trainingConfig['KDE']) \
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
                            random.seed(source_speciesA)
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
        
            if len(speciesA_df)>0 and e_type in ['MEX']: 
                print('Dumping speciesA_df into disk after %s' %e_type)
                speciesA_df = speciesA_df[speciesA_df['hasTrue']>0].reset_index(drop=True)
                speciesA_df = speciesA_df.drop(['hasTrue'],axis=1)
                joblib.dump(speciesA_df, speciesA_f1)
    
    return speciesA_df

def getFBS(speciesA_df,source_speciesA,target_speciesA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig):
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
    fbs_f = "data/tmp/orthologNetwork/s%s_t%s_fbs.pkl" %(source_speciesA,target_speciesA)
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

    return speciesA_df

def getOutput(speciesA_df,source_speciesA,target_speciesA,goldStandardInfo,params):
    visualize_LLR,min_n_sample,evidences_redundancy,method_redundancy,size_redundancy,alpha,\
        evidences_orthology,min_ppv,visualize_PPV,gs_order,ev_order,sp_order,gs_directed,ev_directed = params
    
    speciesA_goldstandard_df,speciesA_goldstandard_list,speciesA_gs_order,speciesA_gs_missing = goldStandardInfo
    speciesA_name = Species.objects.get(Q(tax_id=source_speciesA)).species_name
    speciesA_name_compact = '%s.%s' %(speciesA_name.split()[0][0],speciesA_name.split()[1])
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

    # Extract uniprot_id
    proteome = Proteome.objects.get(Q(species__tax_id=source_speciesA))
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
    write_network_to_tsv(speciesA_df,source_speciesA,target_speciesA,direction_col,gs_order,ev_order,sp_order)

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
        num_protein =  Protein.objects.filter(Q(proteome=proteome)).count()
        sorted_species.append((s,num_protein))
    
    # Sort the list by proteome size
    sorted_species = sorted(sorted_species, key=lambda x: x[1])
    # sorted_species = [(speciesA,num_proteins) for speciesA,num_proteins in sorted_species if speciesA in ['9606']]

    #######################################################
    ## The file instanceConfig['instance']['species_relations'] can be generated with networkAnalysis/getTransferredNetworks.py
    ## Transferred networks target<--source species
    orthologNetwork_df = pd.read_csv(instanceConfig['instance']['species_relations'], sep="\t", dtype=str)
    orthologNetwork_df = orthologNetwork_df[(orthologNetwork_df['FC_species']=='1') & (orthologNetwork_df['spB'].isin([tup[0] for tup in sorted_species[0:8]]))].reset_index(drop=True)
    #######################################################

    ## For testing, adding spA,spA to train spA without spA --> #TEST
    orthologNetwork_df['spB'] = orthologNetwork_df['spA']
    orthologNetwork_df['spA_name'] = orthologNetwork_df['spB_name']
    for i,row in orthologNetwork_df.iterrows():
        target_speciesA,source_speciesA = row['spA'],row['spB']

        species_time_start = time.perf_counter()
        target_speciesA_name = Species.objects.get(Q(tax_id=target_speciesA)).species_name
        target_speciesA_name_compact = '%s.%s' %(target_speciesA_name.split()[0][0],target_speciesA_name.split()[1])
        source_speciesA_name = Species.objects.get(Q(tax_id=source_speciesA)).species_name
        source_speciesA_name_compact = '%s.%s' %(source_speciesA_name.split()[0][0],source_speciesA_name.split()[1])
        print_frame('%s w/ no %s' %(source_speciesA_name_compact,target_speciesA_name_compact))

        # tmp_goldstandard_list = [gt for gt in goldstandard_list if gt.species.tax_id!=target_speciesA]
        ## To assess the quality of networks --> #TEST
        tmp_goldstandard_list = [gt for gt in goldstandard_list]
        tmp_evidence_list = [et for et in evidence_list if et.species.tax_id!=target_speciesA]

        ## Gold standard info
        species_GS_start = time.perf_counter()
        goldStandardInfo = getGoldStandardInfo(source_speciesA,tmp_goldstandard_list,params,parallel)
        species_GS_elapsed = (time.perf_counter() - species_GS_start)
        print('### Extracted GS for s:%s-t:%s in %5.1f secs' %(source_speciesA,target_speciesA,species_GS_elapsed))

        ## If network exists in TSV format, skip all computations
        network_f = 'website/static/website/networks/FunCoup6.0/FC6.0_s%s_t%s_compact.gz' %(source_speciesA,target_speciesA)
        if os.path.exists(network_f): 
            print('### Network already exists for %s' %target_speciesA)
            continue

        fbs_f = "data/tmp/orthologNetwork/s%s_t%s_fbs.pkl" %(source_speciesA,target_speciesA)
        speciesA_df = []
        if os.path.exists(fbs_f):
            species_FBS_start = time.perf_counter()
            print('### 2: FBS')
            speciesA_df = joblib.load(fbs_f)
            species_FBS_elapsed = (time.perf_counter() - species_FBS_start)
            print('### Extracted FBS for s:%s-t:%s in %5.1f secs' %(source_speciesA,target_speciesA,species_FBS_elapsed))
        else:
            ## LLR integrated
            species_LLR_start = time.perf_counter()
            print('### 1: LLR')
            speciesA_df = getLLR(source_speciesA,target_speciesA,goldStandardInfo,tmp_evidence_list,params,instanceConfig,goldStandardConfig,trainingConfig)
            species_LLR_elapsed = (time.perf_counter() - species_LLR_start)
            print('### Extracted LLR for s:%s-t:%s in %5.1f secs' %(source_speciesA,target_speciesA,species_LLR_elapsed))
            
            # Extract FBS
            species_FBS_start = time.perf_counter()
            print('### 2: FBS')
            speciesA_df = getFBS(speciesA_df,source_speciesA,target_speciesA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig)
            species_FBS_elapsed = (time.perf_counter() - species_FBS_start)
            print('### Extracted FBS for s:%s-t:%s in %5.1f secs' %(source_speciesA,target_speciesA,species_FBS_elapsed))

        # Extract PPV
        species_PPV_start = time.perf_counter()
        print('### 3: PPV')
        speciesA_df = getPPV(speciesA_df,source_speciesA,goldStandardInfo,params,instanceConfig,goldStandardConfig,trainingConfig)
        species_PPV_elapsed = (time.perf_counter() - species_PPV_start)
        print('### Extracted ppv for s:%s-t:%s in %5.1f secs' %(source_speciesA,target_speciesA,species_PPV_elapsed))

        # Output networks
        species_output_start = time.perf_counter()
        print('### 4: Output')
        # goldStandardInfo[-1] = speciesA_gs_missing
        speciesA_df = getOutput(speciesA_df,source_speciesA,target_speciesA,goldStandardInfo,params)
        species_output_elapsed = (time.perf_counter() - species_output_start)
        print('### Output network for s:%s-t:%s in %5.1f secs' %(source_speciesA,target_speciesA,species_output_elapsed))

        species_time_elapsed = (time.perf_counter() - species_time_start)
        print('### Extracted network for s:%s-t:%s in %5.1f secs' %(source_speciesA,target_speciesA,species_time_elapsed))

        del speciesA_df

    # Reset the warning filter to its default behavior (optional)
    warnings.resetwarnings()

if __name__ == '__main__':
    if os.path.exists('configFiles/exampleGenerateOrthologNetworks'):
        with open('configFiles/exampleGenerateOrthologNetworks/instanceConfig.yml') as f:
            instanceConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateOrthologNetworks/goldStandardConfig.yml') as f:
            goldStandardConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateOrthologNetworks/proteomeConfig.yml') as f:
            proteomeConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateOrthologNetworks/evidenceConfig.yml') as f:
            evidenceConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateOrthologNetworks/trainingConfig.yml') as f:
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

    speciesA = '44689'

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
