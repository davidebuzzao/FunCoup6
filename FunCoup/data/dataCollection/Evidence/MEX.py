import time
import os
import io
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import GEOparse
import ftplib
import requests
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from contextlib import closing
from io import StringIO
from joblib import Parallel, delayed

def getMEX_EBI(gse_name):
    ## Download RAW counts (or tpms) from EBI ftps
    if not os.path.isfile('data/tmp/'+gse_name+'.tsv'):
        url = 'http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/'+gse_name+'/'+gse_name+'-raw-counts.tsv.undecorated'
        myfile = requests.get(url)
        if myfile.status_code==200:
            open('data/tmp/'+gse_name+'.tsv', "wb").write(myfile.content)
        else:
            # TPM comes in a weird format
            # url = 'http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/'+gse_name+'/'+gse_name+'-tpms.tsv'
            # myfile = requests.get(url)
            # if myfile.status_code==200:
            #     open('data/tmp/'+gse_name+'.tsv', "wb").write(myfile.content)
            # else:
            print("Cant download file from: "+url)
            return None

    ## Read expression data into data frame
    data = pd.read_csv('data/tmp/'+gse_name+'.tsv', sep='\t', comment='#')
    return data

def compute_mapping(gse,gpl_dict):
    # Extract platform/technology
    gpl = gse.metadata['platform_id'][0]
    if gpl in gpl_dict:
        mapping = gse.gpls[gpl].table
        print("Imported: " + gpl)
        mapping = mapping.iloc[:,[0,gpl_dict[gpl]]]
        mapping = mapping.dropna()

        # Solve ambigous mapping probes, e.g. 1007_s_at --> 780 /// 100616237 into:
        # 1007_s_at --> 780
        # 1007_s_at --> 100616237
        new_mapping = []
        for i,row in mapping.iterrows():
            id = row[0]
            row = str(row[1]).split('///')
            for gene in row:
               if gpl == 'GPL6864':
                   #Oryza sativa subsp. japonica
                   gene = gene.split('|')[0]
               new_mapping.append([id,gene])

        new_mapping = pd.DataFrame(new_mapping, columns=['id','gene'])
        return new_mapping,gpl
    else:
        print("####### Could not find: " + gpl + " ###############")
        return None,gpl

def getMEX_GEO(gse_name, gpl_dict, exp='ma'):
    data = None
    mapping = None
    gpl = ''
    if exp == 'ma':
        try:
            if os.path.isfile("data/tmp/" + gse_name + "_family.soft.gz"):
                # Load file
                gse = GEOparse.get_GEO(filepath="data/tmp/" + gse_name + "_family.soft.gz",silent=True)
            else:
                # Download
                gse = GEOparse.get_GEO(gse_name,destdir='data/tmp/',silent=True)
            data = gse.pivot_samples('VALUE')
            mapping,gpl = compute_mapping(gse,gpl_dict)
        except:
            return None,None,None
    elif exp == 'rseq':
        if any(gse_name in f for f in os.listdir("data/tmp/")):
            file_name = [f for f in os.listdir("data/tmp/") if gse_name in f]
            if len(file_name)>1:
                print('Ambiguities for %s' %gse_name)
                return None,None,None
        else:
            try:
                annot_number=gse_name.split("SE")[1][0:-3]
                ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
                ftp.login()
                ftp.cwd("/geo/series/gse"+annot_number+"nnn/"+gse_name+"/suppl/")

                directory_listing = io.StringIO()
                ftp.retrlines("LIST", directory_listing.write)

                files = [line.split()[-1] for line in directory_listing.getvalue().splitlines()]
                ftp.quit()
            except:
                print('FTP error for %s' %gse_name)
                return None,None,None
            else:
                if len(files)>1:
                    print('Possible ambiguities for %s (%s)' %(gse_name,','.join(files)))
                for file_name in files:
                    print(file_name)
                    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/gse"+ annot_number+"nnn/"+gse_name+"/suppl/"+file_name
                    req = requests.get(url)
                    with open('data/tmp/'+file_name, "wb") as filout:
                        filout.write(req.content)

        if '.tar' in file_name: return None,None,None
        if '.xlsx' in file_name:
            data = pd.read_excel('data/tmp/'+file_name)
            if data.index.dtype == 'int' or data.index.dtype == 'int64':
                data.index = data.iloc[:,0]
        elif '.csv' in file_name:
            data = pd.read_csv('data/tmp/'+file_name, index_col=0)
            data.index = data.index.str.strip("'")
        else:
            for F in ['.txt','.tsv','.matrix','.count','.fpkm','.tpm']:
                if F in file_name:
                    data = pd.read_csv('data/tmp/'+file_name, sep="\t", index_col=0)
                    data.index = data.index.str.strip("'")
                    break
                else: return None,None,None

        gpl = ''
        mapping = pd.DataFrame(list(data.index))
        mapping.columns = ['gene']

    return data,mapping,gpl

def quality_check(data):
    ## Select only rows with more than 80% numerical values (or exp>0)
    data_zeros_nan = data.isin([0, np.nan]).sum(axis=1)/data.shape[1]
    data = data[data_zeros_nan < 0.2]
    ## Remove duplicated columns
    data = data.loc[~data.duplicated(),:]
    return data

def excludeHomologLinks(data,tax_id,eval_cutoff):
    formatted_eval_cutoff = '{:.0e}'.format(eval_cutoff)
    diamond_out = 'data/tmp/Diamond_eval%s/diamond%s-%s' %(formatted_eval_cutoff,tax_id,tax_id)
    if os.path.exists(diamond_out): 
        diamond_df = pd.read_csv(diamond_out, sep='\t', header=None)
        diamond_df.columns = ['proteinA','proteinB']
        data = pd.merge(data,diamond_df, indicator=True, how='outer', on=['proteinA','proteinB'])\
                .query('_merge=="left_only"')\
                .drop('_merge', axis=1)\
                .reset_index(drop=True)
    return data

def compute_correlation(data):
    ## Transpose data
    data = data.T

    ## Compute correlations for pairs expressed in min 3 samples,
    ## then extract upper matrix only and exclude diagonal.
    ## As the indices of data matrix are sorted in decreasing order, the stack()
    ## function will return a data frame with sorted proteins (i.e. A>B).
    r_df = data.corr(method="spearman",min_periods=3,numeric_only=True)\
        .where(np.triu(np.ones((data.shape[1],data.shape[1])), k=1).astype(bool))\
            .stack(dropna=True)\
                .reset_index()
    r_df.columns = ['proteinA','proteinB','score']
    return r_df

def get_MEX(gse_name,tax_id,evidenceConfig,proteomeConfig,gpl_list,database,exp,eval_cutoff=0.001):
    ## Extract species name
    species_name = Species.objects.filter(Q(tax_id=tax_id))[0].species_name
    ## Extract known proteins for species under study
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=tax_id))[0]
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    available_protein_df.columns = ['gene','protein_id']

    print('Dataset: %s, %s, %s' %(gse_name,species_name,tax_id))
    evidence_in_db = Evidence.objects.filter(Q(type='MEX') & Q(species__tax_id=tax_id) & Q(version=gse_name) & Q(scoringMethod=evidenceConfig['scoring_method']))
    if evidence_in_db.exists():
        species_in_db = Species.objects.get(tax_id = tax_id)
        evidence_in_db = Evidence(version=gse_name, type="MEX", species=species_in_db, scoringMethod=evidenceConfig['scoring_method'])
        return
    if database == 'EBI':
        data = getMEX_EBI(gse_name)
        if data is None:
            return
        ## Rename data frame columns
        if 'Gene Name' in data.columns:
            data.drop(columns='Gene Name', inplace=True)
        l = list(data.columns)
        l[0] = 'gene'
        data.columns = l

        ## Filter out protein names that are not in the database
        mapping = pd.DataFrame(list(data['gene']))
        mapping.columns = ['gene']
        pre_map = set(mapping['gene'])
        mapping = pd.merge(mapping,available_protein_df, on='gene', how='inner')
        post_map = set(mapping['gene'])
        excluded_due_to_mapping = len(pre_map - post_map)

        ## Remove expression data for non-mapping proteins
        data = pd.merge(data, mapping, on='gene')
        ## Rename the index of expression data frame with protein ids
        data.index = data['gene']
        data.drop(columns='gene', inplace=True)

    elif database == 'GEO':
        gpl_dict = dict([gpl.split()[0],int(gpl.split()[1])] for gpl in gpl_list)
        data,mapping,gpl = getMEX_GEO(gse_name,gpl_dict,exp)
        if data is None or mapping is None: return

        pre_map = set(mapping['gene'])
        mapping = pd.merge(mapping,available_protein_df, on='gene', how='inner')
        post_map = set(mapping['gene'])
        excluded_due_to_mapping = len(pre_map - post_map)

        if exp == 'ma':
            mapping = mapping.drop_duplicates().reset_index(drop=True)
            data['id'] = data.index
            data = pd.merge(data, mapping, on='id')
            data.drop(columns='id',inplace=True)
            data.drop(columns='gene',inplace=True)
        elif exp == 'rseq':
            data['gene'] = data.index
            data = pd.merge(data, mapping, on='gene')
            data.drop(columns='gene',inplace=True)

    ## Average redundant probes/transcripts
    data = data.groupby('protein_id').mean(numeric_only=True)
    if len(data) > 0:
        ## Remove poor quality (or uninformative) genes/samples
        tot_nr_genes_1 = len(data)
        data = quality_check(data)
        tot_nr_genes_2 = len(data)
        print("Extracting MEX for %i genes.\nExcluded:\nMapping: %i\nQuality: %i" \
                %(tot_nr_genes_2,excluded_due_to_mapping,tot_nr_genes_1-tot_nr_genes_2))
    else:
        print("Data %s did not pass the quality test!" %gse_name)

    ## Remove poor quality (or uninformative) genes/samples
    if len(data) > 0:
        data = data.sort_index(ascending=False).rename_axis(None).rename_axis(None, axis=0)
        evidenceLink = compute_correlation(data)
        evidenceLink = excludeHomologLinks(evidenceLink,tax_id,eval_cutoff)
        
        ## Subset for abs(correlation)>0.1
        evidenceLink = evidenceLink[np.abs(evidenceLink['score'])>0.1]
        
        evidenceLink['sign'] = np.where(evidenceLink['score'].values < 0, '-', '')
        abs_scores = np.abs(evidenceLink['score']).values.reshape(-1, 1)
        min_score, max_score = str(np.round(np.min(abs_scores), decimals=4)), str(np.round(np.max(abs_scores), decimals=4))
        scaler = MinMaxScaler(feature_range=(0, 1))
        scaled_abs_scores = scaler.fit_transform(abs_scores)
        evidenceLink['score'] = np.round(scaled_abs_scores, decimals=4)

        # Preparing the dataframe to be added to the database
        species_in_db = Species.objects.get(tax_id = tax_id)
        evidence_in_db = Evidence(type="MEX", species=species_in_db, scoreRange='%s|%s' %(min_score,max_score), version=gse_name, scoringMethod=evidenceConfig['scoring_method'])
        evidence_in_db.save()

        ## We shuffle the data, so that we can sample some of them easily taking first 1M later
        evidenceLink = evidenceLink.sample(frac=1).reset_index(drop=True)
        evidenceLink.insert(0,'evidence',evidence_in_db.id)
        evidenceLink['metadata'] = gse_name
        
        startTime = time.time()
        mem_csv = StringIO()
        evidenceLink.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        # Writing the csv to the evidenceLink table
        print("Writing "+str(len(evidenceLink.index))+" MEX links.")
        doWrite=True
        numberOfAttempts=0
        while doWrite==True and numberOfAttempts<=10:
            numberOfAttempts+=1
            try:
                with closing(mem_csv) as csv_io:
                    MEX.objects.from_csv(csv_io)
                doWrite=False
            except:
                print("Unexpected error when attempting to write MEX to DB. Trying again, attempts: "+str(numberOfAttempts))
                doWrite=True
        #with closing(mem_csv) as csv_io:
        #    ProteinLinkScore.objects.from_csv(csv_io)
        endTime = time.time()
        print('Time: %.2f' %(endTime-startTime))
    else:
        print("No R from MEX_" + gse_name)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='''
    blablabla
    ''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-g", "--gse",
                        help="Name of config file to use",
                        type=str)
    parser.add_argument("-t", "--taxid",
                        help="Name of config file to use",
                        type=str)
    args = parser.parse_args()

    get_MEX([args.gse],args.taxid,database='EBI')
