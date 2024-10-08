import re
import requests
import os 
import pybedtools
import matplotlib.pyplot as plt
import yaml
from yaml.loader import SafeLoader
from sklearn.preprocessing import MinMaxScaler,PowerTransformer
from scipy.stats import combine_pvalues
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from tqdm import tqdm
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO

from auxiliary_functions import *

'''
1. ENCODE API --> extract list of BED narrowPeak files
2. Genome release 38 from gencodegenes
3. Process each file, extract: TF,target,Pvalue,AccessionFile
4. Max out/average pvalue for redundant links. Average? Use Fisher's method.

BED narrowPeak file:
This format is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED6+4 format.
* chrom - Name of the chromosome (or contig, scaffold, etc.).
* chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
* chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
* name - Name given to a region (preferably unique). Use "." if no name is assigned.
* score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were "'0"' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
* strand - +/- to denote strand or orientation (whenever applicable). Use "." if no orientation is assigned.
* signalValue - Measurement of overall (usually, average) enrichment for the region.
* pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
* qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
* peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
'''

def plot_scoredistribution(score,file_accession):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5,5))  # Adjust figsize as needed
    axes.hist(score, bins=20, facecolor='#003049', label='Raw Scores', alpha=1, density=True)
    
    # Add title and axis labels
    axes.set_title(file_accession)
    axes.set_xlabel('Raw score')
    axes.set_ylabel('Density')

    # Adjust spacing between subplots
    plt.tight_layout()
    # Display the plot
    plt.savefig('data/tmp/encode/plots/%s' %file_accession)
    plt.close()

## STEP 1
def get_ENCODE(genome_version,species_name,encode_pipeline,fileout):
    # Define the ENCODE REST API base URL
    encode_api_base_url = "https://www.encodeproject.org/search/"
    if not os.path.exists(fileout):
        # Force return from the server in JSON format
        headers = {'accept': 'application/json'}

        # Define the ENCODE search URL
        # https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=TF+ChIP-seq&assembly=GRCh38&assay_title=total+RNA-seq&target.investigated_as=transcription+factor&files.file_type=bed+narrowPeak&limit=all&frame=object

        query_params = {
            'type': 'Experiment',
            'status': 'released',
            'assay_title': 'TF ChIP-seq',  # You can use a list for multiple values
            'assembly': genome_version,
            'target.investigated_as': 'transcription factor',
            'files.file_type': 'bed narrowPeak',
            'limit': 'all',
            'frame': 'object',
        }

        # Make a GET request to the ENCODE REST API with the query parameters
        response = requests.get(encode_api_base_url, params=query_params, headers=headers)
        # Extract the JSON response as a Python dictionary
        chipseq_data = response.json()
        print('Number of ChIP-Seq experiments: %s' %chipseq_data['total'])
        chips_seq_list = []
        for i in range(len(chipseq_data['@graph'])):
            if species_name in chipseq_data['@graph'][i]['target']:
                tf = chipseq_data['@graph'][i]['target'].split('/targets/')[1].split('-%s' %species_name)[0]
            else: continue
            experiment_accession = chipseq_data['@graph'][i]['accession']
            experiment_type = chipseq_data['@graph'][i]['assay_title']
            chips_seq_list.append([experiment_accession, tf, experiment_type])
        
        chipseq_df = pd.DataFrame(chips_seq_list)
        chipseq_df.columns = ['experiment_accession','tf','experiment_type']

        # Define the query parameters to filter the files
        # Why IDR?
        #https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html
        query_params = {
            'type': 'File',
            'assay_title': 'TF ChIP-seq',
            'file_format': 'bed',
            'output_type': 'IDR thresholded peaks',
            'assembly': genome_version,
            'limit': 'all',
            'frame': 'object',
        }

        # Check the value of lab_pipeline
        if encode_pipeline:
            query_params['lab.title'] = 'ENCODE Processing Pipeline'

        # Make a GET request to the ENCODE REST API with the query parameters
        response = requests.get(encode_api_base_url, params=query_params, headers=headers)
        # Extract the JSON response as a Python dictionary
        bed_data = response.json()
        print('Number of BED narrowPeak files: %s' %bed_data['total'])
        # Extract experiment_accession, target, assay_title and specific formatted file
        bed_data_list = []
        for i in range(len(bed_data['@graph'])):
            file_accession = bed_data['@graph'][i]['accession']
            experiment_accession = bed_data['@graph'][i]['dataset'].split('/')[-2]
            output_type = bed_data['@graph'][i]['output_type']
            try:
                preferred_default = bed_data['@graph'][i]['preferred_default']
            except: 
                preferred_default = False
            lab_pipeline = bed_data['@graph'][i]['lab'].split('/')[-2]
            bed_data_list.append([file_accession,experiment_accession,lab_pipeline,output_type,preferred_default])

        bed_data_df = pd.DataFrame(bed_data_list)
        bed_data_df.columns = ['file_accession','experiment_accession','lab_pipeline','output_type','preferred_default']
        
        encode_df = pd.merge(chipseq_df, bed_data_df, on='experiment_accession', how='inner')
        encode_df.to_csv(fileout,sep='\t', index=False)
    else:
        encode_df = pd.read_csv(fileout, sep='\t')
    return encode_df 
    
def get_annotation(annotation_url):
    # downloadAndUnzipFile(evidenceConfig['url'])
    downloadFile(annotation_url,'data/tmp/encode/'+annotation_url.split('/')[-1])

def extract_gene_info(row,mode='gene_name'):
    # Define a function to extract gene_id or gene_name
    if mode=='gene_id':
        match = re.search(r'gene_id "(.*?)"', row)
    elif mode=='gene_name':
        match = re.search(r'gene_name "(.*?)"', row)
        
    if match:
        return match.group(1) or match.group(2)
    return None

def extract_annotation(available_protein_df,annotation_url,annotation_file,species_chromosomes):
    if not os.path.exists(annotation_file):
        get_annotation(annotation_url)

        # Parse the GTF file and extract gene annotations in BED format
        annotation_df = pd.read_csv('data/tmp/encode/'+annotation_url.split('/')[-1], usecols=[0,2,3,4,6,8], header=None, comment='#', sep='\t')
        annotation_df.columns = ['chrom', 'type', 'start', 'end', 'strand','name']
        annotation_df = annotation_df[(annotation_df['type']=='gene') & annotation_df['name'].str.contains("protein_coding")].reset_index(drop=True)
        annotation_df.drop('type',axis=1,inplace=True)
        
        ## gene_id
        annotation_df['gene_id'] = annotation_df['name'].apply(lambda row: extract_gene_info(row,'gene_id')).str.lower()
        ## gene_name
        annotation_df['gene_name'] = annotation_df['name'].apply(lambda row: extract_gene_info(row,'gene_name')).str.lower()

        annotation_df = annotation_df.drop_duplicates() # Drop duplicates
        annotation_df.drop('name',axis=1,inplace=True)
        
        ## To achieve largest coverage
        gene_id_l = len(annotation_df[annotation_df['gene_id'].isin(available_protein_df['mapped_id'])].reset_index(drop=True))
        gene_name_l = len(annotation_df[annotation_df['gene_name'].isin(available_protein_df['mapped_id'])].reset_index(drop=True))
        if gene_id_l>gene_name_l:
            annotation_df = annotation_df[annotation_df['gene_id'].isin(available_protein_df['mapped_id'])].reset_index(drop=True)
            annotation_df.drop('gene_name',axis=1,inplace=True)
        else:
            annotation_df = annotation_df[annotation_df['gene_name'].isin(available_protein_df['mapped_id'])].reset_index(drop=True)
            annotation_df.drop('gene_id',axis=1,inplace=True)
        annotation_df.columns = ['chrom', 'start', 'end', 'strand', 'name']
        # set(annotation_df['chrom'].to_list())
        annotation_df = annotation_df[annotation_df['chrom'].isin(species_chromosomes)]

        # Check if 'chrom' column doesn't contain 'chr'
        mask = ~annotation_df['chrom'].str.contains('chr')
        # Add 'chr' to the 'chrom' column where the condition is True
        annotation_df.loc[mask, 'chrom'] = 'chr' + annotation_df.loc[mask, 'chrom']

        annotation_df.to_csv(annotation_file, sep='\t', header=False, index=False)

## STEP 2
def convert_logPvalue(logPvalue):
    return 10**(-logPvalue)

def annotate_peaks(encode_row,annotation_file):
    file_accession = encode_row['file_accession']#.iloc[0]
    if not os.path.exists('data/tmp/encode/%s.bed.gz' %file_accession): 
        print('%s #FILE' %file_accession)
        return None
    try:
        # Attempt to create BedTool objects
        narrowPeak_file = pybedtools.BedTool('data/tmp/encode/%s.bed.gz' %file_accession)
        gene_annotations_file = pybedtools.BedTool(annotation_file)
        # Perform the intersection
        narrowPeak_with_gene_names = narrowPeak_file.intersect(gene_annotations_file, wa=True, wb=True)
    except Exception as e:
        # Handle exceptions raised during intersection
        # print(f"Error during intersection: {str(e)}")
        print('%s #BEDtools' %file_accession)
        return None
    else:
        entry = list(narrowPeak_with_gene_names)[0]
        entry_col = str(entry).strip().split('\t')
        # Check if the number of column is larger than 15
        if len(entry_col)!=15:
            print('%s #COL' %file_accession)
            return None
        
        gene_names = []
        for entry in narrowPeak_with_gene_names:
            # entry_col = str(entry).strip().split('\t')
            # print(entry)
            # idr = entry.fields[4]
            gene_name = entry.fields[14]  # Assuming the 'name' field is at index 11
            enrichment_score = float(entry.fields[6])
            # pvalue = float(entry.fields[8])
            gene_names.append([gene_name, enrichment_score])

        if len(gene_names)<2: return None
        # Create a DataFrame from the extracted BED data
        links_df = pd.DataFrame(gene_names, columns=['proteinB', 'score'])
        # print('%s: %f' %(file_accession,np.min(links_df['score'])))
        # Check if the DataFrame is empty or contains only NaN values
        # if links_df.empty or (links_df['score'] == -1).all():
        #     return None  # Return None for empty DataFrames or those with NaN values

        tf = encode_row['tf']#.iloc[0]
        links_df['proteinA'] = tf
        links_df['proteinA'] = links_df['proteinA'].str.lower()
        links_df['proteinB'] = links_df['proteinB'].str.lower()

        links_df['metadata'] = file_accession
        links_df = links_df[['proteinA', 'proteinB', 'score', 'metadata']]
        
        # Apply the function to the 'score' column
        # links_df['score'] = links_df['score'].apply(convert_logPvalue)

        scores = links_df['score'].values.reshape(-1, 1)
        power_transformer = PowerTransformer(method='box-cox')
        power_scores = power_transformer.fit_transform(scores)

        # Min-Max scaling
        minmax_scaler = MinMaxScaler(feature_range=(0, 1))
        minmax_scores = minmax_scaler.fit_transform(power_scores)
        # plot_scoredistribution(minmax_scores,'%s.png' %tax_id)
        links_df['score'] = np.round(minmax_scores,decimals=4)
        # links_df = links_df[(links_df['pvalue']>0) & (convert_logPvalue(links_df['pvalue'])<0.05)]
        # print('%s: %i' %(file_accession,links_df.shape[0]))
        return links_df

def annotate_peaks_parallel(encode_chunk,annotation_file):
    # encode_chunk = encode_chunk.apply(lambda row: annotate_peaks(row,annotation_file), axis=1)
    # Create an empty list to store non-empty DataFrames
    result_dfs = []
    for _, row in encode_chunk.iterrows():
        df = annotate_peaks(row, annotation_file)
        if df is not None:
            result_dfs.append(df)
    
    # Concatenate the non-empty DataFrames
    if result_dfs:
        concatenated_df = pd.concat(result_dfs, ignore_index=True)
        return concatenated_df
    else:
        return pd.DataFrame()

## STEP 3
# Function to transform rows based on conditions
def transform_row(row):
    if row['proteinA_id'] > row['proteinB_id']:
        return row
    else:
        new_row = row.copy()
        new_row['proteinA_id'], new_row['proteinB_id'] = row['proteinB_id'], row['proteinA_id']
        if new_row['direction'] == 2:
            new_row['direction'] = 3
        return new_row

def grg_parallel(grg_chunk):
    grg_chunk = grg_chunk.swifter.apply(lambda row: transform_row(row), axis=1)
    # grg_chunk = grg_chunk.apply(lambda row: transform_row(row), axis=1)
    return grg_chunk

def calculate_combined_pvalue(group):
    score = combine_pvalues(group['score'], method='fisher')[1]
    pA = group['proteinA_id'].iloc[0]
    pB = group['proteinB_id'].iloc[0]
    metadata = ','.join(group['metadata'].unique().tolist())
    new_row = [pA, pB, score, metadata]
    return new_row

def calculateGRGScore_ForSpecies(tax_id,species_name,species_chromosomes,encode_pipeline,annotation_url,genome_version,evidenceConfig,proteomeConfig):
    species_in_db = Species.objects.get(tax_id=tax_id)
    # Getting mappings from the DB for the specified proteome
    proteome = Proteome.objects.filter(Q(version=proteomeConfig['genome']['version']),Q(species=species_in_db))[0]
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id')).drop_duplicates()
    available_protein_df['mapped_id'] = available_protein_df['mapped_id'].str.lower()

    ## To download encode_files.txt go to 
    # "https://www.encodeproject.org/metadata/?control_type%21=%2A&status=released&perturbed=false&assay_title=TF+ChIP-seq&type=Experiment&files.processed=true"
    # encode_files_df = pd.read_csv('data/tmp/encode/bed_files.txt',header=None,skiprows=1)
    # encode_files_df = encode_files_df[0].apply(lambda row: row.split('/')[-1])
    encode_files = [f.split('.bed.gz')[0] for f in os.listdir('data/tmp/encode/') if '.bed.gz' in f]
    encode_metadata = get_ENCODE(genome_version,species_name,encode_pipeline,fileout='data/tmp/encode/ENCODE_%s.tsv' %tax_id)
    encode_metadata = encode_metadata[encode_metadata['file_accession'].isin(encode_files)].reset_index(drop=True)

    ## Use only default=True
    if  evidenceConfig['use_preferred_default']:
        encode_metadata = encode_metadata[encode_metadata['preferred_default']==True].reset_index(drop=True)
    # encode_metadata.groupby('tf')['experiment_accession'].count().sort_values(ascending=False)
    annotation_file = 'data/tmp/encode/%s.bed' %annotation_url.split('/')[-1]
    extract_annotation(available_protein_df,annotation_url,annotation_file,species_chromosomes)

    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    # parallel = Parallel(n_jobs=cpus,prefer="threads")
    chunk_indices = calculateChunkIndex(len(encode_metadata),cpus)
    # chunk_indices = [(0,2),(2,4),(4,6),(6,8),(8,10)]
    encode_df = pd.concat(parallel(delayed(annotate_peaks_parallel)(encode_metadata.iloc[chunk[0]:chunk[1],:], annotation_file) for chunk in chunk_indices), ignore_index=True)
    # encode_df.drop('pvalue',axis=1,inplace=True)
    print(encode_df)

    available_protein_df.columns = ['proteinA','protein_id']
    encode_df = pd.merge(encode_df, available_protein_df, on='proteinA', how='inner')
    encode_df.drop('proteinA',axis=1,inplace=True)
    encode_df.columns = ['proteinB','score','metadata','proteinA_id']
    available_protein_df.columns = ['proteinB','protein_id']
    encode_df = pd.merge(encode_df, available_protein_df, on='proteinB', how='inner')
    encode_df.drop('proteinB',axis=1,inplace=True)
    encode_df.columns = ['score','metadata','proteinA_id','proteinB_id']
    del available_protein_df
    
    encode_df = encode_df.drop_duplicates().reset_index(drop=True)
    
    ## MAX OUT 
    # Group by 'proteinA_id' and 'proteinB_id', keep the row with max 'score'
    # encode_df = encode_df.groupby(['proteinA_id', 'proteinB_id'])['score'].idxmax().reset_index()
    # Find the index of the row with the maximum 'score' within each group
    # Define a boolean mask for the rows with the maximum score in each group
    mask = encode_df.groupby(['proteinA_id', 'proteinB_id'])['score'].transform(max) == encode_df['score']
    encode_df = encode_df[mask].reset_index(drop=True)
    ## AVERAGE
    # encode_df = encode_df.groupby(['proteinA_id', 'proteinB_id']).agg({'score': 'mean', 'metadata': ','.join}).reset_index()
    
    ## COMBINED PVALUE w/ Fisher method
    # Sort and group by 'proteinA_id' and 'proteinB_id', calculating the combined p-value for each group
    # groupedEncode = encode_df.groupby(['proteinA_id', 'proteinB_id'])
    # with tqdm_joblib(tqdm(desc='combinePvalues', total=len(groupedEncode))) as progress_bar:
    #     encode_df = pd.DataFrame(parallel(delayed(calculate_combined_pvalue)(group) for _,group in groupedEncode))
    # encode_df.columns = ['proteinA_id', 'proteinB_id', 'score', 'metadata']
    print(encode_df)

    ## Direction
    chunk_indices = calculateChunkIndex(len(encode_df),cpus)
    encode_df['direction'] = 2
    encode_df = pd.concat(parallel(delayed(grg_parallel)(encode_df.iloc[chunk[0]:chunk[1],:]) for chunk in chunk_indices))\
                .reset_index(drop=True)
    
    print(encode_df)

    if len(encode_df)>0:
        evidence_in_db = Evidence(type="GRG", species=species_in_db, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidence_in_db.save()

        encode_df = encode_df.sample(frac=1).reset_index(drop=True)
        encode_df.insert(0,'evidence',evidence_in_db.id)

        return encode_df
        
def get_GRGLinks(evidenceConfig,proteomeConfig):
    # grep '.bed.gz' encode_files.txt > bed_files.txt
    # cat bed_files.txt | parallel -j+0 curl -O -J -L {}
    for tax_id,species_name,species_chromosomes,genome_version,encode_pipeline,annotation_url in zip(evidenceConfig['species'],evidenceConfig['species_name'],evidenceConfig['species_chromosomes'],evidenceConfig['genome_version'],evidenceConfig['encode_pipeline'],evidenceConfig['url_annotation']):
        print('Working with species: %s' %tax_id)
        links_grg_df = calculateGRGScore_ForSpecies(tax_id,species_name,species_chromosomes,encode_pipeline,annotation_url,genome_version,evidenceConfig,proteomeConfig)
        if len(links_grg_df)==0: continue
        mem_csv = StringIO()
        links_grg_df.to_csv(mem_csv, index=False)
        mem_csv.seek(0)
        print("Writing %i GRG links to DB for %s" %(len(links_grg_df),tax_id))
        # Writing the csv to the evidenceLink table
        with closing(mem_csv) as csv_io:
            GRG.objects.from_csv(csv_io)


if __name__ == '__main__':
    if os.path.exists('configFiles/exampleGenerateEvidenceSmall'):
        with open('configFiles/exampleGenerateEvidenceSmall/instanceConfig.yml') as f:
            instanceConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/goldStandardConfig.yml') as f:
            goldStandardConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/proteomeConfig.yml') as f:
            proteomeConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/evidenceConfig.yml') as f:
            evidenceConfig = yaml.load(f, Loader=SafeLoader)
        with open('configFiles/exampleGenerateEvidenceSmall/trainingConfig.yml') as f:
            trainingConfig = yaml.load(f, Loader=SafeLoader)
    print(instanceConfig['instance']['species'])
    evidenceConfig = evidenceConfig['GRG']
    tax_id=evidenceConfig['species'][0]
    species_name=evidenceConfig['species_name'][0]
    species_chromosomes = evidenceConfig['species_chromosomes'][0]
    genome_version=evidenceConfig['genome_version'][0]
    encode_pipeline=evidenceConfig['encode_pipeline'][0]
    annotation_url=evidenceConfig['url_annotation'][0]
