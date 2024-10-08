import os
from sklearn.preprocessing import MinMaxScaler,PowerTransformer
import yaml
from yaml.loader import SafeLoader
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from contextlib import closing
from io import StringIO
from tqdm import tqdm
from auxiliary_functions import *

from data.dataCollection.Genomes.getGenomes import *
from data.dataCollection.Evidence.PIN import *

def getPINScoresForPathogen(iRefIndex_df,host_id,pathogen_id,evidenceConfig,proteomeConfig):
    ## Parallel config
    cpus = os.cpu_count()
    parallel = Parallel(n_jobs=cpus)
    # parallel = Parallel(n_jobs=cpus,prefer="threads")
    
    species_in_db = Species.objects.get(tax_id=pathogen_id)
    print('Extracting PIN for ' + pathogen_id)
    set(iRefIndex_df['taxA'].to_list())
    iRefIndex_df = iRefIndex_df.drop('taxA',axis=1).drop('taxB',axis=1)

    # Collecting uniProt identifier OR complex id (for complexes). Excluding rows without uniProt identifiers
    iRefIndex_df['proteinA']=iRefIndex_df['idA']+"|"+iRefIndex_df['idA_alt1']# Joining all potentail identifier locations for idA
    iRefIndex_df = iRefIndex_df[(iRefIndex_df.proteinA.str.contains("uniprotkb:") | iRefIndex_df.proteinA.str.contains("complex:"))] # Keeping rows having uniprot OR complex in idA
    iRefIndex_df['proteinA'] = iRefIndex_df['proteinA'].apply(lambda x: x.split("uniprotkb:")[1].split("|")[0] if 'complex:' not in x else x) # Extracting the uniprotID from the joined column (except for complexes)
    iRefIndex_df.proteinA = np.where(iRefIndex_df.idA.str.contains('complex:'), iRefIndex_df.idA, iRefIndex_df.proteinA) # Setting proteinA to complex ID for complexes
    iRefIndex_df['proteinB']=iRefIndex_df['idB']+"|"+iRefIndex_df['idB_alt1'] # Joining all potential identifier locations for idB
    iRefIndex_df = iRefIndex_df[iRefIndex_df.proteinB.str.contains("uniprotkb:")] # Keeping wors having a uniprot id in idB
    iRefIndex_df['proteinB'] = iRefIndex_df['proteinB'].apply(lambda x: x.split("uniprotkb:")[1].split("|")[0]) # Extracting the uniprotID for proteinB
    iRefIndex_df = iRefIndex_df.drop(['idA','idA_alt1','idB','idB_alt1'],axis=1) # Drop cvolumns that are not needed anymore

    iRefIndex_df['experimentalMethod'] = iRefIndex_df['experimentalMethod'].apply(lambda x: x.split("(")[0].replace("psi-mi:","").replace("\"","")) # Clean up experimental method
    iRefIndex_df.experimentalMethod = np.where(iRefIndex_df.experimentalMethod == "-", "MI:0000", iRefIndex_df.experimentalMethod) # If experimental method is -, setting it to term for unspecified instead
    iRefIndex_df['pmid'] = iRefIndex_df['pmid'].apply(lambda x: x.split("|")[-1]) # cleaning up pubmed id
    iRefIndex_df = iRefIndex_df[iRefIndex_df.pmid != "-"] # Only using rows where a pubmed-id is specified
    ## Unreviewd research starts with pubmed:8888 https://wiki.thebiogrid.org/doku.php/covid:unpublished
    iRefIndex_df = iRefIndex_df[~iRefIndex_df['pmid'].str.contains('pubmed:8888')]

    # Extracting complexes into pairs --> TODO parallelize?
    iRefIndex_df = iRefIndex_df.drop_duplicates() # Drop duplicates
    complexDf = iRefIndex_df[iRefIndex_df.proteinA.str.contains("complex:")] # Extracting complex rows from the dataframe
    iRefIndex_df = iRefIndex_df[~iRefIndex_df.proteinA.str.contains("complex:")] # Removing complex rows from the pairwise dataframe
    groupedComplexes = complexDf.groupby(['proteinA','pmid','experimentalMethod']) # Group complexes by complexID-pmid-experimental method

    with tqdm_joblib(tqdm(desc='disentangle_complex', total=len(groupedComplexes))) as progress_bar:
        complex_df = pd.concat(parallel(delayed(disentangle_complex_parallel)(name,group) for name,group in groupedComplexes))
    complex_df.columns = ['experimentalMethod','pmid','proteinA','proteinB']
    iRefIndex_df = pd.concat([iRefIndex_df, complex_df]).reset_index(drop=True)

    # Getting mappings from the DB for the specified proteome
    proteome_covid = Proteome.objects.get(Q(version=proteomeConfig['genome']['version']),Q(species=species_in_db))
    covid_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome_covid)).values('mapped_id', 'protein_id')).drop_duplicates()
    covid_protein_df['tax_id'] = pathogen_id
    ## Host organisms
    species_human = Species.objects.get(tax_id=host_id)
    proteome_human = Proteome.objects.get(Q(version=proteomeConfig['genome']['version']),Q(species=species_human))
    human_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome_human)).values('mapped_id', 'protein_id')).drop_duplicates()
    human_protein_df['tax_id'] = host_id

    available_protein_df = pd.concat([covid_protein_df,human_protein_df],axis=0)
    available_protein_df.columns = ['proteinA','protein_id','taxA_id']

    iRefIndex_df = pd.merge(iRefIndex_df, available_protein_df, on='proteinA', how='inner')
    iRefIndex_df.drop('proteinA',axis=1,inplace=True)
    iRefIndex_df.columns = ['experimentalMethod','pmid','proteinB','proteinA_id','taxA_id']
    available_protein_df.columns = ['proteinB','protein_id','taxB_id']
    iRefIndex_df = pd.merge(iRefIndex_df, available_protein_df, on='proteinB', how='inner')
    iRefIndex_df.drop('proteinB',axis=1,inplace=True)
    iRefIndex_df.columns = ['experimentalMethod','pmid','proteinA_id','taxA_id','proteinB_id','taxB_id']
    del available_protein_df

    ## Delete links between human genes, if any
    ## TODO ---> CHECK THIS OUT
    iRefIndex_df = iRefIndex_df[~((iRefIndex_df['taxA_id'] == host_id) & (iRefIndex_df['taxB_id'] == host_id))]
    iRefIndex_df.drop(['taxA_id','taxB_id'],axis=1,inplace=True)

    # Put IDs in the correct order
    iRefIndex_df['proteinA_id'] = iRefIndex_df['proteinA_id'].astype(int)
    iRefIndex_df['proteinB_id'] = iRefIndex_df['proteinB_id'].astype(int)
    iRefIndex_df.proteinA_id, iRefIndex_df.proteinB_id = np.where(iRefIndex_df.proteinA_id < iRefIndex_df.proteinB_id, [iRefIndex_df.proteinB_id, iRefIndex_df.proteinA_id], [iRefIndex_df.proteinA_id, iRefIndex_df.proteinB_id])
    iRefIndex_df = iRefIndex_df.drop_duplicates() # Drop duplicates
    iRefIndex_df = iRefIndex_df[iRefIndex_df['proteinA_id']!=iRefIndex_df['proteinB_id']]
    
    ## Measure gene frequency
    ## How many distinct links? Remove duplicated lines on ['pmid','proteinA_id','proteinB_id'] and extract frequencies
    print('# Measuring gene and pmid frequencies')
    noduplicates_iRefIndex_df = iRefIndex_df.drop_duplicates(subset=['pmid', 'proteinA_id', 'proteinB_id'])
    genes = noduplicates_iRefIndex_df['proteinA_id'].to_list() + noduplicates_iRefIndex_df['proteinB_id'].to_list()
    gene_dict = compute_frequency(genes)
    pmids = noduplicates_iRefIndex_df['pmid'].to_list()
    pmid_dict = compute_frequency(pmids)

    iRefIndex_df.drop('experimentalMethod',axis=1,inplace=True)
    iRefIndex_df.columns = ['metadata','proteinA_id','proteinB_id']
    iRefIndex_df = iRefIndex_df.groupby(['proteinA_id','proteinB_id'])['metadata'].apply(lambda x: ','.join(x)).reset_index()
    
    # Scoring PIN
    print('# Scoring PIN')
    meta_dict = {'proteinA_id': 'int', 'proteinB_id': 'int', 'metadata':'object', 'score': 'float'}
    iRefIndex_ddf = dd.from_pandas(iRefIndex_df, npartitions=288)  # Adjust the number of partitions as needed
    result_ddf = iRefIndex_ddf.apply(lambda row: calculatePINscore(row,gene_dict,pmid_dict), axis=1, meta=meta_dict)
    with dask.config.set(scheduler='processes'):  # single-threaded, processes or threads
        iRefIndex_df = result_ddf.compute() 
    
    ## To check for SPIKE-ACE2 interaction
    # iRefIndex_df[(iRefIndex_df['proteinA_id']==351481) & (iRefIndex_df['proteinB_id']==334772)]#.sort_values('score')
    # iRefIndex_df[(iRefIndex_df['proteinA_id']==351481)].sort_values('score')
    
    if len(iRefIndex_df)>0:
        scores = iRefIndex_df['score'].values.reshape(-1, 1)
        # Power Transformation https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.PowerTransformer.html#sklearn.preprocessing.PowerTransformer
        # Apply a power transform featurewise to make data more Gaussian-like.
        power_scaler = PowerTransformer(method='box-cox')
        power_scores = power_scaler.fit_transform(scores)
        # Min-Max scaling
        minmax_scaler = MinMaxScaler(feature_range=(0, 1))
        minmax_power_scores = minmax_scaler.fit_transform(power_scores)
        iRefIndex_df['score'] = np.round(minmax_power_scores,decimals=4)

        evidence_in_db = Evidence(type="PIN", species=species_in_db, version=evidenceConfig['version'], scoringMethod=evidenceConfig['scoring_method'])
        evidence_in_db.save()

        iRefIndex_df = iRefIndex_df.sample(frac=1).reset_index(drop=True)
        iRefIndex_df.insert(0,'evidence',evidence_in_db.id)

        return iRefIndex_df

def getPINLinks4Pathogen(pathogenToCollect,evidenceConfig,proteomeConfig):
    # TODO: This will work with iRefIndex files and for COVID, to change for different data sources
    for hostPathogen_pair in pathogenToCollect:
        host_id = hostPathogen[0]
        pathogen_id = hostPathogen[1]
        pathogen_index = evidenceConfig['pathogens'].index(pathogen_id)
        
        # ebiIdentifier = 'UP000464024'
        ### ADD SPECIES
        ebiIdentifier = readREADME([pathogen_id])
        ### ADD PROTEOME
        #s, ebiIdentifiers, proteomeConfig,mappingsOfInterest  = pathogen_id, ebiIdentifier, proteomeConfig, proteomeConfig['genome']['mappingsOfInterest']
        getProteomesAndMappings(pathogen_id, ebiIdentifier, proteomeConfig, proteomeConfig['genome']['mappingsOfInterest'])
        
        ### ADD EVIDENCE 
        # Collecting PIN links from the iRefIndex file. Using complexes as well as pairwise interactions. Caclulates a score for each link
        # based on the number of papers, the number of assays, and the number of members in each assay supporting the link.
        iRefIndexFile = "data/tmp/"+evidenceConfig['filePathogen'][pathogen_index]
        iRefIndex_df = pd.read_csv(iRefIndexFile, sep='\t',usecols=[0,1,2,3,6,8,9,10], header=None)
        print('Dataset ' + iRefIndexFile + ' open!')
        iRefIndex_df.columns = ['idA','idB','idA_alt1','idB_alt1','experimentalMethod','pmid','taxA','taxB']
        iRefIndex_df.taxA = np.where(iRefIndex_df.taxA == "-", iRefIndex_df.taxB, iRefIndex_df.taxA) # For complexes, taxa is -, then set it to the same as taxb
        # iRefIndex_df = iRefIndex_df[iRefIndex_df['taxA'] == iRefIndex_df['taxB']].reset_index(drop=True) # Checking that the two interacting proteins belong to the same species
        
        # Define a regular expression pattern to match the tax ID
        pattern = r'taxid:(\d+)'
        # Use str.extract to extract the tax ID using the defined pattern
        iRefIndex_df['taxA'] = iRefIndex_df['taxA'].str.extract(pattern)
        iRefIndex_df['taxB'] = iRefIndex_df['taxB'].str.extract(pattern)
        # Convert the extracted tax ID to integer type
        iRefIndex_df['taxA'] = iRefIndex_df['taxA'].astype(str)
        iRefIndex_df['taxB'] = iRefIndex_df['taxB'].astype(str)
        
        iRefIndex_df = iRefIndex_df[iRefIndex_df['taxA'].isin(hostPathogen_pair) & iRefIndex_df['taxB'].isin(hostPathogen_pair)]

        links_pin_df = getPINScoresForPathogen(iRefIndex_df,host_id,pathogen_id,evidenceConfig,proteomeConfig)
        # links_pin_df.sort_values('score')
        if len(links_pin_df)>0:
            mem_csv = StringIO()
            links_pin_df.to_csv(mem_csv, index=False)
            mem_csv.seek(0)
            print("Writing %i PIN links to DB for %s" %(len(links_pin_df),pathogen_id))
            # Writing the csv to the evidenceLink table
            with closing(mem_csv) as csv_io:
                PIN.objects.from_csv(csv_io)

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

    hostPathogen = instanceConfig['instance']['hostPathogen'][0]
    host_id = hostPathogen[0]
    pathogen_id = hostPathogen[1]
    evidenceConfig = evidenceConfig['PIN']

