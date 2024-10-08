import os, time, sys
from ..externalAlgorithms.anubix import ANUBIX
from ..externalAlgorithms.ease import EASE
import subprocess
import uuid
from statsmodels.stats.multitest import multipletests
import pandas as pd
from data.models import *
from django.db.models import Q

def getAnubix(query_set, sp, cutoff_ppv, fwer_cutoff, enrichedPathways, pathwayMap):
    # Get enrichments using ANUBIX

    resultsDf = ANUBIX(query_set, sp, ppv_cutoff=cutoff_ppv, num_iterations=2000, degree_bin_size=150, optimizeParams=False)
    for index, row in resultsDf.iterrows():
        if float(row["FWER"])<fwer_cutoff:
            if row["pathway_id"] not in enrichedPathways: # Putting an entry for this pathway in the dict
                enrichedPathways[row["pathway_id"]]={
                    'pathway_id':row["pathway_id"],
                    'pathway_name':pathwayMap[row["pathway_id"]],
                    'anubix':"Undefined",
                    'binox':"Undefined",
                    'ease':"Undefined",
                    'min_fwer':float(row["FWER"])
                    # 'min_fdr':float(row["FDR"])
                }
            else:
                if float(row["FWER"])<enrichedPathways[row["pathway_id"]]['min_fwer']: # Checking for lowest FDR
                    enrichedPathways[row["pathway_id"]]['min_fwer']=float(row["FWER"])
            enrichedPathways[row["pathway_id"]]['anubix']={
                'fwer':float(row["FWER"]),
                'fdr':float(row["FDR"]),
                'pval':float(row["Pvalue"]),
                'observed_crosstalk':row["Xtalk"],
                'expected_crosstalk':round(float(row["expectedCrosstalk"]),2)
            }
    return enrichedPathways

def getBinox(query_set, sp, cutoff_ppv, fwer_cutoff, enrichedPathways, pathwayMap):
    # Get enrichments using BINOX

    # Defining precomputed files where data is available (constructed using the network analysis step in funcoup.py)
    binox="website/static/website/binox/BinoX"
    randFile="website/static/website/binox/randNet_"+sp+".csv"
    pathwayFile="website/static/website/binox/pathways_"+sp
    uniqueFileID=str(uuid.uuid4())
    outputFile="website/static/website/binox/binox_"+uniqueFileID
    subnetFile="website/static/website/binox/subnet_"+uniqueFileID

    # First, removing all temp files from the binox dir that are older than 1 day!
    # Doing this in order to not bloat the folder, as files are not properly deleted if a user leaves the page before computations are done!!!
    now = time.time()
    for f in os.listdir("website/static/website/binox/"):
        if (f.startswith("subnet_") or f.startswith("binox_")) and os.stat("website/static/website/binox/"+f).st_mtime < now - 86400:
            os.remove("website/static/website/binox/"+f)

    if not os.path.exists(randFile):
        print("Could not find precomputed random network")
        return enrichedPathways

    # Writing a temporary file for all subnetwork genes
    f = open(subnetFile, "w")
    for g in query_set:
        f.write(str(g)+"\t"+"subnet"+"\n")
    f.close()

    # Run BINOX via commandline
    binox_call = [binox, '-c', str(cutoff_ppv), '-g', "1", '-p', "large", '-r', randFile, '-a', pathwayFile, '-b', subnetFile, '-o', outputFile]
    response = subprocess.run(binox_call, stdout=subprocess.PIPE)
    
    # Read the resulting file
    binoxOut_df = pd.read_csv(outputFile, sep='\t', dtype=str, usecols=[1,3,4,5,14,13])
    p_vals = [float(i) for i in binoxOut_df["#4:p.value"].fillna(1).tolist()]

    # Get fwer
    if len(p_vals)>0:
        _, fwer, _, _ = multipletests(p_vals, method='bonferroni')
        binoxOut_df['fwer'] = fwer

    # Iterate over all resulting enrichments and organize info
    for index, row in binoxOut_df.iterrows():
        if float(row["#5:FDR"])<fwer_cutoff: 
            if row["#6:relationType"]=="+":
                if row["#2:NameGroupA"] not in enrichedPathways:
                    enrichedPathways[row["#2:NameGroupA"]]={
                        'pathway_id':row["#2:NameGroupA"],
                        'pathway_name':pathwayMap[row["#2:NameGroupA"]],
                        'anubix':"Undefined",
                        'binox':"Undefined",
                        'ease':"Undefined",
                        'min_fdr':float(row["#5:FDR"])
                    }
                else:
                    if float(row["#5:FDR"])<enrichedPathways[row["#2:NameGroupA"]]['min_fdr']:
                        enrichedPathways[row["#2:NameGroupA"]]['min_fdr']=float(row["#5:FDR"])
                enrichedPathways[row["#2:NameGroupA"]]['binox']={
                    'fwer':float(row["fwer"]),
                    'fdr':float(row["#5:FDR"]),
                    'pval':float(row["#4:p.value"]),
                    'observed_crosstalk':row["#15:sharedOriglinks"],
                    'expected_crosstalk':round(float(row["#14:expectedLinks"]),2)
                }
    
    # Removing temporary files
    os.remove(subnetFile)
    os.remove(outputFile)

    return enrichedPathways

def getEase(query_set, sp, fwer_cutoff, enrichedPathways, pathwayMap):
    # Get enrichments using EASE (as in DAVID, sort of)

    query_set = [str(x) for x in query_set]
    # Using the same pathway file as for binox to speed things up and not have to read the entire table
    pathwayFile="website/static/website/binox/pathways_"+sp
    genesets_df=pd.read_csv(pathwayFile, dtype=str, sep="\t")
    genesets_dict = genesets_df.groupby('pathway_id')['protein_id'].apply(list).to_dict()
    genome_coverage = Stats.objects.filter(Q(tax_id=sp))[0].genome_size
    
    # Running EASE
    ease_results = EASE(query_set, genesets_dict, PT=genome_coverage)

    # Putting info into datastructure
    for index, row in ease_results.iterrows():
        if float(row["FWER"])<fwer_cutoff:
            if row["pathway_id"] not in enrichedPathways:
                enrichedPathways[row["pathway_id"]]={
                    'pathway_id':row["pathway_id"],
                    'pathway_name':pathwayMap[row["pathway_id"]],
                    'anubix':"Undefined",
                    'binox':"Undefined",
                    'ease':"Undefined",
                    'min_fwer':float(row["FWER"])
                    # 'min_fdr':float(row["FDR"])
                }
            else:
                if float(row["FWER"])<enrichedPathways[row["pathway_id"]]['min_fwer']:
                    enrichedPathways[row["pathway_id"]]['min_fwer']=float(row["FWER"])
            enrichedPathways[row["pathway_id"]]['ease']={
                'fwer':float(row["FWER"]),
                'fdr':float(row["FDR"]),
                'pval':float(row["Pvalue"]),
                'overlap':row["Overlap"]
            }

    return enrichedPathways