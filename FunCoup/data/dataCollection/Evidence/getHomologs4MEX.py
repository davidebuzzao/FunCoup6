import gzip
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from joblib import Parallel, delayed

def gzipFile(filename):
    # Unzipping file to path
    with gzip.open(filename+".gz", "r") as gzip_ref:
       gzip_content = gzip_ref.read()
    with open(filename, 'wb') as fout:
        fout.write(gzip_content)

def runDiamond(tax_id,baseDir,diamond_out,cpus,eval_cutoff=0.001):
    # Make database:
    fasta = 'data/tmp/%s.fasta' %(tax_id)
    if os.path.exists(fasta+'.gz'): gzipFile(fasta)
    database = '%s/%s.fasta.db' %(baseDir,tax_id)
    os.system('%s/diamond makedb --in %s -d %s --quiet' %(baseDir,fasta,database))

    # Run diamond
    print('Running DIAMOND Blastp for %s' %tax_id)
    
    # If using bitscore cutoff
    # os.system('%s/diamond blastp -d %s.dmnd -q %s -o %s -f 6 qseqid sseqid bitscore --very-sensitive --quiet --comp-based-stats 0 --min-score %i --max-hsps 0 -p %i' %(baseDir,database,fasta,diamond_out,bs_cutoff,cpus))
    
    # If using E-value cutoff
    os.system('%s/diamond blastp -d %s.dmnd -q %s -o %s -f 6 qseqid sseqid evalue bitscore --very-sensitive --quiet --comp-based-stats 0 --evalue %f --max-hsps 0 -p %i' %(baseDir,database,fasta,diamond_out,eval_cutoff,cpus))

def cleanDiamondOutput(tax_id,diamond_out,proteomeConfig,eval_cutoff=0.001):
    ## Extract known proteins for tax_id under study
    proteome = Proteome.objects.get(Q(version=proteomeConfig['genome']['version']),Q(species__tax_id=tax_id))
    available_protein_df = pd.DataFrame.from_records(IdMapping.objects.filter(Q(protein__proteome=proteome)).values('mapped_id', 'protein_id'))
    available_protein_df.columns = ['geneA','proteinA']

    diamond_df = pd.read_csv(diamond_out, sep='\t', header=None)
    diamond_df.columns = ['geneA','geneB','evalue','bitscore']
    diamond_df['geneA'] = diamond_df['geneA'].apply(lambda x: x.split('|')[1]).astype(str)
    diamond_df['geneB'] = diamond_df['geneB'].apply(lambda x: x.split('|')[1]).astype(str)
    diamond_df = diamond_df[diamond_df['geneA'] != diamond_df['geneB']].reset_index(drop=True)

    # Extra check
    # diamond_df = diamond_df[diamond_df['bitscore']>=bs_cutoff].reset_index(drop=True)
    diamond_df = diamond_df[diamond_df['evalue']<=eval_cutoff].reset_index(drop=True)
    diamond_df.drop('evalue',axis=1,inplace=True)
    diamond_df.drop('bitscore',axis=1,inplace=True)

    ## Filter out protein names that are not in the database
    diamond_df = pd.merge(diamond_df,available_protein_df, on='geneA', how='inner')
    available_protein_df.columns = ['geneB','proteinB']
    diamond_df = pd.merge(diamond_df,available_protein_df, on='geneB', how='inner')
    diamond_df = diamond_df.drop(['geneA','geneB'],axis=1).reset_index(drop=True)

    diamond_df.proteinA, diamond_df.proteinB = np.where(diamond_df.proteinA < diamond_df.proteinB, [diamond_df.proteinB, diamond_df.proteinA], [diamond_df.proteinA, diamond_df.proteinB])
    if len(diamond_df)>0:
        diamond_df.to_csv(diamond_out, sep="\t", header=False, index=False)
    else:
        # print('No proteins pair found at BitScore>%i' %bs_cutoff)
        formatted_eval_cutoff = '{:.0e}'.format(eval_cutoff)
        print('No proteins pair found at E-value<%s' %formatted_eval_cutoff)

# def compareHomologs_BSvsEval(tx,eval_cutoff=0.001,bs_cutoff=100):
#     formatted_eval_cutoff = '{:.0e}'.format(eval_cutoff)
#     baseDir_eval = 'data/tmp/Diamond_eval%s/' %formatted_eval_cutoff
#     baseDir_bs = 'data/tmp/Diamond_bs%i/' %bs_cutoff

#     diamond_out = '%s/diamond%s-%s' %(baseDir_bs,tx,tx)
#     diamond_bs_df = pd.read_csv(diamond_out, sep='\t', header=None)

#     diamond_out = '%s/diamond%s-%s' %(baseDir_eval,tx,tx)
#     diamond_eval_df = pd.read_csv(diamond_out, sep='\t', header=None)

#     return [tx,len(diamond_bs_df),len(diamond_eval_df)]

def getHomologs(tax_id,proteomeConfig,eval_cutoff=0.001):
    # baseDir = 'data/tmp/Diamond_bs%i/' %bs_cutoff
    formatted_eval_cutoff = '{:.0e}'.format(eval_cutoff)
    baseDir = 'data/tmp/Diamond_eval%s/' %formatted_eval_cutoff
    if not os.path.exists(baseDir): 
        os.makedirs(baseDir)
        os.system('wget http://github.com/bbuchfink/diamond/releases/download/v2.0.8/diamond-linux64.tar.gz')
        os.system('tar xzf diamond-linux64.tar.gz --directory=%s' %baseDir)
        os.system('rm diamond-linux64.tar.gz')
    
    cpus = os.cpu_count()
    ## 
    for tx in tax_id:
        diamond_out = '%s/diamond%s-%s' %(baseDir,tx,tx)
        if not os.path.exists(diamond_out): 
            runDiamond(tx,baseDir,diamond_out,cpus,eval_cutoff=eval_cutoff)
            cleanDiamondOutput(tx,diamond_out,proteomeConfig,eval_cutoff=eval_cutoff)
        else:
            print('DIAMOND Blastp for %s exists' %tx)

    ## We want to compare the number of homologs using bitscore>100 and Evalue<1e-3
    plot_comparison = False
    if plot_comparison:
        diamond_comparison = []
        for tx in tax_id:
            diamond_comparison.append(compareHomologs_BSvsEval(tx))
    
        if len(diamond_comparison)>0:
            diamond_comparison_df = pd.DataFrame(diamond_comparison)
            diamond_comparison_df.columns = ['Species','BitScore>100','Evalue<1e-3']

            # Set up the figure and axis
            fig, ax = plt.subplots()
            # Define bar width and positions
            bar_width = 0.35
            x = range(len(diamond_comparison_df))

            # Plot the bars
            bar1 = ax.bar(x, diamond_comparison_df['BitScore>100'], width=bar_width, label='BitScore>100')
            bar2 = ax.bar([pos + bar_width for pos in x], diamond_comparison_df['Evalue<1e-3'], width=bar_width, label='Evalue<1e-3')

            # Set the x-axis labels and tick positions
            ax.set_xticks([pos + bar_width/2 for pos in x])
            ax.set_xticklabels(diamond_comparison_df['Species'], rotation=45)  # Rotate labels by 60 degrees

            # Add labels and title
            ax.set_ylabel('Nr. of Homologs')
            ax.set_xlabel('Species')

            # Add a legend
            ax.legend()

            # Set the figure size in inches (10cm width and 5cm height)
            fig.set_size_inches(10, 5)  # 1 inch = 2.54 cm

            # Save the figure as a PNG image
            plt.savefig('data/stats/diamondHomologs_BSvsEVAL.png', dpi=300, bbox_inches='tight')  # Adjust dpi as needed

            # Show the plot (optional)
            # plt.show()

