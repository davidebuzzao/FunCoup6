
import os
import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
import pandas as pd
import yaml
from yaml.loader import SafeLoader

def getNetworkStats():
    print("Initiating retrieval of network statistics to fill the stats table")
    if os.path.exists('configFiles/websiteConfig.yml'):
        with open('configFiles/websiteConfig.yml') as f:
            webConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your websiteConfig.yml does not exist")
        exit()
    instanceName=webConfig['instanceName']
    if os.path.exists('configFiles/'+instanceName):
        if os.path.exists('configFiles/'+instanceName+'/instanceConfig.yml'):
            with open('configFiles/'+instanceName+'/instanceConfig.yml') as f:
                instanceConfig = yaml.load(f, Loader=SafeLoader)
        else:
            print("Your config file does not exist")
            exit()

    ## Evidence-based species + Host-pathogen species
    origin_type = 'Evidence'
    fcSpecies = instanceConfig['instance']['species']
    pathogenA = [id[1] for id in instanceConfig['instance']['hostPathogen']]
    for sp in fcSpecies + pathogenA:
        print("Getting stats for "+str(sp))
        numLinks=Network.objects.filter((Q(proteinA__proteome__species__tax_id=sp))).count()
        print("Links: "+str(numLinks))
        if numLinks>0:
            ## species
            sp_name = Species.objects.get(Q(tax_id=sp)).species_name
            sp_common_name = Species.objects.get(Q(tax_id=sp)).common_name
            ## proteome
            proteomeSize=Protein.objects.filter(Q(proteome__species__tax_id=sp)).count()
            ## network
            print("Proteome size: "+str(proteomeSize))
            protA=list(Network.objects.filter((Q(proteinA__proteome__species__tax_id=sp))).order_by().values_list('proteinA', flat=True).distinct())
            protB=list(Network.objects.filter((Q(proteinA__proteome__species__tax_id=sp))).order_by().values_list('proteinB', flat=True).distinct())
            allProts = set(protA+protB)
            numProts=len(allProts)
            print("Nodes: "+str(numProts))
            existingStat=Stats.objects.filter(Q(instance_name=instanceName) & Q(tax_id=sp))
            if existingStat.exists():
                existingStat.update(num_nodes=numProts, num_links=numLinks, genome_size=proteomeSize)
            else:
                stat = Stats(instance_name=instanceName, tax_id=sp, species_name=sp_name, common_name=sp_common_name, num_nodes=numProts, num_links=numLinks, genome_size=proteomeSize, origin=origin_type)
                stat.save()
    
    ## Orthology-based species
    origin_type = 'Orthology'
    transferNetwork_f = 'data/tmp/orthologNetwork/NCBI&InParanoidSpAB.tsv'
    if os.path.isfile(transferNetwork_f):
        species_df = pd.read_csv(transferNetwork_f, sep="\t", dtype=str)
    else:
        print('NCBI&InParanoidSpAB.tsv is missing!!!')
        return

    speciesNameDict = dict(zip(species_df['spA'], species_df['spA_name']))
    species_df = species_df[species_df['FC_species']=='0']
    closest = dict(zip(species_df['spA'], species_df['spB']))
    for species in closest:
        fcSpecies = closest[species]
        fcSpeciesName_name_clean = speciesNameDict[fcSpecies].split(' (')[0]
        fcSpeciesName = '%s.%s' %(fcSpeciesName_name_clean.split()[0][0],fcSpeciesName_name_clean.split()[-1])
        if fcSpeciesName=='O.japonica': fcSpeciesName = 'O.sativa'
        speciesName = speciesNameDict[species].split()[0]+"_"+speciesNameDict[species].split()[1]
        transferred_network_compact_f = "website/static/website/networks/FunCoup6.0/FC6.0_TRANSFERRED_"+speciesName+"("+species+")"+"_from_"+fcSpeciesName+"_compact.gz"        
        # try:
        #     pd.read_csv(transferred_network_compact_f, sep='\t', usecols=[0,1])
        # except:
        #     print(transferred_network_compact_f)
        print(transferred_network_compact_f)
        network_df = pd.read_csv(transferred_network_compact_f, sep='\t', usecols=[0,1])
        network_df.columns = ['gene1','gene2']
        numLinks=len(network_df)
        print("Links: "+str(numLinks))
        if numLinks>0:
            speciesName = speciesNameDict[species].split()[0]+" "+speciesNameDict[species].split()[1]

            ## network
            protA=network_df['gene1'].to_list()
            protB=network_df['gene2'].to_list()
            allProts = set(protA+protB)
            numProts=len(allProts)
            print("Nodes: "+str(numProts))
            proteomeSize = numProts # we don't know the size of the proteome at this point
            existingStat=Stats.objects.filter(Q(instance_name=instanceName) & Q(tax_id=species))
            if existingStat.exists():
                existingStat.update(species_name=speciesName, common_name=speciesName, num_nodes=numProts, num_links=numLinks, genome_size=proteomeSize)
            else:
                stat = Stats(instance_name=instanceName, tax_id=species, species_name=speciesName, common_name=speciesName, num_nodes=numProts, num_links=numLinks, genome_size=proteomeSize, origin=origin_type)
                stat.save()