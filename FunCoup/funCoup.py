import os
import django
from django.conf import settings
import argparse
import yaml
from yaml.loader import SafeLoader
from auxiliary_functions import *

# Data Collection modules for retrieving genomes, orthologs, gold standards, and evidence
from data.dataCollection.Genomes.getGenomes import getGenome
from data.dataCollection.Genomes.getOrthologs import getOrthologs
from data.dataCollection.GoldStandards.getGoldStandards import getGoldStandards
from data.dataCollection.Evidence.getEvidence import getEvidence

# Network analysis modules for transferring networks and getting pathway, tissue, and statistical data
from data.networkAnalysis.getTransferredNetworks import transferNetworks
from data.networkAnalysis.getKeggPathwayInformation import getKeggPathwayInformation
from data.networkAnalysis.getTissueInformation import getTissueInformation
from data.networkAnalysis.precomputeAnubix import precomputeAnubix
from data.networkAnalysis.precomputeBinox import precomputeBinox
from data.networkAnalysis.getNetworkStats import getNetworkStats

# Data Training modules for training networks and host-pathogen networks
from data.dataTraining.getNetwork import trainNetwork
from data.dataTraining.getHostPathogenNetworks import trainHostPathogenNetwork

# Data Query modules for computing and deleting link score statistics
from data.dataQuery.CountEntry import computeLinksScoreStats
from data.dataQuery.DeleteEntry import deleteLinksScoreStats

# Django setup to initialize the project's settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()

# Function to generate a new instance by fetching genomes, orthologs, gold standards, and evidence
def generateInstance(instanceConfig,goldStandardConfig,proteomeConfig,evidenceConfig):
    getGenome(instanceConfig['instance']['species'],proteomeConfig)
    getOrthologs(instanceConfig['instance']['species'],proteomeConfig)
    getGoldStandards(instanceConfig,goldStandardConfig,proteomeConfig)
    getEvidence(instanceConfig,evidenceConfig,proteomeConfig)

# Function to train networks for the generated instance, including host-pathogen networks
def trainInstance(instanceConfig,goldStandardConfig,proteomeConfig,evidenceConfig,trainingConfig):
    trainNetwork(instanceConfig,goldStandardConfig,evidenceConfig,proteomeConfig,trainingConfig)
    trainHostPathogenNetwork(instanceConfig,goldStandardConfig,evidenceConfig,proteomeConfig,trainingConfig)

# Function to perform network analysis by transferring networks and retrieving pathway and tissue information
def generateNetworkAnalysis():
    transferNetworks()
    getKeggPathwayInformation()
    getTissueInformation()
    getNetworkStats()
    precomputeAnubix()
    precomputeBinox()

# Function to compute and delete link score statistics in the database
def generateStats(instanceConfig,goldStandardConfig,proteomeConfig,evidenceConfig):
    computeLinksScoreStats(instanceConfig,goldStandardConfig,proteomeConfig,evidenceConfig,opt='count_goldstandard')
    deleteLinksScoreStats(instanceConfig,goldStandardConfig,proteomeConfig,evidenceConfig,opt='0links')

# Main function that parses command-line arguments and runs the appropriate task
def main():
    parser = argparse.ArgumentParser(description='''
    FunCoup - Functional Coupling
    ''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-g", "--generateInstance",
                        help="generate an instance by fetching data or using existing data in the DB",
                        default=False, action='store_true')
    parser.add_argument("-t", "--train",
                        help="Train an instance",
                        default=False, action='store_true')
    parser.add_argument("-s", "--generateStats",
                        help="Generate some statistics of the database",
                        default=False, action='store_true')
    parser.add_argument("-a", "--generateNetworkAnalysis",
                        help="Generate all network analysis, including transferred networks, fetching pathway and tissue information",
                        default=False, action='store_true')
    parser.add_argument("-c", "--configDir",
                        help="Name of config directory to use",
                        type=str)
    args = parser.parse_args()


    # Validation of arguments and existence of configuration files
    if not (args.generateInstance or args.train or args.generateStats or args.generateNetworkAnalysis):
        print("Please specify a task: \n-g to generate an instance, \n-t to train an instance, \n-n to construct a network from an instance")
        exit()

    if args.configDir is None:
        print("Please specify a config directory. Config files should be YAML files located under data/configFiles/. See examples in exampleGenerate/ for instanceConfig.yml, evidenceConfig.yml, goldStandardConfig.yml, and proteomeConfig.yml.")
        exit()

    # Load configuration files
    config_path = 'configFiles/' + args.configDir
    if os.path.exists(config_path):
        if os.path.exists(config_path + '/instanceConfig.yml'):
            with open(config_path + '/instanceConfig.yml') as f:
                instanceConfig = yaml.load(f, Loader=SafeLoader)
        else:
            print("Your config file does not exist")
            exit()
    else:
        print("Your config directory does not exist")
        exit()

    if os.path.exists(config_path + '/goldStandardConfig.yml'):
        with open(config_path + '/goldStandardConfig.yml') as f:
            goldStandardConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your goldStandardConfig.yml does not exist")
        exit()

    if os.path.exists(config_path + '/proteomeConfig.yml'):
        with open(config_path + '/proteomeConfig.yml') as f:
            proteomeConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your proteomeConfig.yml does not exist")
        exit()

    if os.path.exists(config_path + '/evidenceConfig.yml'):
        with open(config_path + '/evidenceConfig.yml') as f:
            evidenceConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your evidenceConfig.yml does not exist")
        exit()

    if os.path.exists(config_path + '/trainingConfig.yml'):
        with open(config_path + '/trainingConfig.yml') as f:
            trainingConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your trainingConfig.yml does not exist")
        exit()

    # Execution of the selected task
    if args.generateInstance and not args.train:
        print("Starting instance generation")
        generateInstance(instanceConfig, goldStandardConfig, proteomeConfig, evidenceConfig)
    elif args.train and not args.generateInstance:
        print("Starting training")
        trainInstance(instanceConfig, goldStandardConfig, proteomeConfig, evidenceConfig, trainingConfig)
    elif args.generateStats:
        print("Querying the database")
        generateStats(instanceConfig, goldStandardConfig, proteomeConfig, evidenceConfig)
    elif args.generateNetworkAnalysis:
        print("Starting generation of network analysis, including network transfer and pathway/tissue information")
        generateNetworkAnalysis()
    else:
        print("Please select only one run option")

if __name__ == '__main__':
    main()