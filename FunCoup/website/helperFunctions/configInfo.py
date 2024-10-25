import os
import yaml
from yaml.loader import SafeLoader

def getConfig():
    # Getting information from all configuration files defined as the instance in the websiteConfig
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
    else:
        print("Your config dir does not exist")
        exit()
    if os.path.exists('configFiles/'+instanceName+'/goldStandardConfig.yml'):
        with open('configFiles/'+instanceName+'/goldStandardConfig.yml') as f:
            goldStandardConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your goldStandardConfig.yml does not exist")
        exit()
    if os.path.exists('configFiles/'+instanceName+'/proteomeConfig.yml'):
        with open('configFiles/'+instanceName+'/proteomeConfig.yml') as f:
            proteomeConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your proteomeConfig.yml does not exist")
        exit()
    if os.path.exists('configFiles/'+instanceName+'/evidenceConfig.yml'):
        with open('configFiles/'+instanceName+'/evidenceConfig.yml') as f:
            evidenceConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your evidenceConfig.yml does not exist")
        exit()
    if os.path.exists('configFiles/'+instanceName+'/trainingConfig.yml'):
        with open('configFiles/'+instanceName+'/trainingConfig.yml') as f:
            trainingConfig = yaml.load(f, Loader=SafeLoader)
    else:
        print("Your evidenceConfig.yml does not exist")
        exit()
    instanceConfigurations={'webConfig':webConfig,
                    'instanceConfig':instanceConfig,
                    'proteomeConfig':proteomeConfig,
                    'goldStandardConfig':goldStandardConfig,
                    'evidenceConfig':evidenceConfig,
                    'trainingConfig':trainingConfig}
    return instanceConfigurations