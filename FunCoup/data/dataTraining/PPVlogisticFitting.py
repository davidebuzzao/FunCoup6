import os
import time
import pandas as pd
from scipy.optimize import curve_fit
import random
import numpy as np
import yaml
from yaml.loader import SafeLoader
import matplotlib.pyplot as plt
from sklearn.metrics import precision_score

from joblib import Parallel, delayed
from tqdm import tqdm

import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps
from django.conf import settings

from data.dataTraining.ExtractOrthologLinks import *
from data.dataTraining.NegativeGoldStandard import *
from auxiliary_functions import *

# Define a function to compute precision for a given threshold
def computePPV(fbs_threshold, fbs_sorted, y_class_sorted):
    y_pred = (fbs_sorted > fbs_threshold).astype(int)
    ppv = precision_score(y_class_sorted, y_pred,zero_division=1)
    return ppv

def logistic(x, a, b, c):
    return a / (1 + np.exp(-b * (x - c)))

def fit_logistic(df, x='fbs_max', y='PPV'):
    # Extract data for fitting
    xdata = df[x]
    ydata = df[y]
    # Fit the exponential function with bounds
    try: params, covariance = curve_fit(logistic, xdata, ydata, maxfev=10000)
    except:
        print('No fit!!!')
        params = None
    else: return params

def formatLogistic(params):
    a, b, c = params
    # Format the string representation of the logistic function
    equation = f"{np.round(a, 3)} / (1 + exp(-{np.round(b, 3)} * (x - {np.round(c, 3)})))"
    # Return the reconstructed logistic function equation
    return equation

def write_logPPVtoDatabase(goldstandard_t,params,motivation=None):
    if motivation: expPFC = LogisticPPV(goldStandard=goldstandard_t, function=','.join([str(i) for i in params]), error=motivation)
    else: expPFC = LogisticPPV(goldStandard=goldstandard_t, function=','.join([str(i) for i in params]))
    expPFC.save()

def logisticPPV(positive_goldStandardLink,negative_goldStandardLink,goldstandard_t,speciesA,g_type,trainingConfig,saveToDatabase=True):
    VisualizePPV = trainingConfig['VisualizePPV']
    boot,boot_tot = trainingConfig['Bootstrapping']
    min_r2 = trainingConfig['R2']
    lpos,lneg = len(positive_goldStandardLink),len(negative_goldStandardLink)

    # Compute Precision at different thresholds
    n_thr = 1000
    # Parallelize the loop
    cpus = os.cpu_count()
    num_cores = min(cpus,n_thr)  # Use all available cores, adjust as needed
    parallel = Parallel(n_jobs=num_cores)
    # parallel = Parallel(n_jobs=num_cores,prefer="threads")
    if boot: 
        pfc_df = pd.DataFrame()
        for nseed in range(1,boot_tot+1):
            np.random.seed(int(nseed))
            tmp_negative_goldStandardLink = random.sample(negative_goldStandardLink,np.min([lpos,lneg]))
            lneg = len(tmp_negative_goldStandardLink)
            y_class = np.concatenate((np.ones(len(positive_goldStandardLink)), np.zeros(len(tmp_negative_goldStandardLink))))
            # Sort the indices based on predicted probabilities (fbs)
            fbs = np.array(positive_goldStandardLink+tmp_negative_goldStandardLink)
            indices = np.argsort(fbs)[::-1]

            # Sort y_class and fbs based on the sorted indices
            y_class_sorted = y_class[indices]
            fbs_sorted = fbs[indices]
            # fbs_thresholds = np.linspace(0, np.max(fbs_sorted), n_thr)
            fbs_thresholds = np.linspace(np.min(fbs_sorted), np.max(fbs_sorted), n_thr)

            with tqdm_joblib(tqdm(desc='%s %i/%i' %(g_type,nseed,boot_tot), total=n_thr)) as progress_bar:
                ppv_values = parallel(delayed(computePPV)(threshold, fbs_sorted, y_class_sorted) for threshold in fbs_thresholds[:])
            tmp_pfc_df = pd.DataFrame.from_dict({'FBS': fbs_thresholds, 'PPV': ppv_values})
            # tmp_pfc_df = tmp_pfc_df.iloc[:-1]
            pfc_df = pd.concat([pfc_df,tmp_pfc_df],axis=0)
    else:
        lneg = len(negative_goldStandardLink)
        y_class = np.concatenate((np.ones(len(positive_goldStandardLink)), np.zeros(len(negative_goldStandardLink))))
        fbs = np.array(positive_goldStandardLink+negative_goldStandardLink)
        indices = np.argsort(fbs)[::-1]

        # Sort y_class and fbs based on the sorted indices
        y_class_sorted = y_class[indices]
        fbs_sorted = fbs[indices]
        # fbs_thresholds = np.linspace(0, np.max(fbs_sorted), n_thr)
        fbs_thresholds = np.linspace(np.min(fbs_sorted), np.max(fbs_sorted), n_thr)

        with tqdm_joblib(tqdm(desc='%s' %(g_type), total=n_thr)) as progress_bar:
            ppv_values = parallel(delayed(computePPV)(threshold, fbs_sorted, y_class_sorted) for threshold in fbs_thresholds[:])
        pfc_df = pd.DataFrame.from_dict({'FBS': fbs_thresholds, 'PPV': ppv_values})
        # pfc_df = pfc_df.iloc[:-1]
        
    pfc_df = pfc_df.sort_values(by='PPV').reset_index(drop=True)

    ###### SIGMOIDAL FUNCTION
    logistic_params = fit_logistic(pfc_df, x='FBS', y='PPV')
    if logistic_params is None: return

    pfc_df['pPPV'] = logistic(pfc_df['FBS'], *logistic_params)

    pfc_df['residuals'] = pfc_df['PPV'] - pfc_df['pPPV']
    ss_res = np.sum(pfc_df['residuals']**2)
    ss_tot = np.sum((pfc_df['PPV'] - np.mean(pfc_df['PPV']))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R2: %f' %(np.round(r_squared,2)))

    if VisualizePPV and saveToDatabase:
        fbs = np.array(positive_goldStandardLink+negative_goldStandardLink)
        fbs_thresholds = np.linspace(np.min(fbs), np.max(fbs), n_thr)
        print('4: visualizing PPV')
        ## VISUALIZATION
        species_objects = Species.objects.get(Q(tax_id=speciesA))
        speciesA_name = species_objects.species_name
        speciesA_name_compact = '%s.%s' %(speciesA_name.split()[0][0],speciesA_name.split()[1])
        g_type = goldstandard_t.type

        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))  # Adjust figsize as needed
        fig.suptitle('%s (%s) - %s (%i+, %i-)' %(speciesA_name_compact,speciesA,g_type,lpos,lneg), fontsize=14)

        # Scatter plot of pos_kde and neg_kde after interpolation of common x-range
        x_min, x_max = np.min(pfc_df['FBS']), np.max(pfc_df['FBS'])
        axes[0].hist(positive_goldStandardLink, bins=50, facecolor='#003049', label='Positive Scores', alpha=0.5, density=True, stacked=True)
        axes[0].hist(negative_goldStandardLink, bins=50, facecolor='#c1121f', label='Negative Scores', alpha=0.5, density=True, stacked=True)
        axes[0].set_xlabel('FBS_%s' %g_type)
        axes[0].set_ylabel('Density of links')
        axes[0].set_xlim(x_min, x_max)
        axes[0].legend()


        # Scatter plot of common x-range
        x_min, x_max = np.min(pfc_df['FBS']), np.max(pfc_df['FBS'])
        axes[1].scatter(pfc_df['FBS'], pfc_df['PPV'], color='#a7c957', s=10, alpha=0.5, label='Data Points')
        axes[1].plot(fbs_thresholds, logistic(fbs_thresholds, *logistic_params), label='Logistic fit R2=%.2f' %r_squared, color='#3a86ff')
        axes[1].set_xlabel('FBS_%s' %g_type)
        axes[1].set_ylabel('PPV_%s' %g_type)
        axes[1].set_xlim(x_min, x_max)
        axes[1].legend()
        axes[1].set_title(formatLogistic(logistic_params), fontsize=10)


        # Create scatter plot of the data points
        fbs_to_pfc = np.round(np.clip(logistic(fbs,*logistic_params),0.5,1),3)
        cdf_values, bin_edges, _ = axes[2].hist(fbs_to_pfc, bins=1000, color='skyblue', edgecolor='black', cumulative=True, density=False, alpha=0)
        # Plot the complement of the CDF
        complement_cdf_values = len(fbs_to_pfc) - cdf_values
        ## To align the center of the bins with the points --> bin_edges[1:]
        axes[2].bar(bin_edges[1:], complement_cdf_values, width=bin_edges[1] - bin_edges[0], color='#e63946', edgecolor='#a7c957')
        # Set the y-axis to log scale
        axes[2].set_yscale('log')
        axes[2].grid(True, linestyle='--', alpha=0.5)
        axes[2].set_xlabel('PPV_%s' %g_type)
        axes[2].set_ylabel('Nr. of links')
        axes[2].set_xlim(0.5,1)
        axes[2].set_xticks(np.arange(0.5, 1.01, 0.1))


        # Adjust spacing between subplots
        plt.tight_layout()
        plt.savefig('data/tmp/PPV/%s/%s_%s-PPV.png' %(speciesA,speciesA,g_type))
        plt.close()
    
    if saveToDatabase: 
        if r_squared<min_r2: write_logPPVtoDatabase(goldstandard_t,['Nan'],motivation='No good Logistic fit at R2≥%.2f' %min_r2)
        else: write_logPPVtoDatabase(goldstandard_t,logistic_params)
    return ','.join([str(i) for i in logistic_params])