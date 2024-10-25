import os
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

import django
from django.conf import settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FunCoup.settings')
django.setup()
from data.models import *
from django.db.models import Q
from django.apps import apps

from data.dataTraining.NegativeGoldStandard import *
from data.dataTraining.ExtractOrthologLinks import *
from auxiliary_functions import *

# Define a function to remove outliers using the IQR method
def remove_outliers(scores, lower_percentile, upper_percentile):
    # Calculate the lower and upper values based on percentiles
    lower_value = np.percentile(scores, lower_percentile)
    upper_value = np.percentile(scores, upper_percentile)

    # Filter the data to keep values within the specified percentiles
    return [s for s in scores if lower_value <= s <= upper_value]

def formatPolynomial(coefficients):
    # Get the degree of the polynomial
    degree = len(coefficients) - 1
    # Format the string representation of the equation with aligned exponents
    equation = ""
    for i in range(degree, -1, -1):
        coef = round(coefficients[degree - i],1)
        if coef != 0:
            if i == degree:
                equation += f"{coef}x^{i}"
            elif i == 0:
                equation += f" + {coef}"
            else:
                equation += f" + {coef}x^{i}"
    # Return the reconstructed polynomial equation
    return equation

def write_pLLRToDatabase(evidence_t,goldstandard_t,coefficients,training_params=None,motivation=None):
    if motivation:
        pLLR = PolynomialLLR(evidence=evidence_t, goldStandard=goldstandard_t, function=','.join([str(i) for i in coefficients]), error=motivation)
    else:
        x_min,x_max,pbw,nbw = training_params
        pLLR = PolynomialLLR(evidence=evidence_t, goldStandard=goldstandard_t, function=','.join([str(i) for i in coefficients]), scoreRange='%s|%s' %(str(x_min),str(x_max)), bw='%s|%s' %(str(pbw),str(nbw)))
    pLLR.save()

def polynomialLLR(e_id,positive_goldStandardLink,goldstandard_t,negative_goldStandardLink,trainingConfig,saveToDatabase=True,score_cutoff=None,test=False):

    ### KDE config
    min_n_sample, max_n_sample = trainingConfig['nSample']
    lower_percentile, upper_percentile = trainingConfig['outliersPercentiles']
    bw_method = trainingConfig['Bandwidth']
    min_r2 = trainingConfig['R2']
    poly_degree_min, poly_degree_max = trainingConfig['poly_degree']
    visualize_LLR = trainingConfig['VisualizeLLR']

    evidence_t = Evidence.objects.get(Q(id=e_id))
    species_evidence = evidence_t.species.tax_id
    species_golstandard = goldstandard_t.species.tax_id

    ####################
    ## Positive Gold Standard links
    if positive_goldStandardLink is None or len(positive_goldStandardLink)<min_n_sample:
        if saveToDatabase: write_pLLRToDatabase(evidence_t,goldstandard_t,['Nan'],motivation='|PGS|<%i' %min_n_sample)
        return

    positive_goldStandardLink = remove_outliers(positive_goldStandardLink,lower_percentile,upper_percentile)

    if positive_goldStandardLink is None or len(positive_goldStandardLink)<min_n_sample:
        if saveToDatabase: write_pLLRToDatabase(evidence_t,goldstandard_t,['Nan'],motivation='|PGS|<%i' %min_n_sample)
        print('Too few PGS!')
        return

    if len(positive_goldStandardLink)>max_n_sample:
        random.seed(int(species_golstandard))
        positive_goldStandardLink = random.sample(positive_goldStandardLink, max_n_sample)

    ####################
    ## Negative Gold Standard links
    if negative_goldStandardLink is None or len(negative_goldStandardLink)<min_n_sample:
        if saveToDatabase: write_pLLRToDatabase(evidence_t,goldstandard_t,['Nan'],motivation='|NGS|<%i' %min_n_sample)
        print('Too few NGS!')
        return

    negative_goldStandardLink = remove_outliers(negative_goldStandardLink,lower_percentile,upper_percentile)

    if negative_goldStandardLink is None or len(negative_goldStandardLink)<min_n_sample:
        if saveToDatabase: write_pLLRToDatabase(evidence_t,goldstandard_t,['Nan'],motivation='|NGS|<%i' %min_n_sample)
        print('Too few NGS!')
        return
    if len(negative_goldStandardLink)>max_n_sample:
        random.seed(int(species_golstandard))
        negative_goldStandardLink = random.sample(negative_goldStandardLink, max_n_sample)
        
    ## Save min,max values
    x_min = np.max([np.min(positive_goldStandardLink), np.min(negative_goldStandardLink)])
    x_max = np.min([np.max(positive_goldStandardLink), np.max(negative_goldStandardLink)])
    ## To minimize the number of peaks, we exclude min,max values from training
    negative_goldStandardLink = [x for x in negative_goldStandardLink if x_min < x < x_max]
    positive_goldStandardLink = [x for x in positive_goldStandardLink if x_min < x < x_max]

    ####################
    ## KDE training
    kde_neg = gaussian_kde(negative_goldStandardLink)
    kde_pos = gaussian_kde(positive_goldStandardLink)
    # Obtain the bandwidth estimated with Silverman's rule of thumb
    if bw_method == 'bw_fixed':
        pbw,nbw = 0.1,0.1
    elif bw_method == 'bw_max':
        bw = np.max([kde_neg.silverman_factor(),kde_pos.silverman_factor()])
        pbw,nbw = bw,bw
    elif bw_method == 'bw_min':
        bw = np.min([kde_neg.silverman_factor(),kde_pos.silverman_factor()])
        pbw,nbw = bw,bw
    elif bw_method == 'bw_dynamic':
        pbw = kde_pos.silverman_factor()
        nbw = kde_neg.silverman_factor()

    kde_pos = gaussian_kde(positive_goldStandardLink,bw_method=pbw)
    kde_neg = gaussian_kde(negative_goldStandardLink,bw_method=nbw)

    x_min = max(np.min(positive_goldStandardLink), np.min(negative_goldStandardLink))
    x_max = min(np.max(positive_goldStandardLink), np.max(negative_goldStandardLink))
    positive_raw_score = list(filter(lambda x: x_min < x < x_max, positive_goldStandardLink))
    negative_raw_score = list(filter(lambda x: x_min < x < x_max, negative_goldStandardLink))

    # Calculate the min/max likelihood value (excluding NaNs and infinite values)
    # Replace NaNs and infinite values with the min/max value
    x_common = np.linspace(x_min, x_max, num=1000) 
    pos_kde_interpolated = kde_pos(x_common)
    neg_kde_interpolated = kde_neg(x_common)

    min_value = np.nanmin(pos_kde_interpolated)
    max_value = np.nanmax(pos_kde_interpolated)
    pos_kde_interpolated = np.nan_to_num(pos_kde_interpolated, nan=min_value, posinf=max_value, neginf=min_value)
    
    min_value = np.nanmin(neg_kde_interpolated)
    max_value = np.nanmax(neg_kde_interpolated)
    neg_kde_interpolated = np.nan_to_num(neg_kde_interpolated, nan=min_value, posinf=max_value, neginf=min_value)

    try:
        ratio = np.log(pos_kde_interpolated / neg_kde_interpolated)
    except ValueError: 
        if saveToDatabase: write_pLLRToDatabase(evidence_t,goldstandard_t,['Nan'],motivation='Log-transform likelihoods fail')
        return
    else:
        LLR = np.column_stack((x_common, ratio))

        ### Optimize fit increasing polynomial degree
        r_squared = 0
        # Choose the desired degree for the polynomial model
        tmp_degree = poly_degree_min
        while r_squared<min_r2 and tmp_degree<poly_degree_max:
            # Fit a polynomial model to the LLR points
            coefficients = np.polyfit(LLR[:,0], LLR[:,1], tmp_degree)
            # Evaluate the fitted polynomial model at the x_plot points
            y_pred = np.polyval(coefficients, LLR[:,0])
            residuals = LLR[:,1] - y_pred
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((LLR[:,1] - np.mean(LLR[:,1]))**2)
            r_squared = 1 - (ss_res / ss_tot)
            tmp_degree +=1
        
        # print('max R2: %.2f' %r_squared)
        flag = True
        if r_squared<min_r2 or np.any(np.isnan(coefficients)): 
            if saveToDatabase: write_pLLRToDatabase(evidence_t,goldstandard_t,['Nan'],motivation='No good Polynomial fit w/ d<%i at R2≥%.2f' %(poly_degree_max,min_r2))
            flag = False
            print('No good fitting!')
        
        ## When training was successful
        if saveToDatabase and flag: write_pLLRToDatabase(evidence_t,goldstandard_t,coefficients,[x_min,x_max,pbw,nbw])
        
        if visualize_LLR:  
            ## Plot in data/tmp/LLR RawScore+KDE and LLR+polyFitting 
            positive_raw_score = list(filter(lambda x: x_min < x < x_max, positive_goldStandardLink))
            negative_raw_score = list(filter(lambda x: x_min < x < x_max, negative_goldStandardLink))

            ## VISUALIZATION
            fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))  # Adjust figsize as needed
            fig.suptitle('%s_%s - %s_%s (%i+,%i-)' %(evidence_t.type,species_evidence,goldstandard_t.type,species_golstandard,len(positive_goldStandardLink),len(negative_goldStandardLink)), fontsize=14)

            axes[0].hist(positive_raw_score, bins=50, facecolor='#003049', label='Positive Scores', alpha=0.5, density=True)
            axes[0].hist(negative_raw_score, bins=50, facecolor='#c1121f', label='Negative Scores', alpha=0.5, density=True)

            # Plotting the KDE curve
            # Scatter plot of pos_kde and neg_kde after interpolation of common x-range
            axes[0].scatter(x_common, pos_kde_interpolated, color='#003049', label='Positive KDE (%.2f)' %pbw, s=10)
            axes[0].scatter(x_common, neg_kde_interpolated, color='#c1121f', label='Negative KDE (%.2f)' %nbw, s=10)
            axes[0].set_xlabel('Raw Score')
            axes[0].set_ylabel('Likelihood')
            axes[0].set_xlim(x_min, x_max)
            axes[0].legend()

            # Create scatter plot of the data points
            axes[1].scatter(LLR[:,0], LLR[:,1], color='#a7c957', label='Data Points',s=10)
            # Plot the fitted model
            x_range = np.linspace(x_min, x_max, 1000)
            y_pred = np.polyval(coefficients, x_range)

            axes[1].plot(x_range, y_pred, color='#2c6e49', label='%s R2=%.2f'%('Fitted Model',r_squared))
            # axes[1].axvline(xper x=x_max, color='grey', linestyle='dashed', label='x = max')
            axes[1].set_xlabel('Raw Score')
            axes[1].set_ylabel('LLR')
            axes[1].set_xlim(x_min, x_max)
            axes[1].legend()
            axes[1].set_title(formatPolynomial(coefficients), fontsize=10)

            # Adjust spacing between subplots
            plt.tight_layout()
            evidence_extra = ''
            if evidence_t.type=='MEX':
                evidence_extra = '(%s)' %evidence_t.version
                if score_cutoff is not None:
                    evidence_extra = '(%s_%s)' %(evidence_t.version,score_cutoff)
            elif evidence_t.type=='PEX':
                evidence_extra = '(%s)' %evidence_t.scoringMethod
            if test: 
                print('###This is test!###')
                plt.savefig('data/tmp/LLR/%s/TEST_%s%s%s-%s(%s)_MinMax(%i,%i)_%s.png' %(species_golstandard,species_evidence,evidence_t.type,evidence_extra,goldstandard_t.type,species_golstandard,min_n_sample,max_n_sample,bw_method))
                plt.close()
            elif saveToDatabase: 
                plt.savefig('data/tmp/LLR/%s/%s%s%s-%s(%s)_MinMax(%i,%i)_%s.png' %(species_golstandard,species_evidence,evidence_t.type,evidence_extra,goldstandard_t.type,species_golstandard,min_n_sample,max_n_sample,bw_method))
                plt.close()
            else: return coefficients
            