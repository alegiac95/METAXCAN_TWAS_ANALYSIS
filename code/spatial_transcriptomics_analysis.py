import warnings
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from netneurotools import freesurfer, stats

import seaborn as sns
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None  # default='warn'


def generate_permutations(data, cort_index = 0, n_permutations=1_000, left_only = True):
    ''' Function to generate permutations for spin test
    :param np.ndarray data: imaging data to permute and spin.
    The order of the imaging is first cortical and then subcortical.
    :param int n_permutations: number of spins to generate.
    cort_index: last index of cortex elements'''
    # initialisations
    np.random.seed(2024)
    permuted = np.zeros((data.shape[0], n_permutations))
    cortical = data[:cort_index]
    # Cortical resample
    # Annotation file for the Desikan-Killiany atlas in fs5
    # Get the parcel centroids of the Desikan-Killiany atlas
    parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
                lhannot="atlas-DK_fsa5_lh_aparc.annot", # change to file path if necessary
                rhannot="atlas-DK_fsa5_rh_aparc.annot", # change to file path if necessary
                version="fsaverage5",
                surf="sphere",
                method="surface")
    if left_only:
        # Mask the results to have only the left hemisphere
        left_hemi_mask = parcel_hemi == 0
    parcel_centroids, parcel_hemi = (
                  parcel_centroids[left_hemi_mask],
                  parcel_hemi[left_hemi_mask],
              )
            # Get the spin samples
    spins = stats.gen_spinsamples(
                parcel_centroids, parcel_hemi,
                n_rotate=n_permutations,
                method="vasa",
                seed=2024
            )
    cort_permuted = np.array(cortical[spins]).reshape(cort_index, n_permutations)
    permuted[:cort_index, :] = cort_permuted
    return permuted

def compute_correlation_and_pvalue(data, icv_column, mean_weighted_corr):
    """
    :data: permuted gene data
    :icv_column: imaging data
    :mean_wighted_corr: mean correlation of the original data
    
    Returns:
    mean_weighted_corr: mean correlation of the permuted data
    p_value_M: p-value for the mean correlation
    tmp: array of the correlations with null permutations (useful for plotting later)
    """
    correlation_weighted_average = []

    for i in range(1, data.shape[1]):
        curr_corr = spearmanr(icv_column, data[:, i])
        correlation_weighted_average.append(curr_corr.correlation)

    # p-value for weighted mean
    tmp = np.asarray(correlation_weighted_average)
    idx = np.where((tmp > mean_weighted_corr.correlation) if mean_weighted_corr.correlation > 0 else (tmp < mean_weighted_corr.correlation))[0]
    p_value_M = len(idx) / 1000
   
    ## Generate the plot with correlation of the null distribution and the mean correlation
    ax, fig = plt.subplots(figsize=(7, 5))
    sns.histplot(correlation_weighted_average, kde=True, color='blue', bins=30, ax=ax)
    ax.axvline(x=mean_weighted_corr.correlation, color='red', linestyle='dashed', linewidth=2)

    return mean_weighted_corr.correlation, p_value_M, tmp



if __name__ == '__main__':
    
    base_path = Path("/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/")
    
    if not data_path.exists():
        raise NotADirectoryError("Data path does not exist")
    
    diseases = ["ADHD", "AN", "ASD", "BD", "MDD", "OCD", "SCZ"]
    thresholds = [10, 5, 1, 0.5, 0.1]
    
    for disesease in diseases:
        disease_path = base_path / disease
        if not disease_path.exists():
            raise NotADirectoryError(f"{disease} path does not exist")
        
        for threshold in thresholds:
            data_file = f"{disease}_TPRS_thr_{threshold}.tsv"
            
            data_path = disease_path / data_file
            if not data_path.exists():
                raise FileNotFoundError(f"{data_file} does not exist")
            
            print(f"Analysis {data_file}")
            
            
        