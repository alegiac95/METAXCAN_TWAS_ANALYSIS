import warnings
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr, zscore
from pathlib import Path

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from netneurotools import freesurfer, stats


pd.options.mode.chained_assignment = None  # default='warn'


def generate_permutations(data, cort_index = 0, n_permutations=1_000, left_only = True):
    ''' Function to generate permutations for spin test
    :param np.ndarray data: imaging data to permute and spin.
    The order of the imaging is first cortical and then subcortical.
    :param int n_permutations: number of spins to generate.
    cort_index: last index of cortex elements'''
    # initialisations
    np.random.seed(2024)
    data_path = "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/DK-files/"
    permuted = np.zeros((data.shape[0], n_permutations))
    cortical = data[:cort_index]
    # Cortical resample
    # Annotation file for the Desikan-Killiany atlas in fs5
    # Get the parcel centroids of the Desikan-Killiany atlas
    parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
                lhannot=f"{data_path}atlas-DK_fsa5_lh_aparc.annot", # change to file path if necessary
                rhannot=f"{data_path}atlas-DK_fsa5_rh_aparc.annot", # change to file path if necessary
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
    p_value_M = len(idx) / data.shape[1]

    return mean_weighted_corr.correlation, p_value_M, tmp



if __name__ == '__main__':
    
    base_path = Path("/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/")
    
    # ENIGMA imaging data
    enigma_path = Path("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/ENIGMA/ENIGMA_structural_organised.xlsx")
    
    if not enigma_path.exists():
        raise NotADirectoryError
    
    cortical_data = pd.read_excel(enigma_path, sheet_name = "Thickness", index_col = 0)
    subcortical_data = pd.read_excel(enigma_path, sheet_name = "Subcortical", index_col = 0)


    
    if not base_path.exists():
        raise NotADirectoryError("Data path does not exist")
    
    diseases = ["ADHD", "AN", "ASD", "BD", "MDD", "OCD", "SCZ"]
    thresholds = [10, 5, 1]
    
    for disease in diseases:
        disease_path = base_path / disease
        
        cortical_data_to_save = pd.DataFrame(columns= ["threshold", "corr", "p"])
        subcortical_data_to_save = pd.DataFrame(columns= ["threshold", "corr", "p"])
        
        if not disease_path.exists():
            raise NotADirectoryError(f"{disease} path does not exist")
        
        cort_data = cortical_data[disease].dropna()
        subc_data = subcortical_data[disease].dropna()
        
        for threshold in thresholds:
            data_file = f"{disease}_literature_TPRS_thr_{threshold}.tsv"
            
            data_path = disease_path / data_file
            
            if not data_path.exists():
                raise FileNotFoundError(f"{data_file} does not exist")
            
            print("")
            print(f"Analysing {data_file}....")
            
            gene_data = pd.read_csv(data_path, sep='\t')
            # Cortical data
            cort_gene_data = gene_data[(gene_data['hemisphere'] == 'L') & (gene_data['structure'] == 'cortex')]
            permuted_genes = generate_permutations(cort_gene_data.weighted_avg.to_numpy(), cort_index = 34, n_permutations = 1000, left_only = True)
            orig_corr = spearmanr(cort_data.to_numpy(),cort_gene_data.weighted_avg.to_numpy())
            corr, p, tmp = compute_correlation_and_pvalue(permuted_genes, cort_data, orig_corr)
            
            # data to save to file
            row_data = {'threshold': threshold,'corr': corr,'p': p}
            cortical_data_to_save = cortical_data_to_save.append(row_data, ignore_index=True)
            # save correlations for plot 
            np.savetxt(disease_path / f"{disease}_literature_cortical_corr_{threshold}.tsv", tmp, delimiter='\t')
            
            # Subcortical data
            subc_gene_data = gene_data[(gene_data['hemisphere'] == 'L') & (gene_data['structure'] == 'subcortex/brainstem')]
            choice_len = len(subc_gene_data.weighted_avg.to_numpy())
            sub_permuted = subc_gene_data.weighted_avg.to_numpy()[np.random.choice(choice_len, size=(choice_len, 1000))]
            
            orig_corr = spearmanr(subc_data.to_numpy(), subc_gene_data.weighted_avg.to_numpy())
            corr, p, tmp = compute_correlation_and_pvalue(sub_permuted, subc_data.to_numpy(), orig_corr)
            
            # data to save to file
            row_data = {'threshold': threshold,'corr': corr,'p': p}
            subcortical_data_to_save = subcortical_data_to_save.append(row_data, ignore_index=True)
            # save correlations for plot
            np.savetxt(disease_path / f"{disease}_literature_subcortical_corr_{threshold}.tsv", tmp, delimiter='\t')
            
        print(cortical_data_to_save)
        # Save the data to file
        cortical_data_to_save.to_csv(disease_path / f"{disease}_literature_cortical_correlations_results.tsv", sep='\t', index=False)
        print(subcortical_data_to_save)
        subcortical_data_to_save.to_csv(disease_path / f"{disease}_literature_subcortical_correlations_results.tsv", sep='\t', index=False)