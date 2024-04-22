The `code` folder contains the following scripts:

- `*_TWAS.sh`: file to run both individual GTEx S-PrediXcan and S-MultiXcan models. Each file runs the analysis for a different disorder (e.g., `ASD_TWAS.sh` runs only for Autism Spectrum Disorder (ASD)). _Note_: **Before** running the `MDD_TWAS.sh` on the summary stats from ENIGMA, it is necessary to run the `convert_mdd.py` script to convert the `logOR` to `OR` in order to avoid errors in the TWAS analysis.
- `tissue_correlation.R`: generates overlap plots for each dirorder of the overlap of predicted GTEx genes from S-PrediXcan. 
- `gene_plots.R`: generates volcano plots of S-MultiXcan predicted genes for each disorder, highlighting the top 10% of genes, and additionally naming the top 3% of genes additionally.
- `PTRS_estimation.R`: Calcualtes the Poligenic Transcriptomic Risk Score (PTRS) of genes for each disorder at different thresholds.
- `spatial_trasncriptomic_analysis.py`: Runs spatial permutation correlation analysis between the calculated PTRS and summary stats from ENIGMA. Saves a set of `.csv` and `.tsv` files used later to plot the results more easily in `plot_correlations.R`.