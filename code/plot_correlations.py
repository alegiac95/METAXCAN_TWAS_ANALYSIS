## Generate the plot with correlation of the null distribution and the mean correlation
def plot_correlation_hist(correlation_weighted_average, mean_weighted_corr, save_path):
    """
    Plot histogram of the correlation of the null distribution and the mean correlation.
    
    INPUTS
    ------
    correlation_weighted_average: array of the correlations with null permutations
    mean_weighted_corr: mean correlation of the original data
    save_path: path to save the plot
    
    OUTPUTS
    -------
    Plot saved in the save_path
    """
    ax, fig = plt.subplots(figsize=(7, 5))
    sns.histplot(correlation_weighted_average, kde=True, color='blue', bins=30, ax=ax)
    ax.axvline(x=mean_weighted_corr.correlation, color='red', linestyle='dashed', linewidth=2)
    
    plt.savefig(save_path, dpi=300)