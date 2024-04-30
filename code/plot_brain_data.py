import matplotlib.pyplot as plt
# from enigmatoolbox.plotting import plot_cortical
# from enigmatoolbox.plotting import plot_subcortical
import numpy as np
from pathlib import Path
import neuromaps 
from surfplot import Plot
import nibabel as nib
from neuromaps.datasets import fetch_fslr
from brainspace.datasets import load_conte69




def plot_brains(brain_data: np.ndarray, save_path: str = None, cbar: str = "RdBu", *args, **kwargs) -> None:
    
    # save_path = Path(save_path) if isinstance(save_path, str) else save_path
    # if not save_path.exists():
    #     raise ValueError(f"The path {save_path} does not exist!")
    # if not save_path.is_dir() or save_path.is_file():
    #     raise NotADirectoryError
    # data_shape = brain_data.shape
    surfaces = fetch_fslr()
    #lh, rh = surfaces['inflated']
    lh, rh = load_conte69()
    #lh = "/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/DK-files/atlas-DK_fsa5_lh_aparc.annot"
    p = Plot(lh, size=(400, 200), zoom=1.8, views = ['lateral'], layout = "column")
    fig = p.build()
    fig.show()
    fig.savefig("/Users/alessiogiacomel/Desktop/t_brain.png")
    pass


if __name__ == '__main__':
    
    plot_brains(np.zeros(2))