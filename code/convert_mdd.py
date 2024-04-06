# :author: Alessio Giacomel
# :date: 2024-03-04
# :version: 1.0
# :disease: MDD
# This script is used to convert the MDD GWAS data to work with the TWAS Pipeline. The original data contains LogOR and standard error of the LogOR while the pipeline necessitates of the OR and the standard error of the OR. The script will convert the data and save it in the same folder as the original data.

import pandas as pd
import numpy as np
from pathlib import Path

# File path
file_path = Path('/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/GWAS/MDD/PGC_UKB_depression_genome-wide.txt')
assert file_path.exists(), "File not found!"

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(file_path, delimiter=' ')
print(df.columns)

# Calculate OR from logOR
df['OR'] = np.exp(df['LogOR'])

# Calculate SE of OR from stdErrLogOR
df['SE'] = df['OR'] * df['StdErrLogOR']


# Save the updated DataFrame back to the same file
df.to_csv('/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/GWAS/MDD/PGC_UKB_depression_genome-wide_updated.txt', sep=' ', index=False)

print("Data updated and saved successfully!") 