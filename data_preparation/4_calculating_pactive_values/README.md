### Calculating P(active) values

#### information
In this final step of data preparation, we calculate the enrichment for all building blocks that passed the conformer generation process. We demonstrate how we define and calculate the metric in the notebook provided in this directory.

In this directory we have:

#### scripts
- `4_calculating_pactive_values.ipynb`: a Jupyter notebook demonstrating how to calculate enrichment values for each building block

#### files
input:
- `total_compounds`: list of all cleaned active and inactive compounds from the datasets
- `bb1_info_row.csv`, `bb1_info_col.csv`: SMILES of compound at each index of `bb1_info.npy` in rows and columns, respectively
- `bb2_info_row.csv`, `bb2_info_col.csv`: SMILES of compound at each index of `bb2_info.npy` in rows and columns, respectively
- `bb3_info_row.csv`, `bb3_info_col.csv`: SMILES of compound at each index of `bb3_info.npy` in rows and columns, respectively


output:
- `bb1_pactive.csv`: all SMILES that passed conformer generation and their respective enrichment values in position 1
- `bb2_pactive.csv`: all SMILES that passed conformer generation and their respective enrichment values in position 2
- `bb3_pactive.csv`: all SMILES that passed conformer generation and their respective enrichment values in position 3
