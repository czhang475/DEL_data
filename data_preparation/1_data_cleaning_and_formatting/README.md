### Data cleaning and formatting


#### information
Steps for dataset cleaning will vary from dataset to dataset. In the example script, we demonstrate the specific steps we undertook to go from raw experimental data to the information we needed to begin our building block level analysis.

In this directory, we have:

#### scripts
- `1_data_cleaning_and_formatting.ipynb`: a Jupyter notebook demonstrating procedures to remove undesired features within a dataset and deprotect chemical structures via SMIRKS

#### files
input:
- `del_hits.csv`: initial dataset of hits compounds sent from Anagenex
- `del_inactives.csv`: initial dataset of inactive compounds sent from Anagenex
(these two files are on Google Drive)

output:
- `bb1_list.csv`: list of unique building blocks at position 1 of the library
- `bb2_list.csv`: list of unique building blocks at position 2 of the library
- `bb3_list.csv`: list of unique building blocks at position 3 of the library
- `total_compounds.csv`: list of all cleaned active and inactive compounds from the datasets
