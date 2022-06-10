### Running conformer generation on computing cluster

#### information
This folder contains code to generate conformers for compounds. The exact details of this procedure will vary depending on the specifications of each cluster. The workflow shown here is meant to provide general guidelines.

In this directory, we have:

#### scripts
- `2_generating_conformers.ipynb`: a Jupyter notebook demonstrating the conformer generation process
- `gen_conf.py`: a Python script to generating conformers
- `run_gen_conf.sh`: a sample job submission script (this will be different depending on the cluster)
- `tools.py`: a Python file containing useful premade functions

#### files

input:
- `bb1_list.csv`: list of unique building blocks at position 1 of the library
- `bb2_list.csv`: list of unique building blocks at position 2 of the library
- `bb3_list.csv`: list of unique building blocks at position 3 of the library

output:
- `bb1_info.oeb`, `bb2_info.oeb`, `bb3_info.oeb`: OpenEye binary file with coordinates of generated conformers for each compound
- `bb1_info.log`, `bb2_info.log`, `bb3_info.log`: file of warning messages later used to parse compounds with unspecified stereochemistry
- `bb1_info.pkl`, `bb2_info.pkl`, `bb3_info.pkl`: dictionary containing the indices of compounds and their specified enumerated stereoisomers

example:
- `sample_compounds.oeb`, `sample_compounds.log`, `sample_compounds.pkl`: sample compounds taken as the first 5 compounds from `bb1_info` to demonstrate the conformer generation process in `2_generating_conformers.ipynb`

