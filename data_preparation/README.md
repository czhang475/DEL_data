## Workflow for DEL project (last updated 6/7/22)

### Data Generation

#### `1_data_cleaning_and_formatting`

- Takes unformatted input data and cleans results to prep for analysis
- Outputs a .csv file of combined hits and inactive compounds, all of which pass different filters

#### `2_generating_conformers`

- Generates conformers for compounds provided in the .csv files as SMILES
- Outputs OpenEye binary (.oeb) files containing up to 200 conformers for each compound. In cases a compound has stereocenters but no specified stereochemistry, all stereoisomers are enumerated and conformers are generated for these new compounds as well

#### `3_calculating_3D_tanimoto_scores`

- Calculates 3D Tanimoto similarity scores for all compounds in a given reference and test set. The reference set corresponds to rows of the similarity matrix and the test set corresponds to columns.
- Outputs a NumPy array of the Tanimoto combo score between all compounds in the reference set and all compounds in the test set. The .csv files indicdate which compounds are on each index of the similarity matrix.
 
#### `4_calculating_pactive_values`

- Outputs the enrichment value (P(active)) for the building blocks at each position

#### output: 
- Contains the .csv and .npy files needed to begin building block analysis

#### rdkit_env_spec_file.txt
- File with all package versions and dependencies to reproduce the local environment for analysis
- Build the environment using the command `conda create --name rdkit_env --file rdkit_env_spec_file.txt`

#### oepython_env_spec_file.txt
- File with all package versions and dependencies to reproduce the cluster environment for analysis
- Build the environment using the command `conda create --name oepython --file oepython_env_spec_file.txt`

