### Calculating 3D Tanimoto scores via FastROCS

#### information
We use GPUs on our computing cluster to run FastROCS and find the greatest shape + color overlap scores for all compounds in the reference set to all compounds in the test set. 

In this directory, we provide the following files to reproduce our workflow:

#### scripts
- `calc_3D_sim.py`: Python script containing a modified version of the [BestShapeOverlayMultiConfQuery.py](https://docs.eyesopen.com/toolkits/python/fastrocstk/examples/example_bestshapeoverlaymulticonfquery.html#section-example-fastrocs-bestshapeoverlaymulticonfquery) function from OpenEye documentation
- `clean_3D_sim_matrix.py`: Python script to select only the best stereoisomer for each compound as the representative structure in the final all-by-all similarity matrix
- `run_3D.sh`: job submission script
- `tools.py`: Python file containing useful helper functions (same as the one in `2_generating_conformers`)

#### files
Input:
- `bb1_info.oeb`, `bb2_info.oeb`, `bb3_info.oeb`: OpenEye binary file with coordinates of generated conformers for each compound
- `bb1_info.log`, `bb2_info.log`, `bb3_info.log`: file of warning messages later used to parse compounds with unspecified stereochemistry
- `bb1_info.pkl`, `bb2_info.pkl`, `bb3_info.pkl`: dictionary containing the indices of compounds and their specified enumerated stereoisomers

Output:
- `bb1_info.npy`: all-by-all similarity matrix of Tanimoto combo scores for building blocks in position 1
- `bb2_info.npy`: all-by-all similarity matrix of Tanimoto combo scores for building blocks in position 2
- `bb3_info.npy`: all-by-all similarity matrix of Tanimoto combo scores for building blocks in position 3


- `bb1_info_row.csv`, `bb1_info_col.csv`: SMILES of compound at each index of `bb1_info.npy` in rows and columns, respectively
- `bb2_info_row.csv`, `bb2_info_col.csv`: SMILES of compound at each index of `bb2_info.npy` in rows and columns, respectively
- `bb3_info_row.csv`, `bb3_info_col.csv`: SMILES of compound at each index of `bb3_info.npy` in rows and columns, respectively
(for a symmetric matrix, these two files will be identical)
