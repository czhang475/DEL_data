#! /usr/bin/bash
#SBATCH --job-name=generate_conformers
#SBATCH --partition=titanx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2gb
#SBATCH --gres=gpu:titan:1

#Load in appropriate working environment
source activate oepython

# Initialize variables
python gen_conf.py --infile "../files/sample_compounds.csv"



