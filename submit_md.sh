#!/bin/bash
#SBATCH --job-name=evopred_baseline
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --output=baseline_%j.log

module load miniforge
conda activate evopred

python run_md_baseline.py 6F86.pdb