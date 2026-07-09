#!/bin/bash
#SBATCH --job-name=kivi_mutants
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-5
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --mem=8G

module load miniforge
source activate evopred
module load cuda/11.8

MUTANTS=("K89S" "K89N" "K89D" "K89Y" "K89E")
MUTANT=${MUTANTS[$((SLURM_ARRAY_TASK_ID-1))]}

echo "Starting pipeline for mutant: ${MUTANT} on host $(hostname) with GPU allocation."

cat << EOF > tleap_${MUTANT}.in
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ligand_gaff.frcmod
LIG = loadmol2 ligand_gaff.mol2
PROT = loadpdb ${MUTANT}_heavy.pdb
COMPLEX = combine {PROT LIG}
solvatebox COMPLEX TIP3PBOX 10.0
addions COMPLEX Na+ 0
addions COMPLEX Cl- 0
saveamberparm COMPLEX ${MUTANT}_complex.prmtop ${MUTANT}_complex.inpcrd
quit
EOF

tleap -f tleap_${MUTANT}.in
rm tleap_${MUTANT}.in

python run_md_mutant.py ${MUTANT}_complex.prmtop ${MUTANT}_complex.inpcrd ${MUTANT}_trajectory.xtc