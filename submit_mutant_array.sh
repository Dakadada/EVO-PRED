
```bash
#submit_mutant_array.sh
#!/bin/bash
#SBATCH --job-name=kivi_mutants
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-5
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --mem=8G

# Load required cluster environments
module load miniforge
source activate evopred
module load cuda/11.8

# Map the array task index to our mutant target array
MUTANTS=("K89S" "K89D" "H85G" "F90E" "K89Y")
MUTANT=${MUTANTS}

echo "Starting pipeline for mutant: ${MUTANT} on host $(hostname) with GPU allocation."

# 1. Build Amber topology for this specific mutant in a local temporary file
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

# Clean up input build file
rm tleap_${MUTANT}.in

# 2. Run the production OpenMM molecular dynamics simulation on GPU
python run_md_mutant.py ${MUTANT}_complex.prmtop ${MUTANT}_complex.inpcrd ${MUTANT}_trajectory.xtc