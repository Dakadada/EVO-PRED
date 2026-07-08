import os
import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.toolbox import mutate_residue
import mdtraj as md

# Initialize PyRosetta
pyrosetta.init("-mute all")
scorefxn = get_fa_scorefxn()

# Define selected mutations: (Residue_Index, Target_AA, Name)
selected_mutations =

fixed_pdb = "6F86_fixed.pdb"

for res_num, target_aa, mut_name in selected_mutations:
    print(f"Generating mutant structure: {mut_name}...")
    
    # Load a fresh copy of the cleaned wild-type pose
    pose = pose_from_pdb(fixed_pdb)
    
    # Perform in silico mutagenesis and pack local sidechains within 5.0 Angstroms
    mutate_residue(pose, res_num, target_aa, pack_radius=5.0, pack_scorefxn=scorefxn)
    
    # Dump temporary Rosetta structure
    temp_pdb = f"{mut_name}_temp.pdb"
    pose.dump_pdb(temp_pdb)
    
    # Clean up structure: Use MDTraj to strip ALL hydrogens
    # This bypasses the pyparsing-dependent string select using a clean python list comprehension
    traj = md.load(temp_pdb)
    heavy_atoms = [atom.index for atom in traj.topology.atoms if atom.element.symbol!= 'H']
    heavy_traj = traj.atom_slice(heavy_atoms)
    
    output_pdb = f"{mut_name}_heavy.pdb"
    heavy_traj.save(output_pdb)
    
    # Remove temporary file
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)
        
    print(f"Successfully generated and sanitized: {output_pdb}")

print("All mutated structural templates generated.")