import os
import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.toolbox import mutate_residue
import mdtraj as md

# 1. Initialize PyRosetta
print("Initializing PyRosetta...")
pyrosetta.init("-mute all")

wt_pdb = "6F86_fixed.pdb"
if not os.path.exists(wt_pdb):
    if os.path.exists("6F86.pdb"):
        wt_pdb = "6F86.pdb"
    else:
        raise FileNotFoundError("Could not locate 6F86_fixed.pdb or 6F86.pdb. Ensure you are in the EVO-PRED directory.")

print(f"Loading wild-type template: {wt_pdb}")
pose = pose_from_pdb(wt_pdb)
scorefxn = get_fa_scorefxn()

# Define the five target mutations: (Target_AA, Variant_Name)
mutations =

for mut_aa, name in mutations:
    print(f"\n--- Generating {name} ---")
    
    # Clone the clean wild-type pose
    mut_pose = pose.clone()
    
    # Mutate Lysine at position 89 to the target amino acid and repack within 5.0 Angstroms
    print(f"Applying mutation K89 -> {mut_aa} and repacking local side chains...")
    mutate_residue(mut_pose, 89, mut_aa, pack_radius=5.0, pack_scorefxn=scorefxn)
    
    # Dump a temporary structure
    temp_pdb = f"temp_{name}.pdb"
    mut_pose.dump_pdb(temp_pdb)
    
    # Use MDTraj to load the temp file and strip all hydrogens
    print("Stripping non-compliant hydrogens for AMBER compatibility...")
    traj = md.load(temp_pdb)
    
    # Bypass the pyparsing string selector bug using Pythonic list comprehensions
    heavy_atoms = [atom.index for atom in traj.topology.atoms if atom.element.symbol!= 'H']
    heavy_traj = traj.atom_slice(heavy_atoms)
    
    # Save the polished output
    output_pdb = f"6F86_{name}_heavy.pdb"
    heavy_traj.save(output_pdb)
    print(f"Successfully saved: {output_pdb}")
    
    # Remove the un-sanitized intermediate file
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)

print("\nUnified structural generation complete. All heavy-atom mutants are ready.")