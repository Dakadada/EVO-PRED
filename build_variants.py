import os
import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.toolbox import mutate_residue
import mdtraj as md

print("Initializing PyRosetta...")
pyrosetta.init("-mute all")

wt_pdb = "6F86_fixed.pdb"
if not os.path.exists(wt_pdb):
    if os.path.exists("6F86.pdb"):
        wt_pdb = "6F86.pdb"
    else:
        raise FileNotFoundError("Could not locate 6F86_fixed.pdb or 6F86.pdb.")

print(f"Loading wild-type template: {wt_pdb}")
pose = pose_from_pdb(wt_pdb)
scorefxn = get_fa_scorefxn()

mutations = [("S", "K89S"), ("N", "K89N"), ("D", "K89D"), ("Y", "K89Y"), ("E", "K89E")]

for mut_aa, name in mutations:
    print(f"\n--- Generating {name} ---")
    mut_pose = pose.clone()
    
    print(f"Applying mutation K89 -> {mut_aa} and repacking local side chains...")
    mutate_residue(mut_pose, 89, mut_aa, pack_radius=5.0, pack_scorefxn=scorefxn)
    
    temp_pdb = f"temp_{name}.pdb"
    mut_pose.dump_pdb(temp_pdb)
    
    print("Stripping non-compliant hydrogens for AMBER compatibility...")
    traj = md.load(temp_pdb)
    
    heavy_atoms = [atom.index for atom in traj.topology.atoms if atom.element.symbol != 'H']
    heavy_traj = traj.atom_slice(heavy_atoms)
    
    output_pdb = f"{name}_heavy.pdb"
    heavy_traj.save(output_pdb)
    print(f"Successfully saved: {output_pdb}")
    
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)

print("\nUnified structural generation complete. All heavy-atom mutants are ready.")