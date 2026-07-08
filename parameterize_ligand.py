import subprocess
import mdtraj as md

raw_pdb = "6F86.pdb"
ligand_out = "ligand.pdb"

traj = md.load(raw_pdb)
topology = traj.topology

# Bypass the pyparsing string evaluation bug using native list comprehensions
ligand_atoms = []
for residue in topology.residues:
    if not residue.is_water and not residue.is_protein and not residue.is_nucleic:
        if residue.n_atoms > 10:  # Excludes monoatomic ions and small solvent fragments
            ligand_atoms.extend([atom.index for atom in residue.atoms])

if len(ligand_atoms) == 0:
    print("Error: Ligand not found in topology.")
    exit(1)

ligand_traj = traj.atom_slice(ligand_atoms)
ligand_traj.save(ligand_out)
print(f"Ligand isolated to {ligand_out}")

cmd = [
    "antechamber", "-i", ligand_out, "-fi", "pdb", 
    "-o", "ligand_gaff.mol2", "-fo", "mol2", 
    "-c", "bcc", "-s", "2", "-nc", "0", "-rn", "LIG"
]

print("Running Antechamber (this may take several minutes)...")
subprocess.run(cmd, check=True)

cmd_parmchk = ["parmchk2", "-i", "ligand_gaff.mol2", "-f", "mol2", "-o", "ligand_gaff.frcmod"]
subprocess.run(cmd_parmchk, check=True)

print("Ligand parameterization complete. Generated ligand_gaff.mol2 and ligand_gaff.frcmod")