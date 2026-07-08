import sys
import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.toolbox import mutate_residue
import pandas as pd

pyrosetta.init("-mute all")

pdb_file = "6F86_fixed.pdb"
pose = pose_from_pdb(pdb_file)
scorefxn = get_fa_scorefxn()

wt_score = scorefxn(pose)
print(f"Wild-type Ref2015 Score: {wt_score:.2f}")

# Define target residue range (e.g., sequence indices 80 to 90 for initial test sweep)
target_residues = range(80, 91)
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

results = []

for res_num in target_residues:
    wt_aa = pose.sequence()[res_num - 1]
    for mut_aa in amino_acids:
        if wt_aa == mut_aa:
            continue
            
        mut_pose = pose.clone()
        mutate_residue(mut_pose, res_num, mut_aa, pack_radius=5.0, pack_scorefxn=scorefxn)
        
        mut_score = scorefxn(mut_pose)
        ddg = mut_score - wt_score
        
        results.append({'Residue_Index': res_num, 'WT_AA': wt_aa, 'Mut_AA': mut_aa, 'Total_Score': mut_score, 'ddG': ddg})
        print(f"Mut {wt_aa}{res_num}{mut_aa} | ddG: {ddg:.2f}")

df = pd.DataFrame(results)
df.sort_values(by='ddG', inplace=True)
df.to_csv("6F86_ddg_predictions.csv", index=False)
print("Pipeline complete. Data exported to 6F86_ddg_predictions.csv")