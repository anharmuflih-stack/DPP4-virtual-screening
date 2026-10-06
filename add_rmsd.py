import pandas as pd
import os

df = pd.read_csv('results/02_docking/multi_target_screening_results.csv')
rmsd_lb = []
rmsd_ub = []

for idx, row in df.iterrows():
    protein = row['Protein']
    chembl = row['ChEMBL_ID']
    log_file = f"results/02_docking/{chembl}_{protein}.log"
    lb, ub = 0.0, 0.0
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            for line in f:
                if line.strip().startswith('1 '):
                    parts = line.split()
                    if len(parts) >= 4:
                        lb = float(parts[2])
                        ub = float(parts[3])
                    break
    rmsd_lb.append(lb)
    rmsd_ub.append(ub)

df['RMSD_lb'] = rmsd_lb
df['RMSD_ub'] = rmsd_ub

df.to_csv('results/02_docking/multi_target_screening_results.csv', index=False)

# Update individual files
proteins = ['3G0B', '5T4B']
for p in proteins:
    p_df = df[df['Protein'] == p][['ChEMBL_ID', 'Binding_Affinity_kcal_mol', 'RMSD_lb', 'RMSD_ub']].copy()
    p_df = p_df.sort_values(by='Binding_Affinity_kcal_mol')
    p_df.to_csv(f'results/02_docking/screening_results_{p}.csv', index=False)

print("RMSD berhasil ditambahkan!")
