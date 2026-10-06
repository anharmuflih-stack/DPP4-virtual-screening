import pandas as pd
import requests
import os
import time

print("Membaca hasil screening...")
res_df = pd.read_csv('results/02_docking/multi_target_screening_results.csv')

# Ambil chembl_id unik
chembl_ids = res_df['ChEMBL_ID'].unique()
names_dict = {}

print("Mengambil nama senyawa dari ChEMBL API...")
for cid in chembl_ids:
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            pref_name = data.get('pref_name')
            if pref_name:
                names_dict[cid] = pref_name
            else:
                names_dict[cid] = "Unnamed/No Pref Name"
        else:
            names_dict[cid] = "API Error"
    except Exception as e:
        names_dict[cid] = "Request Error"
    time.sleep(0.1) # Be nice to the API

res_df['Nama_Senyawa'] = res_df['ChEMBL_ID'].map(names_dict)

# Reorder & rename columns
final_df = res_df[['ChEMBL_ID', 'Nama_Senyawa', 'Protein', 'Binding_Affinity_kcal_mol', 'RMSD_lb', 'RMSD_ub']].copy()
final_df.rename(columns={
    'Protein': 'Protein_Target',
    'Binding_Affinity_kcal_mol': 'Affinity_Score',
    'RMSD_lb': 'RMSD_Lower_Bound',
    'RMSD_ub': 'RMSD_Upper_Bound'
}, inplace=True)

# Sort by Affinity
final_df = final_df.sort_values(by='Affinity_Score')

out_path = 'results/02_docking/Final_Docking_Report.csv'
final_df.to_csv(out_path, index=False)

print(f"File laporan akhir berhasil disimpan di: {out_path}")
