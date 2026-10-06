import requests
import pandas as pd
import time
import os

def fetch_dpp4_activities():
    print("Memulai pengambilan data bioaktivitas untuk target DPP-4 (CHEMBL284)...")
    base_url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    
    # Target DPP-4
    target_chembl_id = "CHEMBL284"
    
    # Parameter filter: IC50, Kd, Ki
    params = {
        'target_chembl_id': target_chembl_id,
        'standard_type__in': 'IC50,Ki,Kd',
        'limit': 1000,
        'offset': 0
    }
    
    activities = []
    
    while True:
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            page_activities = data.get('activities', [])
            activities.extend(page_activities)
            
            print(f"Mengambil {len(page_activities)} data (Total: {len(activities)})")
            
            meta = data.get('page_meta', {})
            if meta.get('next'):
                params['offset'] += params['limit']
                time.sleep(1) # Jeda untuk menghormati API rate limit
            else:
                break
                
        except Exception as e:
            print(f"Terjadi kesalahan saat mengambil data: {e}")
            break
            
    print(f"Berhasil mengambil total {len(activities)} data aktivitas.")
    
    # Mengubah ke format DataFrame pandas untuk disimpan sebagai CSV
    df = pd.DataFrame(activities)
    
    # Filter kolom-kolom yang relevan agar file tidak terlalu besar
    columns_to_keep = [
        'activity_id', 'molecule_chembl_id', 'standard_type', 'standard_relation', 
        'standard_value', 'standard_units', 'pchembl_value', 'target_chembl_id', 
        'target_pref_name', 'canonical_smiles'
    ]
    
    # Hanya simpan kolom yang ada
    existing_columns = [col for col in columns_to_keep if col in df.columns]
    df_filtered = df[existing_columns]
    
    # Menentukan path output relatif terhadap lokasi skrip ini
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    output_dir = os.path.join(project_root, 'data', '01_raw')
    output_file = os.path.join(output_dir, 'dpp4_activities_raw.csv')
    
    # Menyimpan ke file CSV
    df_filtered.to_csv(output_file, index=False)
    
    print(f"Data mentah berhasil disimpan di {output_file}")
    print(f"Total baris data: {len(df_filtered)}")
    print(f"Total senyawa unik (ChEMBL ID): {df_filtered['molecule_chembl_id'].nunique()}")
    
if __name__ == "__main__":
    fetch_dpp4_activities()
