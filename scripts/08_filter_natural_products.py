import pandas as pd
import requests
import time
import os

def filter_natural_products():
    print("Memulai proses identifikasi Senyawa Bahan Alam (Natural Products) dari ChEMBL...")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    input_data_path = os.path.join(project_root, 'data', '02_processed', 'dpp4_lipinski_passed.csv')
    output_data_path = os.path.join(project_root, 'data', '02_processed', 'dpp4_natural_products.csv')
    
    if not os.path.exists(input_data_path):
        print(f"File data tidak ditemukan di {input_data_path}")
        return
        
    df = pd.read_csv(input_data_path)
    print(f"Total kandidat senyawa awal: {len(df)}")
    
    # Mengambil daftar ID ChEMBL unik
    chembl_ids = df['molecule_chembl_id'].unique().tolist()
    
    # ChEMBL API membatasi jumlah ID yang bisa dicek sekaligus. Kita akan membaginya (chunking) per 100 ID
    chunk_size = 100
    natural_products = set()
    
    base_url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
    
    print(f"Akan memeriksa {len(chembl_ids)} senyawa. Mohon tunggu, proses ini akan memakan waktu beberapa menit...")
    
    for i in range(0, len(chembl_ids), chunk_size):
        chunk = chembl_ids[i:i+chunk_size]
        
        # Format parameter pencarian: molecule_chembl_id__in=ID1,ID2,ID3...
        ids_str = ",".join(chunk)
        params = {
            'molecule_chembl_id__in': ids_str,
            'limit': chunk_size
        }
        
        try:
            response = requests.get(base_url, params=params)
            if response.status_code == 200:
                data = response.json()
                molecules = data.get('molecules', [])
                
                for mol in molecules:
                    # Mengecek tanda 'natural_product'
                    if mol.get('natural_product') == 1 or mol.get('natural_product') == True:
                        natural_products.add(mol.get('molecule_chembl_id'))
                        
            time.sleep(0.5) # Menghindari batas rate-limit API
            
            # Menampilkan progress setiap 500 senyawa
            if (i + chunk_size) % 500 == 0:
                print(f"Progress: {i + chunk_size} senyawa diperiksa...")
                
        except Exception as e:
            print(f"Gagal memeriksa chunk pada indeks {i}: {e}")
            continue
            
    print(f"Selesai! Ditemukan {len(natural_products)} senyawa bahan alam dari ChEMBL.")
    
    # Filter dataset utama hanya untuk yang ada di set natural_products
    df_np = df[df['molecule_chembl_id'].isin(natural_products)].copy()
    
    # Urutkan berdasarkan pIC50 tertinggi
    df_np = df_np.sort_values(by='pIC50', ascending=False)
    
    # Menyimpan hasil
    df_np.to_csv(output_data_path, index=False)
    print(f"Data Bahan Alam berhasil disimpan di: {output_data_path}")
    
    if len(df_np) > 0:
        print("\n--- Senyawa Bahan Alam Terbaik (Hits) ---")
        for i, (_, row) in enumerate(df_np.head(5).iterrows()):
            print(f"{i+1}. ChEMBL ID: {row['molecule_chembl_id']} | pIC50: {row['pIC50']:.2f} | MW: {row['MW']:.2f}")
    else:
        print("Sayangnya tidak ditemukan senyawa dengan anotasi 'Natural Product' yang mematuhi filter sebelumnya.")
        print("Catatan: Banyak senyawa bahan alam di ChEMBL yang tidak diberi anotasi secara eksplisit. Anda mungkin butuh dataset eksternal seperti COCONUT database.")

if __name__ == "__main__":
    filter_natural_products()
