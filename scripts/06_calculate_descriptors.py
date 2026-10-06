import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
import warnings
warnings.filterwarnings("ignore")

def calculate_lipinski(smiles):
    """
    Menghitung 4 deskriptor aturan Lipinski: 
    Berat Molekul (MW), LogP, Donor Ikatan H (HBD), dan Akseptor Ikatan H (HBA).
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            return mw, logp, hbd, hba
        else:
            return None, None, None, None
    except:
        return None, None, None, None

def process_data():
    print("Memulai kurasi data dan perhitungan deskriptor...")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    raw_data_path = os.path.join(project_root, 'data', '01_raw', 'dpp4_activities_raw.csv')
    processed_data_path = os.path.join(project_root, 'data', '02_processed', 'dpp4_activities_processed.csv')
    
    if not os.path.exists(raw_data_path):
        print(f"File data mentah tidak ditemukan di {raw_data_path}")
        return
        
    df = pd.read_csv(raw_data_path)
    print(f"Data awal dimuat: {len(df)} baris")
    
    # 1. Menghapus data yang tidak memiliki SMILES atau nilai standard_value yang valid
    df = df.dropna(subset=['canonical_smiles', 'standard_value'])
    print(f"Setelah menghapus data tanpa SMILES/nilai: {len(df)} baris")
    
    # 2. Hanya menggunakan pengujian IC50 sebagai patokan awal (opsional)
    df = df[df['standard_type'] == 'IC50']
    
    # Memastikan standard_value dalam format numerik positif
    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
    df = df[df['standard_value'] > 0]
    
    # 3. Menghitung pIC50 (-log10(IC50 dalam Molar)) untuk kemudahan pemodelan
    # Asumsi unit IC50 dari ChEMBL adalah nM (nanoMolar), sehingga dikali 10^-9
    def convert_to_pic50(value):
        try:
            molar = value * 1e-9 # Konversi nM ke Molar
            return -np.log10(molar)
        except:
            return None
            
    df['pIC50'] = df['standard_value'].apply(convert_to_pic50)
    df = df.dropna(subset=['pIC50'])
    print(f"Setelah menghitung pIC50 dan membersihkan tipe pengujian: {len(df)} baris unik")
    
    # 4. Menghapus duplikasi molekul (menyimpan nilai terbaik/rata-rata bisa dilakukan, di sini menyimpan yang pertama)
    df = df.drop_duplicates(subset=['canonical_smiles'])
    print(f"Senyawa unik setelah hapus duplikat: {len(df)}")
    
    # 5. Menghitung deskriptor molekuler (Hukum Lima Lipinski)
    print("Sedang menghitung deskriptor molekuler dengan RDKit...")
    descriptors = df['canonical_smiles'].apply(calculate_lipinski)
    
    # Menambahkan hasil deskriptor sebagai kolom baru
    df['MW'] = [d[0] for d in descriptors]
    df['LogP'] = [d[1] for d in descriptors]
    df['HBD'] = [d[2] for d in descriptors]
    df['HBA'] = [d[3] for d in descriptors]
    
    # Hapus senyawa yang gagal diurai oleh RDKit
    df = df.dropna(subset=['MW'])
    
    # 6. Menyimpan data hasil kurasi
    df.to_csv(processed_data_path, index=False)
    print(f"Data bersih dan ber-deskriptor berhasil disimpan di:\n{processed_data_path}")
    print(f"Total dataset akhir: {len(df)} senyawa")

if __name__ == "__main__":
    process_data()
