import pandas as pd
import os

def filter_lipinski():
    print("Memulai penyaringan berdasarkan Aturan Lima Lipinski (Lipinski's Rule of Five)...")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    processed_data_path = os.path.join(project_root, 'data', '02_processed', 'dpp4_activities_processed.csv')
    output_data_path = os.path.join(project_root, 'data', '02_processed', 'dpp4_lipinski_passed.csv')
    
    if not os.path.exists(processed_data_path):
        print(f"File data tidak ditemukan di {processed_data_path}")
        return
        
    df = pd.read_csv(processed_data_path)
    print(f"Total senyawa sebelum penyaringan: {len(df)}")
    
    # Aturan Lipinski:
    # 1. Berat Molekul (MW) <= 500 Dalton
    # 2. Koefisien partisi (LogP) <= 5
    # 3. Donor Ikatan Hidrogen (HBD) <= 5
    # 4. Akseptor Ikatan Hidrogen (HBA) <= 10
    
    # Fungsi untuk mengevaluasi apakah sebuah molekul melanggar aturan
    def check_lipinski_violations(row):
        violations = 0
        if row['MW'] > 500:
            violations += 1
        if row['LogP'] > 5:
            violations += 1
        if row['HBD'] > 5:
            violations += 1
        if row['HBA'] > 10:
            violations += 1
        return violations
        
    # Menghitung jumlah pelanggaran untuk setiap senyawa
    df['Lipinski_Violations'] = df.apply(check_lipinski_violations, axis=1)
    
    # Aturan tradisional memperbolehkan maksimal 1 pelanggaran (violations <= 1)
    # Namun untuk hasil terbaik, kita cari yang sepenuhnya mematuhi aturan (violations == 0)
    df_passed = df[df['Lipinski_Violations'] == 0].copy()
    
    print(f"Senyawa yang memenuhi SEMUA Aturan Lipinski: {len(df_passed)}")
    
    # Kita juga bisa menyaring berdasarkan aktivitas. Misalnya kita hanya mengambil senyawa
    # yang tergolong 'Aktif' (IC50 < 1000 nM atau pIC50 > 6.0)
    df_active_passed = df_passed[df_passed['pIC50'] > 6.0].copy()
    
    # Mengurutkan dari yang paling aktif (pIC50 tertinggi) ke terendah
    df_active_passed = df_active_passed.sort_values(by='pIC50', ascending=False)
    
    print(f"Senyawa aktif (pIC50 > 6.0) yang mematuhi Lipinski: {len(df_active_passed)}")
    
    # Menyimpan hasil akhir
    df_active_passed.to_csv(output_data_path, index=False)
    
    print(f"\nData berhasil disimpan di: {output_data_path}")
    print("\n--- 5 Senyawa Kandidat Terbaik (Top 5 Hits) ---")
    for i, (_, row) in enumerate(df_active_passed.head(5).iterrows()):
        print(f"{i+1}. ChEMBL ID: {row['molecule_chembl_id']} | pIC50: {row['pIC50']:.2f} | MW: {row['MW']:.2f}")

if __name__ == "__main__":
    filter_lipinski()
