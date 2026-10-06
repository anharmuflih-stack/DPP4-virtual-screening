import os
import requests

def download_dpp4_structure():
    print("Memulai pengunduhan struktur 3D protein DPP-4 dari PDB...")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    pdb_dir = os.path.join(project_root, 'data', '03_pdb')
    
    # Membuat folder pdb jika belum ada
    if not os.path.exists(pdb_dir):
        os.makedirs(pdb_dir)
        print(f"Folder {pdb_dir} berhasil dibuat.")
        
    # Daftar PDB ID untuk target DPP-4 (2ONC dihapus karena tidak valid)
    pdb_ids = ["3G0B", "5T4B"]
    
    for pdb_id in pdb_ids:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        output_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
        
        print(f"Mengunduh PDB ID: {pdb_id}...")
        
        try:
            response = requests.get(pdb_url)
            response.raise_for_status()
            
            with open(output_path, 'wb') as f:
                f.write(response.content)
                
            print(f"Berhasil mengunduh! File disimpan di: {output_path}")
        except Exception as e:
            print(f"Gagal mengunduh file PDB {pdb_id}: {e}")
            
    print("\nTahap selanjutnya untuk AutoDock Vina:")
    print("1. Buka file PDB ini (misalnya menggunakan PyMOL, AutoDockTools, atau Discovery Studio).")
    print("2. Hapus molekul air (HOH) dan ligan bawaan (inhibitor).")
    print("3. Tambahkan atom hidrogen polar dan muatan (Kollman charges).")
    print("4. Simpan makromolekul tersebut sebagai file berformat .pdbqt")

if __name__ == "__main__":
    download_dpp4_structure()
