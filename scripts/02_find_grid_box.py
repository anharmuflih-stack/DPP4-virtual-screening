import os
import numpy as np

def calculate_grid_box():
    print("Menghitung pusat Grid Box berdasarkan Ligan Bawaan (Native Ligand)...")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    pdb_dir = os.path.join(project_root, 'data', '03_pdb')
    
    pdb_files = ['2ONC.pdb', '3G0B.pdb', '5T4B.pdb']
    
    for pdb_file in pdb_files:
        path = os.path.join(pdb_dir, pdb_file)
        if not os.path.exists(path):
            print(f"File {pdb_file} tidak ditemukan.")
            continue
            
        print(f"\n--- Menganalisis {pdb_file} ---")
        
        # Mengumpulkan koordinat HETATM (Ligan)
        # Akan mengabaikan air (HOH) dan ion sederhana (seperti CL, NA, dll.)
        ligand_atoms = []
        ligand_names = set()
        
        with open(path, 'r') as f:
            for line in f:
                if line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    chain_id = line[21:22].strip()
                    
                    # Hanya fokus pada rantai 'A' agar tidak menggabungkan ligan dari rantai lain (jika proteinnya dimer/tetramer)
                    if chain_id == 'A':
                        # Abaikan pelarut dan molekul kristalisasi umum
                        if res_name not in ['HOH', 'GOL', 'EDO', 'SO4', 'PO4', 'CL', 'NA', 'NAG', 'DMS', 'ACT']:
                            try:
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                ligand_atoms.append([x, y, z])
                                ligand_names.add(res_name)
                            except ValueError:
                                pass
                            
        if len(ligand_atoms) == 0:
            print("Tidak ditemukan ligan bawaan yang valid!")
            continue
            
        # Menghitung Center of Mass (pusat ligan)
        coords = np.array(ligand_atoms)
        center_x, center_y, center_z = np.mean(coords, axis=0)
        
        # Menghitung seberapa luas ligannya untuk panduan Grid Size
        min_coords = np.min(coords, axis=0)
        max_coords = np.max(coords, axis=0)
        dimensions = max_coords - min_coords
        
        # Biasanya Grid Size dilebihkan beberapa Angstrom dari ukuran ligan asli (buffer)
        buffer_size = 10.0 # 10 Angstrom extra space
        size_x = dimensions[0] + buffer_size
        size_y = dimensions[1] + buffer_size
        size_z = dimensions[2] + buffer_size
        
        print(f"Ligan yang ditemukan: {', '.join(ligand_names)}")
        print(">> Koordinat Pusat (Center):")
        print(f"   center_x = {center_x:.3f}")
        print(f"   center_y = {center_y:.3f}")
        print(f"   center_z = {center_z:.3f}")
        
        print(">> Rekomendasi Ukuran Kotak (Size):")
        print(f"   size_x = {size_x:.1f}")
        print(f"   size_y = {size_y:.1f}")
        print(f"   size_z = {size_z:.1f}")

if __name__ == "__main__":
    calculate_grid_box()
