import os
import subprocess

def prepare_receptors():
    print("=== Persiapan Reseptor (PDB -> PDBQT) ===")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    pdb_dir = os.path.join(project_root, 'data', '03_pdb')
    receptor_out_dir = os.path.join(project_root, 'data', '05_receptors_pdbqt')
    
    if not os.path.exists(receptor_out_dir):
        os.makedirs(receptor_out_dir)
        
    pdb_files = ['3G0B.pdb', '5T4B.pdb']
    
    for pdb_name in pdb_files:
        pdb_path = os.path.join(pdb_dir, pdb_name)
        if not os.path.exists(pdb_path):
            continue
            
        print(f"Membersihkan protein {pdb_name} (Hanya mengambil Rantai A)...")
        clean_pdb_path = os.path.join(pdb_dir, f"{pdb_name.split('.')[0]}_clean.pdb")
        
        # Membersihkan PDB: Hanya menyimpan baris ATOM untuk rantai A
        with open(pdb_path, 'r') as fin, open(clean_pdb_path, 'w') as fout:
            for line in fin:
                if line.startswith("ATOM") and line[21:22] == 'A':
                    fout.write(line)
        
        # Konversi menggunakan Meeko mk_prepare_receptor.py via command line
        out_pdbqt = os.path.join(receptor_out_dir, f"{pdb_name.split('.')[0]}.pdbqt")
        print(f"Mengonversi {clean_pdb_path} ke PDBQT...")
        
        try:
            cmd = f"mk_prepare_receptor.py -i {clean_pdb_path} -o {out_pdbqt}"
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"Berhasil! Tersimpan di {out_pdbqt}")
        except subprocess.CalledProcessError as e:
            print(f"Gagal mengonversi protein {pdb_name}")

if __name__ == "__main__":
    prepare_receptors()
