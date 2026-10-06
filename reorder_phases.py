import os

base = '/Users/aanmuf/drug_discovery/dpp4_project/scripts/'

# Move everything to temp names first to avoid conflicts
os.rename(base+'01_fetch_chembl_data.py', base+'temp_fetch.py')
os.rename(base+'02_calculate_descriptors.py', base+'temp_calc.py')
os.rename(base+'03_filter_lipinski.py', base+'temp_lip.py')
os.rename(base+'04_filter_natural_products.py', base+'temp_nat.py')
os.rename(base+'05_download_protein.py', base+'01_download_protein.py')
os.rename(base+'06_find_grid_box.py', base+'02_find_grid_box.py')
os.rename(base+'08_validate_docking.py', base+'04_validate_docking.py')
os.rename(base+'09_run_virtual_screening.py', base+'10_run_virtual_screening.py')

os.rename(base+'temp_fetch.py', base+'05_fetch_chembl_data.py')
os.rename(base+'temp_calc.py', base+'06_calculate_descriptors.py')
os.rename(base+'temp_lip.py', base+'07_filter_lipinski.py')
os.rename(base+'temp_nat.py', base+'08_filter_natural_products.py')

# Read 07 and split it
with open(base+'07_prepare_docking_files.py', 'r') as f:
    content = f.read()

# We know the content. We can just create 03 and 09 directly
rec_content = """import os
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
"""

lig_content = """import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

def prepare_ligands():
    print("=== Persiapan Ligan (SMILES -> PDBQT) ===")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    input_csv = os.path.join(project_root, 'data', '02_processed', 'dpp4_natural_products.csv')
    ligand_out_dir = os.path.join(project_root, 'data', '04_ligands_pdbqt')
    
    if not os.path.exists(ligand_out_dir):
        os.makedirs(ligand_out_dir)
        
    df = pd.read_csv(input_csv)
    
    success_count = 0
    for idx, row in df.iterrows():
        chembl_id = row['molecule_chembl_id']
        smiles = row['canonical_smiles']
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: continue
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            preparator = MoleculePreparation()
            preparator.prepare(mol)
            pdbqt_string = preparator.write_pdbqt_string()
            out_path = os.path.join(ligand_out_dir, f"{chembl_id}.pdbqt")
            with open(out_path, 'w') as f:
                f.write(pdbqt_string)
            success_count += 1
        except Exception as e:
            pass
            
    print(f"Selesai! {success_count}/{len(df)} ligan berhasil dikonversi ke format PDBQT.")

if __name__ == "__main__":
    prepare_ligands()
"""

with open(base+'03_prepare_receptors.py', 'w') as f:
    f.write(rec_content)

with open(base+'09_prepare_ligands.py', 'w') as f:
    f.write(lig_content)

os.remove(base+'07_prepare_docking_files.py')

print("Selesai mengurutkan!")
