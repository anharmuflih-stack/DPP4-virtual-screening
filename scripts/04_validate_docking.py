import os
import subprocess
import prody
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
import pandas as pd

def extract_native_ligand(pdb_file, resname, out_pdb):
    """Mengekstrak ligan bawaan dari file PDB kristal."""
    pdb = prody.parsePDB(pdb_file)
    # Pilih ligan dari Rantai A saja
    ligand = pdb.select(f'resname {resname} and chain A')
    
    if ligand is not None:
        # Jika ada beberapa residu (misal 2 molekul di rantai A), ambil residu pertama saja
        resnums = set(ligand.getResnums())
        first_resnum = min(resnums)
        single_ligand = ligand.select(f'resnum {first_resnum}')
        
        prody.writePDB(out_pdb, single_ligand)
        return True
    return False

import math

def calculate_rmsd_from_pdb(ref_pdbqt, docked_pdbqt):
    """Menghitung RMSD manual (Heavy-Atom) antara hasil docking dan referensi PDBQT."""
    try:
        def get_heavy_atoms(filepath):
            coords = []
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('MODEL 2'):
                        break
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # Cek elemen (kolom 77-78 atau nama atom kolom 13-16)
                        atom_name = line[12:16].strip()
                        # Lewati Hidrogen
                        if atom_name.startswith('H'):
                            continue
                        
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
            return coords

        ref_coords = get_heavy_atoms(ref_pdbqt)
        dock_coords = get_heavy_atoms(docked_pdbqt)
        
        if not ref_coords or not dock_coords or len(ref_coords) != len(dock_coords):
            print(f"Jumlah atom berat tidak cocok: ref={len(ref_coords)}, dock={len(dock_coords)}")
            return None
            
        sum_sq = 0.0
        for (rx, ry, rz), (dx, dy, dz) in zip(ref_coords, dock_coords):
            sum_sq += (rx - dx)**2 + (ry - dy)**2 + (rz - dz)**2
            
        rmsd = math.sqrt(sum_sq / len(ref_coords))
        return rmsd
    except Exception as e:
        print(f"Error kalkulasi RMSD: {e}")
        return None

def run_validation():
    print("=== Memulai Validasi Redocking (Menghitung RMSD) ===")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    pdb_dir = os.path.join(project_root, 'data', '03_pdb')
    receptor_dir = os.path.join(project_root, 'data', '05_receptors_pdbqt')
    ligand_dir = os.path.join(project_root, 'data', '04_ligands_pdbqt')
    results_dir = os.path.join(project_root, 'results', '01_validation')
    
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        
    targets = {
        '3G0B': {'ligand': 'T22', 'center': (42.05, 34.29, 14.62), 'size': (18, 18, 18), 'smiles': 'NC1CC(F)(F)CC(N2CC3(F)CCC3C2)C1'}, # Dummy
        '5T4B': {'ligand': '75N', 'center': (37.61, 49.92, 40.33), 'size': (18, 16, 23), 'smiles': 'O=C1CN(C2CCNCC2)C(=O)N1c3cccc(Cl)c3'} # Dummy
    }
    
    vina_exec = os.path.join(script_dir, 'vina')
    results = []
    
    for protein, params in targets.items():
        print(f"\\n>> Validasi Protein {protein} dengan Native Ligand {params['ligand']}")
        
        ref_pdb = os.path.join(results_dir, f"{protein}_native.pdb")
        extract_native_ligand(os.path.join(pdb_dir, f"{protein}.pdb"), params['ligand'], ref_pdb)
        
        # Konversi PDB ke PDBQT menggunakan meeko secara langsung
        # (Idealnya meeko bisa membaca RDKit mol dari PDB)
        try:
            mol = Chem.MolFromPDBFile(ref_pdb)
            if mol is not None:
                mol = Chem.AddHs(mol)
                preparator = MoleculePreparation()
                preparator.prepare(mol)
                ligand_pdbqt = os.path.join(results_dir, f"{protein}_native.pdbqt")
                with open(ligand_pdbqt, 'w') as f:
                    f.write(preparator.write_pdbqt_string())
            else:
                print(f"Gagal memproses ligan {params['ligand']}. RDKit tidak bisa membaca PDB-nya.")
                continue
                
        except Exception as e:
            print(f"Gagal konversi ke PDBQT: {e}")
            continue
            
        # Docking
        out_pdbqt = os.path.join(results_dir, f"{protein}_native_docked.pdbqt")
        c_x, c_y, c_z = params['center']
        s_x, s_y, s_z = params['size']
        receptor_path = os.path.join(receptor_dir, f"{protein}.pdbqt")
        
        cmd = [
            vina_exec,
            '--receptor', receptor_path,
            '--ligand', ligand_pdbqt,
            '--center_x', str(c_x),
            '--center_y', str(c_y),
            '--center_z', str(c_z),
            '--size_x', str(s_x),
            '--size_y', str(s_y),
            '--size_z', str(s_z),
            '--out', out_pdbqt,
            '--exhaustiveness', '8'
        ]
        
        print("Sedang melakukan Redocking (Vina)...")
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Hitung RMSD
        print("Menghitung RMSD...")
        rmsd = calculate_rmsd_from_pdb(ligand_pdbqt, out_pdbqt)
        
        if rmsd is not None:
            results.append({
                'Protein': protein,
                'Native_Ligand': params['ligand'],
                'RMSD_Angstrom': round(rmsd, 3),
                'Status': 'VALID' if rmsd <= 2.0 else 'TIDAK VALID'
            })
            print(f"Hasil RMSD {protein}: {rmsd:.3f} Å")
        else:
            print("Gagal menghitung RMSD.")
            
    df_results = pd.DataFrame(results)
    df_results.to_csv(os.path.join(results_dir, 'validation_rmsd.csv'), index=False)
    print("\\n=== HASIL VALIDASI REDOCKING ===")
    print(df_results)

if __name__ == "__main__":
    run_validation()
