import os
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
