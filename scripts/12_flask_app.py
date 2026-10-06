from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import os
import subprocess
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from meeko import MoleculePreparation
import base64

app = Flask(__name__)
CORS(app)

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
df_raw = pd.read_csv(os.path.join(project_root, 'data', '01_raw', 'dpp4_activities_raw.csv'))

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/dock', methods=['POST'])
def dock():
    data = request.json
    smiles = data.get('smiles')
    protein = data.get('protein')
    
    if not smiles or not protein:
        return jsonify({'error': 'SMILES and protein target are required'}), 400
        
    try:
        # RDKit Analysis
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return jsonify({'error': 'Invalid SMILES string'}), 400
            
        canonical = Chem.MolToSmiles(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        # Match ChEMBL ID
        chembl_id = "Novel Compound"
        ic50 = "Unknown"
        pref_name = "Unnamed"
        
        match = df_raw[df_raw['canonical_smiles'] == canonical]
        if not match.empty:
            chembl_id = match.iloc[0]['molecule_chembl_id']
            ic50 = f"{match.iloc[0]['standard_value']} {match.iloc[0]['standard_units']}"
            try:
                resp = requests.get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json", timeout=3)
                if resp.status_code == 200:
                    pref_name = resp.json().get('pref_name') or "Unnamed"
            except:
                pass
                
        # Docking
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_ligand = preparator.write_pdbqt_string()
        
        temp_ligand = os.path.join(project_root, 'results', 'flask_ligand.pdbqt')
        temp_out = os.path.join(project_root, 'results', 'flask_out.pdbqt')
        with open(temp_ligand, 'w') as f: f.write(pdbqt_ligand)
        
        vina_exec = os.path.join(project_root, 'scripts', 'vina')
        receptor_path = os.path.join(project_root, 'data', '05_receptors_pdbqt', f"{protein}.pdbqt")
        
        targets = {
            '3G0B': {'center': (42.05, 34.29, 14.62), 'size': (18, 18, 18)},
            '5T4B': {'center': (37.61, 49.92, 40.33), 'size': (18, 16, 23)}
        }
        params = targets[protein]
        
        cmd = [
            vina_exec, '--receptor', receptor_path, '--ligand', temp_ligand,
            '--center_x', str(params['center'][0]), '--center_y', str(params['center'][1]), '--center_z', str(params['center'][2]),
            '--size_x', str(params['size'][0]), '--size_y', str(params['size'][1]), '--size_z', str(params['size'][2]),
            '--out', temp_out, '--exhaustiveness', '1'
        ]
        
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        
        poses = []
        with open(temp_out, 'r') as f:
            out_ligand_data = f.read()
            
        with open(os.path.join(project_root, 'data', '03_pdb', f"{protein}_clean.pdb"), 'r') as f:
            rec_data = f.read()
            
        # Mocking parsing log from stdout (just parse from pdbqt remark lines)
        # Actually Vina outputs Affinity in the PDBQT REMARK VINA RESULT
        for line in out_ligand_data.split('\n'):
            if line.startswith('REMARK VINA RESULT:'):
                parts = line.split()
                if len(parts) >= 4:
                    poses.append({
                        'Affinity': float(parts[3]),
                        'RMSD_lb': float(parts[4]),
                        'RMSD_ub': float(parts[5])
                    })
                    
        return jsonify({
            'success': True,
            'molecule': {
                'smiles': canonical,
                'formula': formula,
                'chembl_id': chembl_id,
                'name': pref_name,
                'ic50': ic50,
                'mw': round(mw, 2),
                'logp': round(logp, 2),
                'hbd': hbd,
                'hba': hba
            },
            'docking': poses,
            '3d': {
                'receptor': rec_data,
                'ligand': out_ligand_data
            }
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(port=8080)
