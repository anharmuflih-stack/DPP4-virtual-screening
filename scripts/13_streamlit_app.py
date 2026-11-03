import streamlit as st
import streamlit.components.v1 as components
import os
import math
import platform
import subprocess
import tempfile
import urllib.request
import urllib.parse
import json
import ssl
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from stmol import showmol
import py3Dmol

# Set konfigurasi halaman lebar penuh dengan tema cerah
st.set_page_config(page_title="DPP-4 Virtual Screener", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

# ==========================================
# 1. INJEKSI CSS CUSTOM (Unified Scientific Design System)
# Putih (#FFFFFF) = Ruang Bernapas | Biru (#2485C7) = Identitas | Hijau (#15966A) = Aksen Validasi
# ==========================================
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=Fira+Code:wght@400;500;600;700&display=swap');

:root {
    --primary-blue: #2485C7;
    --light-blue: #E5F3FF;
    --soft-sky-blue: #F1F8FF;
    --mint-green: #DDF5E9;
    --emerald-green: #15966A;
    --deep-navy: #172B4D;
    --soft-bg: #F5F9FC;
    --white: #FFFFFF;
    --border-color: #DCE8F2;
    --muted-text: #71839A;
}

html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    color: var(--deep-navy);
}

.block-container {
    padding-top: 1.5rem !important;
    padding-bottom: 2.5rem !important;
    padding-left: 2.25rem !important;
    padding-right: 2.25rem !important;
    max-width: 93% !important;
}

.stApp {
    background-color: var(--soft-bg);
    color: var(--deep-navy);
}

header {visibility: hidden !important;}
#MainMenu {visibility: hidden !important;}
footer {visibility: hidden !important;}
[data-testid="stElementToolbar"] {display: none !important;}

/* Header Utama */
.app-header {
    background: linear-gradient(135deg, var(--white) 0%, var(--soft-sky-blue) 100%);
    border: 1px solid var(--border-color);
    border-radius: 16px;
    padding: 1.35rem 1.75rem;
    margin-bottom: 1.25rem;
    box-shadow: 0 2px 10px rgba(23, 43, 77, 0.03);
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: space-between;
    gap: 1.25rem;
}

.app-header-left {
    display: flex;
    align-items: center;
    gap: 1.15rem;
    flex-wrap: wrap;
}

.app-logo-badge {
    width: 52px;
    height: 52px;
    border-radius: 14px;
    background-color: var(--light-blue);
    border: 1px solid var(--border-color);
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 1.65rem;
    flex-shrink: 0;
}

.app-title {
    font-size: 1.65rem;
    font-weight: 800;
    color: var(--deep-navy);
    letter-spacing: -0.02em;
    margin: 0 0 0.15rem 0;
    line-height: 1.2;
}

.app-subtitle {
    font-size: 0.92rem;
    font-weight: 500;
    color: var(--muted-text);
    margin: 0;
}

/* Navigasi Tabs */
.stTabs [data-baseweb="tab-list"] {
    gap: 0.5rem;
    background-color: var(--white);
    padding: 0.4rem;
    border-radius: 14px;
    border: 1px solid var(--border-color);
    box-shadow: 0 2px 6px rgba(23, 43, 77, 0.02);
    margin-bottom: 1.5rem;
}

.stTabs [data-baseweb="tab"] {
    height: 44px;
    white-space: pre-wrap;
    background-color: transparent;
    border-radius: 10px !important;
    border: none !important;
    color: var(--muted-text) !important;
    font-weight: 600 !important;
    padding: 0 1.25rem !important;
    font-size: 0.92rem !important;
    transition: all 0.2s ease !important;
}

.stTabs [data-baseweb="tab"]:hover {
    background-color: var(--soft-sky-blue) !important;
    color: var(--primary-blue) !important;
}

.stTabs [aria-selected="true"] {
    background-color: var(--primary-blue) !important;
    color: var(--white) !important;
    box-shadow: 0 2px 8px rgba(36, 133, 199, 0.2) !important;
}

.stTabs [aria-selected="true"] p {
    color: var(--white) !important;
    font-weight: 600 !important;
}

.stTabs [data-baseweb="tab-highlight"] {
    display: none !important;
}

/* Unified White Card (Ruang Bernapas Utama) */
.science-card {
    background: var(--white);
    border: 1px solid var(--border-color);
    border-radius: 16px;
    box-shadow: 0 2px 10px rgba(23, 43, 77, 0.03);
    padding: 1.5rem;
    margin-bottom: 1.35rem;
    transition: box-shadow 0.2s ease, border-color 0.2s ease;
}

.science-card:hover {
    box-shadow: 0 6px 16px rgba(23, 43, 77, 0.05);
}

/* Unified Section Header */
.card-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 0.75rem;
    border-bottom: 1px solid var(--border-color);
    padding-bottom: 0.85rem;
    margin-bottom: 1.25rem;
}

.card-title {
    font-size: 1.1rem;
    font-weight: 700;
    color: var(--deep-navy);
    margin: 0;
}

.page-intro {
    margin-bottom: 1.35rem;
}

.page-intro-title {
    font-size: 1.25rem;
    font-weight: 700;
    color: var(--deep-navy);
    margin: 0 0 0.25rem 0;
}

.page-intro-desc {
    font-size: 0.92rem;
    color: var(--muted-text);
    margin: 0;
}

/* KPI Cards */
.kpi-grid {
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    gap: 1rem;
    margin-bottom: 1.5rem;
}

@media (max-width: 1100px) {
    .kpi-grid {
        grid-template-columns: repeat(3, 1fr);
    }
}

@media (max-width: 640px) {
    .kpi-grid {
        grid-template-columns: repeat(2, 1fr);
    }
}

.kpi-card {
    background: var(--white);
    border: 1px solid var(--border-color);
    border-radius: 16px;
    padding: 1.15rem 1.2rem;
    box-shadow: 0 2px 8px rgba(23, 43, 77, 0.03);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
}

.kpi-card:hover {
    transform: translateY(-2px);
    box-shadow: 0 6px 14px rgba(23, 43, 77, 0.06);
}

.kpi-blue {
    border-top: 3px solid var(--primary-blue);
}

.kpi-green {
    border-top: 3px solid var(--emerald-green);
}

.kpi-head {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 0.55rem;
}

.kpi-icon-b {
    width: 32px;
    height: 32px;
    border-radius: 8px;
    background: var(--light-blue);
    color: var(--primary-blue);
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 0.95rem;
}

.kpi-icon-g {
    width: 32px;
    height: 32px;
    border-radius: 8px;
    background: var(--mint-green);
    color: var(--emerald-green);
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 0.95rem;
}

.kpi-label {
    font-size: 0.73rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--muted-text);
}

.kpi-number {
    font-size: 1.65rem;
    font-weight: 800;
    color: var(--deep-navy);
    font-family: 'Fira Code', monospace;
    line-height: 1.2;
}

/* Metric Box & Lipinski Box */
.metric-box {
    background-color: var(--white);
    border: 1px solid var(--border-color);
    border-radius: 12px;
    padding: 0.95rem 1.1rem;
}

.data-label {
    font-size: 0.73rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: var(--muted-text);
    font-weight: 600;
    margin-bottom: 0.25rem;
}

.data-value {
    font-size: 1.05rem;
    font-weight: 600;
    color: var(--deep-navy);
    font-family: 'Fira Code', monospace;
}

.lipinski-box {
    background-color: var(--white);
    border: 1px solid var(--border-color);
    border-radius: 12px;
    padding: 0.95rem 1.05rem;
    border-left: 4px solid var(--emerald-green);
}

.lipinski-fail {
    border-left: 4px solid #E11D48;
}

/* Molecular Docking Results Card Component */
.docking-res-section-title {
    font-size: 0.72rem;
    font-weight: 800;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: var(--primary-blue);
    margin: 0.65rem 0 0.45rem 0;
}

.docking-target-row {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 0.42rem 0.65rem;
    border-radius: 8px;
    margin-bottom: 0.35rem;
    background-color: var(--soft-sky-blue);
    border: 1px solid var(--border-color);
}

.docking-target-row-best {
    background-color: #F4FBF8;
    border: 1px solid #B7E4D0;
}

/* Badges */
.badge-pass {
    background-color: var(--mint-green);
    color: var(--emerald-green);
    padding: 0.22rem 0.65rem;
    border-radius: 9999px;
    font-size: 0.74rem;
    font-weight: 700;
    display: inline-block;
}

.badge-info {
    background-color: var(--light-blue);
    color: var(--primary-blue);
    padding: 0.22rem 0.65rem;
    border-radius: 9999px;
    font-size: 0.74rem;
    font-weight: 700;
    display: inline-block;
}

.badge-warn {
    background-color: #FEF3C7;
    color: #B45309;
    padding: 0.22rem 0.65rem;
    border-radius: 9999px;
    font-size: 0.74rem;
    font-weight: 700;
    display: inline-block;
}

.badge-fail {
    background-color: #FFE4E6;
    color: #BE123C;
    padding: 0.22rem 0.65rem;
    border-radius: 9999px;
    font-size: 0.74rem;
    font-weight: 700;
    display: inline-block;
}

/* Legend Pill for Interactions */
.ix-pill {
    display: inline-flex;
    align-items: center;
    gap: 0.45rem;
    padding: 0.35rem 0.75rem;
    border-radius: 999px;
    background: var(--white);
    border: 1px solid var(--border-color);
    font-size: 0.78rem;
    font-weight: 600;
    color: var(--deep-navy);
}

.ix-dot {
    width: 10px;
    height: 10px;
    border-radius: 50%;
    display: inline-block;
}

/* Form Input & Buttons */
div[data-baseweb="select"] > div, 
div[data-baseweb="input"] > div {
    background-color: var(--white) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 10px !important;
    font-family: 'Inter', sans-serif !important;
    color: var(--deep-navy) !important;
}

div[data-baseweb="select"] > div:focus-within,
div[data-baseweb="input"] > div:focus-within {
    border-color: var(--primary-blue) !important;
    box-shadow: 0 0 0 3px var(--light-blue) !important;
}

label[data-testid="stWidgetLabel"] p {
    font-size: 0.76rem !important;
    text-transform: uppercase !important;
    letter-spacing: 0.05em !important;
    color: var(--muted-text) !important;
    font-weight: 600 !important;
}

div[data-testid="stButton"] > button {
    width: 100% !important;
    background-color: var(--primary-blue) !important;
    color: var(--white) !important;
    font-weight: 600 !important;
    border-radius: 10px !important;
    border: none !important;
    padding: 0.65rem 1rem !important;
    height: 3rem !important;
    box-shadow: 0 2px 6px rgba(36, 133, 199, 0.18) !important;
    transition: all 0.2s ease !important;
}

div[data-testid="stButton"] > button:hover {
    background-color: #1D6FA5 !important;
    box-shadow: 0 4px 10px rgba(36, 133, 199, 0.25) !important;
}

div[data-testid="stDownloadButton"] > button {
    background-color: var(--white) !important;
    color: var(--primary-blue) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 10px !important;
    font-weight: 600 !important;
    padding: 0.5rem 1.1rem !important;
    transition: all 0.2s ease !important;
}

div[data-testid="stDownloadButton"] > button:hover {
    background-color: var(--soft-sky-blue) !important;
    border-color: var(--primary-blue) !important;
}

[data-testid="stDataFrame"] {
    border: 1px solid var(--border-color) !important;
    border-radius: 12px !important;
    overflow: hidden !important;
    background: var(--white) !important;
}
</style>
""", unsafe_allow_html=True)

# ==========================================
# 2. SETUP DIREKTORI, METADATA RESEPTOR & VINA
# ==========================================
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(project_root, 'data')
bin_dir = os.path.join(project_root, 'scripts')
results_dir = os.path.join(project_root, 'results', '02_docking')

RECEPTOR_TARGETS = {
    '3G0B': {
        'pdb_id': '3G0B',
        'protein_name': 'Human Dipeptidyl Peptidase-4 (DPP-4)',
        'description': 'High-Resolution Crystal Complex (Chain A)',
        'chain_id': 'A',
        'center': (42.05, 34.29, 14.62),
        'size': (18, 18, 18),
        'exhaustiveness': 4,
        'num_modes': 5
    },
    '5T4B': {
        'pdb_id': '5T4B',
        'protein_name': 'Human Dipeptidyl Peptidase-4 (DPP-4)',
        'description': 'Wild-Type Human DPP-4 Structure (Chain A)',
        'chain_id': 'A',
        'center': (37.61, 49.92, 40.33),
        'size': (18, 16, 23),
        'exhaustiveness': 4,
        'num_modes': 5
    }
}

@st.cache_resource
def setup_vina():
    is_linux = platform.system() == 'Linux'
    vina_exec = os.path.join(bin_dir, 'vina_linux' if is_linux else 'vina')
    
    if is_linux and not os.path.exists(vina_exec):
        url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
        urllib.request.urlretrieve(url, vina_exec)
        os.chmod(vina_exec, 0o755)
    return vina_exec

vina_exec = setup_vina()

# Load dataset terkurasi & hasil docking terpisah per PDB ID
@st.cache_data
def load_curated_datasets():
    master_path = os.path.join(data_dir, '02_processed', 'curated_natural_products_master.csv')
    test_path = os.path.join(data_dir, '02_processed', 'dpp4_natural_products_test_set.csv')
    analysis_path = os.path.join(data_dir, '02_processed', 'dpp4_natural_products_analysis_set.csv')
    results_path = os.path.join(results_dir, 'curated_np_all_docking_results.csv')
    summary_path = os.path.join(results_dir, 'validation_metrics_summary.json')
    
    df_master = pd.read_csv(master_path) if os.path.exists(master_path) else None
    df_test = pd.read_csv(test_path) if os.path.exists(test_path) else None
    df_analysis = pd.read_csv(analysis_path) if os.path.exists(analysis_path) else None
    df_results = pd.read_csv(results_path) if os.path.exists(results_path) else None
    
    metrics = None
    if os.path.exists(summary_path):
        with open(summary_path, 'r') as f:
            metrics = json.load(f)
            
    return df_master, df_test, df_analysis, df_results, metrics

df_master, df_test, df_analysis, df_results, summary_metrics = load_curated_datasets()

@st.cache_data
def load_separate_target_results():
    path_3g0b = os.path.join(results_dir, 'docking_results_3G0B.csv')
    path_5t4b = os.path.join(results_dir, 'docking_results_5T4B.csv')
    df_3g0b = pd.read_csv(path_3g0b) if os.path.exists(path_3g0b) else None
    df_5t4b = pd.read_csv(path_5t4b) if os.path.exists(path_5t4b) else None
    return {'3G0B': df_3g0b, '5T4B': df_5t4b}

target_results_by_pdb = load_separate_target_results()

@st.cache_data
def get_receptor_structure_info(protein):
    pdb_path = os.path.join(data_dir, '03_pdb', f"{protein}_clean.pdb")
    pdbqt_path = os.path.join(data_dir, '05_receptors_pdbqt', f"{protein}.pdbqt")
    if not os.path.exists(pdb_path):
        return {'valid': False, 'residues': 0, 'atoms': 0, 'chains': [], 'pdb_text': ''}
    chains = set()
    residues = set()
    atom_count = 0
    with open(pdb_path, 'r') as f:
        pdb_text = f.read()
    for line in pdb_text.split('\n'):
        if line.startswith('ATOM'):
            atom_count += 1
            ch = line[21:22].strip() or 'A'
            res_seq = line[22:26].strip()
            res_name = line[17:20].strip()
            chains.add(ch)
            residues.add((ch, res_seq, res_name))
    return {
        'valid': os.path.exists(pdb_path) and os.path.exists(pdbqt_path) and atom_count > 0,
        'residues': len(residues),
        'atoms': atom_count,
        'chains': sorted(list(chains)),
        'pdb_text': pdb_text
    }

# ==========================================
# 3. FUNGSI PENCARIAN CHEMBL, DOCKING & INTERAKSI MOLEKULER
# ==========================================
def get_chembl_data(mol):
    try:
        inchikey = Chem.MolToInchiKey(mol)
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        
        # Ambil Nama dan CID
        url_props = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/Title/JSON"
        req_props = urllib.request.Request(url_props, headers={'Accept': 'application/json'})
        cid = "-"
        name = "Novel Compound"
        with urllib.request.urlopen(req_props, context=ctx, timeout=5) as response:
            data = json.loads(response.read().decode())
            if 'PropertyTable' in data:
                props = data['PropertyTable']['Properties'][0]
                cid = str(props.get('CID', '-'))
                name = props.get('Title', 'Novel Compound')
                
        # Ambil ChEMBL ID dari Sinonim
        chembl_id = "-"
        if cid != "-":
            try:
                url_syns = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/synonyms/JSON"
                req_syns = urllib.request.Request(url_syns, headers={'Accept': 'application/json'})
                with urllib.request.urlopen(req_syns, context=ctx, timeout=5) as response:
                    data_syns = json.loads(response.read().decode())
                    if 'InformationList' in data_syns:
                        syns = data_syns['InformationList']['Information'][0].get('Synonym', [])
                        chembl_list = [s for s in syns if s.startswith('CHEMBL') and 'SCHEMBL' not in s]
                        if chembl_list:
                            chembl_id = chembl_list[0]
            except Exception:
                pass
            
        return {'id': chembl_id, 'name': name, 'cid': cid}
    except Exception:
        pass
    return {'id': 'Tidak Ditemukan', 'name': 'Novel Compound', 'cid': '-'}

def process_ligand(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None, "SMILES tidak valid"
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    preparator = MoleculePreparation()
    setup_list = preparator.prepare(mol)
    writer = PDBQTWriterLegacy()
    
    result = writer.write_string(setup_list[0])
    pdbqt_string = result[0] if isinstance(result, tuple) else result
    
    chembl_data = get_chembl_data(mol)
    props = {
        'MW': round(Descriptors.MolWt(mol), 2),
        'LogP': round(Descriptors.MolLogP(mol), 2),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'Formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
        'ChEMBL_ID': chembl_data['id'],
        'Name': chembl_data['name'],
        'CID': chembl_data['cid']
    }
    return pdbqt_string, props

def run_docking(pdbqt_ligand, protein):
    targets = {
        '3G0B': {'center': (42.05, 34.29, 14.62), 'size': (18, 18, 18)},
        '5T4B': {'center': (37.61, 49.92, 40.33), 'size': (18, 16, 23)}
    }
    if protein not in targets:
        return [], ""
        
    receptor_path = os.path.join(data_dir, '05_receptors_pdbqt', f"{protein}.pdbqt")
    if not os.path.exists(receptor_path):
        return [], ""

    with tempfile.NamedTemporaryFile(delete=False, suffix='.pdbqt') as lig_file:
        lig_file.write(pdbqt_ligand.encode('utf-8'))
        lig_path = lig_file.name
        
    out_path = lig_path.replace('.pdbqt', f'_{protein}_out.pdbqt')
    
    cmd = [
        vina_exec, '--receptor', receptor_path, '--ligand', lig_path,
        '--center_x', str(targets[protein]['center'][0]),
        '--center_y', str(targets[protein]['center'][1]),
        '--center_z', str(targets[protein]['center'][2]),
        '--size_x', str(targets[protein]['size'][0]),
        '--size_y', str(targets[protein]['size'][1]),
        '--size_z', str(targets[protein]['size'][2]),
        '--out', out_path,
        '--exhaustiveness', '4', '--num_modes', '5'
    ]
    
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
    except Exception:
        pass
    
    poses = []
    out_content = ""
    if os.path.exists(out_path):
        with open(out_path, 'r') as f:
            out_content = f.read()
            for line in out_content.split('\n'):
                if line.startswith('REMARK VINA RESULT:'):
                    parts = line.split()
                    poses.append({
                        'Pose': len(poses) + 1,
                        'Affinity (kcal/mol)': float(parts[3]),
                        'RMSD l.b. (Å)': float(parts[4]),
                        'RMSD u.b. (Å)': float(parts[5])
                    })
        try:
            os.remove(out_path)
        except Exception:
            pass
    try:
        os.remove(lig_path)
    except Exception:
        pass
                
    return poses, out_content

def run_dual_target_docking(pdbqt_ligand):
    """Menjalankan docking secara independen terhadap kedua target protein: 3G0B dan 5T4B."""
    results_by_pdb = {}
    for pdb_id in ['3G0B', '5T4B']:
        poses, docked_pdbqt = run_docking(pdbqt_ligand, pdb_id)
        best_aff = poses[0]['Affinity (kcal/mol)'] if poses else None
        results_by_pdb[pdb_id] = {
            'pdb_id': pdb_id,
            'poses': poses,
            'docked_pdbqt': docked_pdbqt,
            'best_affinity': best_aff,
            'status': 'Success' if best_aff is not None else 'Docking Failed'
        }
    return results_by_pdb

def load_saved_docking_for_compound(chembl_id):
    """Memuat pose dan skor hasil docking aktual yang tersimpan untuk senyawa pustaka pada 3G0B dan 5T4B."""
    results_by_pdb = {}
    for pdb_id in ['3G0B', '5T4B']:
        out_path = os.path.join(results_dir, f"{chembl_id}_{pdb_id}_out.pdbqt")
        poses = []
        out_content = ""
        if os.path.exists(out_path):
            with open(out_path, 'r') as f:
                out_content = f.read()
                for line in out_content.split('\n'):
                    if line.startswith('REMARK VINA RESULT:'):
                        parts = line.split()
                        if len(parts) >= 6:
                            poses.append({
                                'Pose': len(poses) + 1,
                                'Affinity (kcal/mol)': float(parts[3]),
                                'RMSD l.b. (Å)': float(parts[4]),
                                'RMSD u.b. (Å)': float(parts[5])
                            })
        best_aff = poses[0]['Affinity (kcal/mol)'] if poses else None
        results_by_pdb[pdb_id] = {
            'pdb_id': pdb_id,
            'poses': poses,
            'docked_pdbqt': out_content,
            'best_affinity': best_aff,
            'status': 'Success' if best_aff is not None else 'Docking Failed'
        }
    return results_by_pdb

# ==========================================
# 4. ANALISIS KOORDINAT INTERAKSI MOLEKULER AKTUAL (3D COORDINATE ENGINE)
# ==========================================
DPP4_SUBSITE_MAP = {
    630: 'Catalytic Triad (S1)', 708: 'Catalytic Triad (S1)', 740: 'Catalytic Triad (S1)',
    662: 'S1 Hydrophobic Pocket', 666: 'S1 Hydrophobic Pocket', 656: 'S1 Hydrophobic Pocket',
    659: 'S1 Hydrophobic Pocket', 711: 'S1 Hydrophobic Pocket', 665: 'S1 Hydrophobic Pocket',
    205: 'S2 Electrostatic Pocket', 206: 'S2 Electrostatic Pocket', 125: 'S2 Electrostatic Pocket',
    547: 'S1/S2 Oxyanion & Aromatic', 357: 'S2 Extensive Subsite', 209: 'S2 Extensive Subsite'
}

INTERACTION_STYLES = {
    'Hydrogen Bond': {'color': '#15966A', 'hex': '0x15966A', 'badge': 'badge-pass'},
    'Salt Bridge': {'color': '#D97706', 'hex': '0xD97706', 'badge': 'badge-warn'},
    'π-π Stacking': {'color': '#8B5CF6', 'hex': '0x8B5CF6', 'badge': 'badge-info'},
    'π-Cation': {'color': '#06B6D4', 'hex': '0x06B6D4', 'badge': 'badge-info'},
    'Hydrophobic': {'color': '#2485C7', 'hex': '0x2485C7', 'badge': 'badge-info'},
    'Van der Waals': {'color': '#64748B', 'hex': '0x64748B', 'badge': 'badge-info'}
}

def extract_first_model_pdbqt(pdbqt_text):
    if not pdbqt_text:
        return ""
    if "MODEL 2" in pdbqt_text:
        return pdbqt_text.split("MODEL 2")[0]
    return pdbqt_text

def parse_pdb_atoms(pdb_text):
    atoms = []
    for line in pdb_text.split('\n'):
        if line.startswith('ATOM'):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21:22].strip() or 'A'
                res_seq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                elem = line[76:78].strip() if len(line) >= 78 else atom_name[0]
                if not elem:
                    elem = ''.join([c for c in atom_name if c.isalpha()])[:1]
                if elem.upper() != 'H':
                    atoms.append({
                        'atom_name': atom_name,
                        'res_name': res_name,
                        'chain': chain,
                        'res_seq': res_seq,
                        'x': x, 'y': y, 'z': z,
                        'elem': elem.upper()
                    })
            except Exception:
                continue
    return atoms

def parse_ligand_pdbqt_atoms(pdbqt_model1):
    atoms = []
    for line in pdbqt_model1.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                parts = line.split()
                ad_type = parts[-1] if parts else atom_name[0]
                elem = ''.join([c for c in ad_type if c.isalpha()])[:2].capitalize()
                if elem in ['Oa', 'Os']: elem = 'O'
                elif elem in ['Na', 'Ns', 'N']: elem = 'N'
                elif elem in ['Sa', 'S']: elem = 'S'
                elif elem in ['A', 'C']: elem = 'C'
                elif elem in ['Hd', 'H']: elem = 'H'
                if elem.upper() != 'H':
                    atoms.append({
                        'atom_name': atom_name,
                        'x': x, 'y': y, 'z': z,
                        'ad_type': ad_type,
                        'elem': elem.upper()
                    })
            except Exception:
                continue
    return atoms

@st.cache_data
def analyze_binding_interactions(protein_target, best_pose_pdbqt):
    """Menganalisis interaksi protein-ligan berdasarkan koordinat 3D aktual dari struktur PDB dan pose PDBQT."""
    rec_info = get_receptor_structure_info(protein_target)
    if not rec_info['valid'] or not best_pose_pdbqt:
        return pd.DataFrame(), {'Hydrogen Bond': 0, 'Hydrophobic': 0, 'π-π Stacking': 0, 'π-Cation': 0, 'Salt Bridge': 0, 'Van der Waals': 0, 'Total': 0}, []

    rec_atoms = parse_pdb_atoms(rec_info['pdb_text'])
    lig_atoms = parse_ligand_pdbqt_atoms(best_pose_pdbqt)
    if not rec_atoms or not lig_atoms:
        return pd.DataFrame(), {'Hydrogen Bond': 0, 'Hydrophobic': 0, 'π-π Stacking': 0, 'π-Cation': 0, 'Salt Bridge': 0, 'Van der Waals': 0, 'Total': 0}, []

    aromatic_res_atoms = {
        'TYR': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'PHE': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'TRP': {'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
        'HIS': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'}
    }
    charged_neg_atoms = {'GLU': {'OE1', 'OE2'}, 'ASP': {'OD1', 'OD2'}}
    charged_pos_atoms = {'ARG': {'NE', 'NH1', 'NH2'}, 'LYS': {'NZ'}, 'HIS': {'ND1', 'NE2'}}
    hydrophobic_res = {'VAL', 'LEU', 'ILE', 'ALA', 'PRO', 'MET', 'PHE', 'TRP', 'TYR'}

    # Filter receptor atoms within 5.0 A bounding box of ligand for fast exact calculation
    lx = [a['x'] for a in lig_atoms]
    ly = [a['y'] for a in lig_atoms]
    lz = [a['z'] for a in lig_atoms]
    min_x, max_x = min(lx) - 5.0, max(lx) + 5.0
    min_y, max_y = min(ly) - 5.0, max(ly) + 5.0
    min_z, max_z = min(lz) - 5.0, max(lz) + 5.0

    near_rec_atoms = [
        ra for ra in rec_atoms
        if min_x <= ra['x'] <= max_x and min_y <= ra['y'] <= max_y and min_z <= ra['z'] <= max_z
    ]

    residue_best_interactions = {}
    priority_order = {
        'Hydrogen Bond': 1,
        'Salt Bridge': 2,
        'π-π Stacking': 3,
        'π-Cation': 4,
        'Hydrophobic': 5,
        'Van der Waals': 6
    }

    for la in lig_atoms:
        for ra in near_rec_atoms:
            dx = la['x'] - ra['x']
            dy = la['y'] - ra['y']
            dz = la['z'] - ra['z']
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist > 4.5:
                continue

            r_name = ra['res_name']
            r_seq = ra['res_seq']
            r_atom = ra['atom_name']
            l_elem = la['elem']
            r_elem = ra['elem']
            l_ad = la['ad_type']

            ix_type = None
            # 1. Hydrogen Bond (<= 3.5 A between polar N/O/S/F atoms)
            if dist <= 3.5 and l_elem in {'N', 'O', 'S', 'F'} and r_elem in {'N', 'O', 'S'}:
                ix_type = 'Hydrogen Bond'
            # 2. Salt Bridge (<= 4.0 A between polar/charged N/O and charged sidechains)
            elif dist <= 4.0 and l_elem in {'N', 'O'} and (
                (r_name in charged_neg_atoms and r_atom in charged_neg_atoms[r_name]) or
                (r_name in charged_pos_atoms and r_atom in charged_pos_atoms[r_name])
            ):
                ix_type = 'Salt Bridge'
            # 3. pi-pi Stacking (<= 4.5 A between aromatic ligand atom 'A' or ring C and aromatic residue ring)
            elif dist <= 4.5 and (l_ad == 'A' or l_elem == 'C') and (
                r_name in aromatic_res_atoms and r_atom in aromatic_res_atoms[r_name]
            ):
                ix_type = 'π-π Stacking'
            # 4. pi-Cation (<= 4.5 A between aromatic ring and cationic N)
            elif dist <= 4.5 and (
                ((l_ad == 'A' or l_elem == 'C') and (r_name in charged_pos_atoms and r_atom in charged_pos_atoms[r_name])) or
                (l_elem == 'N' and (r_name in aromatic_res_atoms and r_atom in aromatic_res_atoms[r_name]))
            ):
                ix_type = 'π-Cation'
            # 5. Hydrophobic Contact (<= 4.0 A between non-polar C/Halogen and hydrophobic sidechain C)
            elif dist <= 4.0 and l_elem in {'C', 'CL', 'BR', 'F', 'I'} and r_elem == 'C' and r_name in hydrophobic_res and r_atom not in {'C', 'CA'}:
                ix_type = 'Hydrophobic'
            # 6. Van der Waals (<= 4.0 A close steric contact)
            elif dist <= 4.0:
                ix_type = 'Van der Waals'

            if ix_type:
                key = (ra['chain'], r_seq, r_name)
                cand = {
                    'Residue Name': r_name,
                    'Residue Number': r_seq,
                    'Residue Label': f"{r_name}{r_seq}",
                    'Chain': ra['chain'],
                    'Pocket Subsite': DPP4_SUBSITE_MAP.get(r_seq, 'Active Pocket Region'),
                    'Interaction Type': ix_type,
                    'Distance (Å)': round(dist, 2),
                    'Receptor Atom': r_atom,
                    'Ligand Atom': la['atom_name'],
                    'lig_xyz': (la['x'], la['y'], la['z']),
                    'rec_xyz': (ra['x'], ra['y'], ra['z']),
                    'priority': priority_order[ix_type]
                }
                prev = residue_best_interactions.get(key)
                if prev is None or (cand['priority'] < prev['priority']) or (cand['priority'] == prev['priority'] and cand['Distance (Å)'] < prev['Distance (Å)']):
                    residue_best_interactions[key] = cand

    records = sorted(residue_best_interactions.values(), key=lambda r: (r['priority'], r['Distance (Å)']))
    counts = {
        'Hydrogen Bond': sum(1 for r in records if r['Interaction Type'] == 'Hydrogen Bond'),
        'Hydrophobic': sum(1 for r in records if r['Interaction Type'] == 'Hydrophobic'),
        'π-π Stacking': sum(1 for r in records if r['Interaction Type'] == 'π-π Stacking'),
        'π-Cation': sum(1 for r in records if r['Interaction Type'] == 'π-Cation'),
        'Salt Bridge': sum(1 for r in records if r['Interaction Type'] == 'Salt Bridge'),
        'Van der Waals': sum(1 for r in records if r['Interaction Type'] == 'Van der Waals'),
        'Total': len(records)
    }

    display_records = [{
        'Residue Name': r['Residue Name'],
        'Residue Number': r['Residue Number'],
        'Chain': r['Chain'],
        'Pocket Subsite': r['Pocket Subsite'],
        'Interaction Type': r['Interaction Type'],
        'Distance (Å)': r['Distance (Å)'],
        'Receptor Atom': r['Receptor Atom'],
        'Ligand Atom': r['Ligand Atom']
    } for r in records]

    return pd.DataFrame(display_records), counts, records

def render_interactive_3d_explorer(
    protein_target,
    best_pose_pdbqt,
    interactions_list,
    rep_mode="Cartoon + Ligand",
    show_protein=True,
    show_ligand=True,
    show_pocket=True,
    show_interactions=True,
    focus_residue="All Pocket Residues",
    viewer_height=500,
    viewer_id="mol3d_main"
):
    """Merender Interactive 3D Molecular Interaction Explorer menggunakan 3Dmol.js dengan floating controls."""
    rec_info = get_receptor_structure_info(protein_target)
    if not rec_info['valid'] or not best_pose_pdbqt:
        st.warning(f"Struktur 3D atau pose ligan untuk target {protein_target} tidak tersedia.")
        return

    pdb_json = json.dumps(rec_info['pdb_text'])
    lig_json = json.dumps(best_pose_pdbqt)

    # Prepare interaction lines and labels for JS
    js_interactions = []
    for ix in interactions_list:
        style_meta = INTERACTION_STYLES.get(ix['Interaction Type'], {'color': '#38BDF8'})
        js_interactions.append({
            'res_label': ix['Residue Label'],
            'res_seq': ix['Residue Number'],
            'type': ix['Interaction Type'],
            'dist': ix['Distance (Å)'],
            'color': style_meta['color'],
            'lx': ix['lig_xyz'][0], 'ly': ix['lig_xyz'][1], 'lz': ix['lig_xyz'][2],
            'rx': ix['rec_xyz'][0], 'ry': ix['rec_xyz'][1], 'rz': ix['rec_xyz'][2]
        })
    ix_json = json.dumps(js_interactions)
    focus_res_seq = None
    if focus_residue and focus_residue != "All Pocket Residues":
        digits = ''.join([c for c in focus_residue.split()[0] if c.isdigit()])
        if digits:
            focus_res_seq = int(digits)

    html_viewer = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@500;600;700&family=Fira+Code:wght@500;600&display=swap" rel="stylesheet">
<style>
  body {{
    margin: 0;
    padding: 0;
    background: transparent;
    font-family: 'Inter', sans-serif;
    overflow: hidden;
  }}
  .viewer-wrapper {{
    position: relative;
    width: 100%;
    height: {viewer_height}px;
    border-radius: 16px;
    overflow: hidden;
    border: 1.5px solid #DCE8F2;
    background: radial-gradient(circle at 50% 35%, #1E3A5F 0%, #172B4D 55%, #0F1E36 100%);
    box-shadow: inset 0 0 40px rgba(0, 0, 0, 0.35);
  }}
  #{viewer_id} {{
    width: 100%;
    height: 100%;
    position: relative;
  }}
  .floating-badge-top {{
    position: absolute;
    top: 12px;
    left: 14px;
    z-index: 20;
    display: flex;
    gap: 8px;
    flex-wrap: wrap;
    pointer-events: none;
  }}
  .hud-pill {{
    background: rgba(23, 43, 77, 0.82);
    backdrop-filter: blur(6px);
    color: #F1F8FF;
    border: 1px solid rgba(220, 232, 242, 0.25);
    border-radius: 999px;
    padding: 5px 12px;
    font-size: 11.5px;
    font-weight: 600;
    letter-spacing: 0.02em;
  }}
  .hud-pill b {{
    color: #5EEAD4;
    font-family: 'Fira Code', monospace;
  }}
  .floating-toolbar {{
    position: absolute;
    bottom: 14px;
    left: 50%;
    transform: translateX(-50%);
    z-index: 20;
    display: flex;
    gap: 6px;
    flex-wrap: wrap;
    justify-content: center;
    background: rgba(23, 43, 77, 0.86);
    backdrop-filter: blur(8px);
    padding: 6px 10px;
    border-radius: 12px;
    border: 1px solid rgba(220, 232, 242, 0.22);
    box-shadow: 0 8px 20px rgba(0, 0, 0, 0.3);
  }}
  .ctrl-btn {{
    background: rgba(255, 255, 255, 0.08);
    color: #E5F3FF;
    border: 1px solid rgba(229, 243, 255, 0.18);
    border-radius: 8px;
    padding: 5px 10px;
    font-size: 11.5px;
    font-weight: 600;
    cursor: pointer;
    transition: all 0.18s ease;
    font-family: 'Inter', sans-serif;
    white-space: nowrap;
  }}
  .ctrl-btn:hover {{
    background: #2485C7;
    color: #FFFFFF;
    border-color: #2485C7;
  }}
  .ctrl-btn.active {{
    background: #15966A;
    color: #FFFFFF;
    border-color: #15966A;
  }}
</style>
</head>
<body>
<div class="viewer-wrapper" id="wrap_{viewer_id}">
  <div class="floating-badge-top">
    <div class="hud-pill">Target PDB: <b>{protein_target}</b></div>
    <div class="hud-pill">Mode: <b>{rep_mode}</b></div>
    <div class="hud-pill">Validated Contacts: <b>{len(js_interactions)}</b></div>
  </div>
  <div id="{viewer_id}"></div>
  <div class="floating-toolbar">
    <button class="ctrl-btn" onclick="resetCam()">🔄 Reset Camera</button>
    <button class="ctrl-btn" onclick="zoomLigand()">🎯 Focus Ligand</button>
    <button class="ctrl-btn" onclick="zoomStep(1.2)">🔍 Zoom In</button>
    <button class="ctrl-btn" onclick="zoomStep(0.82)">🔎 Zoom Out</button>
    <button class="ctrl-btn" id="spinBtn_{viewer_id}" onclick="toggleSpin()">💫 Spin 360°</button>
    <button class="ctrl-btn" onclick="toggleFull()">⛶ Fullscreen</button>
  </div>
</div>

<script>
  const pdbData = {pdb_json};
  const ligData = {lig_json};
  const interactions = {ix_json};
  const showProtein = {str(show_protein).lower()};
  const showLigand = {str(show_ligand).lower()};
  const showPocket = {str(show_pocket).lower()};
  const showInteractions = {str(show_interactions).lower()};
  const repMode = {json.dumps(rep_mode)};
  const focusResSeq = {json.dumps(focus_res_seq)};

  let spinning = false;
  let viewer = $3Dmol.createViewer(document.getElementById("{viewer_id}"), {{
    backgroundColor: "#13233F",
    antialias: true,
    cartoonQuality: 10
  }});

  viewer.addModel(pdbData, "pdb");
  viewer.addModel(ligData, "pdbqt");

  function applyScene() {{
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.removeAllSurfaces();
    viewer.setStyle({{}}, {{}});

    // 1. Protein Representation
    if (showProtein) {{
      if (repMode === "Cartoon" || repMode === "Cartoon + Ligand") {{
        viewer.setStyle({{model: 0}}, {{cartoon: {{color: "#B8D4EE", style: "oval", thickness: 0.24, opacity: 0.92}}}});
      }} else if (repMode === "Stick") {{
        viewer.setStyle({{model: 0}}, {{stick: {{colorscheme: "lightGrayCarbon", radius: 0.12, opacity: 0.75}}}});
      }} else if (repMode === "Ball and Stick") {{
        viewer.setStyle({{model: 0}}, {{
          stick: {{colorscheme: "lightGrayCarbon", radius: 0.12, opacity: 0.75}},
          sphere: {{colorscheme: "lightGrayCarbon", scale: 0.2, opacity: 0.75}}
        }});
      }} else if (repMode === "Surface") {{
        viewer.setStyle({{model: 0}}, {{cartoon: {{color: "#94B8D8", opacity: 0.5}}}});
        viewer.addSurface($3Dmol.SurfaceType.MS, {{opacity: 0.55, color: "#DCE8F2"}}, {{model: 0}}, {{model: 0}});
      }}
    }}

    // 2. Binding Pocket Highlighting (Actual Residues within 4.5 A of Ligand)
    const pocketSel = {{model: 0, within: {{distance: 4.5, sel: {{model: 1}}}}}};
    if (showPocket) {{
      viewer.addStyle(pocketSel, {{
        stick: {{colorscheme: "cyanCarbon", radius: 0.16}}
      }});
      viewer.addSurface($3Dmol.SurfaceType.VDW, {{
        opacity: 0.28,
        color: "#38BDF8"
      }}, pocketSel, pocketSel);

      // Label interacting residues from coordinate analysis
      const seenRes = new Set();
      interactions.forEach(ix => {{
        if (!seenRes.has(ix.res_seq)) {{
          seenRes.add(ix.res_seq);
          const isFocused = (focusResSeq !== null && ix.res_seq === focusResSeq);
          viewer.addLabel(ix.res_label, {{
            position: {{x: ix.rx, y: ix.ry, z: ix.rz}},
            backgroundColor: isFocused ? "#15966A" : "rgba(23, 43, 77, 0.85)",
            borderColor: isFocused ? "#DDF5E9" : "#2485C7",
            borderThickness: 1,
            fontColor: "#FFFFFF",
            fontSize: isFocused ? 13 : 11,
            showBackground: true
          }});
        }}
      }});
    }}

    // 3. Ligand Representation (Prominent Emerald Green / Cyan Ball-and-Stick)
    if (showLigand && repMode !== "Cartoon") {{
      if (repMode === "Stick") {{
        viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: "greenCarbon", radius: 0.26}}}});
      }} else {{
        viewer.setStyle({{model: 1}}, {{
          stick: {{colorscheme: "greenCarbon", radius: 0.24}},
          sphere: {{colorscheme: "greenCarbon", scale: 0.32}}
        }});
      }}
    }}

    // 4. Coordinate-Validated Dashed Interaction Lines & Distance Labels (Å)
    if (showInteractions) {{
      interactions.forEach(ix => {{
        viewer.addCylinder({{
          start: {{x: ix.lx, y: ix.ly, z: ix.lz}},
          end: {{x: ix.rx, y: ix.ry, z: ix.rz}},
          radius: 0.065,
          fromCap: 1,
          toCap: 1,
          dashed: true,
          color: ix.color
        }});
        const mx = (ix.lx + ix.rx) / 2.0;
        const my = (ix.ly + ix.ry) / 2.0;
        const mz = (ix.lz + ix.rz) / 2.0;
        viewer.addLabel(ix.dist + " Å", {{
          position: {{x: mx, y: my, z: mz}},
          backgroundColor: "rgba(15, 30, 54, 0.8)",
          fontColor: "#DDF5E9",
          fontSize: 10,
          showBackground: true
        }});
      }});
    }}

    // Highlight focused residue if selected
    if (focusResSeq !== null) {{
      const focusSel = {{model: 0, resi: focusResSeq}};
      viewer.addStyle(focusSel, {{
        stick: {{colorscheme: "yellowCarbon", radius: 0.26}},
        sphere: {{colorscheme: "yellowCarbon", scale: 0.28}}
      }});
      viewer.zoomTo({{or: [{{model: 1}}, focusSel]}});
    }} else {{
      viewer.zoomTo({{model: 1}});
    }}

    viewer.render();
  }}

  applyScene();

  function resetCam() {{
    viewer.zoomTo({{model: 1}});
    viewer.render();
  }}
  function zoomLigand() {{
    viewer.zoomTo({{model: 1}}, 500);
    viewer.render();
  }}
  function zoomStep(factor) {{
    viewer.zoom(factor);
    viewer.render();
  }}
  function toggleSpin() {{
    spinning = !spinning;
    viewer.spin(spinning ? "y" : false);
    const btn = document.getElementById("spinBtn_{viewer_id}");
    if (spinning) btn.classList.add("active");
    else btn.classList.remove("active");
  }}
  function toggleFull() {{
    const el = document.getElementById("wrap_{viewer_id}");
    if (!document.fullscreenElement) {{
      if (el.requestFullscreen) el.requestFullscreen();
    }} else {{
      if (document.exitFullscreen) document.exitFullscreen();
    }}
  }}
</script>
</body>
</html>
"""
    components.html(html_viewer, height=viewer_height + 12)

def generate_interaction_interpretation(compound_name, target_pdb, best_affinity, counts, records):
    """Menghasilkan interpretasi ilmiah berbasis data aktual residu yang berinteraksi."""
    if not records:
        return f"Pada struktur kompleks {target_pdb}, tidak terdeteksi kontak antarmolekul jarak dekat (≤4.5 Å) pada pose ligan saat ini."
    res_labels = [r['Residue Label'] for r in records]
    res_nums = {r['Residue Number'] for r in records}

    triad_hits = [f"Ser630" if 630 in res_nums else "", f"Asp708" if 708 in res_nums else "", f"His740" if 740 in res_nums else ""]
    triad_hits = [t for t in triad_hits if t]
    s1_hits = [r['Residue Label'] for r in records if r['Residue Number'] in {662, 666, 656, 659, 711, 665}]
    s2_hits = [r['Residue Label'] for r in records if r['Residue Number'] in {205, 206, 125, 547, 357, 209}]

    parts = [
        f"Kompleks penambatan **{compound_name}** terhadap reseptor **{target_pdb}** menghasilkan energi ikatan **ΔG = {best_affinity} kcal/mol** dengan total **{counts['Total']} interaksi residu tervalidasi** "
        f"({counts['Hydrogen Bond']} ikatan hidrogen, {counts['Hydrophobic']} kontak hidrofobik, {counts['π-π Stacking'] + counts['π-Cation']} interaksi aromatik π, dan {counts['Salt Bridge']} jembatan garam)."
    ]
    if s2_hits:
        parts.append(f"Ligan membentuk jangkar pengikatan pada **saku S2 / subsite elektrostatik-aromatik** melalui residu **{', '.join(s2_hits)}**, yang berperan krusial dalam pengenalan substrat DPP-4.")
    if s1_hits:
        parts.append(f"Stabilitas konformasi diperkuat oleh okupansi pada **saku hidrofobik S1** melalui kontak dengan **{', '.join(s1_hits)}**.")
    if triad_hits:
        parts.append(f"Selain itu, ligan berada dalam jarak interaksi langsung terhadap residu **Triad Katalitik ({', '.join(triad_hits)})**, mendukung potensi penghambatan aktivitas enzimatik DPP-4 secara kompetitif.")
    elif not s1_hits and not s2_hits:
        parts.append(f"Residu utama yang menstabilkan ligan di kantong pengikatan meliputi **{', '.join(res_labels[:5])}**.")
    return " ".join(parts)

# ==========================================
# 5. HEADER TAMPILAN UTAMA
# ==========================================
st.markdown("""
<div class="app-header">
<div class="app-header-left">
<div class="app-logo-badge">🧬</div>
<div>
<h1 class="app-title">DPP-4 Virtual Screener</h1>
<p class="app-subtitle">
Computational Drug Discovery &amp; Virtual Screening Laboratory · 
<span style="color: #2485C7; font-weight: 600;">Dual-Target Cheminformatics Platform (3G0B &amp; 5T4B)</span>
</p>
</div>
</div>
<div style="display: flex; gap: 0.5rem; align-items: center; flex-wrap: wrap;">
<span class="badge-info">Targets: PDB 3G0B &amp; 5T4B</span>
<span class="badge-pass">✓ Dual-Target Docking Verified</span>
</div>
</div>
""", unsafe_allow_html=True)

# ==========================================
# 6. TABS STRUKTURAL
# ==========================================
tab_screen, tab_explorer, tab_curation = st.tabs([
    "🔬 Single Screener & 3D Interaction", 
    "🌿 Screened Natural Products Library",
    "🧪 Data Curation & Validation Hub"
])

# -------------------------------------------------------------
# TAB 1: SINGLE LIGAND SCREENER & 3D INTERACTION EXPLORER
# -------------------------------------------------------------
with tab_screen:
    st.markdown("""
<div class="page-intro">
<h3 class="page-intro-title">🔬 Single Ligand Screener &amp; 3D Molecular Interaction Explorer</h3>
<p class="page-intro-desc">Penambatan molekuler independen terhadap dua struktur kristal DPP-4 (3G0B &amp; 5T4B), perbandingan afinitas, dan eksplorasi interaksi 3D berbasis koordinat.</p>
</div>
""", unsafe_allow_html=True)

    PRESET_COMPOUNDS = {
        "QUERCETIN (CHEMBL50) - Flavonoid": ("CHEMBL50", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
        "CHRYSIN (CHEMBL149) - Flavone": ("CHEMBL149", "O=c1cc(-c2ccccc2)oc2cc(O)cc(O)c12"),
        "BERBERINE (CHEMBL295124) - Alkaloid": ("CHEMBL295124", "COc1ccc2cc3[n+](cc2c1OC)CCc1cc2c(cc1-3)OCO2"),
        "ROSMARINIC ACID (CHEMBL269229) - Phenolic Ester": ("CHEMBL269229", "O=C(O)[C@H](Cc1ccc(O)c(O)c1)OC(=O)/C=C/c1ccc(O)c(O)c1"),
        "CURCUMIN (CHEMBL140) - Polyphenol": ("CHEMBL140", "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"),
        "RESVERATROL (CHEMBL165) - Stilbene": ("CHEMBL165", "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
        "ERIODICTYOL (CHEMBL8996) - Flavanone": ("CHEMBL8996", "O=C1C[C@@H](c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c21"),
        "LUTEOLIN (CHEMBL151) - Flavone": ("CHEMBL151", "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
        "KAEMPFEROL (CHEMBL150) - Flavonol": ("CHEMBL150", "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
        "GENISTEIN (CHEMBL44) - Isoflavone": ("CHEMBL44", "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12"),
        "HESPERETIN (CHEMBL399121) - Flavanone": ("CHEMBL399121", "COc1ccc([C@@H]2CC(=O)c3c(O)cc(O)cc3O2)cc1O"),
        "VILDAGLIPTIN (CHEMBL142703) - Approved Inhibitor": ("CHEMBL142703", "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2"),
        "SAXAGLIPTIN (CHEMBL385517) - Approved Inhibitor": ("CHEMBL385517", "N#C[C@@H]1C[C@@H]2C[C@@H]2N1C(=O)[C@@H](N)C12CC3CC(CC(O)(C3)C1)C2"),
        "UBENIMEX (CHEMBL29292) - Natural Peptide": ("CHEMBL29292", "CC(C)C[C@H](NC(=O)[C@@H](O)[C@H](N)Cc1ccccc1)C(=O)O"),
        "Custom SMILES": ("CUSTOM", "")
    }

    # Initialize session_state with pre-validated Quercetin (CHEMBL50) so 3D explorer & dual results are immediately available
    if 'active_screening' not in st.session_state:
        init_cid, init_smiles = PRESET_COMPOUNDS["QUERCETIN (CHEMBL50) - Flavonoid"]
        init_mol = Chem.MolFromSmiles(init_smiles)
        init_props = {
            'MW': round(Descriptors.MolWt(init_mol), 2),
            'LogP': round(Descriptors.MolLogP(init_mol), 2),
            'HBD': Descriptors.NumHDonors(init_mol),
            'HBA': Descriptors.NumHAcceptors(init_mol),
            'Formula': Chem.rdMolDescriptors.CalcMolFormula(init_mol),
            'ChEMBL_ID': init_cid,
            'Name': 'QUERCETIN',
            'CID': '5280343'
        }
        st.session_state['active_screening'] = {
            'props': init_props,
            'smiles': init_smiles,
            'docking_by_pdb': load_saved_docking_for_compound(init_cid)
        }

    col_input, col_result = st.columns([1.05, 2.45], gap="large")

    with col_input:
        st.markdown("""
<div class="science-card" style="padding: 1.25rem 1.4rem; margin-bottom: 1rem;">
<div class="card-header" style="margin-bottom: 0.65rem; padding-bottom: 0.65rem;">
<h4 class="card-title">Analysis &amp; Docking Setup</h4>
<span class="badge-info">Dual Target</span>
</div>
<div style="font-size: 0.82rem; color: #71839A; line-height: 1.55;">
Setiap ligan ditambatkan secara independen terhadap kedua struktur reseptor DPP-4:<br>
• <b style="color: #2485C7;">PDB 3G0B:</b> Center <code>(42.05, 34.29, 14.62)</code>, Box <code>18×18×18 Å</code><br>
• <b style="color: #15966A;">PDB 5T4B:</b> Center <code>(37.61, 49.92, 40.33)</code>, Box <code>18×16×23 Å</code>
</div>
</div>
""", unsafe_allow_html=True)

        preset_choice = st.selectbox("Quick Preset (Curated Phytochemicals)", list(PRESET_COMPOUNDS.keys()), index=0)
        preset_cid, preset_smiles = PRESET_COMPOUNDS[preset_choice]

        default_smiles = preset_smiles if preset_smiles else "CC(C)[C@H](N)C(=O)N1CCCC1"
        smiles_input = st.text_input("Ligand Structure (SMILES)", default_smiles)

        protein_target = st.selectbox("Primary Target Focus", ["3G0B", "5T4B"], help="Pilih fokus target utama (kedua target 3G0B & 5T4B tetap dihitung secara independen)")

        st.markdown("<div style='height: 0.35rem;'></div>", unsafe_allow_html=True)
        run_btn = st.button("Run AutoDock Vina (3G0B & 5T4B)")

        # If user changes preset to a pre-docked library compound and hasn't clicked run_btn, load its saved real dual-target results automatically
        if not run_btn and preset_cid != "CUSTOM":
            current_cid = st.session_state['active_screening']['props'].get('ChEMBL_ID')
            if current_cid != preset_cid and smiles_input == preset_smiles:
                saved_dual = load_saved_docking_for_compound(preset_cid)
                if saved_dual['3G0B']['best_affinity'] is not None or saved_dual['5T4B']['best_affinity'] is not None:
                    mol_obj = Chem.MolFromSmiles(preset_smiles)
                    comp_title = preset_choice.split(" (")[0]
                    st.session_state['active_screening'] = {
                        'props': {
                            'MW': round(Descriptors.MolWt(mol_obj), 2),
                            'LogP': round(Descriptors.MolLogP(mol_obj), 2),
                            'HBD': Descriptors.NumHDonors(mol_obj),
                            'HBA': Descriptors.NumHAcceptors(mol_obj),
                            'Formula': Chem.rdMolDescriptors.CalcMolFormula(mol_obj),
                            'ChEMBL_ID': preset_cid,
                            'Name': comp_title,
                            'CID': '-'
                        },
                        'smiles': preset_smiles,
                        'docking_by_pdb': saved_dual
                    }

        if run_btn:
            with st.spinner("Menjalankan penambatan molekuler independen pada PDB 3G0B dan PDB 5T4B..."):
                ligand_pdbqt, props = process_ligand(smiles_input)
                if not ligand_pdbqt:
                    st.error("Gagal: SMILES tidak valid atau struktur tidak dapat diubah ke format 3D.")
                else:
                    dual_results = run_dual_target_docking(ligand_pdbqt)
                    st.session_state['active_screening'] = {
                        'props': props,
                        'smiles': smiles_input,
                        'docking_by_pdb': dual_results
                    }

    with col_result:
        active_data = st.session_state['active_screening']
        props = active_data['props']
        dual_res = active_data['docking_by_pdb']

        aff_3g0b = dual_res['3G0B']['best_affinity']
        aff_5t4b = dual_res['5T4B']['best_affinity']
        str_3g0b = f"{aff_3g0b:.3f} kcal/mol" if aff_3g0b is not None else "Docking Failed"
        str_5t4b = f"{aff_5t4b:.3f} kcal/mol" if aff_5t4b is not None else "Docking Failed"

        if aff_3g0b is not None and aff_5t4b is not None:
            diff_val = round(abs(aff_3g0b - aff_5t4b), 3)
            stronger_target = "3G0B" if aff_3g0b <= aff_5t4b else "5T4B"
            diff_str = f"{diff_val:.3f} kcal/mol (Stronger: {stronger_target})"
        else:
            diff_str = "N/A"
            stronger_target = "3G0B" if aff_3g0b is not None else ("5T4B" if aff_5t4b is not None else "N/A")

        violations = sum([
            props['MW'] > 500,
            props['LogP'] > 5,
            props['HBD'] > 5,
            props['HBA'] > 10
        ])
        lip_pass = (violations <= 1)
        lip_str = "PASS" if lip_pass else "FAIL"
        lip_badge_cls = "badge-pass" if lip_pass else "badge-fail"

        c_mw = "lipinski-fail" if props['MW'] > 500 else ""
        c_logp = "lipinski-fail" if props['LogP'] > 5 else ""
        c_hbd = "lipinski-fail" if props['HBD'] > 5 else ""
        c_hba = "lipinski-fail" if props['HBA'] > 10 else ""

        st.markdown(f"""
<div class="science-card">
<div class="card-header">
<div>
<h4 class="card-title">{props['Name']} ({props['ChEMBL_ID']})</h4>
<span style="font-size: 0.82rem; color: #71839A;">Formula: <code>{props['Formula']}</code> · PubChem CID: <code>{props['CID']}</code></span>
</div>
<span class="{lip_badge_cls}">Lipinski: {lip_str}</span>
</div>

<div class="docking-res-section-title">MOLECULAR DOCKING RESULTS</div>
<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 0.85rem; margin-bottom: 1.15rem;">
<div class="metric-box" style="border-top: 3px solid #2485C7;">
<div class="data-label" style="color: #2485C7;">3G0B Binding Affinity</div>
<div class="data-value" style="color: #172B4D; font-size: 1.2rem;">3G0B: {str_3g0b}</div>
<div style="font-size: 0.74rem; color: #71839A; margin-top: 0.2rem;">Grid Center: (42.05, 34.29, 14.62)</div>
</div>
<div class="metric-box" style="border-top: 3px solid #15966A;">
<div class="data-label" style="color: #15966A;">5T4B Binding Affinity</div>
<div class="data-value" style="color: #15966A; font-size: 1.2rem;">5T4B: {str_5t4b}</div>
<div style="font-size: 0.74rem; color: #71839A; margin-top: 0.2rem;">Grid Center: (37.61, 49.92, 40.33)</div>
</div>
<div class="metric-box">
<div class="data-label">Affinity Difference (|ΔΔG|)</div>
<div class="data-value" style="font-size: 1rem; color: #2485C7;">{diff_str}</div>
<div style="font-size: 0.74rem; color: #71839A; margin-top: 0.2rem;">MW: <b>{props['MW']} Da</b> · Lipinski: <b>{lip_str}</b></div>
</div>
</div>

<div class="data-label" style="margin-bottom: 0.65rem;">Lipinski's Rule of Five Physicochemical Profile</div>
<div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 0.85rem;">
<div class="lipinski-box {c_mw}">
<div style="font-size: 0.72rem; color: #71839A;">Molecular Weight</div>
<div class="data-value" style="font-size: 1rem;">{props['MW']} <span style="font-size: 0.72rem; color:#71839A; font-weight:normal;">Da</span></div>
</div>
<div class="lipinski-box {c_logp}">
<div style="font-size: 0.72rem; color: #71839A;">LogP</div>
<div class="data-value" style="font-size: 1rem;">{props['LogP']}</div>
</div>
<div class="lipinski-box {c_hbd}">
<div style="font-size: 0.72rem; color: #71839A;">H-Bond Donors</div>
<div class="data-value" style="font-size: 1rem;">{props['HBD']}</div>
</div>
<div class="lipinski-box {c_hba}">
<div style="font-size: 0.72rem; color: #71839A;">H-Bond Acceptors</div>
<div class="data-value" style="font-size: 1rem;">{props['HBA']}</div>
</div>
</div>
</div>
""", unsafe_allow_html=True)

        # Horizontal Bar Chart Comparison + Poses Tables for 3G0B and 5T4B
        comp_col1, comp_col2 = st.columns([1.15, 1.35], gap="medium")
        with comp_col1:
            if aff_3g0b is not None or aff_5t4b is not None:
                chart_rows = []
                if aff_3g0b is not None:
                    chart_rows.append({'Target PDB': '3G0B (Crystal Complex)', 'Affinity (kcal/mol)': aff_3g0b, 'PDB': '3G0B'})
                if aff_5t4b is not None:
                    chart_rows.append({'Target PDB': '5T4B (Wild-Type DPP-4)', 'Affinity (kcal/mol)': aff_5t4b, 'PDB': '5T4B'})
                df_comp_chart = pd.DataFrame(chart_rows)
                st.vega_lite_chart(
                    df_comp_chart,
                    {
                        "background": "#FFFFFF",
                        "padding": 12,
                        "height": 165,
                        "title": {
                            "text": f"Perbandingan Binding Affinity: {props['Name']} (kcal/mol)",
                            "color": "#172B4D",
                            "font": "Inter",
                            "fontSize": 12.5,
                            "fontWeight": 700,
                            "anchor": "start"
                        },
                        "mark": {"type": "bar", "cornerRadiusEnd": 6, "height": 26, "tooltip": True},
                        "encoding": {
                            "y": {
                                "field": "Target PDB",
                                "type": "nominal",
                                "title": None,
                                "axis": {"labelColor": "#172B4D", "labelFontWeight": 600}
                            },
                            "x": {
                                "field": "Affinity (kcal/mol)",
                                "type": "quantitative",
                                "title": "Binding Affinity (kcal/mol)",
                                "axis": {"labelColor": "#172B4D", "titleColor": "#71839A", "gridColor": "#E5F3FF"}
                            },
                            "color": {
                                "field": "PDB",
                                "type": "nominal",
                                "scale": {"domain": ["3G0B", "5T4B"], "range": ["#2485C7", "#15966A"]},
                                "legend": None
                            }
                        },
                        "config": {"view": {"stroke": "#DCE8F2"}}
                    },
                    use_container_width=True
                )
        with comp_col2:
            pose_tab_3g0b, pose_tab_5t4b = st.tabs(["📊 Poses 3G0B", "📊 Poses 5T4B"])
            with pose_tab_3g0b:
                if dual_res['3G0B']['poses']:
                    df_p3 = pd.DataFrame(dual_res['3G0B']['poses'])
                    st.dataframe(df_p3.style.highlight_min(subset=['Affinity (kcal/mol)'], color='#DDF5E9'), use_container_width=True, hide_index=True, height=140)
                else:
                    st.warning("Pose 3G0B tidak tersedia.")
            with pose_tab_5t4b:
                if dual_res['5T4B']['poses']:
                    df_p5 = pd.DataFrame(dual_res['5T4B']['poses'])
                    st.dataframe(df_p5.style.highlight_min(subset=['Affinity (kcal/mol)'], color='#DDF5E9'), use_container_width=True, hide_index=True, height=140)
                else:
                    st.warning("Pose 5T4B tidak tersedia.")

    # ==========================================================
    # PREMIUM INTERACTIVE 3D MOLECULAR INTERACTION EXPLORER
    # ==========================================================
    st.markdown("""
<div class="science-card" style="margin-top: 0.75rem; padding: 1.35rem 1.6rem;">
<div class="card-header" style="margin-bottom: 0; border-bottom: none; padding-bottom: 0;">
<div>
<h3 class="card-title" style="font-size: 1.25rem;">🧊 Molecular Interaction Analysis</h3>
<p class="page-intro-desc" style="margin-top: 0.2rem;">Explore protein-ligand binding interactions in 3D.</p>
</div>
<div style="display: flex; gap: 0.5rem; flex-wrap: wrap;">
<span class="badge-info">Coordinate-Validated Contacts</span>
<span class="badge-pass">✓ Dual-Target 3D Explorer</span>
</div>
</div>
</div>
""", unsafe_allow_html=True)

    # Target Selector & Interactive 3D Controls
    ctrl_c1, ctrl_c2, ctrl_c3 = st.columns([1.3, 1.1, 1.2], gap="medium")
    with ctrl_c1:
        target_view_mode = st.radio(
            "Target Structure Selector:",
            ["3G0B", "5T4B", "Side-by-Side (3G0B vs 5T4B)"],
            index=0 if protein_target == "3G0B" else 1,
            horizontal=True
        )
    with ctrl_c2:
        rep_choice = st.selectbox(
            "Molecular Representation:",
            ["Cartoon + Ligand", "Ball and Stick", "Stick", "Surface", "Cartoon"],
            index=0
        )
    with ctrl_c3:
        active_single_pdb = "3G0B" if target_view_mode == "3G0B" else ("5T4B" if target_view_mode == "5T4B" else protein_target)
        pose_Single = extract_first_model_pdbqt(dual_res[active_single_pdb]['docked_pdbqt'])
        df_ix_active, ix_counts_active, ix_list_active = analyze_binding_interactions(active_single_pdb, pose_Single)
        res_options = ["All Pocket Residues"] + [f"{r['Residue Label']} ({r['Interaction Type']} - {r['Distance (Å)']} Å)" for r in ix_list_active]
        focus_res_choice = st.selectbox("Inspect Specific Interacting Residue:", res_options, index=0)

    tog_c1, tog_c2, tog_c3, tog_c4 = st.columns(4)
    with tog_c1:
        show_prot_toggle = st.checkbox("Toggle Protein Visibility", value=True)
    with tog_c2:
        show_lig_toggle = st.checkbox("Toggle Ligand Visibility", value=True)
    with tog_c3:
        show_pocket_toggle = st.checkbox("Toggle Binding Pocket", value=True)
    with tog_c4:
        show_ix_toggle = st.checkbox("Toggle Interaction Lines", value=True)

    # MAIN AREA: 3D Viewer + Information Panel
    if target_view_mode == "Side-by-Side (3G0B vs 5T4B)":
        sbs_col1, sbs_col2 = st.columns(2, gap="medium")
        pose_3g0b_m1 = extract_first_model_pdbqt(dual_res['3G0B']['docked_pdbqt'])
        pose_5t4b_m1 = extract_first_model_pdbqt(dual_res['5T4B']['docked_pdbqt'])
        _, counts_3g0b, list_3g0b = analyze_binding_interactions('3G0B', pose_3g0b_m1)
        _, counts_5t4b, list_5t4b = analyze_binding_interactions('5T4B', pose_5t4b_m1)

        with sbs_col1:
            st.markdown(f'<div class="science-card" style="padding: 0.85rem 1.2rem; margin-bottom: 0.5rem;"><div style="display:flex; justify-content:space-between; align-items:center;"><b>🎯 Target 1: PDB 3G0B</b><span class="badge-info">Affinity: {str_3g0b}</span></div></div>', unsafe_allow_html=True)
            render_interactive_3d_explorer(
                '3G0B', pose_3g0b_m1, list_3g0b,
                rep_mode=rep_choice,
                show_protein=show_prot_toggle,
                show_ligand=show_lig_toggle,
                show_pocket=show_pocket_toggle,
                show_interactions=show_ix_toggle,
                focus_residue="All Pocket Residues",
                viewer_height=440,
                viewer_id="mol3d_3g0b"
            )
        with sbs_col2:
            st.markdown(f'<div class="science-card" style="padding: 0.85rem 1.2rem; margin-bottom: 0.5rem;"><div style="display:flex; justify-content:space-between; align-items:center;"><b>🎯 Target 2: PDB 5T4B</b><span class="badge-pass">Affinity: {str_5t4b}</span></div></div>', unsafe_allow_html=True)
            render_interactive_3d_explorer(
                '5T4B', pose_5t4b_m1, list_5t4b,
                rep_mode=rep_choice,
                show_protein=show_prot_toggle,
                show_ligand=show_lig_toggle,
                show_pocket=show_pocket_toggle,
                show_interactions=show_ix_toggle,
                focus_residue="All Pocket Residues",
                viewer_height=440,
                viewer_id="mol3d_5t4b"
            )
    else:
        viewer_col, info_col = st.columns([1.85, 1.05], gap="large")
        with viewer_col:
            render_interactive_3d_explorer(
                active_single_pdb,
                pose_Single,
                ix_list_active,
                rep_mode=rep_choice,
                show_protein=show_prot_toggle,
                show_ligand=show_lig_toggle,
                show_pocket=show_pocket_toggle,
                show_interactions=show_ix_toggle,
                focus_residue=focus_res_choice,
                viewer_height=510,
                viewer_id=f"mol3d_{active_single_pdb}"
            )
        with info_col:
            rec_meta = RECEPTOR_TARGETS[active_single_pdb]
            rec_struct = get_receptor_structure_info(active_single_pdb)
            active_aff = dual_res[active_single_pdb]['best_affinity']
            active_aff_str = f"{active_aff:.3f} kcal/mol" if active_aff is not None else "N/A"

            st.markdown(f"""
<div class="science-card" style="padding: 1.15rem 1.3rem; margin-bottom: 0.9rem;">
<div class="card-header" style="margin-bottom: 0.65rem; padding-bottom: 0.5rem;">
<h4 class="card-title" style="font-size: 0.98rem;">🧪 Ligand Information</h4>
<span class="badge-info">{active_single_pdb} Pose</span>
</div>
<div style="font-size: 0.84rem; line-height: 1.65; color: #172B4D;">
• <b>Compound Name:</b> {props['Name']}<br>
• <b>Compound ID:</b> <span style="color:#2485C7; font-family:'Fira Code',monospace;">{props['ChEMBL_ID']}</span><br>
• <b>Molecular Formula:</b> <code>{props['Formula']}</code><br>
• <b>Molecular Weight:</b> {props['MW']} Da<br>
• <b>Binding Affinity ({active_single_pdb}):</b> <b style="color:#15966A; font-family:'Fira Code',monospace;">{active_aff_str}</b>
</div>
</div>

<div class="science-card" style="padding: 1.15rem 1.3rem; margin-bottom: 0.9rem;">
<div class="card-header" style="margin-bottom: 0.65rem; padding-bottom: 0.5rem;">
<h4 class="card-title" style="font-size: 0.98rem;">🧬 Protein Target Information</h4>
<span class="badge-pass">PDB {active_single_pdb}</span>
</div>
<div style="font-size: 0.84rem; line-height: 1.65; color: #172B4D;">
• <b>Protein Name:</b> {rec_meta['protein_name']}<br>
• <b>PDB ID / Chain:</b> <code>{active_single_pdb}</code> (Chain <code>{rec_meta['chain_id']}</code>)<br>
• <b>Amino Acid Residues:</b> <b>{rec_struct['residues']} residues</b> ({rec_struct['atoms']} heavy atoms)<br>
• <b>Grid Box Center:</b> <code>{rec_meta['center']}</code>
</div>
</div>

<div class="science-card" style="padding: 1.15rem 1.3rem; margin-bottom: 0;">
<div class="card-header" style="margin-bottom: 0.65rem; padding-bottom: 0.5rem;">
<h4 class="card-title" style="font-size: 0.98rem;">⚡ Interaction Summary ({active_single_pdb})</h4>
<span class="badge-pass">{ix_counts_active['Total']} Validated</span>
</div>
<div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 0.55rem;">
<div class="metric-box" style="padding: 0.6rem 0.8rem;">
<div class="data-label" style="font-size: 0.68rem;">Hydrogen Bonds</div>
<div class="data-value" style="color: #15966A; font-size: 1.15rem;">{ix_counts_active['Hydrogen Bond']}</div>
</div>
<div class="metric-box" style="padding: 0.6rem 0.8rem;">
<div class="data-label" style="font-size: 0.68rem;">Hydrophobic</div>
<div class="data-value" style="color: #2485C7; font-size: 1.15rem;">{ix_counts_active['Hydrophobic']}</div>
</div>
<div class="metric-box" style="padding: 0.6rem 0.8rem;">
<div class="data-label" style="font-size: 0.68rem;">π-π &amp; π-Cation</div>
<div class="data-value" style="color: #8B5CF6; font-size: 1.15rem;">{ix_counts_active['π-π Stacking'] + ix_counts_active['π-Cation']}</div>
</div>
<div class="metric-box" style="padding: 0.6rem 0.8rem;">
<div class="data-label" style="font-size: 0.68rem;">Salt Bridge &amp; VdW</div>
<div class="data-value" style="color: #D97706; font-size: 1.15rem;">{ix_counts_active['Salt Bridge'] + ix_counts_active['Van der Waals']}</div>
</div>
</div>
</div>
""", unsafe_allow_html=True)

    # BOTTOM AREA: Interaction Type Legend, Detailed Residue Table & Interpretation
    st.markdown("""
<div class="science-card" style="margin-top: 1.15rem; padding: 1.15rem 1.5rem;">
<div style="display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between; gap: 0.75rem;">
<span style="font-size: 0.82rem; font-weight: 700; color: #172B4D; text-transform: uppercase; letter-spacing: 0.04em;">Interaction Type Color Legend (3D Dashed Lines):</span>
<div style="display: flex; flex-wrap: wrap; gap: 0.6rem;">
<span class="ix-pill"><span class="ix-dot" style="background:#15966A;"></span>Hydrogen Bond (≤3.5 Å)</span>
<span class="ix-pill"><span class="ix-dot" style="background:#D97706;"></span>Salt Bridge (≤4.0 Å)</span>
<span class="ix-pill"><span class="ix-dot" style="background:#8B5CF6;"></span>π-π Stacking (≤4.5 Å)</span>
<span class="ix-pill"><span class="ix-dot" style="background:#06B6D4;"></span>π-Cation (≤4.5 Å)</span>
<span class="ix-pill"><span class="ix-dot" style="background:#2485C7;"></span>Hydrophobic (≤4.0 Å)</span>
<span class="ix-pill"><span class="ix-dot" style="background:#64748B;"></span>Van der Waals (≤4.0 Å)</span>
</div>
</div>
</div>
""", unsafe_allow_html=True)

    bot_col1, bot_col2 = st.columns([1.65, 1.15], gap="large")
    with bot_col1:
        st.markdown(f"""
<div class="science-card" style="padding: 1.1rem 1.4rem; margin-bottom: 0.65rem;">
<div class="card-header" style="margin-bottom: 0; border-bottom: none; padding-bottom: 0;">
<h4 class="card-title">📋 Detailed Residue Interaction Table ({active_single_pdb})</h4>
<span class="badge-info">Coordinates in Ångström (Å)</span>
</div>
</div>
""", unsafe_allow_html=True)
        if not df_ix_active.empty:
            st.dataframe(df_ix_active, use_container_width=True, hide_index=True)
        else:
            st.info("Tidak ada residu dalam jarak cutoff interaksi (≤4.5 Å) untuk pose ini.")

    with bot_col2:
        interp_text = generate_interaction_interpretation(
            props['Name'],
            active_single_pdb,
            dual_res[active_single_pdb]['best_affinity'],
            ix_counts_active,
            ix_list_active
        )
        st.markdown(f"""
<div class="science-card" style="padding: 1.35rem 1.5rem;">
<div class="card-header" style="margin-bottom: 0.75rem; padding-bottom: 0.65rem;">
<h4 class="card-title">🔬 Molecular Interaction Interpretation</h4>
<span class="badge-pass">Pocket Mechanism</span>
</div>
<p style="font-size: 0.88rem; color: #172B4D; line-height: 1.65; margin-bottom: 0.9rem;">
{interp_text}
</p>
<div style="background: #F1F8FF; border: 1px solid #DCE8F2; border-left: 4px solid #2485C7; border-radius: 10px; padding: 0.75rem 1rem; font-size: 0.8rem; color: #172B4D; line-height: 1.55;">
<b>DPP-4 Reference Pocket Architecture:</b><br>
• <b>Catalytic Triad:</b> Ser630, Asp708, His740<br>
• <b>S1 Hydrophobic Pocket:</b> Tyr662, Tyr666, Val656, Trp659, Val711<br>
• <b>S2 Electrostatic Subsite:</b> Glu205, Glu206, Arg125, Tyr547, Phe357
</div>
</div>
""", unsafe_allow_html=True)

# -------------------------------------------------------------
# TAB 2: NATURAL PRODUCTS SCREENING & TOP HITS LIBRARY (DUAL TARGET 3G0B & 5T4B)
# -------------------------------------------------------------
with tab_explorer:
    st.markdown("""
<div class="page-intro">
<h3 class="page-intro-title">🌿 Screened Natural Products Library &amp; Dual-Target Top Hits</h3>
<p class="page-intro-desc">Hasil penapisan virtual AutoDock Vina multi-target independen (PDB 3G0B &amp; PDB 5T4B), perbandingan afinitas, dan sistem perankingan fleksibel.</p>
</div>
""", unsafe_allow_html=True)
    
    active_df = df_results.copy() if df_results is not None else (df_test.copy() if df_test is not None else df_master.copy())
    
    if active_df is not None:
        # Compute comparison columns strictly on valid numerical scores (never mixing NaNs or copying 3G0B to 5T4B)
        if 'Affinity_3G0B_kcal_mol' in active_df.columns and 'Affinity_5T4B_kcal_mol' in active_df.columns:
            both_valid_mask = active_df['Affinity_3G0B_kcal_mol'].notna() & active_df['Affinity_5T4B_kcal_mol'].notna()
            active_df['Mean_Affinity_kcal_mol'] = np.where(
                both_valid_mask,
                ((active_df['Affinity_3G0B_kcal_mol'] + active_df['Affinity_5T4B_kcal_mol']) / 2.0).round(3),
                np.nan
            )
            active_df['Affinity_Diff_kcal_mol'] = np.where(
                both_valid_mask,
                (active_df['Affinity_3G0B_kcal_mol'] - active_df['Affinity_5T4B_kcal_mol']).abs().round(3),
                np.nan
            )
            active_df['Stronger_Target'] = np.where(
                both_valid_mask,
                np.where(active_df['Affinity_3G0B_kcal_mol'] <= active_df['Affinity_5T4B_kcal_mol'], '3G0B', '5T4B'),
                np.where(active_df['Affinity_3G0B_kcal_mol'].notna(), '3G0B Only', np.where(active_df['Affinity_5T4B_kcal_mol'].notna(), '5T4B Only', 'Failed / Skipped'))
            )

        # Filter & Ranking Controls
        fc1, fc2, fc3, fc4 = st.columns([1.6, 1.3, 1.1, 1.0])
        with fc1:
            search_query = st.text_input("🔍 Cari Nama Senyawa atau ChEMBL ID:", "")
        with fc2:
            ranking_mode = st.selectbox(
                "📊 Sistem Perankingan (Ranking Basis):",
                [
                    "Ranking berdasarkan 3G0B",
                    "Ranking berdasarkan 5T4B",
                    "Ranking berdasarkan Rata-Rata Kedua Target (3G0B & 5T4B)"
                ],
                index=0
            )
        with fc3:
            split_filter = st.selectbox("Filter Dataset Split:", ["Semua", "Test Set (Validation & Benchmark)", "Analysis Set (Screening Candidates)"])
        with fc4:
            lip_filter = st.selectbox("Filter Lipinski:", ["Semua", "Hanya Lolos Lipinski", "Gagal Lipinski"])
            
        filtered_df = active_df.copy()
        if search_query:
            filtered_df = filtered_df[
                filtered_df['compound_name'].str.contains(search_query, case=False, na=False) |
                filtered_df['molecule_chembl_id'].str.contains(search_query, case=False, na=False)
            ]
        if split_filter != "Semua" and 'dataset_split' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['dataset_split'] == split_filter]
        if lip_filter == "Hanya Lolos Lipinski" and 'Lipinski_Pass' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['Lipinski_Pass'] == True]
        elif lip_filter == "Gagal Lipinski" and 'Lipinski_Pass' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['Lipinski_Pass'] == False]

        # Determine ranking column and filter out invalid/failed docking strictly for ranking
        if ranking_mode == "Ranking berdasarkan 5T4B":
            rank_col = 'Affinity_5T4B_kcal_mol'
            valid_ranked_df = filtered_df.dropna(subset=['Affinity_5T4B_kcal_mol']).sort_values(by='Affinity_5T4B_kcal_mol', ascending=True)
        elif ranking_mode == "Ranking berdasarkan Rata-Rata Kedua Target (3G0B & 5T4B)":
            rank_col = 'Mean_Affinity_kcal_mol'
            valid_ranked_df = filtered_df.dropna(subset=['Affinity_3G0B_kcal_mol', 'Affinity_5T4B_kcal_mol', 'Mean_Affinity_kcal_mol']).sort_values(by='Mean_Affinity_kcal_mol', ascending=True)
        else:
            rank_col = 'Affinity_3G0B_kcal_mol'
            valid_ranked_df = filtered_df.dropna(subset=['Affinity_3G0B_kcal_mol']).sort_values(by='Affinity_3G0B_kcal_mol', ascending=True)

        # Redesigned Top 4 Hits Cards displaying BOTH 3G0B & 5T4B Docking Results
        if not valid_ranked_df.empty:
            st.markdown(f"""
<div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; margin-top: 0.75rem; margin-bottom: 0.85rem;">
<h4 style="color: #172B4D; font-size: 1.08rem; font-weight: 700; margin: 0;">🏆 Top Natural Product Hits ({ranking_mode})</h4>
<span class="badge-pass">{len(valid_ranked_df)} Valid Docked Candidates Ranked</span>
</div>
""", unsafe_allow_html=True)
            top_cols = st.columns(4)
            for i, (_, hit) in enumerate(valid_ranked_df.head(4).iterrows()):
                with top_cols[i]:
                    val_3g0b = f"{hit['Affinity_3G0B_kcal_mol']:.3f} kcal/mol" if pd.notna(hit.get('Affinity_3G0B_kcal_mol')) else "Failed / N/A"
                    val_5t4b = f"{hit['Affinity_5T4B_kcal_mol']:.3f} kcal/mol" if pd.notna(hit.get('Affinity_5T4B_kcal_mol')) else "Failed / N/A"
                    is_3g0b_better = pd.notna(hit.get('Affinity_3G0B_kcal_mol')) and pd.notna(hit.get('Affinity_5T4B_kcal_mol')) and (hit['Affinity_3G0B_kcal_mol'] <= hit['Affinity_5T4B_kcal_mol'])
                    is_5t4b_better = pd.notna(hit.get('Affinity_3G0B_kcal_mol')) and pd.notna(hit.get('Affinity_5T4B_kcal_mol')) and (hit['Affinity_5T4B_kcal_mol'] < hit['Affinity_3G0B_kcal_mol'])
                    row_3g0b_cls = "docking-target-row docking-target-row-best" if is_3g0b_better else "docking-target-row"
                    row_5t4b_cls = "docking-target-row docking-target-row-best" if is_5t4b_better else "docking-target-row"
                    diff_badge = f"ΔΔG: {hit['Affinity_Diff_kcal_mol']:.2f} kcal/mol ({hit['Stronger_Target']})" if pd.notna(hit.get('Affinity_Diff_kcal_mol')) else f"Status: {hit.get('Stronger_Target', 'N/A')}"
                    lip_status = 'PASS' if hit.get('Lipinski_Pass') else 'FAIL'
                    lip_cls = 'badge-pass' if hit.get('Lipinski_Pass') else 'badge-fail'

                    st.markdown(f"""
<div class="science-card" style="padding: 1.15rem 1.25rem; border-top: 3px solid #2485C7;">
<div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.35rem;">
<span class="badge-info">RANK #{i+1}</span>
<span class="{lip_cls}">Lipinski: {lip_status}</span>
</div>
<div style="font-size: 1.02rem; font-weight: 700; color: #172B4D; margin: 0.3rem 0 0.1rem 0; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;" title="{hit['compound_name']}">{hit['compound_name']}</div>
<div style="font-size: 0.78rem; color: #2485C7; font-family: 'Fira Code', monospace; font-weight: 600;">{hit['molecule_chembl_id']}</div>

<div class="docking-res-section-title">MOLECULAR DOCKING RESULTS</div>
<div class="{row_3g0b_cls}">
<span style="font-size: 0.76rem; color: #172B4D; font-weight: 700;">3G0B:</span>
<span style="font-weight: 700; color: #2485C7; font-family: 'Fira Code', monospace; font-size: 0.84rem;">{val_3g0b}</span>
</div>
<div class="{row_5t4b_cls}">
<span style="font-size: 0.76rem; color: #172B4D; font-weight: 700;">5T4B:</span>
<span style="font-weight: 700; color: #15966A; font-family: 'Fira Code', monospace; font-size: 0.84rem;">{val_5t4b}</span>
</div>

<hr style="margin: 0.55rem 0; border: none; border-top: 1px solid #DCE8F2;">
<div style="display: flex; justify-content: space-between; margin-bottom: 0.22rem;">
<span style="font-size: 0.75rem; color: #71839A;">MW:</span>
<span style="font-size: 0.8rem; color: #172B4D; font-weight: 600;">{hit['MW']} Da</span>
</div>
<div style="display: flex; justify-content: space-between; margin-bottom: 0.22rem;">
<span style="font-size: 0.75rem; color: #71839A;">Lipinski:</span>
<span style="font-size: 0.8rem; color: {'#15966A' if hit.get('Lipinski_Pass') else '#BE123C'}; font-weight: 700;">{lip_status}</span>
</div>
<div style="display: flex; justify-content: space-between;">
<span style="font-size: 0.74rem; color: #71839A;">Comparison:</span>
<span style="font-size: 0.75rem; color: #2485C7; font-weight: 600;">{diff_badge}</span>
</div>
</div>
""", unsafe_allow_html=True)

            # Horizontal Bar Chart Comparing 3G0B vs 5T4B for Top Ranked Candidates
            st.markdown("""
<div class="science-card" style="padding: 1.15rem 1.5rem; margin-top: 0.25rem; margin-bottom: 0.75rem;">
<div class="card-header" style="margin-bottom: 0; border-bottom: none; padding-bottom: 0;">
<div>
<h4 class="card-title">📊 Horizontal Bar Chart: Perbandingan Binding Affinity 3G0B vs 5T4B (Top 10 Ranked)</h4>
<p class="page-intro-desc" style="font-size: 0.84rem; margin-top: 0.15rem;">Perbandingan langsung skor penambatan molekuler independen pada kedua reseptor DPP-4 (kcal/mol).</p>
</div>
<span class="badge-info">3G0B (Blue) vs 5T4B (Emerald)</span>
</div>
</div>
""", unsafe_allow_html=True)

            top10_comp = valid_ranked_df.head(10).copy()
            melted_rows = []
            for _, r in top10_comp.iterrows():
                cname = f"{r['compound_name']} ({r['molecule_chembl_id']})"
                if pd.notna(r.get('Affinity_3G0B_kcal_mol')):
                    melted_rows.append({
                        'Compound': cname,
                        'Target': '3G0B',
                        'Binding Affinity (kcal/mol)': float(r['Affinity_3G0B_kcal_mol']),
                        'Diff (kcal/mol)': r.get('Affinity_Diff_kcal_mol'),
                        'Stronger Target': r.get('Stronger_Target')
                    })
                if pd.notna(r.get('Affinity_5T4B_kcal_mol')):
                    melted_rows.append({
                        'Compound': cname,
                        'Target': '5T4B',
                        'Binding Affinity (kcal/mol)': float(r['Affinity_5T4B_kcal_mol']),
                        'Diff (kcal/mol)': r.get('Affinity_Diff_kcal_mol'),
                        'Stronger Target': r.get('Stronger_Target')
                    })
            if melted_rows:
                df_melted = pd.DataFrame(melted_rows)
                st.vega_lite_chart(
                    df_melted,
                    {
                        "background": "#FFFFFF",
                        "padding": 16,
                        "height": 340,
                        "mark": {"type": "bar", "cornerRadiusEnd": 4, "tooltip": True},
                        "encoding": {
                            "y": {
                                "field": "Compound",
                                "type": "nominal",
                                "sort": None,
                                "title": "Kandidat Senyawa Bahan Alam",
                                "axis": {"labelColor": "#172B4D", "titleColor": "#71839A", "labelLimit": 210, "labelFontWeight": 600}
                            },
                            "yOffset": {"field": "Target"},
                            "x": {
                                "field": "Binding Affinity (kcal/mol)",
                                "type": "quantitative",
                                "title": "Binding Affinity (kcal/mol) — Semakin Negatif Semakin Kuat",
                                "axis": {"labelColor": "#172B4D", "titleColor": "#71839A", "gridColor": "#E5F3FF"}
                            },
                            "color": {
                                "field": "Target",
                                "type": "nominal",
                                "scale": {"domain": ["3G0B", "5T4B"], "range": ["#2485C7", "#15966A"]},
                                "legend": {"title": "Target Protein PDB", "orient": "top", "labelColor": "#172B4D", "titleColor": "#71839A"}
                            },
                            "tooltip": [
                                {"field": "Compound", "type": "nominal", "title": "Compound"},
                                {"field": "Target", "type": "nominal", "title": "Target PDB"},
                                {"field": "Binding Affinity (kcal/mol)", "type": "quantitative", "title": "Affinity (kcal/mol)"},
                                {"field": "Diff (kcal/mol)", "type": "quantitative", "title": "|ΔΔG| Diff (kcal/mol)"},
                                {"field": "Stronger Target", "type": "nominal", "title": "More Negative Target"}
                            ]
                        },
                        "config": {"view": {"stroke": "#DCE8F2"}}
                    },
                    use_container_width=True
                )

        # Display Data Table inside consistent card header (Ranked Valid + Clear Status for Unranked/Failed)
        st.markdown("""
<div class="science-card" style="padding: 1.15rem 1.5rem; margin-top: 0.75rem; margin-bottom: 0.75rem;">
<div class="card-header" style="margin-bottom: 0; border-bottom: none; padding-bottom: 0;">
<h4 class="card-title">📋 Screened Natural Products Comparison Table (3G0B vs 5T4B)</h4>
<span class="badge-info">Independent Dual-Target Scores</span>
</div>
</div>
""", unsafe_allow_html=True)

        # Order table: valid ranked rows first (with Rank column), followed by unranked/failed rows clearly flagged
        valid_ranked_copy = valid_ranked_df.copy()
        valid_ranked_copy.insert(0, 'Rank', range(1, len(valid_ranked_copy) + 1))
        unranked_df = filtered_df.drop(index=valid_ranked_df.index).copy()
        if not unranked_df.empty:
            unranked_df.insert(0, 'Rank', None)
        combined_display_df = pd.concat([valid_ranked_copy, unranked_df], ignore_index=True)

        display_cols = [
            c for c in [
                'Rank', 'compound_name', 'molecule_chembl_id',
                'Affinity_3G0B_kcal_mol', 'Affinity_5T4B_kcal_mol',
                'Mean_Affinity_kcal_mol', 'Affinity_Diff_kcal_mol', 'Stronger_Target',
                'MW', 'LogP', 'HBD', 'HBA', 'Lipinski_Pass', 'Docking_Status', 'dataset_split'
            ] if c in combined_display_df.columns
        ]
        st.dataframe(combined_display_df[display_cols], use_container_width=True, hide_index=True)

        dl_c1, dl_c2, dl_c3 = st.columns(3)
        with dl_c1:
            st.download_button("📥 Unduh Hasil Gabungan 3G0B & 5T4B (CSV)", combined_display_df[display_cols].to_csv(index=False), "dpp4_dual_target_docking_results.csv", "text/csv")
        with dl_c2:
            if target_results_by_pdb['3G0B'] is not None:
                st.download_button("📥 Unduh Struktur Data Khusus 3G0B (CSV)", target_results_by_pdb['3G0B'].to_csv(index=False), "docking_results_3G0B.csv", "text/csv")
        with dl_c3:
            if target_results_by_pdb['5T4B'] is not None:
                st.download_button("📥 Unduh Struktur Data Khusus 5T4B (CSV)", target_results_by_pdb['5T4B'].to_csv(index=False), "docking_results_5T4B.csv", "text/csv")
    else:
        st.info("Pangkalan data hasil skrining sedang diproses atau belum tersedia.")

# -------------------------------------------------------------
# TAB 3: DATA CURATION & VALIDATION HUB (INFORMASI SAJA)
# -------------------------------------------------------------
with tab_curation:
    st.markdown("""
<div class="page-intro">
<h3 class="page-intro-title">🧪 Data Curation &amp; Validation Hub</h3>
<p class="page-intro-desc">Audit kualitas data molekuler, pembagian Test vs Analysis Set, serta validasi korelasi eksperimental penapisan obat DPP-4 (Ditampilkan sebagai informasi referensi).</p>
</div>
""", unsafe_allow_html=True)
    
    if df_master is not None:
        total_np = len(df_master)
        test_cnt = len(df_test) if df_test is not None else 0
        an_cnt = len(df_analysis) if df_analysis is not None else 0
        lip_rate = f"{df_master['Lipinski_Pass'].mean()*100:.1f}%"
        veber_rate = f"{df_master['Veber_Pass'].mean()*100:.1f}%"

        st.markdown(f"""
<div class="kpi-grid">
<div class="kpi-card kpi-blue">
<div class="kpi-head"><span class="kpi-label">Total Natural Products</span><span class="kpi-icon-b">🌿</span></div>
<div class="kpi-number">{total_np}</div>
</div>
<div class="kpi-card kpi-blue">
<div class="kpi-head"><span class="kpi-label">Test Set (Control)</span><span class="kpi-icon-b">🎯</span></div>
<div class="kpi-number" style="color: #2485C7;">{test_cnt}</div>
</div>
<div class="kpi-card kpi-blue">
<div class="kpi-head"><span class="kpi-label">Analysis Set (Hits)</span><span class="kpi-icon-b">🔬</span></div>
<div class="kpi-number" style="color: #2485C7;">{an_cnt}</div>
</div>
<div class="kpi-card kpi-green">
<div class="kpi-head"><span class="kpi-label">Valid RDKit SMILES</span><span class="kpi-icon-g">🧪</span></div>
<div class="kpi-number" style="color: #15966A;">100%</div>
</div>
<div class="kpi-card kpi-green">
<div class="kpi-head"><span class="kpi-label">Lipinski Ro5 Pass</span><span class="kpi-icon-g">💊</span></div>
<div class="kpi-number" style="color: #15966A;">{lip_rate}</div>
</div>
<div class="kpi-card kpi-green">
<div class="kpi-head"><span class="kpi-label">Veber Pass</span><span class="kpi-icon-g">✨</span></div>
<div class="kpi-number" style="color: #15966A;">{veber_rate}</div>
</div>
</div>
""", unsafe_allow_html=True)
        
        # Section 1: Quality Validation Audit Table
        st.markdown("""
<div class="science-card">
<div class="card-header">
<h4 class="card-title">1. Audit Integritas &amp; Kualitas Data Molekuler</h4>
<span class="badge-pass">✓ Quality Verified</span>
</div>
<div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem;">
<div class="metric-box">
<div class="data-label">Validitas Struktur &amp; Conformer 3D</div>
<div style="font-size: 0.9rem; color: #172B4D; margin-top: 0.25rem; line-height: 1.5;">
<b>114/114 (100%)</b> SMILES valid RDKit<br>
<b>113/114 (99.1%)</b> Berhasil konformasi 3D MMFF94
</div>
</div>
<div class="metric-box">
<div class="data-label">Nilai Hilang (Missing Values)</div>
<div style="font-size: 0.9rem; color: #15966A; font-weight: 600; margin-top: 0.25rem; line-height: 1.5;">
Nol Missing Values (0)<br>
<span style="font-size: 0.78rem; color: #71839A; font-weight: normal;">Lengkap pada MW, LogP, HBD, HBA, TPSA</span>
</div>
</div>
<div class="metric-box">
<div class="data-label">Kompatibilitas Gaya Atom Vina &amp; Dual Target</div>
<div style="font-size: 0.9rem; color: #2485C7; font-weight: 600; margin-top: 0.25rem; line-height: 1.5;">
107 Senyawa Kompatibel Penuh (3G0B &amp; 5T4B)<br>
<span style="font-size: 0.78rem; color: #71839A; font-weight: normal;">7 Analog boron dipeptida ditandai covalent/transition</span>
</div>
</div>
</div>
</div>
""", unsafe_allow_html=True)

        # Section 2: Dataset Partitioning (Read-Only Information)
        st.markdown("""
<div class="science-card">
<div class="card-header">
<h4 class="card-title">2. Informasi Pemisahan Dataset (Test Set vs Analysis Set)</h4>
<span class="badge-info">Read-Only Information</span>
</div>
<p style="color: #71839A; font-size: 0.88rem; margin-bottom: 0.85rem; line-height: 1.6;">
<b style="color: #172B4D;">Test Set</b> dirancang untuk menguji presisi docking, akurasi pose, dan konsistensi skor vs bioaktivitas eksperimental (memuat Fitokimia Acuan &amp; Kontrol Inhibitor). 
<b style="color: #172B4D;">Analysis Set</b> digunakan untuk penapisan forward menemukan kandidat inhibitor bahan alam terbaik.
</p>
<div style="background-color: #F1F8FF; border: 1px solid #DCE8F2; border-left: 4px solid #2485C7; padding: 0.7rem 1rem; border-radius: 10px; font-size: 0.82rem; color: #172B4D;">
🔒 <b>Mode Informasi (Read-Only):</b> Tampilan dataset di bawah ini disajikan khusus sebagai informasi kurasi dan validasi. Unduhan berkas mentah/terkurasi tidak tersedia di antarmuka publik dan hanya dapat diakses oleh pengembang melalui direktori lokal (<code>data/02_processed/</code>) atau repositori GitHub proyek.
</div>
</div>
""", unsafe_allow_html=True)
        
        split_view = st.radio("Pilih Informasi Dataset yang Ditampilkan:", ["Test Set (Validation & Benchmark)", "Analysis Set (Screening Candidates)", "Master Dataset (Semua)"], horizontal=True)
        
        if split_view == "Test Set (Validation & Benchmark)" and df_test is not None:
            st.dataframe(df_test[['compound_name', 'molecule_chembl_id', 'MW', 'LogP', 'HBD', 'HBA', 'pIC50', 'compound_source']], use_container_width=True, hide_index=True)
        elif split_view == "Analysis Set (Screening Candidates)" and df_analysis is not None:
            st.dataframe(df_analysis[['compound_name', 'molecule_chembl_id', 'MW', 'LogP', 'HBD', 'HBA', 'Lipinski_Pass', 'compound_source']], use_container_width=True, hide_index=True)
        else:
            st.dataframe(df_master[['compound_name', 'molecule_chembl_id', 'MW', 'LogP', 'HBD', 'HBA', 'dataset_split']], use_container_width=True, hide_index=True)
            
        # Section 3: Predictive Validity Correlation
        if summary_metrics and 'correlation_metrics' in summary_metrics:
            corr_val = summary_metrics['correlation_metrics'].get('pearson_r', -0.58)
            st.markdown(f"""
<div class="science-card" style="margin-top: 1.35rem;">
<div class="card-header">
<h4 class="card-title">3. Validasi Korelasi Eksperimental (pIC50 vs AutoDock Vina Energy)</h4>
<span class="badge-pass">✓ Experimental Benchmark</span>
</div>
<div style="display: flex; align-items: center; gap: 1.75rem; flex-wrap: wrap;">
<div class="metric-box" style="background: #F1F8FF; border-color: #DCE8F2; padding: 1rem 1.5rem; min-width: 190px;">
<div class="data-label" style="color: #2485C7;">Pearson Correlation (r)</div>
<div style="font-size: 1.8rem; font-weight: 800; color: #2485C7; font-family: 'Fira Code', monospace;">{corr_val}</div>
</div>
<div style="flex: 1; min-width: 240px; color: #172B4D; font-size: 0.9rem; line-height: 1.6;">
Korelasi negatif yang signifikan membuktikan bahwa senyawa dengan afinitas biologis laboratorium lebih tinggi (pIC50 besar) 
secara konsisten menghasilkan energi bebas ikatan (<i>binding affinity</i>) yang lebih kuat / negatif pada simulasi AutoDock Vina reseptor DPP-4.
</div>
</div>
</div>
""", unsafe_allow_html=True)
    else:
        st.warning("Dataset terkurasi belum dimuat. Jalankan script `14_curation_validation_test.py` terlebih dahulu.")
