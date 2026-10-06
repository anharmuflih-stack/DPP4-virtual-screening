import streamlit as st
import os
import platform
import subprocess
import tempfile
import urllib.request
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from stmol import showmol
import py3Dmol

# Set konfigurasi halaman lebar penuh dengan tema cerah
st.set_page_config(page_title="Dirof DPP-4", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

# ==========================================
# 1. INJEKSI CSS CUSTOM (Gaya Clean Lab dari Localhost)
# ==========================================
st.markdown("""
    <style>
        /* Mengubah font utama */
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Fira+Code&display=swap');
        
        html, body, [class*="css"]  {
            font-family: 'Inter', sans-serif;
        }
        
        /* Menghilangkan ruang kosong berlebih di atas halaman */
        .block-container {
            padding-top: 1rem !important;
            padding-left: 2rem !important;
            padding-right: 2rem !important;
            max-width: 90% !important;
        }
        
        /* Tema Background Slate 100 */
        .stApp {
            background-color: #f1f5f9; 
            color: #334155;
        }
        
        /* Menyembunyikan elemen bawaan Streamlit yang mengganggu */
        header {visibility: hidden !important;}
        #MainMenu {visibility: hidden !important;}
        footer {visibility: hidden !important;}
        
        /* Gaya Kartu Sains (Science Card) */
        .science-card {
            background: #ffffff;
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05), 0 2px 4px -1px rgba(0, 0, 0, 0.03);
            padding: 1.5rem;
            margin-bottom: 1.5rem;
        }
        
        /* Mengganti desain Input bawaan Streamlit (Mirip Tailwind) */
        div[data-baseweb="select"] > div, 
        div[data-baseweb="input"] > div {
            background-color: #f8fafc !important;
            border: 1px solid #cbd5e1 !important;
            border-radius: 0.375rem !important;
            font-family: 'Inter', sans-serif !important;
        }
        
        label[data-testid="stWidgetLabel"] p {
            font-size: 0.75rem !important;
            text-transform: uppercase !important;
            letter-spacing: 0.05em !important;
            color: #64748b !important;
            font-weight: 600 !important;
        }
        
        /* Gaya Teks */
        .data-label {
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #64748b;
            font-weight: 600;
            margin-bottom: 0.25rem;
        }
        .data-value {
            font-size: 1.125rem;
            font-weight: 600;
            color: #1e293b;
            font-family: 'Fira Code', monospace;
        }
        
        /* Header Kustom */
        .header-title {
            font-size: 2rem;
            font-weight: 700;
            color: #1e293b;
            border-bottom: 1px solid #cbd5e1;
            padding-bottom: 1rem;
            margin-bottom: 2rem;
            display: flex;
            align-items: center;
        }
        .header-subtitle {
            font-size: 1rem;
            color: #64748b;
            font-weight: 400;
            margin-top: 0.25rem;
        }
        
        /* Metrik Grid HTML */
        .metric-box {
            background-color: #f8fafc;
            border: 1px solid #e2e8f0;
            border-radius: 0.5rem;
            padding: 1rem;
        }
        
        .lipinski-box {
            border-left: 4px solid #10b981; /* Emerald green (Pass) */
            padding-left: 0.75rem;
        }
        .lipinski-fail {
            border-left: 4px solid #f43f5e; /* Rose red (Fail) */
        }
        
        /* Modifikasi Tombol Streamlit secara Ekstrem */
        div[data-testid="stButton"] > button {
            width: 100% !important;
            background-color: #0284c7 !important;
            color: white !important;
            font-weight: 600 !important;
            border-radius: 0.375rem !important;
            border: none !important;
            padding: 0.6rem !important;
            height: 3rem !important;
            box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05) !important;
            transition: all 0.2s !important;
        }
        div[data-testid="stButton"] > button:hover {
            background-color: #0369a1 !important;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1) !important;
        }
        
        div[data-testid="stSpinner"] > div {
            border-color: #0284c7 transparent transparent transparent !important;
        }
    </style>
""", unsafe_allow_html=True)

# ==========================================
# 2. SETUP DIREKTORI & VINA
# ==========================================
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(project_root, 'data')
bin_dir = os.path.join(project_root, 'scripts')

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

# ==========================================
# 3. FUNGSI LOGIKA DOCKING
# ==========================================
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
    
    props = {
        'MW': round(Descriptors.MolWt(mol), 2),
        'LogP': round(Descriptors.MolLogP(mol), 2),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'Formula': Chem.rdMolDescriptors.CalcMolFormula(mol)
    }
    return pdbqt_string, props

def run_docking(pdbqt_ligand, protein):
    targets = {
        '3G0B': {'center': (42.05, 34.29, 14.62), 'size': (18, 18, 18)},
        '5T4B': {'center': (37.61, 49.92, 40.33), 'size': (18, 16, 23)}
    }
    
    with tempfile.NamedTemporaryFile(delete=False, suffix='.pdbqt') as lig_file:
        lig_file.write(pdbqt_ligand.encode('utf-8'))
        lig_path = lig_file.name
        
    out_path = lig_path.replace('.pdbqt', '_out.pdbqt')
    receptor_path = os.path.join(data_dir, '05_receptors_pdbqt', f"{protein}.pdbqt")
    
    cmd = [
        vina_exec, '--receptor', receptor_path, '--ligand', lig_path,
        '--center_x', str(targets[protein]['center'][0]),
        '--center_y', str(targets[protein]['center'][1]),
        '--center_z', str(targets[protein]['center'][2]),
        '--size_x', str(targets[protein]['size'][0]),
        '--size_y', str(targets[protein]['size'][1]),
        '--size_z', str(targets[protein]['size'][2]),
        '--exhaustiveness', '8', '--num_modes', '5'
    ]
    
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    poses = []
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
                
    return poses, out_content

# ==========================================
# 4. TATA LETAK UI (Meniru HTML)
# ==========================================
st.markdown("""
    <div class="header-title">
        🧬 DPP-4 Virtual Screener 
        <div class="header-subtitle">Computational Drug Discovery & Virtual Screening Laboratory</div>
    </div>
""", unsafe_allow_html=True)

# Layout Grid (mirip cols lg:col-span-1 & lg:col-span-3)
col_input, col_result = st.columns([1, 2.5], gap="large")

with col_input:
    st.markdown('<h3 style="border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; color: #0f172a; font-size: 1.1rem; margin-bottom: 1rem;">Analysis Setup</h3>', unsafe_allow_html=True)
    
    protein_target = st.selectbox("Receptor Model", ["3G0B", "5T4B"], help="Pilih resolusi target protein")
    smiles_input = st.text_input("Ligand Structure (SMILES)", "CC(C)[C@H](N)C(=O)N1CCCC1")
    
    st.markdown("<br>", unsafe_allow_html=True)
    run_btn = st.button("Run AutoDock Vina")

with col_result:
    if run_btn:
        with st.spinner("Mengalkulasi Interaksi Molekuler..."):
            ligand_pdbqt, props = process_ligand(smiles_input)
            
            if not ligand_pdbqt:
                st.error("Gagal: SMILES tidak valid.")
            else:
                poses, docked_pdbqt = run_docking(ligand_pdbqt, protein_target)
                best_affinity = poses[0]['Affinity (kcal/mol)'] if poses else 0.0
                
                # Cek batas Lipinski
                c_mw = "lipinski-fail" if props['MW'] > 500 else ""
                c_logp = "lipinski-fail" if props['LogP'] > 5 else ""
                c_hbd = "lipinski-fail" if props['HBD'] > 5 else ""
                c_hba = "lipinski-fail" if props['HBA'] > 10 else ""
                
                # Kartu 1: Identitas & Properti (HTML rata kiri absolut)
                html_card = f"""
<div class="science-card">
<h3 style="border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; color: #0f172a; font-size: 1.1rem; margin-bottom: 1rem;">Molecular Properties & Identity</h3>
<div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem; margin-bottom: 1.5rem;">
<div class="metric-box">
<div class="data-label">Chemical Formula</div>
<div class="data-value">{props['Formula']}</div>
</div>
<div class="metric-box">
<div class="data-label">ChEMBL ID</div>
<div class="data-value" style="color: #0284c7;">-</div>
</div>
<div class="metric-box">
<div class="data-label">Binding Affinity</div>
<div class="data-value" style="color: #f43f5e;">{best_affinity} kcal/mol</div>
</div>
</div>
<div class="data-label" style="margin-bottom: 1rem;">Lipinski's Rule of Five Analysis</div>
<div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 1rem;">
<div class="lipinski-box {c_mw}">
<div style="font-size: 0.7rem; color: #64748b;">Molecular Weight</div>
<div class="data-value" style="font-size: 1rem;">{props['MW']} <span style="font-size: 0.7rem; color:#94a3b8; font-weight:normal;">Da</span></div>
</div>
<div class="lipinski-box {c_logp}">
<div style="font-size: 0.7rem; color: #64748b;">LogP</div>
<div class="data-value" style="font-size: 1rem;">{props['LogP']}</div>
</div>
<div class="lipinski-box {c_hbd}">
<div style="font-size: 0.7rem; color: #64748b;">H-Bond Donors</div>
<div class="data-value" style="font-size: 1rem;">{props['HBD']}</div>
</div>
<div class="lipinski-box {c_hba}">
<div style="font-size: 0.7rem; color: #64748b;">H-Bond Acceptors</div>
<div class="data-value" style="font-size: 1rem;">{props['HBA']}</div>
</div>
</div>
</div>
"""
                st.markdown(html_card, unsafe_allow_html=True)
                
                # Tabel Top Poses
                if poses:
                    st.markdown('<div class="science-card"><h3 style="color: #0f172a; font-size: 1.1rem; margin-bottom: 1rem;">Docking Conformational Poses (Top 5)</h3>', unsafe_allow_html=True)
                    df_poses = pd.DataFrame(poses)
                    st.dataframe(df_poses.style.highlight_min(subset=['Affinity (kcal/mol)'], color='#fca5a5'), use_container_width=True, hide_index=True)
                    st.markdown('</div>', unsafe_allow_html=True)
                
                # Kartu 2: Interaksi 3D (Hanya Title yang dibungkus HTML)
                html_title = """
<div style="margin-top: 1rem; border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; margin-bottom: 1rem;">
<h3 style="color: #0f172a; font-size: 1.1rem; margin: 0;">Molecular Interaction Analysis</h3>
</div>
"""
                st.markdown(html_title, unsafe_allow_html=True)
                
                view = py3Dmol.view(width="100%", height=450)
                view.setBackgroundColor('#1e293b') # Dark slate background for contrast
                
                # Add Protein
                receptor_pdb = os.path.join(data_dir, '03_pdb', f"{protein_target}_clean.pdb")
                with open(receptor_pdb, 'r') as f:
                    view.addModel(f.read(), 'pdb')
                
                # Ekstrak hanya Pose Terbaik (MODEL 1) untuk ditampilkan di 3D Viewer
                best_pose_pdbqt = docked_pdbqt.split("MODEL 2")[0] if "MODEL 2" in docked_pdbqt else docked_pdbqt
                
                # Add Ligand
                view.addModel(best_pose_pdbqt, 'pdbqt')
                
                # Styles (Kontras Ekstrem untuk Kejelasan Maksimal)
                view.setStyle({'model': 0}, {'cartoon': {'color': '#e2e8f0', 'style': 'oval', 'thickness': 0.2}})
                
                # Ligand: Kuning Cerah dan sangat tebal
                view.setStyle({'model': 1}, {'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.25}})
                
                # Interacting Residues (Jarak 3.5 Angstrom, Warna Magenta Menyala, plus Label Tebal)
                interaction_sel = {'model': 0, 'within': {'distance': 3.5, 'sel': {'model': 1}}}
                view.addStyle(interaction_sel, {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.15}})
                view.addResLabels(interaction_sel, {
                    'fontOpacity': 1.0, 
                    'fontSize': 14, 
                    'fontColor': '#ffffff', 
                    'backgroundColor': '#db2777', 
                    'showBackground': True
                })
                
                # Permukaan saku ikat (Binding Pocket Surface) berwarna biru transparan
                view.addSurface(py3Dmol.VDW, {'opacity': 0.35, 'color': '#0ea5e9'}, interaction_sel, interaction_sel)
                
                view.zoomTo({'model': 1})
                showmol(view, height=450, width="100%")
    else:
        # Tampilan kosong di awal
        st.info("👈 Masukkan struktur SMILES dan pilih model protein di panel sebelah kiri untuk memulai penapisan (*screening*).")
