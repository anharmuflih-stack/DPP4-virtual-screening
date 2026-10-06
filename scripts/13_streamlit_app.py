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

st.set_page_config(page_title="Dirof DPP-4", page_icon="🧬", layout="wide")

# ==========================================
# 1. SETUP LINGKUNGAN & VINA
# ==========================================
st.title("🧬 Dirof DPP-4 Virtual Screening")
st.markdown("Platform komputasional penemuan obat untuk reseptor *Dipeptidyl Peptidase-4*.")

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(project_root, 'data')
bin_dir = os.path.join(project_root, 'scripts')

# Fungsi mendownload Vina jika berjalan di Linux (Streamlit Cloud)
@st.cache_resource
def setup_vina():
    is_linux = platform.system() == 'Linux'
    vina_exec = os.path.join(bin_dir, 'vina_linux' if is_linux else 'vina')
    
    if is_linux and not os.path.exists(vina_exec):
        st.info("Mengunduh modul AutoDock Vina untuk server Linux...")
        url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
        urllib.request.urlretrieve(url, vina_exec)
        os.chmod(vina_exec, 0o755)
    return vina_exec

vina_exec = setup_vina()

# ==========================================
# 2. UI SIDEBAR
# ==========================================
st.sidebar.header("Konfigurasi Simulasi")
protein_target = st.sidebar.selectbox("Pilih Target Protein (Receptor)", ["3G0B", "5T4B"])
smiles_input = st.sidebar.text_input("SMILES Senyawa", "CC(C)[C@H](N)C(=O)N1CCCC1")
run_btn = st.sidebar.button("Mulai Docking Vina", type="primary")

targets = {
    '3G0B': {'center': (42.05, 34.29, 14.62), 'size': (18, 18, 18)},
    '5T4B': {'center': (37.61, 49.92, 40.33), 'size': (18, 16, 23)}
}

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
    preparator.prepare(mol)
    writer = PDBQTWriterLegacy()
    pdbqt_string = writer.write_string(preparator.setup[0])
    
    props = {
        'MW': round(Descriptors.MolWt(mol), 2),
        'LogP': round(Descriptors.MolLogP(mol), 2),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol)
    }
    return pdbqt_string, props

def run_docking(pdbqt_ligand, protein):
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
        '--exhaustiveness', '4', '--num_modes', '1'
    ]
    
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    affinity = 0.0
    with open(out_path, 'r') as f:
        out_content = f.read()
        for line in out_content.split('\n'):
            if line.startswith('REMARK VINA RESULT:'):
                affinity = float(line.split()[3])
                break
                
    return affinity, out_content

# ==========================================
# 4. EKSEKUSI & VISUALISASI
# ==========================================
if run_btn:
    with st.spinner("Memproses Ligan dan Menjalankan Docking..."):
        ligand_pdbqt, props = process_ligand(smiles_input)
        
        if not ligand_pdbqt:
            st.error(props) # Menampilkan pesan error SMILES
        else:
            col1, col2 = st.columns([1, 2])
            
            with col1:
                st.subheader("Aturan Lipinski (Ro5)")
                st.metric("Molecular Weight", f"{props['MW']} Da", "≤ 500", delta_color="inverse" if props['MW']>500 else "normal")
                st.metric("LogP (Lipofilisitas)", props['LogP'], "≤ 5", delta_color="inverse" if props['LogP']>5 else "normal")
                st.metric("H-Bond Donors", props['HBD'], "≤ 5", delta_color="inverse" if props['HBD']>5 else "normal")
                st.metric("H-Bond Acceptors", props['HBA'], "≤ 10", delta_color="inverse" if props['HBA']>10 else "normal")
            
            # Jalankan Docking
            affinity, docked_pdbqt = run_docking(ligand_pdbqt, protein_target)
            
            with col2:
                st.success(f"**Docking Berhasil!** Binding Affinity Terbaik: **{affinity} kcal/mol**")
                
                # Visualisasi 3D dengan Py3DMol via stmol
                st.subheader("Interaksi Molekuler 3D")
                view = py3Dmol.view(width=800, height=500)
                view.setBackgroundColor('#1e293b')
                
                # Muat Protein
                receptor_pdb = os.path.join(data_dir, '03_pdb', f"{protein_target}_clean.pdb")
                with open(receptor_pdb, 'r') as f:
                    view.addModel(f.read(), 'pdb')
                view.setStyle({'model': 0}, {'cartoon': {'color': '#cbd5e1', 'style': 'oval', 'opacity': 0.8}})
                
                # Muat Ligan Hasil Docking
                view.addModel(docked_pdbqt, 'pdbqt')
                view.setStyle({'model': 1}, {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.2}})
                
                # Interaksi (Residu 5 Angstrom)
                view.addStyle({'model': 0, 'within': {'distance': 5.0, 'sel': {'model': 1}}}, 
                              {'stick': {'colorscheme': 'whiteCarbon', 'radius': 0.15}})
                
                # Surface Kantong Ikatan
                view.addSurface(py3Dmol.VDW, {'opacity': 0.4, 'color': '#38bdf8'}, 
                                {'model': 0, 'within': {'distance': 5.0, 'sel': {'model': 1}}})
                
                view.zoomTo({'model': 1})
                showmol(view, height=500, width=800)
