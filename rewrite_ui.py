import os

with open('/Users/aanmuf/drug_discovery/dpp4_project/scripts/11_dashboard_ui.py', 'r') as f:
    content = f.read()

# Since replacing the whole content is easier to get right:
new_content = """import streamlit as st
import pandas as pd
import plotly.express as px
import os
import subprocess
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
from meeko import MoleculePreparation

st.set_page_config(page_title="DPP-4 Drug Discovery Dashboard", page_icon="🧬", layout="wide")
st.title("🧬 DPP-4 Inhibitor Virtual Screening Dashboard")
st.markdown("Antarmuka interaktif (UI/UX) untuk memvisualisasikan hasil penapisan virtual senyawa bahan alam dari pangkalan data ChEMBL terhadap reseptor DPP-4.")

@st.cache_data
def load_data():
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    raw_path = os.path.join(project_root, 'data', '01_raw', 'dpp4_activities_raw.csv')
    processed_path = os.path.join(project_root, 'data', '02_processed', 'dpp4_natural_products.csv')
    res_3g0b_path = os.path.join(project_root, 'results', '02_docking', 'Final_Report_3G0B.csv')
    res_5t4b_path = os.path.join(project_root, 'results', '02_docking', 'Final_Report_5T4B.csv')
    
    df_raw = pd.read_csv(raw_path) if os.path.exists(raw_path) else None
    df_proc = pd.read_csv(processed_path) if os.path.exists(processed_path) else None
    df_3g0b = pd.read_csv(res_3g0b_path) if os.path.exists(res_3g0b_path) else None
    df_5t4b = pd.read_csv(res_5t4b_path) if os.path.exists(res_5t4b_path) else None
    return df_raw, df_proc, df_3g0b, df_5t4b

df_raw, df_proc, df_3g0b, df_5t4b = load_data()

tab1, tab2, tab3, tab4 = st.tabs(["🚀 Live Screening", "📊 Hasil Molecular Docking", "🧪 Data ChEMBL", "⚙️ Analisis Deskriptor"])

with tab1:
    st.header("🚀 Uji Senyawa Anda Sendiri (Live Virtual Screening)")
    st.write("Masukkan struktur SMILES senyawa kandidat Anda untuk langsung ditambatkan (di-docking) ke dalam kantong aktif protein DPP-4.")
    
    target_protein = st.selectbox("Pilih Target Protein:", ["3G0B", "5T4B"])
    smiles_input = st.text_input("SMILES Ligan:", "CC(C)[C@H](N)C(=O)N1CCCC1")
    
    if st.button("Jalankan Molecular Docking", type="primary"):
        if smiles_input:
            with st.spinner("Sedang memproses struktur dan menjalankan Vina..."):
                try:
                    mol = Chem.MolFromSmiles(smiles_input)
                    if mol is None:
                        st.error("Format SMILES tidak valid!")
                    else:
                        canonical = Chem.MolToSmiles(mol)
                        formula = rdMolDescriptors.CalcMolFormula(mol)
                        
                        # Pencarian di database lokal untuk metadata tambahan
                        chembl_id = "Tidak terdaftar (Senyawa Baru)"
                        ic50 = "N/A"
                        nama_senyawa = "Unknown"
                        
                        if df_raw is not None:
                            match = df_raw[df_raw['canonical_smiles'] == canonical]
                            if len(match) > 0:
                                chembl_id = match.iloc[0]['molecule_chembl_id']
                                ic50 = f"{match.iloc[0]['standard_value']} {match.iloc[0]['standard_units']}"
                                # Coba ambil nama dari API
                                try:
                                    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
                                    response = requests.get(url, timeout=3)
                                    if response.status_code == 200:
                                        data = response.json()
                                        nama_senyawa = data.get('pref_name') or "Unnamed"
                                except:
                                    pass
                        
                        # Tampilan Identitas Molekul
                        st.subheader("📋 Identitas & Properti Senyawa")
                        col_id1, col_id2 = st.columns([1,2])
                        with col_id1:
                            img = Draw.MolToImage(mol, size=(300, 300))
                            st.image(img, caption="Struktur 2D")
                        with col_id2:
                            st.markdown(f'''
                            - **SMILES**: `{canonical}`
                            - **Rumus Senyawa**: **{formula}**
                            - **ChEMBL ID**: **{chembl_id}**
                            - **Nama Senyawa**: **{nama_senyawa}**
                            - **Nilai IC50 (Eksperimen in-vitro)**: **{ic50}**
                            ''')
                            
                            # Lipinski
                            mw = Descriptors.MolWt(mol)
                            logp = Descriptors.MolLogP(mol)
                            hbd = Descriptors.NumHDonors(mol)
                            hba = Descriptors.NumHAcceptors(mol)
                            
                            st.markdown("---")
                            st.write(f"**Berat Molekul (MW):** {mw:.2f} g/mol")
                            st.write(f"**Lipofilisitas (LogP):** {logp:.2f}")
                            
                        # Docking
                        st.subheader("⚙️ Hasil Molecular Docking")
                        mol = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(mol, randomSeed=42)
                        AllChem.MMFFOptimizeMolecule(mol)
                        
                        preparator = MoleculePreparation()
                        preparator.prepare(mol)
                        pdbqt_string = preparator.write_pdbqt_string()
                        
                        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                        temp_ligand_path = os.path.join(project_root, 'results', 'temp_ligand.pdbqt')
                        temp_out_path = os.path.join(project_root, 'results', 'temp_out.pdbqt')
                        with open(temp_ligand_path, 'w') as f: f.write(pdbqt_string)
                            
                        vina_exec = os.path.join(project_root, 'scripts', 'vina')
                        receptor_path = os.path.join(project_root, 'data', '05_receptors_pdbqt', f"{target_protein}.pdbqt")
                        pdb_receptor_path = os.path.join(project_root, 'data', '03_pdb', f"{target_protein}_clean.pdb")
                        
                        targets = {
                            '3G0B': {'center': (42.05, 34.29, 14.62), 'size': (18, 18, 18)},
                            '5T4B': {'center': (37.61, 49.92, 40.33), 'size': (18, 16, 23)}
                        }
                        params = targets[target_protein]
                        
                        cmd = [
                            vina_exec, '--receptor', receptor_path, '--ligand', temp_ligand_path,
                            '--center_x', str(params['center'][0]), '--center_y', str(params['center'][1]), '--center_z', str(params['center'][2]),
                            '--size_x', str(params['size'][0]), '--size_y', str(params['size'][1]), '--size_z', str(params['size'][2]),
                            '--out', temp_out_path, '--exhaustiveness', '1'
                        ]
                        
                        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        if process.returncode == 0:
                            stdout_text = process.stdout.decode('utf-8')
                            poses = []
                            in_table = False
                            for line in stdout_text.split('\\n'):
                                if line.startswith('---+'):
                                    in_table = True; continue
                                if line.startswith('Writing output'):
                                    in_table = False; break
                                if in_table and line.strip():
                                    parts = line.strip().split()
                                    if len(parts) >= 4:
                                        poses.append({'Mode': parts[0], 'Affinity (kcal/mol)': parts[1], 'RMSD l.b.': parts[2], 'RMSD u.b.': parts[3]})
                            if poses:
                                st.success(f"🎉 Afinitas Terbaik: **{poses[0]['Affinity (kcal/mol)']} kcal/mol** | RMSD Terbaik: **{poses[0]['RMSD l.b.']} Å**")
                                st.table(pd.DataFrame(poses))
                                
                                st.subheader("👁️ Visualisasi 3D Kompleks Protein-Ligan")
                                try:
                                    import py3Dmol
                                    from stmol import showmol
                                    with open(pdb_receptor_path, 'r') as f: rec_data = f.read()
                                    with open(temp_out_path, 'r') as f: lig_data = f.read()
                                        
                                    view = py3Dmol.view(width=800, height=500)
                                    view.addModel(rec_data, 'pdb')
                                    view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})
                                    view.addModel(lig_data, 'pdbqt')
                                    view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})
                                    view.zoomTo({'model': -1})
                                    showmol(view, height=500, width=800)
                                except:
                                    st.warning("Visualisasi 3D membutuhkan py3Dmol & stmol.")
                        else:
                            st.error("Gagal menjalankan Vina.")
                except Exception as e:
                    st.error(f"Error: {e}")
        else:
            st.warning("Mohon masukkan kode SMILES.")

with tab2:
    st.header("Hasil Akhir Penambatan Molekuler (Docking)")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("3G0B")
        if df_3g0b is not None:
            st.dataframe(df_3g0b, use_container_width=True)
    with col2:
        st.subheader("5T4B")
        if df_5t4b is not None:
            st.dataframe(df_5t4b, use_container_width=True)

with tab3:
    st.header("Data Mentah Bioaktivitas ChEMBL")
    if df_raw is not None: st.dataframe(df_raw.head(100), use_container_width=True)

with tab4:
    st.header("Data Molekul Alam Terkurasi")
    if df_proc is not None: st.dataframe(df_proc, use_container_width=True)
"""

with open('/Users/aanmuf/drug_discovery/dpp4_project/scripts/11_dashboard_ui.py', 'w') as f:
    f.write(new_content)

print("Berhasil diperbarui!")
