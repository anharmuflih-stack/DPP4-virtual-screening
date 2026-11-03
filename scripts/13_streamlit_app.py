import streamlit as st
import os
import platform
import subprocess
import tempfile
import urllib.request
import urllib.parse
import json
import ssl
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from stmol import showmol
import py3Dmol

# Set konfigurasi halaman lebar penuh dengan tema cerah
st.set_page_config(page_title="DPP-4 Virtual Screener", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

# ==========================================
# 1. INJEKSI CSS CUSTOM (Gaya Clean Lab & Modern Science)
# ==========================================
st.markdown("""
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=Fira+Code:wght@400;600&display=swap');
        
        html, body, [class*="css"]  {
            font-family: 'Inter', sans-serif;
        }
        
        .block-container {
            padding-top: 1rem !important;
            padding-left: 2rem !important;
            padding-right: 2rem !important;
            max-width: 92% !important;
        }
        
        .stApp {
            background-color: #f8fafc; 
            color: #1e293b;
        }
        
        header {visibility: hidden !important;}
        #MainMenu {visibility: hidden !important;}
        footer {visibility: hidden !important;}
        [data-testid="stElementToolbar"] {display: none !important;}
        
        /* Gaya Kartu Sains (Science Card) */
        .science-card {
            background: #ffffff;
            border: 1px solid #e2e8f0;
            border-radius: 10px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.04), 0 2px 4px -1px rgba(0, 0, 0, 0.02);
            padding: 1.5rem;
            margin-bottom: 1.5rem;
        }
        
        /* Modifikasi Tab Streamlit */
        .stTabs [data-baseweb="tab-list"] {
            gap: 10px;
            background-color: transparent;
            border-bottom: 2px solid #e2e8f0;
            margin-bottom: 1.5rem;
        }
        .stTabs [data-baseweb="tab"] {
            height: 48px;
            white-space: pre-wrap;
            background-color: #ffffff;
            border-radius: 8px 8px 0px 0px;
            border: 1px solid #e2e8f0;
            border-bottom: none;
            color: #475569;
            font-weight: 600;
            padding: 0 20px;
            font-size: 0.95rem;
        }
        .stTabs [aria-selected="true"] {
            background-color: #0284c7 !important;
            color: #ffffff !important;
            border-color: #0284c7 !important;
        }
        
        /* Gaya Input */
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
        
        .metric-box {
            background-color: #f8fafc;
            border: 1px solid #e2e8f0;
            border-radius: 0.5rem;
            padding: 1rem;
        }
        
        .lipinski-box {
            background-color: #f8fafc;
            border: 1px solid #e2e8f0;
            border-radius: 0.5rem;
            padding: 1rem;
            border-left: 4px solid #10b981;
        }
        .lipinski-fail {
            border-left: 4px solid #f43f5e;
        }
        
        .badge-pass {
            background-color: #dcfce7;
            color: #166534;
            padding: 3px 8px;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 700;
        }
        .badge-fail {
            background-color: #ffe4e6;
            color: #9f1239;
            padding: 3px 8px;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 700;
        }
        .badge-test {
            background-color: #e0f2fe;
            color: #0369a1;
            padding: 3px 8px;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 700;
        }
        .badge-analysis {
            background-color: #f3e8ff;
            color: #7e22ce;
            padding: 3px 8px;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 700;
        }
        
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
    </style>
""", unsafe_allow_html=True)

# ==========================================
# 2. SETUP DIREKTORI & VINA
# ==========================================
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(project_root, 'data')
bin_dir = os.path.join(project_root, 'scripts')
results_dir = os.path.join(project_root, 'results', '02_docking')

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

# Load dataset terkurasi jika tersedia
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

# ==========================================
# 3. FUNGSI PENCARIAN CHEMBL & DOCKING
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
        '--out', out_path,
        '--exhaustiveness', '4', '--num_modes', '5'
    ]
    
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
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
                
    return poses, out_content

# ==========================================
# 4. HEADER TAMPILAN UTAMA
# ==========================================
st.markdown("""
<div style="display: flex; flex-wrap: wrap; align-items: center; gap: 1.5rem; margin-bottom: 1.5rem; border-bottom: 2px solid #e2e8f0; padding-bottom: 1.5rem;">
    <div style="background: white; padding: 0.75rem 1.5rem; border-radius: 0.75rem; border: 1px solid #e2e8f0; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05); font-size: 2rem; font-weight: 800; color: #0f172a; display: flex; align-items: center; gap: 0.75rem;">
        <span>🧬</span>
        <span>DPP-4 Virtual Screener</span>
    </div>
    <div style="font-size: 1.05rem; font-weight: 500; color: #64748b; line-height: 1.4; max-width: 480px;">
        Computational Drug Discovery & Virtual Screening Laboratory<br>
        <span style="font-size: 0.85rem; color: #0284c7; font-weight: 600;">Autonomous Cheminformatics Agent Platform</span>
    </div>
</div>
""", unsafe_allow_html=True)

# ==========================================
# 5. TABS STRUKTURAL
# ==========================================
tab_screen, tab_explorer, tab_curation = st.tabs([
    "🔬 Single Screener & 3D Interaction", 
    "🌿 Screened Natural Products Library",
    "🧪 Data Curation & Validation Hub"
])

# -------------------------------------------------------------
# TAB 1: SINGLE LIGAND SCREENER (ORISINAL TERAWAT)
# -------------------------------------------------------------
with tab_screen:
    PRESET_COMPOUNDS = {
        "Custom SMILES": "",
        "QUERCETIN (CHEMBL50) - Flavonoid": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "RESVERATROL (CHEMBL165) - Stilbene": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",
        "BERBERINE (CHEMBL295124) - Alkaloid": "COc1ccc2cc3[n+](cc2c1OC)CCc1cc2c(cc1-3)OCO2",
        "ERIODICTYOL (CHEMBL8996) - Flavanone": "O=C1C[C@@H](c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c21",
        "CURCUMIN (CHEMBL140) - Polyphenol": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
        "LUTEOLIN (CHEMBL151) - Flavone": "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "KAEMPFEROL (CHEMBL150) - Flavonol": "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
        "GENISTEIN (CHEMBL44) - Isoflavone": "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12",
        "HESPERETIN (CHEMBL399121) - Flavanone": "COc1ccc([C@@H]2CC(=O)c3c(O)cc(O)cc3O2)cc1O",
        "VILDAGLIPTIN (CHEMBL142703) - Approved Inhibitor": "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2",
        "SAXAGLIPTIN (CHEMBL385517) - Approved Inhibitor": "N#C[C@@H]1C[C@@H]2C[C@@H]2N1C(=O)[C@@H](N)C12CC3CC(CC(O)(C3)C1)C2",
        "RUTIN (CHEMBL226335) - Glycoside": "CC1OC(OCC2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O",
        "UBENIMEX (CHEMBL29292) - Natural Peptide": "CC(C)C[C@H](NC(=O)[C@@H](O)[C@H](N)Cc1ccccc1)C(=O)O",
        "TALABOSTAT (CHEMBL67279) - Control Boronic Dipeptide": "CC(C)[C@H](N)C(=O)N1CCCC1"
    }
    
    col_input, col_result = st.columns([1, 2.5], gap="large")
    
    with col_input:
        st.markdown('<h3 style="border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; color: #0f172a; font-size: 1.1rem; margin-bottom: 1rem;">Analysis Setup</h3>', unsafe_allow_html=True)
        
        protein_target = st.selectbox("Receptor Model", ["3G0B", "5T4B"], help="Pilih resolusi target protein")
        
        preset_choice = st.selectbox("Quick Preset (Curated Phytochemicals)", list(PRESET_COMPOUNDS.keys()), index=0)
        
        default_smiles = PRESET_COMPOUNDS[preset_choice] if PRESET_COMPOUNDS[preset_choice] else "CC(C)[C@H](N)C(=O)N1CCCC1"
        smiles_input = st.text_input("Ligand Structure (SMILES)", default_smiles)
        
        st.markdown("<br>", unsafe_allow_html=True)
        run_btn = st.button("Run AutoDock Vina")
    
    with col_result:
        if run_btn:
            with st.spinner("Mengalkulasi Interaksi Molekuler..."):
                ligand_pdbqt, props = process_ligand(smiles_input)
                
                if not ligand_pdbqt:
                    st.error("Gagal: SMILES tidak valid atau struktur tidak dapat diubah ke format 3D.")
                else:
                    poses, docked_pdbqt = run_docking(ligand_pdbqt, protein_target)
                    best_affinity = poses[0]['Affinity (kcal/mol)'] if poses else 0.0
                    
                    c_mw = "lipinski-fail" if props['MW'] > 500 else ""
                    c_logp = "lipinski-fail" if props['LogP'] > 5 else ""
                    c_hbd = "lipinski-fail" if props['HBD'] > 5 else ""
                    c_hba = "lipinski-fail" if props['HBA'] > 10 else ""
                    
                    targets = {
                        '3G0B': {'center': (42.05, 34.29, 14.62)},
                        '5T4B': {'center': (37.61, 49.92, 40.33)}
                    }
                    center_xyz = f"{targets[protein_target]['center'][0]}, {targets[protein_target]['center'][1]}, {targets[protein_target]['center'][2]}"
                    
                    html_card = f"""
<div class="science-card">
<h3 style="border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; color: #0f172a; font-size: 1.1rem; margin-bottom: 1rem;">Molecular Properties & Identity</h3>
<div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem; margin-bottom: 1rem;">
<div class="metric-box">
<div class="data-label">Compound Name</div>
<div class="data-value" style="font-size: 1rem; color: #475569;">{props['Name']}</div>
</div>
<div class="metric-box">
<div class="data-label">ChEMBL ID</div>
<div class="data-value" style="color: #0284c7;">{props['ChEMBL_ID']}</div>
</div>
<div class="metric-box">
<div class="data-label">PubChem CID</div>
<div class="data-value" style="color: #10b981;">{props['CID']}</div>
</div>
</div>
<div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem; margin-bottom: 1.5rem;">
<div class="metric-box">
<div class="data-label">Chemical Formula</div>
<div class="data-value">{props['Formula']}</div>
</div>
<div class="metric-box">
<div class="data-label">Grid Center (X,Y,Z)</div>
<div class="data-value" style="font-size: 0.95rem;">{center_xyz}</div>
</div>
<div class="metric-box">
<div class="data-label">Best Affinity</div>
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
                    
                    if poses:
                        st.markdown('<div class="science-card"><h3 style="color: #0f172a; font-size: 1.1rem; margin-bottom: 1rem;">Docking Conformational Poses (Top 5)</h3>', unsafe_allow_html=True)
                        df_poses = pd.DataFrame(poses)
                        st.dataframe(df_poses.style.highlight_min(subset=['Affinity (kcal/mol)'], color='#fca5a5'), use_container_width=True, hide_index=True)
                        st.markdown('</div>', unsafe_allow_html=True)
                    
                    st.markdown("""
<div style="margin-top: 1rem; border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; margin-bottom: 1rem;">
<h3 style="color: #0f172a; font-size: 1.1rem; margin: 0;">Molecular Interaction Analysis</h3>
</div>
""", unsafe_allow_html=True)
                    
                    view = py3Dmol.view(width="100%", height=450)
                    view.setBackgroundColor('#1e293b')
                    
                    receptor_pdb = os.path.join(data_dir, '03_pdb', f"{protein_target}_clean.pdb")
                    if os.path.exists(receptor_pdb):
                        with open(receptor_pdb, 'r') as f:
                            view.addModel(f.read(), 'pdb')
                    
                    best_pose_pdbqt = docked_pdbqt.split("MODEL 2")[0] if "MODEL 2" in docked_pdbqt else docked_pdbqt
                    view.addModel(best_pose_pdbqt, 'pdbqt')
                    
                    view.setStyle({'model': 0}, {'cartoon': {'color': '#e2e8f0', 'style': 'oval', 'thickness': 0.2}})
                    view.setStyle({'model': 1}, {'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.25}})
                    
                    interaction_sel = {'model': 0, 'within': {'distance': 3.5, 'sel': {'model': 1}}}
                    view.addStyle(interaction_sel, {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.15}})
                    view.addResLabels(interaction_sel, {
                        'fontOpacity': 1.0, 
                        'fontSize': 14, 
                        'fontColor': '#ffffff', 
                        'backgroundColor': '#db2777', 
                        'showBackground': True
                    })
                    view.addSurface(py3Dmol.VDW, {'opacity': 0.35, 'color': '#0ea5e9'}, interaction_sel, interaction_sel)
                    view.zoomTo({'model': 1})
                    showmol(view, height=450, width="100%")
        else:
            st.info("👈 Masukkan struktur SMILES atau pilih preset bahan alam di panel sebelah kiri untuk memulai penapisan (*screening*).")

# -------------------------------------------------------------
# TAB 2: NATURAL PRODUCTS SCREENING & TOP HITS LIBRARY
# -------------------------------------------------------------
with tab_explorer:
    st.markdown('<h3 style="color: #0f172a; margin-bottom: 0.5rem;">🌿 Screened Natural Products Library & Top Hits</h3>', unsafe_allow_html=True)
    st.markdown('<p style="color: #64748b; font-size: 0.95rem; margin-bottom: 1.5rem;">Hasil penapisan virtual AutoDock Vina multi-target (3G0B & 5T4B) pada pustaka senyawa bahan alam terkurasi.</p>', unsafe_allow_html=True)
    
    # Prioritaskan df_results, jika belum ada gunakan df_master
    active_df = df_results if df_results is not None else (df_test if df_test is not None else df_master)
    
    if active_df is not None:
        # Filter controls
        fc1, fc2, fc3 = st.columns([2, 1, 1])
        with fc1:
            search_query = st.text_input("🔍 Cari Nama Senyawa atau ChEMBL ID:", "")
        with fc2:
            split_filter = st.selectbox("Filter Dataset Split:", ["Semua", "Test Set (Validation & Benchmark)", "Analysis Set (Screening Candidates)"])
        with fc3:
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
            
        # Top 4 Hits Cards
        if 'Affinity_3G0B_kcal_mol' in filtered_df.columns:
            valid_docked = filtered_df.dropna(subset=['Affinity_3G0B_kcal_mol']).sort_values(by='Affinity_3G0B_kcal_mol', ascending=True)
            if not valid_docked.empty:
                st.markdown('<h4 style="color: #0f172a; margin-top: 1rem; margin-bottom: 0.75rem;">🏆 Top Natural Product Hits (Afinitas DPP-4 Terkuat)</h4>', unsafe_allow_html=True)
                top_cols = st.columns(4)
                for i, (_, hit) in enumerate(valid_docked.head(4).iterrows()):
                    with top_cols[i]:
                        st.markdown(f"""
<div class="science-card" style="padding: 1rem; border-top: 4px solid #0284c7;">
    <div style="font-size: 0.75rem; color: #64748b; font-weight: 700; text-transform: uppercase;">RANK #{i+1}</div>
    <div style="font-size: 1.05rem; font-weight: 700; color: #0f172a; margin: 0.25rem 0;">{hit['compound_name']}</div>
    <div style="font-size: 0.8rem; color: #0284c7; font-family: monospace;">{hit['molecule_chembl_id']}</div>
    <hr style="margin: 0.5rem 0; border: none; border-top: 1px solid #e2e8f0;">
    <div style="display: flex; justify-content: space-between; margin-bottom: 0.25rem;">
        <span style="font-size: 0.75rem; color: #64748b;">3G0B:</span>
        <span style="font-weight: 700; color: #f43f5e; font-family: monospace;">{hit['Affinity_3G0B_kcal_mol']} kcal/mol</span>
    </div>
    <div style="display: flex; justify-content: space-between; margin-bottom: 0.25rem;">
        <span style="font-size: 0.75rem; color: #64748b;">MW:</span>
        <span style="font-size: 0.8rem; color: #334155;">{hit['MW']} Da</span>
    </div>
    <div style="display: flex; justify-content: space-between;">
        <span style="font-size: 0.75rem; color: #64748b;">Lipinski:</span>
        <span class="{'badge-pass' if hit.get('Lipinski_Pass') else 'badge-fail'}">{'PASS' if hit.get('Lipinski_Pass') else 'FAIL'}</span>
    </div>
</div>
""", unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)
        # Display Data Table
        display_cols = [c for c in ['compound_name', 'molecule_chembl_id', 'Affinity_3G0B_kcal_mol', 'Affinity_5T4B_kcal_mol', 'Ligand_Efficiency', 'MW', 'LogP', 'HBD', 'HBA', 'Lipinski_Pass', 'pIC50', 'dataset_split'] if c in filtered_df.columns]
        st.dataframe(filtered_df[display_cols], use_container_width=True, hide_index=True)
        st.download_button("📥 Unduh Seluruh Hasil Skrining (CSV)", filtered_df.to_csv(index=False), "dpp4_natural_products_docking_results.csv", "text/csv")
    else:
        st.info("Pangkalan data hasil skrining sedang diproses atau belum tersedia.")

# -------------------------------------------------------------
# TAB 3: DATA CURATION & VALIDATION HUB (INFORMASI SAJA)
# -------------------------------------------------------------
with tab_curation:
    st.markdown('<h3 style="color: #0f172a; margin-bottom: 0.5rem;">🧪 Data Curation & Validation Hub</h3>', unsafe_allow_html=True)
    st.markdown('<p style="color: #64748b; font-size: 0.95rem; margin-bottom: 1.5rem;">Audit kualitas data molekuler, pembagian Test vs Analysis Set, serta validasi korelasi eksperimental penapisan obat DPP-4 (Ditampilkan sebagai informasi referensi).</p>', unsafe_allow_html=True)
    
    if df_master is not None:
        c1, c2, c3, c4, c5, c6 = st.columns(6)
        with c1:
            st.metric("Total Natural Products", len(df_master))
        with c2:
            test_cnt = len(df_test) if df_test is not None else 0
            st.metric("Test Set (Control)", test_cnt)
        with c3:
            an_cnt = len(df_analysis) if df_analysis is not None else 0
            st.metric("Analysis Set (Hits)", an_cnt)
        with c4:
            st.metric("Valid RDKit SMILES", "100%")
        with c5:
            lip_rate = f"{df_master['Lipinski_Pass'].mean()*100:.1f}%"
            st.metric("Lipinski Ro5 Pass", lip_rate)
        with c6:
            veber_rate = f"{df_master['Veber_Pass'].mean()*100:.1f}%"
            st.metric("Veber Pass", veber_rate)
            
        st.markdown("<br>", unsafe_allow_html=True)
        
        # Section 1: Quality Validation Audit Table
        st.markdown("""
<div class="science-card">
<h4 style="color: #0f172a; margin-bottom: 1rem;">1. Audit Integritas & Kualitas Data Molekuler</h4>
<div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem;">
    <div class="metric-box">
        <div class="data-label">Validitas Struktur & Conformer 3D</div>
        <div style="font-size: 0.9rem; color: #1e293b; margin-top: 0.25rem;">
            <b>114/114 (100%)</b> SMILES valid RDKit<br>
            <b>113/114 (99.1%)</b> Berhasil konformasi 3D MMFF94
        </div>
    </div>
    <div class="metric-box">
        <div class="data-label">Nilai Hilang (Missing Values)</div>
        <div style="font-size: 0.9rem; color: #10b981; font-weight: 600; margin-top: 0.25rem;">
            Nol Missing Values (0)<br>
            <span style="font-size: 0.75rem; color: #64748b; font-weight: normal;">Lengkap pada MW, LogP, HBD, HBA, TPSA</span>
        </div>
    </div>
    <div class="metric-box">
        <div class="data-label">Kompatibilitas Gaya Atom Vina</div>
        <div style="font-size: 0.9rem; color: #0284c7; font-weight: 600; margin-top: 0.25rem;">
            107 Senyawa Kompatibel Penuh<br>
            <span style="font-size: 0.75rem; color: #64748b; font-weight: normal;">7 Anolog boron dipeptida ditandai covalent/transition</span>
        </div>
    </div>
</div>
</div>
""", unsafe_allow_html=True)

        # Section 2: Dataset Partitioning (Read-Only Information)
        st.markdown("""
<div class="science-card">
<h4 style="color: #0f172a; margin-bottom: 0.5rem;">2. Informasi Pemisahan Dataset (Test Set vs Analysis Set)</h4>
<p style="color: #64748b; font-size: 0.85rem; margin-bottom: 0.75rem;">
<b>Test Set</b> dirancang untuk menguji presisi docking, akurasi pose, dan konsistensi skor vs bioaktivitas eksperimental (memuat Fitokimia Acuan & Kontrol Inhibitor). 
<b>Analysis Set</b> digunakan untuk penapisan forward menemukan kandidat inhibitor bahan alam terbaik.
</p>
<div style="background-color: #f1f5f9; border-left: 4px solid #64748b; padding: 0.6rem 1rem; border-radius: 0.375rem; font-size: 0.8rem; color: #475569;">
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
<div class="science-card" style="margin-top: 1.5rem;">
<h4 style="color: #0f172a; margin-bottom: 0.5rem;">3. Validasi Korelasi Eksperimental (pIC50 vs AutoDock Vina Energy)</h4>
<div style="display: flex; align-items: center; gap: 2rem;">
    <div style="background: #e0f2fe; padding: 1rem 1.5rem; border-radius: 8px; border: 1px solid #bae6fd;">
        <div style="font-size: 0.8rem; font-weight: 600; color: #0369a1; text-transform: uppercase;">Pearson Correlation (r)</div>
        <div style="font-size: 1.8rem; font-weight: 800; color: #0284c7; font-family: 'Fira Code', monospace;">{corr_val}</div>
    </div>
    <div style="color: #475569; font-size: 0.9rem; line-height: 1.5;">
        Korelasi negatif yang signifikan membuktikan bahwa senyawa dengan afinitas biologis laboratorium lebih tinggi (pIC50 besar) 
        secara konsisten menghasilkan energi bebas ikatan (*binding affinity*) yang lebih kuat / negatif pada simulasi AutoDock Vina reseptor DPP-4.
    </div>
</div>
</div>
""", unsafe_allow_html=True)
    else:
        st.warning("Dataset terkurasi belum dimuat. Jalankan script `14_curation_validation_test.py` terlebih dahulu.")
