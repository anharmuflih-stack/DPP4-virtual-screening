import os
import glob

base_dir = '/Users/aanmuf/drug_discovery/dpp4_project'

# 1. Rename Scripts
scripts_map = {
    'install_docking_tools.sh': '00_setup_environment.sh',
    'fetch_chembl_data.py': '01_fetch_chembl_data.py',
    'calculate_descriptors.py': '02_calculate_descriptors.py',
    'filter_lipinski.py': '03_filter_lipinski.py',
    'filter_natural_products.py': '04_filter_natural_products.py',
    'download_protein.py': '05_download_protein.py',
    'find_grid_box.py': '06_find_grid_box.py',
    'prepare_docking_files.py': '07_prepare_docking_files.py',
    'validate_docking.py': '08_validate_docking.py',
    'run_vina_all.py': '09_run_virtual_screening.py'
}

for old, new in scripts_map.items():
    old_p = os.path.join(base_dir, 'scripts', old)
    new_p = os.path.join(base_dir, 'scripts', new)
    if os.path.exists(old_p):
        os.rename(old_p, new_p)

if os.path.exists(os.path.join(base_dir, 'scripts', 'run_vina.py')):
    os.remove(os.path.join(base_dir, 'scripts', 'run_vina.py'))

# 2. Rename Data and Results Subfolders
folder_map = {
    'data/raw': 'data/01_raw',
    'data/processed': 'data/02_processed',
    'data/pdb': 'data/03_pdb',
    'data/ligands_pdbqt': 'data/04_ligands_pdbqt',
    'data/receptors_pdbqt': 'data/05_receptors_pdbqt',
    'results/validation': 'results/01_validation',
    'results/docking': 'results/02_docking'
}

for old, new in folder_map.items():
    old_p = os.path.join(base_dir, old)
    new_p = os.path.join(base_dir, new)
    if os.path.exists(old_p):
        os.rename(old_p, new_p)

# 3. Update paths inside all python scripts
path_replacements = {
    "'data', 'raw'": "'data', '01_raw'",
    "'data', 'processed'": "'data', '02_processed'",
    "'data', 'pdb'": "'data', '03_pdb'",
    "'data', 'ligands_pdbqt'": "'data', '04_ligands_pdbqt'",
    "'data', 'receptors_pdbqt'": "'data', '05_receptors_pdbqt'",
    "'results', 'validation'": "'results', '01_validation'",
    "'results', 'docking'": "'results', '02_docking'",
    "data/raw": "data/01_raw",
    "data/processed": "data/02_processed",
    "data/pdb": "data/03_pdb",
    "data/ligands_pdbqt": "data/04_ligands_pdbqt",
    "data/receptors_pdbqt": "data/05_receptors_pdbqt",
    "results/validation": "results/01_validation",
    "results/docking": "results/02_docking",
}

for py_file in glob.glob(os.path.join(base_dir, 'scripts', '*.py')):
    with open(py_file, 'r') as f:
        content = f.read()
    for old, new in path_replacements.items():
        content = content.replace(old, new)
    with open(py_file, 'w') as f:
        f.write(content)

# Update logbook too
logbook_path = os.path.join(base_dir, 'docs', 'logbook.ipynb')
if os.path.exists(logbook_path):
    with open(logbook_path, 'r') as f:
        content = f.read()
    for old, new in path_replacements.items():
        content = content.replace(old, new)
    with open(logbook_path, 'w') as f:
        f.write(content)

print("Penamaan ulang dan pembaruan path selesai!")
