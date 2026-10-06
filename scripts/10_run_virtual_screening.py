import os
import subprocess
import glob
import pandas as pd
import time

def run_virtual_screening_all():
    print("=== Memulai Multi-Target Virtual Screening (3G0B, 5T4B) ===")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    # Path direktori
    receptor_dir = os.path.join(project_root, 'data', '05_receptors_pdbqt')
    ligand_dir = os.path.join(project_root, 'data', '04_ligands_pdbqt')
    results_dir = os.path.join(project_root, 'results', '02_docking')
    
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        
    ligand_files = glob.glob(os.path.join(ligand_dir, '*.pdbqt'))
    # Skip _out files if any exist
    ligand_files = [f for f in ligand_files if not f.endswith('_out.pdbqt')]
    
    print(f"Ditemukan {len(ligand_files)} ligan siap ditambatkan.")
    
    # Parameter Grid Box untuk masing-masing protein
    targets = {
        '3G0B': {'center': (42.05, 34.29, 14.62), 'size': (18, 18, 18)},
        '5T4B': {'center': (37.61, 49.92, 40.33), 'size': (18, 16, 23)}
    }
    
    vina_exec = os.path.join(script_dir, 'vina')
    all_results = []
    
    start_time = time.time()
    
    for protein, params in targets.items():
        receptor_path = os.path.join(receptor_dir, f"{protein}.pdbqt")
        if not os.path.exists(receptor_path):
            print(f"Melewati {protein}: File reseptor tidak ditemukan di {receptor_path}")
            continue
            
        print(f"\n>> Menambatkan ligan ke protein {protein}...")
        c_x, c_y, c_z = params['center']
        s_x, s_y, s_z = params['size']
        
        for i, ligand_path in enumerate(ligand_files):
            ligand_name = os.path.basename(ligand_path).replace('.pdbqt', '')
            out_pdbqt = os.path.join(results_dir, f"{ligand_name}_{protein}_out.pdbqt")
            log_file = os.path.join(results_dir, f"{ligand_name}_{protein}.log")
            
            # Gunakan exhaustiveness 1 untuk screening cepat (demo purpose) 
            cmd = [
                vina_exec,
                '--receptor', receptor_path,
                '--ligand', ligand_path,
                '--center_x', str(c_x),
                '--center_y', str(c_y),
                '--center_z', str(c_z),
                '--size_x', str(s_x),
                '--size_y', str(s_y),
                '--size_z', str(s_z),
                '--out', out_pdbqt,
                '--exhaustiveness', '1' 
            ]
            
            try:
                process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                if process.returncode != 0:
                    continue
                    
                stdout_text = process.stdout.decode('utf-8')
                with open(log_file, 'w') as f:
                    f.write(stdout_text)
                    
                affinity = None
                rmsd_lb = None
                rmsd_ub = None
                for line in stdout_text.split('\n'):
                    if line.strip().startswith('1 '):
                        parts = line.strip().split()
                        if len(parts) >= 4:
                            affinity = float(parts[1])
                            rmsd_lb = float(parts[2])
                            rmsd_ub = float(parts[3])
                        break
                                
                if affinity is not None:
                    all_results.append({
                        'Protein': protein,
                        'ChEMBL_ID': ligand_name,
                        'Binding_Affinity_kcal_mol': affinity,
                        'RMSD_lb': rmsd_lb,
                        'RMSD_ub': rmsd_ub
                    })
                    print(f"   [{i+1}/{len(ligand_files)}] {ligand_name} -> {affinity} kcal/mol")
                    
            except Exception as e:
                pass
                
    if all_results:
        df_results = pd.DataFrame(all_results)
        
        # Simpan hasil komprehensif
        output_csv = os.path.join(results_dir, 'multi_target_screening_results.csv')
        df_results.to_csv(output_csv, index=False)
        
        # Buat tabel pivot untuk perbandingan yang lebih mudah dibaca
        pivot_df = df_results.pivot(index='ChEMBL_ID', columns='Protein', values='Binding_Affinity_kcal_mol')
        # Menghitung rata-rata afinitas di 3 protein
        pivot_df['Average_Affinity'] = pivot_df.mean(axis=1)
        pivot_df = pivot_df.sort_values(by='Average_Affinity')
        
        pivot_csv = os.path.join(results_dir, 'multi_target_summary_pivot.csv')
        pivot_df.to_csv(pivot_csv)
        
        duration = time.time() - start_time
        print("\n" + "="*50)
        print(f"MULTI-TARGET VIRTUAL SCREENING SELESAI ({duration:.1f} detik)!")
        print(f"Hasil detail: {output_csv}")
        print(f"Hasil pivot (Rangkuman): {pivot_csv}")
        print("="*50)
        print("\n--- TOP 3 KANDIDAT OBAT (Rata-Rata Terbaik di 3 Protein) ---")
        print(pivot_df.head(3))

if __name__ == "__main__":
    run_virtual_screening_all()
