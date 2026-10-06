import os
import glob

log_files = glob.glob('results/02_docking/*.log')
log_files.sort()

output_md = '/Users/aanmuf/.gemini/antigravity-ide/brain/91949754-b51f-4d40-8c6b-6926bab5fc39/all_poses_analysis.md'

with open(output_md, 'w') as out:
    out.write("# Analisis Seluruh Pose (Mode 1-9) Docking\n\n")
    out.write("Dokumen ini berisi detail 9 mode konformasi (pose) yang dihasilkan oleh algoritma Genetik Vina untuk setiap ligan dan target protein.\n\n")
    
    for log_path in log_files:
        basename = os.path.basename(log_path)
        if basename.endswith('_out.log'): continue
        name = basename.replace('.log', '')
        parts = name.split('_')
        if len(parts) >= 2:
            chembl = parts[0]
            protein = parts[1]
        else:
            chembl = name
            protein = "Unknown"
            
        out.write(f"## {chembl} - Protein Target {protein}\n")
        out.write("| Mode | Affinity (kcal/mol) | RMSD l.b. | RMSD u.b. |\n")
        out.write("|------|---------------------|-----------|-----------|\n")
        
        with open(log_path, 'r') as f:
            lines = f.readlines()
            
        # find the table
        in_table = False
        for line in lines:
            if line.startswith('---+'):
                in_table = True
                continue
            if line.startswith('Writing output'):
                in_table = False
                break
            if in_table and line.strip():
                parts = line.strip().split()
                if len(parts) >= 4:
                    out.write(f"| {parts[0]} | {parts[1]} | {parts[2]} | {parts[3]} |\n")
                    
        out.write("\n")

print("Artifact created!")
