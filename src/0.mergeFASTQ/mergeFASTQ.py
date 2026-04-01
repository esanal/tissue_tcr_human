import pandas as pd
import os
import ast
import subprocess

data_main_dir = "/home/erdem/Hosts/vacuole1/former-NOBINFBACKUP/human_tissue_tcr/data"
output_dir = "/home/erdem/Hosts/vacuole1/former-NOBINFBACKUP/human_tissue_tcr/data/combined_fastq"

# Generates metadata from the excel with sample names
metadata = pd.read_csv("../../metadata/metadata_merged.tsv", sep="\t")

metadata["gs.id"] = metadata["gs.id"].apply(ast.literal_eval)

# Index the directory once for speed
all_files = os.listdir(data_main_dir)

# Iterate through rows
for index, row in metadata.iterrows():
    # Capture all metadata levels
    indiv = str(row['individual']).replace(" ", "_")
    organ = str(row['organ']).replace(" ", "_")
    subset = str(row['subset']).replace(" ", "_")
    ids = row['gs.id']
    
    # Construct a clean base name for the output
    base_name = f"{indiv}_{organ}_{subset}"
    
    print(f"Row {index}: Processing {base_name} ({len(ids)} IDs)")

    for direction in ['R1', 'R2']:
        matching_files = []
        
        # Find all files matching the IDs for this direction
        for sample_id in ids:
            matches = [os.path.join(data_main_dir, f) for f in all_files 
                       if str(sample_id) in f and direction in f and f.endswith(".fastq.gz")]
            matching_files.extend(matches)
        
        # Deduplicate matches
        matching_files = list(set(matching_files))
        
        if not matching_files:
            print(f"  [!] No {direction} files found for IDs: {ids}")
            continue

        output_file = os.path.join(output_dir, f"{base_name}_{direction}.fastq.gz")

        # --- LOGIC: Concatenate or Copy ---
        if len(matching_files) > 1:
            # Multiple IDs: Concatenate
            inputs = " ".join(f"'{f}'" for f in matching_files) # Quoting paths for safety
            print(f"  [Merge]\n {"\n".join(matching_files)} files\n ->\n {output_file}")
            subprocess.run(f"cat {inputs} > '{output_file}'", shell=True, check=True)
        
        else:
            # Only one ID: Copy
            source_file = matching_files[0]
            print(f"  [Copy] Single file -> {output_file}")
            subprocess.run(f"cp '{source_file}' '{output_file}'", shell=True, check=True)

print("\nProcess complete. Merged files are in:", output_dir)
