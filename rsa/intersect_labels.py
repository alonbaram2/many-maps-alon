def load_label_vertices(label_file):
    with open(label_file, 'r') as f:
        lines = f.readlines()

    # Skip comments and header
    data_lines = [line for line in lines if not line.startswith('#')]
    return set(int(line.strip().split()[0]) for line in data_lines)

def write_label(vertices, output_path):
    with open(output_path, 'w') as f:
        f.write('#!ascii label file\n')
        f.write(f'{len(vertices)}\n')
        for v in sorted(vertices):
            f.write(f'{v} 0.0 0.0 0.0 1.0\n')

# === Edit these paths ===
label1_path = '/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/groupStats/correlation/Wilcoxon/P_diff_distRel_bothMaps_xRun1324_smth5_rh_on_rh_thrsh0p99.label'  # The label you made from the MGH file
label2_path = '/vols/Scratch/mgarvert/ManyMaps/imagingData/FS/fsaverage/label/rh.PALS_B12_Brodmann_labels/rh.Brodmann.24_32.label' # this should usually be lh - indeces are all in lh
output_label = '/vols/Scratch/mgarvert/ManyMaps/imagingData/masks/fsaverage/rh_on_rh.P_diff_distRel_bothMaps_xRun1324_smth5_thrsh0p99_intersect_BA24_32.label'

# === Intersect ===
vertices1 = load_label_vertices(label1_path)
vertices2 = load_label_vertices(label2_path)

intersected = vertices1 & vertices2

write_label(intersected, output_label)

print(f"Saved intersected label with {len(intersected)} vertices to: {output_label}")

print(f"to run: mri_label2label --s fsaverage --srclabel {output_label} --trglabel {output_label} --regmethod surface --hemi lh --outmask {output_label}.mgh")
