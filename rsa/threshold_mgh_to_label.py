import os
import nibabel as nib
import numpy as np

# === Configuration ===
stats_dir = '/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/groupStats/correlation/Wilcoxon'
mgh_path = os.path.join(stats_dir, 'P_diff_distRel_bothMaps_xRun1324_smth5_rh_on_rh.mgh')
output_label = os.path.join(stats_dir, 'P_diff_distRel_bothMaps_xRun1324_smth5_rh_on_rh_thrsh0p99.label')
output_mgh = os.path.join(stats_dir, 'P_diff_distRel_bothMaps_xRun1324_smth5_rh_on_rh_thrsh0p99.mgh')
threshold = 0.99

# === Load .mgh file ===
mgh = nib.load(mgh_path)
data = mgh.get_fdata().squeeze()  # Shape: (n_vertices,) or (n_vertices, 1, 1, 1)
affine = mgh.affine
header = mgh.header

# === Find vertices above threshold ===
vertices = np.where(data > threshold)[0]

# === Write FreeSurfer .label format ===
with open(output_label, 'w') as f:
    f.write('#!ascii label file\n')
    f.write(f'{len(vertices)}\n')
    for v in vertices:
        f.write(f'{v} 0.0 0.0 0.0 1.0\n')

# === Create and save thresholded .mgh file ===
thresholded_data = np.zeros_like(data)
thresholded_data[vertices] = data[vertices]

# Expand dimensions back to original shape for saving
thresholded_data = thresholded_data[:, np.newaxis, np.newaxis]  # (n_vertices, 1, 1)

new_img = nib.MGHImage(thresholded_data, affine=affine, header=header)
nib.save(new_img, output_mgh)

print(f"Saved label with {len(vertices)} vertices to: {output_label}")
print(f"Saved thresholded .mgh file to: {output_mgh}")