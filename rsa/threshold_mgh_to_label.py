import os
import nibabel as nib
import numpy as np
import shutil

# === Configuration ===
stats_dir = '/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/groupStats/correlation/Wilcoxon'
mask_dir = '/vols/Scratch/mgarvert/ManyMaps/imagingData/masks/fsaverage'

mgh_path = os.path.join(stats_dir, 'P_diff_distRel_diffMaps_xRun1324_smth5_rh.mgh')
output_fname = 'rh.diff_distRel_diffMaps_xRun1324_smth5_thrsh0p99'
output_label = os.path.join(stats_dir, f"{output_fname}.label")
output_mgh = os.path.join(stats_dir, f"{output_fname}.mgh")
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

# copy to masks folder
shutil.copy(output_label, mask_dir)
shutil.copy(output_mgh, mask_dir)