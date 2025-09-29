import numpy as np
import nibabel as nib
from nibabel.freesurfer.io import read_label

# Paths
overlay_path = "/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/allSubjStacked/correlation/diff/rh_on_rh/distRel_sameMap_xRun1324_smth5_rh_allSubj.mgh"
label_path = "/vols/Scratch/mgarvert/ManyMaps/imagingData/masks/fsaverage/rh_on_rh.diff_distRel_diffMaps_xRun1324_smth5_thrsh0p95_intersect_BA24_32.label"
output_path = "/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/allSubjStacked/correlation/diff/rh_on_rh/distRel_sameMap_xRun1324_smth5_rh_" \
"subjectMeansInMask_diff_distRel_diffMaps_xRun1324_smth5_thrsh0p95_intersect_BA24_32.txt"

# Load overlay (shape: 163842 x 1 x 1 x 22)
overlay_img = nib.load(overlay_path)
overlay_data = overlay_img.get_fdata().squeeze()  # Resulting shape: 163842 x 22

# Load label (vertex indices)
label_vertices = read_label(label_path)

# Extract values from the mask for each subject and compute the mean
means = []
for subj in range(overlay_data.shape[1]):
    subj_data = overlay_data[:, subj]
    masked_data = subj_data[label_vertices]
    mean_val = np.mean(masked_data)
    means.append(mean_val)

# Save to text file
np.savetxt(output_path, means, fmt="%.6f")
print(f"Saved subject means to: {output_path}")


