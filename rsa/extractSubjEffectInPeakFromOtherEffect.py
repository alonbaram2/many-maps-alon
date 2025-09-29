import os
import numpy as np
import nibabel as nib
from nibabel.freesurfer.io import read_label
import scipy.stats as stats
import matplotlib.pyplot as plt

# Extract subject-specific effects (in overlay_4d) at the peak of a group effect 
# (scalar_overlay) in a mask. Output to .txt with the MNI coords of the peak 
# (Note these coords are from a linear registration so not %100 accurate - use 
# mri_surf2vol for accurate volume-based MNI coords.

# === INPUTS ===

root = '/vols/Scratch/mgarvert/ManyMaps/imagingData'
hemi = 'rh'

mask_path = os.path.join(
    root, 'FS', 'fsaverage', 'label',
    f'{hemi}.PALS_B12_Brodmann_labels',
    f'{hemi}.Brodmann.24_32.label'
)

scalar_overlay_path = os.path.join(
    root, '2ndLevel', 'design_322_fsl_', 'session_diff', '04_rel_dist_switch',
    'rh_on_rh.spmT_rel_dist_switch_session_diff_0001_oppositeSign.mgh'
)

# === LOOP THROUGH SESSIONS AND EXTRACT PEAK VALUES ===

for sess in ['session_1', 'session_2', 'diff']:
    # === Paths ===
    overlay_4d_path = os.path.join(
        root,
        'rsa_alon',
        'allSubjStacked',
        'correlation',
        sess,
        'rh_on_rh',
        'distRel_bothMaps_xRun1324_smth5_rh_allSubj.mgh'
    )
    output_txt_path = os.path.join(
        root,
        'rsa_alon',
        'allSubjStacked',
        'correlation',
        sess,
        'rh_on_rh',
        'distRel_bothMaps_xRun1324_smth5_rh_allSubj_in_peakOf_spmT_rel_dist_switch_session_diff_0001_oppositeSign_inMask_BA24_32.txt'
    )

    # === LOAD DATA ===
    overlay_4d = nib.load(overlay_4d_path).get_fdata().squeeze()
    scalar_overlay = nib.load(scalar_overlay_path).get_fdata().squeeze()

    # === LOAD MASK ===
    if mask_path.endswith('.label'):
        mask_vertices = read_label(mask_path)
        mask = np.zeros_like(scalar_overlay, dtype=bool)
        mask[mask_vertices] = True
    else:
        mask = nib.load(mask_path).get_fdata().squeeze().astype(bool)

    # === FIND PEAK VERTEX ===
    masked_scalar = scalar_overlay.copy()
    masked_scalar[~mask] = -np.inf
    peak_vertex = np.argmax(masked_scalar)

    # === SUBJECT VALUES AT PEAK ===
    peak_values = overlay_4d[peak_vertex, :]

    # === SAVE RESULTS ===
    with open(output_txt_path, 'w') as f:
        for val in peak_values:
            f.write(f"{val:.6f}\n")

    print(f"✅ Output written to: {output_txt_path}")

# === STATS + PLOT ===

# Define base path pattern for reloading saved text files
txt_base = os.path.join(
    root, 'rsa_alon', 'allSubjStacked', 'correlation',
    '{sess}', 'rh_on_rh',
    'distRel_bothMaps_xRun1324_smth5_rh_allSubj_in_peakOf_spmT_rel_dist_switch_session_diff_0001_oppositeSign_inMask_BA24_32.txt'
)

# Load paired values
session1_vals = np.loadtxt(txt_base.format(sess='session_1'))
session2_vals = np.loadtxt(txt_base.format(sess='session_2'))

print(f"peak vertex number: {peak_vertex}")
print(f"\n🔍 mean session 1: {np.mean(session1_vals)}")
print(f"\n🔍 mean session 2: {np.mean(session2_vals)}")

# Perform one-sided paired t-test: session 2 > session 1
t_stat, p_val_two_sided = stats.ttest_rel(session2_vals, session1_vals)
p_val_one_sided = p_val_two_sided / 2 if t_stat > 0 else 1 - p_val_two_sided / 2

print(f"\n🔍 Paired one-sided t-test (session2 > session1):")
print(f"    t = {t_stat:.3f}, p = {p_val_one_sided:.5f}")

# Plot paired values
plt.figure(figsize=(6, 5))
plt.plot([1, 2], [session1_vals, session2_vals], 'o-', alpha=0.6)
plt.xticks([1, 2], ['Session 1', 'Session 2'])
plt.ylabel('Effect size at peak')
plt.title('Paired Values at Peak Vertex\n(session2 > session1)')
plt.grid(True)
plt.tight_layout()
plt.show()