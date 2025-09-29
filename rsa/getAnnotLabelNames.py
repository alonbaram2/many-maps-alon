import nibabel.freesurfer.io as fsio

# Path to your annotation file
annot_path = "/vols/Scratch/mgarvert/ManyMaps/imagingData/FS/fsaverage/label/lh.PALS_B12_Brodmann.annot"

# Load annotation file
labels, ctab, region_names = fsio.read_annot(annot_path)

# Print all region names
print("Regions in annotation file:")
for name in region_names:
    print(name.decode('utf-8'))