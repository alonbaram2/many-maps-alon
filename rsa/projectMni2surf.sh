
root=/vols/Scratch/mgarvert/ManyMaps/imagingData
export SUBJECTS_DIR=$root/FS
inDir=$root/2ndLevel/design_322_fsl_/session_diff/36_rel_dist_stayAndSwitch
fname=spmT_rel_dist_stayAndSwitch_session_diff_0001

mri_vol2surf \
  --mov $inDir/$fname.nii \
  --regheader fsaverage \
  --hemi lh \
  --projfrac 0.5 \
  --interp trilinear \
  --out $inDir/lh.$fname.mgh

mri_vol2surf \
  --mov $inDir/$fname.nii \
  --regheader fsaverage \
  --hemi rh \
  --projfrac 0.5 \
  --interp trilinear \
  --out $inDir/rh_on_rh.$fname.mgh




#####

  mris_apply_reg --src /vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/allSubjStacked/correlation/diff/distRel_sameMap_xRun1324_smth5_rh_allSubj.mgh \
                 --trg /vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/allSubjStacked/correlation/diff/rh_on_rh/distRel_sameMap_xRun1324_smth5_rh_allSubj.mgh \
                 --streg /vols/Scratch/mgarvert/ManyMaps/imagingData/FS/fsaverage/surf/lh.sphere.left_right \
                         /vols/Scratch/mgarvert/ManyMaps/imagingData/FS/fsaverage/surf/rh.sphere.left_right
