sub=$1
sess=$2

root=/vols/Scratch/mgarvert/ManyMaps/imagingData
for run in run_1 run_2 run_3 run_4; do
    cd $root/Subj_${sub}/session_${sess}/${run}/preprocess_noSmooth.feat
    pwd
    applywarp \
    --in=filtered_func_data.nii.gz \
    --ref=reg/standard \
    --warp=reg/highres2standard_warp \
    --premat=reg/example_func2highres.mat \
    --out=filtered_func_MNI_2mm
done;

    