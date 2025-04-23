#!/bin/bash

root=/home/fs0/mgarvert/scratch/ManyMaps/imagingData
for subj in {1..25}; do
    for sess in 1 2; do
        spmDir=$root/Subj_${subj}/session_${sess}/1stLevel/design_401_noSmooth
        betaFiles=$spmDir/beta*.nii
        # arbitrarily use run 1 folder for registration - anyway just using warps from structural image
        regDir=$root/Subj_${subj}/session_${sess}/run_1/preprocess_noSmooth.feat/reg
        mkdir $spmDir/MNI;
        for betaF in $betaFiles; do
            betaNum=${betaF: -8:4}
            jid=`fsl_sub -q short applywarp -i $spmDir/beta_${betaNum}.nii -o $spmDir/MNI/MNI_beta_${betaNum}.nii -r $regDir/standard.nii.gz -w $regDir/highres2standard_warp.nii.gz`
            echo Subj_$subj sess_${session} beta_$betaNum $jid
        done
    done
done