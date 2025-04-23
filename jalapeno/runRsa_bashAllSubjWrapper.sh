#!/bin/bash

# for subj in  {2..25}
    # for sess in 1 2

sub=$1
sess=$2

# origScript=/home/fs0/mgarvert/scratch/ManyMaps/imagingData/scripts/alon/jalapeno/runRsa_jalapenoSubjWrapper.m
origScript=/vols/Scratch/abaram/ManyMaps/jalapenoTmp/runRsa_jalapenoSubjWrapper.m
# tmpOutScriptsDir=/home/fs0/mgarvert/scratch/ManyMaps/imagingData/scripts/alon/jalapeno
tmpOutScriptsDir=/vols/Scratch/abaram/ManyMaps/jalapenoTmp
cat ${origScript} | sed "s/XXsubjIDXX/Subj_${sub}/g" | sed "s/XXsessionIDXX/session_${sess}/g" > ${tmpOutScriptsDir}/Subj_${sub}_session_${sess}.m;
jobID1=`fsl_sub -T 1600 -R 16 -N sj${sub}_ses${sess} matlab -nodisplay -nosplash \< ${tmpOutScriptsDir}/Subj_${sub}_session_${sess}.m`;
echo $jobID1
