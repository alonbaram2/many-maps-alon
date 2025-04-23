#!/bin/bash

# for subj in  {1..25}
    # for sess in 1 2

sub=$1
sess=$2

origScript=/home/fs0/mgarvert/scratch/ManyMaps/imagingData/scripts/alon/decoding_getData.m
tmpOutScriptsDir=/home/fs0/mgarvert/scratch/ManyMaps/imagingData/scripts/alon/jalapeno
cat ${origScript} | sed "s/XXsubjIDXX/Subj_${sub}/g" | sed "s/XXsessionIDXX/session_${sess}/g" > ${tmpOutScriptsDir}/Subj_${sub}_session_${sess}.m;
jobID1=`fsl_sub -N sj${sub}_ses${sess} -q short.q matlab -nodisplay -nosplash \< ${tmpOutScriptsDir}/Subj_${sub}_session_${sess}.m`;
echo $jobID1