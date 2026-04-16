#!/bin/bash

sub=$1
sess=$2

origScript=/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/designs/design_322.m
tmpOutScriptsDir=/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/alon/jalapeno
cat ${origScript} | sed "s/XXsubjIDXX/${sub}/g" | sed "s/XXsessionIDXX/${sess}/g" > ${tmpOutScriptsDir}/Subj_${sub}_session_${sess}.m;
jobID1=`fsl_sub -N sj${sub}_ses${sess} -T 480 matlab -nodisplay -nosplash \< ${tmpOutScriptsDir}/Subj_${sub}_session_${sess}.m`;
echo $jobID1


