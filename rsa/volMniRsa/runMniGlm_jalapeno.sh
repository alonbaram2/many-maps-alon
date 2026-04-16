#!/bin/bash

sub=$1
sess=$2

glmDesignScript=/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/designs/design_307_MNI.m
scriptsDir=/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/alon/rsa/volMniRsa


cat ${glmDesignScript} | sed "s/XXsubjIDXX/${sub}/g" | sed "s/XXsessionIDXX/${sess}/g" > ${scriptsDir}/Subj_${sub}_session_${sess}.m;
jobID1=`fsl_sub -N s${sub}s${sess}_glm -T 480 matlab -nodisplay -nosplash \< ${scriptsDir}/Subj_${sub}_session_${sess}.m`;
echo $jobID1


