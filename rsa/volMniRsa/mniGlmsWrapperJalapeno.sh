scriptsDir=/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/alon/rsa/volMniRsa
glmDesignScript=/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/designs/design_407_MNI.m
module load MATLAB
module load fsl
for sub in {2..25}; do
    for sess in 1 2; do
        jobID0=`fsl_sub -N s${sub}s${sess}_warp -T 180 $scriptsDir/movePreproc2Mni.sh $sub $sess`;
        echo $sub $sess $jobID0
        cat ${glmDesignScript} | sed "s/XXsubjIDXX/${sub}/g" | sed "s/XXsessionIDXX/${sess}/g" > ${scriptsDir}/Subj_${sub}_session_${sess}.m;   
        jobID1=`fsl_sub -N s${sub}s${sess}_glm -T 180 -j $jobID0 matlab -nodisplay -nosplash \< ${scriptsDir}/Subj_${sub}_session_${sess}.m`;
        echo $sub $sess $jobID1
    done
done;