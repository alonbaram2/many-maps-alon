
#!/bin/bash

sub=$1
glm=$2
roi=$3
queue=$4
iBlock=$5

scriptsDir=/home/fs0/abaram/scripts/NN2020/utilities
cat ${scriptsDir}/clusterWrapper.m | sed "s/sub-00/${sub}/g" | sed "s/GLM1/${glm}/g" | sed "s/surf/${roi}/g" | sed "s/999/${iBlock}/g"> ${scriptsDir}/${sub}_${glm}_${roi}_${iBlock}.m;
jobID1=`fsl_sub -N ${sub}${iBlock}${glm}${roi}${queue} -q ${queue}.q matlab -nodisplay -nosplash \< ${scriptsDir}/${sub}_${glm}_${roi}_${iBlock}.m`;
echo $jobID1




# for sub in sub-01 sub-02 sub-03 sub-04 sub-05 sub-06 sub-07 sub-08 sub-09 sub-10 sub-11 sub-12 sub-13 sub-14 sub-15 sub-16 sub-17 sub-18 sub-19 sub-20 sub-21 sub-22 sub-23 sub-24 sub-25 sub-26 sub-27 sub-28; do /home/fs0/abaram/scripts/NN2020/utilities/glmCluster.sh
# for sub in sub-01 sub-02 sub-03 sub-04 sub-05 sub-06 sub-07 sub-11 sub-12 sub-13 sub-14 sub-15 sub-16 sub-17 sub-18 sub-19 sub-20 sub-21 sub-22 sub-23 sub-24 sub-25 sub-26 sub-27 sub-28; do /home/fs0/abaram/scripts/NN2020/utilities/glmCluster.sh 
