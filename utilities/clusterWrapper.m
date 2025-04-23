rootData    = '/vols/Scratch/abaram/NN2020';

rootScripts = '/home/fs0/abaram/scripts/NN2020';
addpath(genpath(rootScripts));
spmPath = '/vols/Scratch/abaram/MATLAB/spm12';
addpath(genpath(spmPath));

sub='sub-00';
glm='GLM1';
roi='surf';
fwhm = 5;
iBlock = 999;
% assuming already ran getEvs.m for current subject

% % run SPM first-level GLM
% runSpmGlmEst(rootData,sub,'confounds')
% runRSA(rootData,sub,'GLM2','xRun','surf',fwhm)

% 
% subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09',...
%     'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
%     'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};
% 
% subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07',...
%     'sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
%     'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};



setenv('SUBJECTS_DIR','/vols/Scratch/abaram/NN2020/freesurfer-subjects/')
setenv('FREESURFER_HOME','/opt/fmrib/FreeSurfer_releases/6.0')


% runSpmGlmEst_MNI(rootData,sub,glm)
% runRSA(rootData,sub,glm,'xRun','lAccumbens20',fwhm)
% 

% runSpmGlmEst(rootData,sub,glm)

% 


spm_get_defaults('stats.fmri.hrf',[6,16,1,1,6,0,32])

% run SPM first-level GLM
runSpmGlmEst(rootData,sub,glm)

% calculate contrasts
runSpmContrast(rootData,sub,glm)

% warp contrasts to standard space and smooth
smoothAndWarpContrasts(rootData,sub,glm,fwhm)

% runRSA(rootData,sub,glm,'xRun','surf',fwhm)



% runRSA_withinCluster(rootData,sub,'GLM2','corrWithinBlock',roi,5,iBlock)
