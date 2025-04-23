function run_rsa(root, subj, distType,session)

% distType (str): correlation / euclidean (inputed to pdist)
% session(str): session_1 or session_2

nConditions = 34; % 17 objects, 2 maps
nElemFinal = nConditions*(nConditions-1)/2;
nRuns=4;

glm = 'design_307_fsl__noSmooth';
dataRdmDir= fullfile(root,'rsa_alon',subj,'dataRdms', distType, session);
L = load(fullfile(root,'rsa_alon',subj,'searchlight','WhBr_surf_r10_v100'));
spmFile = fullfile(root,subj,session,'1stLevel',glm,'SPM.mat'); % directory for SPM.mat file
load(spmFile);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Running searchlight analysis ' subj, ' ' session])
disp(['Data RDMs will be stored here: ',dataRdmDir])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if ~exist('dataRdmDir','dir')
    mkdir (dataRdmDir)
end
cd(dataRdmDir);

regInds4RSA = 1:nConditions;
condition = false(1,size(SPM.xX.xKXs.X,2));
for iRun = 1:nRuns
     condition(SPM.Sess(iRun).col(regInds4RSA)) = true;
end
% number of elements in the exapnded RDM, split into runs and
% before averaging the distance between conditions i,j in runs k,l by "folding"
% the matrix.
% Note that here we find both within and across-run distances, but we will
% only use the across runs distances (this is an unefficient way of
% computing this but it's easier to undrstand what's going on).
nCondSplitRunsBeforeCollapse = sum(condition); % 34*4 = 136
nDistsSplitRunsBeforeCollapse = nCondSplitRunsBeforeCollapse*(nCondSplitRunsBeforeCollapse-1)/2;

% file name for RDM
outFiles={};
for k=1:nDistsSplitRunsBeforeCollapse % k is 4th dim of RDM NIFTI.
    outFiles{end+1}=sprintf('dist_%s.nii,%d',distType,k);
end

        
% run the analysis
analysisName = distType;
rsa.runSearchlight(L,SPM.xY.VY,outFiles,analysisName,@calcRDMs,'optionalParams',{SPM,condition,distType});
