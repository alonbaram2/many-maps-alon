
scriptsDir = '/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/alon';
addpath(genpath(scriptsDir));
spmPath = '/vols/Scratch/abaram/MATLAB/spm12';
addpath(spmPath);

root = '/vols/Scratch/mgarvert/ManyMaps/imagingData';
% subj = 'XXsubjIDXX';
% sess = 'XXsessionIDXX';
subs = {'1','2','3','4','5','6','7','8','9','11','12','13','14',...
    '16','19','20','21','22','23','24','25'};
for iSub=1:length(subs)
    subj= ['Subj_' subs{iSub}]
    for iSess = 1:2
        sess = ['session_' int2str(iSess)]
        spmDirMni = fullfile(root,subj,sess,'1stLevel','design_401_noSmooth','MNI');

        nTrials = 136;
        nRegsPerRun = 144; % including nuisance regressors - choiceTrials, buttonPress, 6*motionParams
        nRuns = 4;

        % load mask
        mask =  'vmPFC_alon_2mm'%'juelich_V4_thr20';
        mask_file = fullfile(root,'masks',[mask, '.nii']);
        mask_nii = spm_vol(mask_file);
        mask_data = spm_read_vols(mask_nii);
        % Find the indices of non-zero voxels in the mask
        [mask_i, mask_j, mask_k] = ind2sub(size(mask_data), find(mask_data ~= 0));
        nVox = numel(mask_i);

        decodingDataDir = fullfile(root,'decoding','data',mask);
        mkdir(decodingDataDir)

        D = nan(nVox,nTrials,nRuns);

        for iRun=1:nRuns
            for iTrial = 1:nTrials
                betaNum = (iRun-1)*nRegsPerRun + iTrial;
                betaNumStr = sprintf('%0*d', 4, betaNum);
                betaFile = fullfile(spmDirMni,['MNI_beta_' betaNumStr '.nii']);
                beta_nii = spm_vol(betaFile);
                beta_data = spm_read_vols(beta_nii);
                D(:,iTrial,iRun) = beta_data(mask_data ~= 0);
            end
        end

        % get rid of voxls that have NaNs in any run
        [I,J,K] = find(isnan(D));
        disp(sprintf('num NaN vox: %d', length(I)))
        D(I,:,:) = [];

        save(fullfile(decodingDataDir,[subj '_' sess]),'D');
        %
    end
end
