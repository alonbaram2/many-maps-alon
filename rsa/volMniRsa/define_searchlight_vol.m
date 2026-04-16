function define_searchlight_vol(rootData,subj,roi)
masksDir = fullfile(rootData,'masks');

% Define searchlights within a volumetric ROI. This is done in MNI space,
% as otherwise it's difficult later to match between grey matter of
% different subjects. 

% mask to find searchlights in. 
Vmask = spm_vol(fullfile(masksDir,[roi '.nii'])); % mask
Vmask.data = spm_read_vols(Vmask);
Vmask.sphere=[10, 100]; %[maxRadius, nVoxInSphere]

% File name of searchlight matlab structure
fname = [roi '.mat']; % 
slDir= fullfile(rootData,'rsa_alon',subj,'searchlight')
if ~exist(slDir,'dir') % create a folder for 1st-level results
    mkdir(slDir);
end

L = rsa.defineSearchlight({Vmask},Vmask);
save(fullfile(slDir,fname),'-struct','L'); 