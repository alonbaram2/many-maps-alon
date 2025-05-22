function groupLevelPermTests(rootData,clusterThresh,nPerm,maskName,sl,pathIn)

% sl: either 'surf' or 'vol'. If it's a surface map already projected to
% MNI volume, it should be 'vol'. 

% rootData=/vols/Scratch/mgarvert/ManyMaps/imagingData;
% pathIn = '/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/allSubjStacked/correlation/diff/distRel_diffMaps_xRun1324_smth5_MNI.nii'; 

% make sure palm is installed an in path.
palmDir = '/vols/Scratch/abaram/MATLAB/PALM/';
addpath(genpath(palmDir));


% mask to run the permutation tests in and perform multiple comparisons in.
% mask is intentionally of left hemi (lh) as the right hemi (rh) was
% previously registered to the left, so it now has the lh indeces. This
% happened inside searchlightDefinitionSurfaceWrapper.m


if strcmp(sl,'surf')
    % outDir = fullfile(rootData,'rsa_alon','groupStats','correlation','perm',maskName);
    % mkdir(outDir);
    % maskFile = fullfile(rootData,'masks','fsaverage',['lh.' maskName '.mgh']);
    % surface = 'pial';
    % surfDir  = fullfile(rootData,'freesurfer-subjects','fsaverage','surf');
    % surfFile = fullfile([surfDir '/lh.' surface]); % this is intentionally always lh because rh data was flipped to be with lh indeces.
    % 
    % if ~exist(maskFile,'file')
    %     maskLabelFName = [maskFile(1:end-4) '.label'];
    %     % change format of mask from .label to .mgh as that's what PALM
    %     % wants
    %     system (['mri_label2label --s fsaverage --srclabel ' maskLabelFName ...
    %         ' --trglabel ' maskLabelFName ' --regmethod surface --hemi lh --outmask ' maskFile]);
    % end
    % for iAnalysis = 1:length(fname) % better to parallelise this
    %     outFile = fullfile(outDir,[fname(iAnalysis).name(1:end-4) '_nPerm' nPerm '_clstrTh' clusterThresh]);        
    %     outFile = strrep(outFile,'.','p'); % switch decimal point to 'p', as PALM doesn't like points in filenames.        
    %     str = ['palm -i ' fname(iAnalysis).name ' -s ' surfFile ...
    %         ' -n ' nPerm ' -o ' outFile ' -ise -save1-p -m ' maskFile];
    %     if ~strcmp(clusterThresh,'None')
    %         str = [str ' -C ' clusterThresh ' -Cstat mass' ];
    %     end
    %     eval(str)
    % end    
else % volumetric
    outDir = fullfile(rootData,'rsa_alon','groupStats','correlation','perm',maskName);
    mkdir(outDir);
    maskFile = fullfile(rootData,'masks',[maskName '.nii']);
    outFile = fullfile(outDir,[pathIn(1:end-4) '_nPerm' nPerm '_clstrTh' clusterThresh]);
    outFile = strrep(outFile,'.','p'); % switch decimal point to 'p', as PALM doesn't like points in filenames.
    str = ['palm -i ' pathIn ' -n ' nPerm ' -o ' outFile ' -ise -save1-p -m ' maskFile];
    if ~strcmp(clusterThresh,'None')
        str = [str ' -C ' clusterThresh ' -Cstat mass' ];
    end
    eval(str)
end


