function groupLevelPermTests(rootData,clusterThresh,nPerm,maskName,sl,pathIn)

% sl: either 'surf' or 'vol'. If it's a surface map already projected to
% MNI volume, it should be 'vol'. 

% rootData=/vols/Scratch/mgarvert/ManyMaps/imagingData;
% pathIn = '/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/allSubjStacked/correlation/diff/distRel_diffMaps_xRun1324_smth5_MNI.nii'; 

% pathIn = '/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/allSubjStacked/correlation/diff/distRel_sameMap_xRun1324_smth5_rh_allSubj.mgh';
% maskName for surface blobs: e.g. rh.diff_distRel_diffMaps_xRun1324_smth5_thrsh0p95_intersect_BA24_32

% make sure palm is installed an in path.
palmDir = '/vols/Scratch/abaram/MATLAB/PALM/';
addpath(genpath(palmDir));


% mask to run the permutation tests in and perform multiple comparisons in.
% mask is intentionally of left hemi (lh) as the right hemi (rh) was
% previously registered to the left, so it now has the lh indeces. This
% happened inside searchlightDefinitionSurfaceWrapper.m

[~, pathInFileName, ~] = fileparts(pathIn);


if strcmp(sl,'surf')
    hemi = input('hemisphere? ''lh'' or ''rh''');
    outDir = fullfile(rootData,'rsa_alon','groupStats','correlation','perm',maskName);
    mkdir(outDir);
    maskFile = fullfile(rootData,'masks','fsaverage',[maskName '.mgh']);
    surface = 'pial';
    surfDir  = fullfile(rootData,'FS','fsaverage','surf');
    surfFile = fullfile([surfDir '/' hemi '.' surface]); % this is intentionally always lh because rh data was flipped to be with lh indeces.

    if ~exist(maskFile,'file')
        maskLabelFName = [maskFile(1:end-4) '.label'];
        % change format of mask from .label to .mgh as that's what PALM
        % wants. If you have trouble running this on the cluster (because
        % it can't find the command mri_label2label) just copy the string
        % and run it in terminal
            system (['mri_label2label --s fsaverage --srclabel ' maskLabelFName ...
                ' --trglabel ' maskLabelFName ' --regmethod surface --hemi ' hemi ' --outmask ' maskFile]);
    end
    outFile = fullfile(outDir,[pathInFileName '_nPerm' nPerm '_clstrTh' clusterThresh]);
%    outFile = strrep(outFile,'.','p'); % switch decimal point to 'p', as PALM doesn't like points in filenames.
    str = ['palm -i ' pathIn ' -s ' surfFile ...
        ' -n ' nPerm ' -o ' outFile ' -ise -save1-p -m ' maskFile];
    if ~strcmp(clusterThresh,'None')
        str = [str ' -C ' clusterThresh ' -Cstat extent' ];
    end
        eval(str)
        
else % volumetric
    outDir = fullfile(rootData,'rsa_alon','groupStats','correlation','perm',maskName);
    mkdir(outDir);
    maskFile = fullfile(rootData,'masks',[maskName '.nii']);
    outFile = fullfile(outDir,[pathInFileName '_nPerm' nPerm '_clstrTh' clusterThresh]);
    outFile = strrep(outFile,'.','p'); % switch decimal point to 'p', as PALM doesn't like points in filenames.
    str = ['palm -i ' pathIn ' -n ' nPerm ' -o ' outFile ' -ise -save1-p -m ' maskFile];
    if ~strcmp(clusterThresh,'None')
        str = [str ' -C ' clusterThresh ' -Cstat extent' ];
    end
    eval(str)
end


