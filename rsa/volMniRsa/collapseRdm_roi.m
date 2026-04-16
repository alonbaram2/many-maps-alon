function collapseRdm_roi(rootData,subj,session, roi)

fsl_dir = getenv('FSLDIR');
if isempty(fsl_dir)
    setenv('FSLDIR','/cvmfs/fsl.fmrib.ox.ac.uk/el9/fsl/6.0.7.16');
end

nRun = 4;
nElem = 34;
allRunsSquareMat = zeros(nRun*nElem);


%%  load the expanded (nRuns*nElem x nRuns*Elem, in squareform vec format) RDM
rdmDir = fullfile(rootData,'rsa_alon',subj,'dataRdms', roi, session);
dataFile = fullfile(rdmDir,['dist_correlation.nii']);
outputFile = fullfile(rdmDir,['dist_correlation_xRunCollapsed.mat']);
collapsedRDM = nan(nElem*(nElem-1)/2 + nElem,1);

% get NIFTI header info
V = niftiinfo(dataFile);
% read NIFTI array
dist_nii = niftiread(dataFile);

mask = any(~isnan(dist_nii), 4);   % 3D logical mask
[X,Y,Z] = ind2sub(size(mask), find(mask));
x = X(1); y = Y(1); z = Z(1); % it doesn't matter which voxel from the mask we're choosing. 
rdm_dist = squeeze(dist_nii(x,y,z,:)); % 
clear dist_nii

%% get the indeces
nRunPairs = (nRun)*(nRun-1)/2;  %6
runPairsToUse = 1:6; %  all run pairs
indsVec = false(nRunPairs,nElem*nRun*(nElem*nRun-1)/2); % indeces for each pair of runs
runPairCounter=0;
for iRun=1:nRun
    for jRun=1:nRun
        if iRun>jRun %lower block triangle (not including same run)
            runPairCounter = runPairCounter +1;
            M = allRunsSquareMat;
            M((iRun-1)*nElem+1:iRun*nElem,(jRun-1)*nElem+1:jRun*nElem)=1;
            % make symmetric
            M((jRun-1)*nElem+1:jRun*nElem,(iRun-1)*nElem+1:iRun*nElem)=1;
            indsVec(runPairCounter,:) = logical(squareform(M)); % a vector of length nElem*4*(nElem*4-1)/2, with nElem x nElem 1s in it
        end
    end
end

%% Organise data and average over run pairs
Xmat= single(zeros(nElem,nElem,length(runPairsToUse)));
for iRunPair = 1:2
    runPair = runPairsToUse(iRunPair);
    Xvec = rdm_dist(indsVec(runPair,:)); % the non-symetrical RDM for the pair, flattened (nElem*nElem)
    % Change to a square matrix, instead of flatened matrix.
    Xmat(:,:,iRunPair) = reshape(Xvec,nElem,nElem);
end

% average over the runPairsToUse
Xmat = mean(Xmat,3, 'omitnan'); 


%% Average lower and upper triangles - i.e. average the two distances
% (d1,d2) between each couple of conditions: d1 is between condition i in run 
% 1 and condition j in run 2 and d2 is between conditino i in run 2 and
% condition j in run 1). This will end up a symmetric matrix. 

% Xmat = (Xmat + permute(Xmat,[1,2,3,5,4])) / 2;
Xmat = mean(cat(3, Xmat, permute(Xmat,[2,1])), 3, 'omitnan');

% Extract diagonal elements (distance between a condition to itself - across 
% runs). We will later append these elements to the flattned version of the 
% symmetric RDM of the off-diagonal elements
indDiag = 1:nElem+1:nElem*nElem;
diagDist = Xmat(indDiag);
% zero diagonal (to make it a distance matrix which can be easily
% flattened).
Xmat = reshape(Xmat,1,nElem*nElem);
Xmat(indDiag) = 0;
Xmat = reshape(Xmat,nElem,nElem);

% For each voxel, flatten the (off-diagonals) RDM using squareform. This is quite
% unefficient - there's probably a vectorised way of doing this. 
xRunsDist = single(nan(nElem*(nElem-1)/2,1));

xRunsDist(:) = squareform(reshape((Xmat(:,:)),nElem,nElem));

% append the diagonal elemnts at the end of the flattened off-diagonal
% vector. 
vec2save = cat(1,xRunsDist,diagDist');
if (size(vec2save,1) ~= (nElem*(nElem-1)/2 + nElem))
    warning('something is wrong with the dims')
end
% save
save(outputFile,'vec2save');

% delete full RDM (because it's huge)
delete(dataFile);