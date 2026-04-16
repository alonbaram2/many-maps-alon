subjects = {'Subj_1','Subj_2','Subj_3','Subj_4','Subj_5','Subj_6','Subj_7','Subj_8',...
    'Subj_9','Subj_10','Subj_11','Subj_12','Subj_13','Subj_14','Subj_16',...
    'Subj_19','Subj_20','Subj_21','Subj_22','Subj_23','Subj_24','Subj_25' };
sessions = {'session_1','session_2'};
models = {'distRel_bothMaps', 'distIrrel_sameMap'};
hemis = {'rh','lh'};
sl='WhBr_surf_r10_v100';
rootData    = '/vols/Scratch/mgarvert/ManyMaps/imagingData/';
withDiagFlag = false;

% get neural data

fname_example = ['distRel_bothMaps_xRun_smth5_rh.mgh'];
example_subjDir = fullfile(rootData,'rsa_alon','Subj_1','statistics',sl,'session_1');
analyses_fnames = dir(fullfile(example_subjDir,'*.mgh'));
example_subjFile = MRIread(fullfile(example_subjDir,analyses_fnames(1).name));
mkdir(fullfile(rootData,'rsa_alon','groupStats',sl,'Wilcoxon'));
medianData = nan(length(subjects),length(sessions),length(models),length(hemis));
medianData_1324 = nan(length(subjects),length(sessions),length(models),length(hemis));
corrTimeAndModelRDMs = nan(length(subjects),length(sessions));

for iSub= 1:length(subjects)
    subjNo = erase(subjects{iSub},'Subj_');
    for iSess=1:2
        for iHemi = 1:2
            for iModel=1:length(models)
                % get neural data
                fname = [models{iModel} '_xRun_smth5_' hemis{iHemi} '.mgh'];
                subjFile = MRIread(fullfile(rootData,'rsa_alon',subjects{iSub},'statistics',sl,sessions{iSess},fname));
                sessData = subjFile.vol;
                medianData(iSub,iSess,iModel,iHemi) = median(subjFile.vol);

                % get neural data - only 1-3,2-4 runs comparisons in RSA
                fname_1324 = [models{iModel} '_xRun1324_smth5_' hemis{iHemi} '.mgh'];
                subjFile_1324 = MRIread(fullfile(rootData,'rsa_alon',subjects{iSub},'statistics',sl,sessions{iSess},fname_1324));
                medianData_1324(iSub,iSess,iModel,iHemi) = median(subjFile_1324.vol);
            end
        end
        % get beh data
        behFile = strtrim(ls(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/Subj_',subjNo,'/*_session_',num2str(iSess),'/data_', subjNo,'_',num2str(iSess),'.mat']));
        load(behFile); % into "data"
        for iRun=1:4
            runs{iRun} = data.scan{iRun}.seq(1,:)';
        end
        [avgDist_dir, meanSym, countsSym, probBoundary_sym] = adjacency_measures(runs, 34, 10);
        chosenMeasure = probBoundary_sym;
        chosenMeasureVec = vectoriseRDM(chosenMeasure);
        chosenMeasureVecWithDiag = vertcat(chosenMeasureVec,zeros(34,1));
        [modelRdmWithDiag , modelRdmNoDiag] = getModelRdm(models{iModel},false,subjects{iSub},sessions{iSess});
        if withDiagFlag
            corrTimeAndModelRDMs(iSub,iSess)= corr(chosenMeasureVecWithDiag,modelRdmWithDiag,'Type','Spearman','Rows','complete');
        else
            corrTimeAndModelRDMs(iSub,iSess)= corr(chosenMeasureVec,modelRdmNoDiag,'Type','Spearman','Rows','complete');
        end
    end
end


corrTimeAndModelRDMs = repmat(corrTimeAndModelRDMs,[1 1 2 2]);

%%
figure; 
subplot(1,2,1)
hold on
x = medianData(:);
y = corrTimeAndModelRDMs(:);
[r,p] = corr(x,y);

% ----- Plot -----
scatter(x, y, 50, 'filled', 'MarkerFaceAlpha', 0.6)
lsline
hold on

% Clean look
box off
set(gca,'FontSize',14,'LineWidth',1.2)
xlabel('Median of whole-brain histogram of \tau, all runs','FontSize',16)
ylabel({'\rho between model RDM and', 'run boundary co-occurence probability'},'FontSize',16)
title('Correlation across subjects','FontSize',16)
xlim([min(x)-0.001,max(x)+0.001])

% Text annotation
txt = sprintf('r = %.3f\np = %.4f', r, p);
text(0.1, 0.95, txt, 'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontSize',14, ...
    'BackgroundColor','white', ...
    'EdgeColor','black');

%%
subplot(1,2,2)
hold on
x = medianData_1324(:);
y = corrTimeAndModelRDMs(:);
[r,p] = corr(x,y);

% ----- Plot -----
scatter(x, y, 50, 'filled', 'MarkerFaceAlpha', 0.6)
lsline
hold on

% Clean look
box off
set(gca,'FontSize',14,'LineWidth',1.2)
xlabel('Median of whole-brain histogram of \tau, runs 1-3,2-4','FontSize',16)
ylabel({'\rho between model RDM and', 'run boundary co-occurence probability'},'FontSize',16)
title('Correlation across subjects','FontSize',16)
xlim([min(x)-0.001,max(x)+0.001])

% Text annotation
txt = sprintf('r = %.3f\np = %.4f', r, p);
text(0.1, 0.95, txt, 'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontSize',14, ...
    'BackgroundColor','white', ...
    'EdgeColor','black');

%%

function [avgDist_dir, meanSym, countsSym, probBoundary_sym] = adjacency_measures(runs, nStim, K)

% Inputs:
% runs: 1x4 cell each length L (e.g. 136) OR 4xL matrix
% nStim: number of stimuli (34)
% K: window size for boundary co-occurrence (e.g. 5)


if nargin<2, nStim=34; end
if nargin<3, K=5; end



% normalize runs -> 4 x L
if iscell(runs)
    L = numel(runs{1});
    R = zeros(4,L);
    for r=1:4, R(r,:) = runs{r}(:)'; end
else
    R = runs;
    if size(R,1)==size(runs,2) && size(R,2)==size(runs,1)
        R = R';
    end
    L = size(R,2);
end

% --- directional avg distances (i in run r -> j in run r+1) ---
sumDist = zeros(nStim);
counts  = zeros(nStim);
for r = 1:3
    cur = R(r,:);
    nxt = R(r+1,:);
    % positions per stimulus
    posCur = cell(nStim,1); posNext = cell(nStim,1);
    for s=1:nStim
        posCur{s} = find(cur==s);
        posNext{s}= find(nxt==s);
    end
    for i=1:nStim
        a = posCur{i};
        if isempty(a), continue; end
        a_col = a(:);
        for j=1:nStim
            b = posNext{j};
            if isempty(b), continue; end
            D = (b(:).' - a_col) + L;      % distances for i->j in this pair of runs
            sumDist(i,j) = sumDist(i,j) + sum(D(:));
            counts(i,j)  = counts(i,j)  + numel(D);
        end
    end
end
avgDist_dir = nan(nStim);
mask = counts>0;
avgDist_dir(mask) = sumDist(mask) ./ counts(mask);

% --- pooled symmetric mean ---
sumMat = sumDist; sumMat(~(counts>0)) = 0;
symSum = sumMat + sumMat.';
symCounts = counts + counts.';
meanSym = nan(nStim);
mMask = symCounts>0;
meanSym(mMask) = symSum(mMask) ./ symCounts(mMask);
countsSym = symCounts;

% --- boundary co-occurrence probability (X in last K of run r AND Y in first K of r+1)
probBoundary = nan(nStim);
% count occurrences across the 3 adjacent pairs
coCounts = zeros(nStim);
denom = 3; % total adjacent run pairs
for r = 1:3
    cur = R(r,:);
    nxt = R(r+1,:);
    lastIdx = (L-K+1):L;
    firstIdx= 1:K;
    for i=1:nStim
        % was i present in last K of run r?
        if any(cur(lastIdx)==i)
            for j=1:nStim
                if any(nxt(firstIdx)==j)
                    coCounts(i,j) = coCounts(i,j) + 1;
                end
            end
        end
    end
end
probBoundary = coCounts / denom; % proportion of adjacent-run pairs showing X-last & Y-first
probBoundary_sym = (probBoundary + probBoundary.') / 2;
end


function v = vectoriseRDM(RDM)
if ~any(any(isnan(RDM))) % if there are nans, matlab will treat a symmetrical matrix as asymmetric!
    if any(any(RDM ~= RDM')) || size(RDM,1) ~= size(RDM,2)
        error('RDM must be square and symmetric')
    end
end
if any(diag(RDM))
    warning('non-zero elements on diagonal of RDM!')
end

v=RDM(tril(true(length(RDM)),-1));

end