function [vecRdmXRun, vecRdmWithin] = getModelRdm(analysisName,plotFlag,subj,session)

nElem = 34;
scanFile = ls('/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/Subj_1/*_session_1/data_1_1.mat');

if exist(fullfile(['/vols/Scratch/mgarvert/ManyMaps/imagingData/rsa_alon/' subj]),'dir') % exclude subj 15,17,18 which don't have 2 sessions
    load(strtrim(ls(fullfile('/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles',subj,['*_' session],['data_' extractAfter(subj,'_') '_' extractAfter(session,'_') '.mat']))))
else
    error('no RDMs for this subject, probably only has 1 session')
end

% in all RDMs, 0 is similar and 1 is dissimilar, nan is ignore.  
switch analysisName

    case {'identity_bothMaps','identity_noSameCond','identity_diffMaps'}
        clear ix ix1 ix2 order
        for i = 1:17
            ix1(i,data.stimuli(2,:) == data.stimuli(1,i)) = 1;
            ix2(i,data.stimuli(1,:) == data.stimuli(2,i)) = 1;
            order(i) = find(data.stimuli(2,:) == data.stimuli(1,i));
            order2(i) = find(data.stimuli(1,:) == data.stimuli(2,i));
        end
        ix1 = 1-ix1;   ix2 = 1-ix2;
        switch analysisName
            case 'identity_bothMaps'
                RDM = [1-eye(17) ix1; ix2 1-eye(17)] ;
            case 'identity_noSameCond'
                RDM = [1-eye(17) ix1; ix2 1-eye(17)] ;
                % nan diagonal
                RDM(logical(eye(length(RDM)))) = nan(nElem,1);
            case 'identity_diffMaps'
                RDM = [nan(17) ix1; ix2 nan(17)] ;
        end

    case 'context'
        RDM = 1-[ones(17) zeros(17); zeros(17) ones(17)];
    case 'context_noSameCond'
        RDM = 1-[ones(17) zeros(17); zeros(17) ones(17)];
        % nan diag
        RDM(logical(eye(length(RDM))))  = nan(nElem,1);
    case 'position_bothMaps'
        RDM = 1-[eye(17) eye(17); eye(17) eye(17)];
    case 'position_noSameCond'
        RDM = 1-[eye(17) eye(17); eye(17) eye(17)];
        % nan diag
        RDM(logical(eye(length(RDM)))) = nan(nElem,1);
    case 'position_diffMaps'
        RDM = 1-[nan(17) eye(17); eye(17) nan(17)];

    case 'distRel_bothMaps'
        RDM = [data.map{1,1} data.map{1,1}; data.map{2,2} data.map{2,2}];
    case 'distRel_noSameCond'
        RDM = [data.map{1,1} data.map{1,1}; data.map{2,2} data.map{2,2}];
        RDM(logical(eye(length(RDM)))) = nan(nElem,1);
    case 'distRel_sameMap'
        RDM = [data.map{1,1} NaN(17); NaN(17) data.map{2,2}];
    case 'distRel_diffMaps'
        RDM = [NaN(17) data.map{1,1}; data.map{1,1} NaN(17)];

    case 'distIrrel_sameMap'
        RDM  = [data.map{1,2} NaN(17); NaN(17) data.map{2,1}];

    % case 'arena_rel'
    %     try
    %         arena{1} = load(fullfile(arenadir,['similarityJudgementData_Subj_',num2str(subj),'_map_1_trial1.mat']));
    %         arena{2} = load(fullfile(arenadir,['similarityJudgementData_Subj_',num2str(subj),'_map_2_trial1.mat']));
    %         RDM = [squareform(arena{1}.estimate_RDM_ltv) NaN(17); NaN(17) squareform(arena{2}.estimate_RDM_ltv)];
    %     catch
    %         RDM = zeros(nElem);
    %     end
    % case 'arena_irrel'
    %     try
    %         arena{1} = load(fullfile(arenadir,['similarityJudgementData_Subj_',num2str(subj),'_map_1_trial1.mat']));
    %         arena{2} = load(fullfile(arenadir,['similarityJudgementData_Subj_',num2str(subj),'_map_2_trial1.mat']));
    %         RDM  = [squareform(arena{2}.estimate_RDM_ltv) NaN(17); NaN(17) squareform(arena{1}.estimate_RDM_ltv)];
    %     catch
    %         RDM = zeros(nElem);
    %     end

end

diagon = diag(RDM);

vecRdmWithin = vectoriseRDM(RDM); % this will not take the diagonal
vecRdmXRun = vertcat(vectoriseRDM(RDM),diagon);





if plotFlag

    figure;
    addpath('/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/alon/utilities/') 
    subplot(1,2,1)
    imagescwithnan(RDM,hot,[0 1 1])
    hold on 
    title('xRun')
    ax = gca();
    set(gca, 'FontSize', 14);    
    axis square
    
    subplot(1,2,2)
    RDM_noDiag=RDM;
    RDM_noDiag(logical(eye(length(RDM)))) = nan;
    imagescwithnan(RDM_noDiag,hot,[0 1 1])
    title('withinRun')
    ax = gca();
    set(gca, 'FontSize', 14);
    % ax.XTick=1:length(labels);
    % ax.YTick=1:length(labels);
    % ax.XTickLabels = labels;
    % ax.YTickLabels = labels;
    % ax.XTickLabelRotation = 45;
    
    axis square
    sgtitle(analysisName,'FontSize', 18)
end

function RDM = makeRdmSymmetric(RDM)
% make symmetric
for i=1:length(RDM)
    for j=i:length(RDM)
        if i<j
            RDM(j,i)=RDM(i,j);
        end
    end
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

