function stats1324_roi(rootData, subjects, sessions, rois, model)

%%
% model = 'relIrrel_sameMap';
% subjects = {'Subj_1','Subj_2','Subj_3','Subj_4','Subj_5','Subj_6','Subj_7','Subj_8',...
%     'Subj_9','Subj_10','Subj_11','Subj_12','Subj_13','Subj_14','Subj_16',...
%     'Subj_19','Subj_20','Subj_21','Subj_22','Subj_23','Subj_24','Subj_25' };
% 
% sessions = {'session_1','session_2'};
% roi = '322_both_07_left_parahippoc_2p5_mask';


%% plotting configs

fs = 12;   % <-- global font size
rel_color   = [109 199 183] / 255;   % #6DC7B7
irrel_color = [27 119 96]  / 255;    % #1B7760


%%

if (strcmp(model,'relIrrel_sameMap') | strcmp(model,'relIrrel_sameMapNoSameCond'))

    for iRoi=1:2
        % betas for a single model (with more than one regressor, like rel+irrel_sameMap)
        if strcmp(model,'relIrrel_sameMap')
            analyses.names = {'distRel_sameMap','distIrrel_sameMap'};
        elseif strcmp(model,'relIrrel_sameMapNoSameCond')
            analyses.names = {'distRel_sameMapNoSameCond','distIrrel_sameMapNoSameCond'};
        end
        betas = nan(length(subjects),length(sessions),length(analyses.names));

        for iSub = 1:length(subjects)
            for iSess = 1:length(sessions)
                subj = subjects{iSub};
                session = sessions{iSess};

                for iAn = 1:length(analyses.names)
                    [analyses.RDMs.xRun{iAn}, analyses.RDMs.within{iAn}] = getModelRdm(analyses.names{iAn},false,subj,session);
                end
                rdmDir = fullfile(rootData,'rsa_alon',subj,'dataRdms', rois{iRoi}, session);
                rdmMatFile1 = fullfile(rdmDir,'dist_correlation_xRun1324Collapsed.mat');
                rdmMatFile2 = fullfile(rdmDir,[rois{iRoi} '_ROI_RDM_correlation_collapsed1324.mat']);

                rdmNiiFile = fullfile(rdmDir,'dist_correlation_xRun1324Collapsed.nii');
                if exist(rdmMatFile1,'file')
                    load (rdmMatFile1,'vec2save')
                    rdm = vec2save;
                elseif exist(rdmMatFile2,'file')
                    load (rdmMatFile2,'vec2save')
                    rdm = vec2save;
                elseif (exist(rdmNiiFile,'file')) % read RDM from one of the voxels - they're all the same
                    rdmNii = niftiread(rdmNiiFile);
                    mask = any(~isnan(rdmNii), 4);   % 3D logical mask
                    [X,Y,Z] = ind2sub(size(mask), find(mask));
                    x = X(1); y = Y(1); z = Z(1); % it doesn't matter which voxel from the mask we're choosing.
                    rdm = squeeze(rdmNii(x,y,z,:)); %
                else
                    collapseRdm1324_roi(rootData,subj,session, rois{iRoi})
                end

                % linear model - distRel_sameMap, distIrrel_sameMap or
                % distRel_sameMapNoSameCond, distIrrel_sameMapNoSameCond

                X1 = analyses.RDMs.xRun{1}; %distRel
                X2 = analyses.RDMs.xRun{2}; % distIrrel
                y = rdm;
                tbl = table(X1, X2, y);
                % Fit model
                mdl = fitlm(tbl, 'y ~ X1 + X2');
                betas(iSub,iSess, iRoi,1) = mdl.Coefficients.Estimate(2); % X1
                betas(iSub,iSess, iRoi,2) = mdl.Coefficients.Estimate(3); % X2
            end
        end
    end
    

    % %
    % % [1 1] contrast across rel and irrel, and also average over
    % sessions, difference between ROIs
    [h,p,ci,stats] = ttest(mean(mean(betas(:,:,1,:),4),2),mean(mean(betas(:,:,2,:),4),2));
    fprintf('%s\n%s - %s\ntwo-sided p = %.3f\ndf = %d\nt = %.3f\n ',model,rois{1},rois{2},p,stats.df,stats.tstat)
    %%



end