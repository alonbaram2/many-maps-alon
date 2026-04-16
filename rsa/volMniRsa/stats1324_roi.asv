function stats1324_roi(rootData, subjects, sessions, roi, model)

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
            rdmDir = fullfile(rootData,'rsa_alon',subj,'dataRdms', roi, session);
            rdmMatFile1 = fullfile(rdmDir,'dist_correlation_xRun1324Collapsed.mat');
            rdmMatFile2 = fullfile(rdmDir,[roi '_ROI_RDM_correlation_collapsed1324.mat']);

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
                collapseRdm1324_roi(rootData,subj,session, roi)
            end

            % linear model - distRel_sameMap, distIrrel_sameMap or
            % distRel_sameMapNoSameCond, distIrrel_sameMapNoSameCond

            X1 = analyses.RDMs.xRun{1}; %distRel
            X2 = analyses.RDMs.xRun{2}; % distIrrel
            y = rdm;
            tbl = table(X1, X2, y);
            % Fit model
            mdl = fitlm(tbl, 'y ~ X1 + X2');
            betas(iSub,iSess,1) = mdl.Coefficients.Estimate(2); % X1
            betas(iSub,iSess,2) = mdl.Coefficients.Estimate(3); % X2
        end
    end
    
    %% save the betas to plot nicely in python
    save(fullfile(rootData,'FigsForCurrentBiologyProofs',['betas_' roi '_' model '.mat']),'betas')


    %% plotting the betas

    figure('Color','w');
    x =[betas(:,1,1), betas(:,1,2),betas(:,2,1), betas(:,2,2)];

    boxplot(x, ...
        'Labels', {'Relevant','Irrelevant','Relevant','Irrelevant'}, ...
        'Symbol','k');   % outliers as black +
    set(gca, ...
        'TickDir','out', ...
        'Box','off', ...
        'LineWidth',1.2, ...
        'FontSize',fs);

    ylabel('betas','FontSize',fs);
    title('RSA betas','FontSize',fs);

    % ---- Color the boxes ----
    h = findobj(gca,'Tag','Box'); %order is according to plotting -h(4) is session 1, relevant

    patch(get(h(4),'XData'), get(h(4),'YData'), rel_color, ...
        'FaceAlpha',0.6, 'EdgeColor','None');

    patch(get(h(3),'XData'), get(h(3),'YData'), irrel_color, ...
        'FaceAlpha',0.6, 'EdgeColor','None');

    patch(get(h(2),'XData'), get(h(2),'YData'), rel_color, ...
        'FaceAlpha',0.6, 'EdgeColor','None');

    patch(get(h(1),'XData'), get(h(1),'YData'), irrel_color, ...
        'FaceAlpha',0.6, 'EdgeColor','None');    
    hold on
    % ---- Connect paired observations ----
    for i = 1:length(subjects)
        plot([1 2], ...
            [betas(i,1,1), betas(i,1,2)], ...
            'Color', [0.6 0.6 0.6], 'LineWidth',1);
        plot([3 4], ...
            [betas(i,2,1), betas(i,2,2)], ...
            'Color', [0.6 0.6 0.6], 'LineWidth',1);
    end

    % ---- Overlay individual subject points ----
    hold on
    scatter(ones(length(subjects)), betas(:,1,1), 40, ...
        'MarkerFaceColor', rel_color, ...
        'MarkerEdgeColor',rel_color);

    scatter(2*ones(length(subjects)), betas(:,1,2), 40, ...
        'MarkerFaceColor', irrel_color, ...
        'MarkerEdgeColor',irrel_color);

    scatter(3*ones(length(subjects)), betas(:,2,1), 40, ...
        'MarkerFaceColor', rel_color, ...
        'MarkerEdgeColor',rel_color);

    scatter(4*ones(length(subjects)), betas(:,2,2), 40, ...
        'MarkerFaceColor', irrel_color, ...
        'MarkerEdgeColor',irrel_color);    




    %%
    %-----------------------------------------------------------------------
    % From here it's chatGPT code to run ANOVA on the betas

    % Inputs: betas (nSubjects x nSessions x nModels) e.g. 22 x 2 x 2
    % Assumptions: session dimension ordering = [1 2] (S1, S2)
    %              model dimension ordering   = [1 2] (M1, M2)
    %-----------------------------------------------------------------------
    if ~exist('betas','var')
        error('This script expects a variable named ''betas'' in the workspace.');
    end

    [nSubjects, nSessions, nModels] = size(betas);
    if ~isequal(nSessions,2) || ~isequal(nModels,2)
        error('This script assumes exactly 2 sessions and 2 models (size should be n x 2 x 2).');
    end

    % Make a table with columns: M1_S1, M1_S2, M2_S1, M2_S2
    % Column order chosen to match earlier examples
    M1_S1 = squeeze(betas(:,1,1)); % n x 1
    M1_S2 = squeeze(betas(:,2,1));
    M2_S1 = squeeze(betas(:,1,2));
    M2_S2 = squeeze(betas(:,2,2));

    T = table(M1_S1, M1_S2, M2_S1, M2_S2);

    % Optionally add subject ID column (not required by fitrm but useful)
    T.Subject = (1:nSubjects)';

    % Re-order columns so Subject is first (not necessary but clean)
    T = movevars(T,'Subject','Before',1);

    % Define within-subject factors for fitrm
    % rows correspond to the repeated measures variables in the order of table columns:
    % here the repeated measures variables are the 4 cells: M1_S1, M1_S2, M2_S1, M2_S2
    Within = table( ...
        {'M1'; 'M1'; 'M2'; 'M2'}, ...   % Model factor labels
        {'S1'; 'S2'; 'S1'; 'S2'}, ...   % Session factor labels
        'VariableNames', {'Model','Session'});

    % Fit the repeated-measures model and run ANOVA
    % The dependent variables in the formula are the 4 repeated columns
    rm = fitrm(T,'M1_S1-M2_S2 ~ 1','WithinDesign',Within);

    % Run repeated-measures ANOVA for Model * Session
    ranovatbl = ranova(rm,'WithinModel','Model*Session');

    fprintf('\n===== Repeated-measures ANOVA (2x2 Model x Session) =====\n');
    disp(ranovatbl);

    % Extract and print means and SEM for each cell
    means = [mean(M1_S1) mean(M1_S2) mean(M2_S1) mean(M2_S2)];
    sems  = [std(M1_S1)/sqrt(nSubjects) std(M1_S2)/sqrt(nSubjects) ...
        std(M2_S1)/sqrt(nSubjects) std(M2_S2)/sqrt(nSubjects)];

    cellNames = {'M1\_S1','M1\_S2','M2\_S1','M2\_S2'};
    fprintf('\nCell means (SEM):\n');
    for k=1:4
        fprintf(' %s: %.4f (%.4f)\n', cellNames{k}, means(k), sems(k));
    end

    % 1) Effect of Session (averaged over models): paired t-test S1 vs S2
    % session1 per subject = mean over models at session 1
    session1 = squeeze(mean(betas(:,1,:),3)); % n x 1
    session2 = squeeze(mean(betas(:,2,:),3)); % n x 1

    [~,p_sess,ci_sess,stats_sess] = ttest(session1, session2);
    diff_sess = session1 - session2;
    cohend_sess = mean(diff_sess)./std(diff_sess); % paired Cohen's d (difference / SDdiff)

    fprintf('\n===== Session effect (collapsed over models): S1 vs S2 =====\n');
    fprintf('t(%d)=%.3f, p=%.4f, mean diff=%.4f, Cohen''s d (paired)=%.3f\n', ...
        stats_sess.df, stats_sess.tstat, p_sess, mean(diff_sess), cohend_sess);

    % 2) Effect of Model (averaged over sessions): paired t-test M1 vs M2
    % model1 per subject = mean over sessions at model 1
    model1 = squeeze(mean(betas(:,:,1),2)); % n x 1
    model2 = squeeze(mean(betas(:,:,2),2)); % n x 1

    [~,p_mod,ci_mod,stats_mod] = ttest(model1, model2);
    diff_mod = model1 - model2;
    cohend_mod = mean(diff_mod)./std(diff_mod);

    fprintf('\n===== Model effect (collapsed over sessions): M1 vs M2 =====\n');
    fprintf('t(%d)=%.3f, p=%.4f, mean diff=%.4f, Cohen''s d (paired)=%.3f\n', ...
        stats_mod.df, stats_mod.tstat, p_mod, mean(diff_mod), cohend_mod);

    % 3) Grand mean across sessions & models (one-sample test vs 0)
    % grand mean per subject: average over sessions and models
    grandmean = squeeze(mean(mean(betas,2),3)); % n x 1

    [~,p_grand,ci_grand,stats_grand] = ttest(grandmean);
    cohend_grand = mean(grandmean)./std(grandmean); % one-sample Cohen's d (mean/SD)

    fprintf('\n===== Grand mean (collapsed over models and sessions) vs 0 =====\n');
    fprintf('t(%d)=%.3f, p=%.4f, mean=%.4f, Cohen''s d (one-sample)=%.3f\n', ...
        stats_grand.df, stats_grand.tstat, p_grand, mean(grandmean), cohend_grand);

    % Optional: one-sample t-tests for each cell (vs zero)
    fprintf('\n===== One-sample t-tests for each cell (vs 0) =====\n');
    for s = 1:2
        for m = 1:2
            cellvec = squeeze(betas(:,s,m));
            [~,p_cell,ci_cell,stats_cell] = ttest(cellvec);
            d_cell = mean(cellvec)./std(cellvec);
            fprintf(' %s: t(%d)=%.3f, p=%.4f, mean=%.4f, d=%.3f\n', ...
                sprintf('M%d_S%d',m,s), stats_cell.df, stats_cell.tstat, p_cell, mean(cellvec), d_cell);
        end
    end

    % Summary outputs as a structure (optional, useful for programmatic use)
    results.ranova = ranovatbl;
    results.means = means;
    results.sems = sems;
    results.session = struct('tstat', stats_sess.tstat, 'df', stats_sess.df, 'p', p_sess, 'cohend', cohend_sess);
    results.model   = struct('tstat', stats_mod.tstat, 'df', stats_mod.df, 'p', p_mod, 'cohend', cohend_mod);
    results.grand   = struct('tstat', stats_grand.tstat, 'df', stats_grand.df, 'p', p_grand, 'cohend', cohend_grand);
    results.cell = table(cellNames', means', sems', 'VariableNames', {'Cell','Mean','SEM'});

    % Return results struct to workspace
    assignin('base','ANOVA_results',results);
    fprintf('\nResults saved to variable ''ANOVA_results'' in base workspace.\n');

    % %
    % % [1 1] contrast across rel and irrel, and also average over sessions
    [h,p,ci,stats] = ttest(mean(mean(betas,3),2), 0, 'Tail', 'right');
    fprintf('%s\n%s\none-sided p = %.3f\ndf = %d\nt = %.3f\n ',model,roi,p,stats.df,stats.tstat)
    %%






elseif strcmp(model, 'distRel_diffMaps')
    %%
%    roi = '322_diff_04_mPFCbigger_2p5';
    % roi = '322_diff_04_mPFCsmaller_2p5';
     % roi = '322_both_04_mPFC_2p5_mask';
  
    % to test in mPFC effects.     
    % when only have a single regressor, better to use Kendall's Tau than a
    % linear regression
    tau = nan(length(subjects),length(sessions));

    for iSub = 1:length(subjects)
        for iSess = 1:length(sessions)
            subj = subjects{iSub};
            session = sessions{iSess};

            
            [modelRdm, ~] = getModelRdm(model,false,subj,session);
            
            rdmDir = fullfile(rootData,'rsa_alon',subj,'dataRdms', roi, session);
            rdmMatFile = fullfile(rdmDir,'dist_correlation_xRun1324Collapsed.mat');
            rdmMatFile2 = fullfile(rdmDir,[roi '_ROI_RDM_correlation_collapsed1324.mat']);            
            rdmNiiFile = fullfile(rdmDir,'dist_correlation_xRun1324Collapsed.nii');
            if exist(rdmMatFile,'file')
                load (rdmMatFile,'vec2save')
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
                collapseRdm1324_roi(rootData,subj,session, roi)
            end

            
            tau(iSub,iSess) = rsa.stat.rankCorr_Kendall_taua(rdm, modelRdm);

        end
    end

    
    % test if the tau correcoefs of session 2 are bigger than session 1
    [h,p,ci,stats] = ttest(tau(:,2),tau(:,1), 'Tail', 'right');
    fprintf('sess 2 > sess 1: %s\n%s\n one sided p = %.3f\nt_21=%.3f',model,roi,p,stats.tstat)
    

    %% save taus of to plot in python
    save(fullfile(rootData,'FigsForCurrentBiologyProofs',['tau_' roi '_' model '.mat']),'tau')

end