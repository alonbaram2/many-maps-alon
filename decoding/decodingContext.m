clear

root = '/vols/Scratch/mgarvert/ManyMaps/imagingData';

% currently excluding Subj_10 becasue in runs 3 and 4 has less than 136 trials
subs = {'1','2','3','4','5','6','7','8','9','11','12','13','14',... 
    '16','19','20','21','22','23','24','25'};


corr_irrelScore_dIrrel = nan(length(subs),2);
accuracy = nan(length(subs),2,4);


for iSub =1:length(subs)
    subj = ['Subj_' subs{iSub}];
    for iSess=1:2
        sess= ['session_' int2str(iSess)];
        spmDir = fullfile(root,subj,sess,'1stLevel','design_401_noSmooth');
        mask =  'MNI_bothCons_thrAbs0p5_X_322_both_07_bilateral_phc_1p5'; 
        decodingDataDir = fullfile(root,'decoding','data',mask);

        % Load your data, into variable D
        load(fullfile(decodingDataDir,[subj '_' sess '.mat']));
        [nVoxels, nTrials, nRuns] = size(D);

        % load labels (loading variable 'data')
        labels = nan(nTrials,nRuns);
        stayTrials = false(nTrials,nRuns); % logical array
        load(strtrim(ls(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/' subj '/*_' sess,'/data_',subj(6:end),'_',sess(end),'.mat'])))

        for iRun=1:nRuns
            objTr = data.scan{iRun}; % object trials

            % get Rel and Irrel distances
            ix = 2:length(objTr.map); % ignore first trial because it doesn't have a preceding object
            curr_obj = objTr.seq(1,ix);
            prev_obj = objTr.seq(1,ix-1);
            curr_obj(curr_obj>17) = curr_obj(curr_obj>17) -17; % subtract 17 if context 2
            prev_obj(prev_obj>17) = prev_obj(prev_obj>17) -17; % subtract 17 if context 2
            for i = 1:length(curr_obj) % go through all the objects
                m = objTr.map(i);
                dRel(i)   = data.map{m,m}(curr_obj(i),prev_obj(i)); % compute the distance between the two stimuli on the irrelevant$
                dIrrel(i)   = data.map{m,mod(m,2)+1}(curr_obj(i),prev_obj(i));                   % compute the distance between the two stimuli on the irrelev$
            end
            % add back the first trial so that the trial indeces
            % match
            dRel = [nan dRel];
            dIrrel = [nan dIrrel];

            % get labels - contexts
            labels(:,iRun) = objTr.map';

            % get stay trials
            stayTrials(:,iRun) =  [false diff(objTr.map)==0];
        end

        clear data

        % Assuming you have your data stored in 'D' and labels in 'labels'

        %%

        scoreIrrelStayTrials = []; % will add to this from every run, because diff number of stay trials in each run
        dIrrelStayTrials = [];

        % Cross-validation over runs
        for iRun = 1:nRuns
            % Create training and testing sets
            training_data = (reshape(D(:,:,[1:iRun-1, iRun+1:end]),[nVoxels, nTrials * (nRuns-1)]))'; % nTrials x nVox 
            % training_data_zscored = (training_data - mean(training_data)) ./ std(training_data); % zscore the features (voxels)
            training_labels = reshape(labels(:, [1:iRun-1, iRun+1:end]), [], 1);

            testing_data = (squeeze(D(:,:,iRun)))'; % nTrials x nVox
            % testing_data_zscored = (testing_data - mean(training_data)) ./ std(training_data); % use mean and std from training data to zscore
            testing_labels = labels(:, iRun);

            % % Specify regularization strength (lambda). A smaller value indicates stronger regularization.
            % lambda = 1;


            % Train the classifier (e.g., linear SVM)
            % currently not zscoring

            if strcmp (mask,'MNI_bothCons_thrAbs0p5_X_harvardoxford_HPC_bilateral')
                lambda = 2636.65090;
                classifier = fitclinear(training_data, training_labels,'Learner', 'logistic', 'Regularization', 'ridge','Solver','lbfgs','Lambda', lambda); % transpose so that it's nObservations x nPredictors, i.e. nTrials x nVox
                [predicted_labels, score] = predict(classifier, testing_data);
            elseif strcmp(mask,'MNI_bothCons_thrAbs0p5_X_322_both_07_bilateral_phc_1p5')
                lambda = 0.04833;
                [B, FitInfo] = lassoglm(training_data, training_labels-1, 'binomial', 'Alpha', 0.5,'Lambda', lambda); % alpha=0.5 corresponds to same weighting of lasso and ridge penalty
                linear_predictor = FitInfo.Intercept + testing_data * B;
                % Apply the logistic transformation
                score = 1 ./ (1 + exp(-(FitInfo.Intercept + testing_data * B)));
                predicted_labels = score >= 0.5;  % Binary decision rule
                predicted_labels = predicted_labels + 1;
            end
            accuracy(iSub,iSess,iRun) = sum(predicted_labels == testing_labels) ./ length(predicted_labels);
            % get probabilities for the relevant map
            scoreRel = score(:,1);
            scoreRel(testing_labels==2) = 1-scoreRel(testing_labels==2);
            scoreIrrelStayTrials = vertcat(scoreIrrelStayTrials, 1-scoreRel(stayTrials(:,iRun)));
            dIrrelStayTrials = vertcat(dIrrelStayTrials, dIrrel(stayTrials(:,iRun))');
        end

        corr_irrelScore_dIrrel (iSub, iSess) = atanh(corr(scoreIrrelStayTrials,dIrrelStayTrials));
    end
end

disp(['mask: ' mask ', numVox: ' int2str(nVoxels)])

% run ttest. the expected relationship is negative: higher decoding score
% for irrelevant map when the distance in irrlevant map is small. 

[~, p,~,S] = ttest(corr_irrelScore_dIrrel(:));
disp(sprintf('both sessions. mean decoding accuracy %.2f; ttest corr dIrrel with decoding score irrel in stay trials: p = %.2f, t = %.2f',mean(mean(mean(accuracy(:,:,:)))),p,S.tstat))
[~, p,~,S] = ttest(corr_irrelScore_dIrrel(:,1));
disp(sprintf('session 1. mean decoding accuracy %.2f; ttest corr dIrrel with decoding score irrel in stay trials: p = %.2f, t = %.2f',mean(mean(squeeze(accuracy(:,1,:)))),p,S.tstat))
[~, p,~,S] = ttest(corr_irrelScore_dIrrel(:,2));
disp(sprintf('session 2. mean decoding accuracy %.2f; ttest corr dIrrel with decoding score irrel in stay trials: p = %.2f, t = %.2f',mean(mean(squeeze(accuracy(:,2,:)))),p,S.tstat))
