clear 
close all

spmPath = '/vols/Scratch/abaram/MATLAB/spm12';
addpath(spmPath);

root='/home/fs0/mgarvert/scratch/ManyMaps/imagingData';

subs = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14',...
    '16','19','20','21','22','23','24','25'};

for s = 1:length(subs)
    subStr = ['Subj_' subs{s}];
    for session = 1:2
        sessionStr = ['session_' num2str(session)];
        load(strtrim(ls(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/Subj_',subs{s},'/*_session_',num2str(session),'/data_',subs{s},'_',num2str(session),'.mat'])))
        for r = 1:4
            runStr = ['run_' num2str(r)];

            d = data.scan{r}.objDiff;
            objTr = data.scan{r}; % object trials

            % get irrelevant distance, indexed in object trials

            % Exclude first trial: doesn't have a preceding object
            ix = 2:length(objTr.map); % ignore first trial because it doesn't have a preceding object
            curr_obj = data.scan{r}.seq(1,ix);
            prev_obj = data.scan{r}.seq(1,ix-1);
            curr_obj(curr_obj>17) = curr_obj(curr_obj>17) -17; % subtract 17 if context 2
            prev_obj(prev_obj>17) = prev_obj(prev_obj>17) -17; % subtract 17 if context 2
            for i = 1:length(curr_obj) % go through all the objects
                m = objTr.map(i);
                dRel(i)   = data.map{m,m}(curr_obj(i),prev_obj(i)); % compute the distance between the two stimuli on the irrelevant$
                dIrrel(i)   = data.map{m,mod(m,2)+1}(curr_obj(i),prev_obj(i));                   % compute the distance between the two stimuli on the irrelev$
            end
            % add back the first trial so that the trial indeces match
            dRel = [nan dRel];
            dIrrel = [nan dIrrel];

            % get switch trials, indexed in Obj trials
            switchTrials = find([false objTr.map(2:end) ~= objTr.map(1:end-1)]);
            stayTrials = find([false objTr.map(2:end) == objTr.map(1:end-1)]);

            stayWithPrecedingSwitchInds = stayTrials(ismember(stayTrials-1,switchTrials));
            stayWithoutPrecedingSwitchInds = setdiff(stayTrials,stayWithPrecedingSwitchInds);
            switchWithSuccedingStayInds = switchTrials(ismember(switchTrials+1,stayTrials)); % same as stayWithPrecedingSwitch-1
            irrelDistInValidStayTrials = dIrrel(stayWithPrecedingSwitchInds);

            % load onsets of all objects
            onsetsAllObj = load(fullfile(root,subStr,sessionStr,'behaviour','EVs',[subStr '_' sessionStr '_' runStr],'allObjects.txt'));
        
            % get onsets of valid stay trials
            onsetsStayWithPrecedingSwitch = onsetsAllObj(stayWithPrecedingSwitchInds,:);
            stayWithPrecedingSwitch_irrelDist = onsetsStayWithPrecedingSwitch;
            stayWithPrecedingSwitch_irrelDist(:,3) = irrelDistInValidStayTrials;

            onsetsStayWithoutPrecedingSwitch = onsetsAllObj(stayWithoutPrecedingSwitchInds,:);

            % get all non valid stay trials onsets
            % get parameteric modulator of mPFC signal in preceding switch
            % trials
            
            % get mask
            mask =  'spmT_switch_stay_session_2_0001_mask_1p5_mPFC'; %switchMinusStay from design_322
            mask_file = fullfile(root,'masks',[mask, '.nii']);
            mask_nii = spm_vol(mask_file);
            mask_data = spm_read_vols(mask_nii);

            % single-trial betas spmDir
            spmDirMni = fullfile(root,subStr,sessionStr,'1stLevel','design_401_noSmooth','MNI');
            nRegsPerRun = 144; % including nuisance regressors - choiceTrials, buttonPress, 6*motionParams
            stayWithPrecedingSwitch_switchMPfc = stayWithPrecedingSwitch_irrelDist; % get the onsets and the 3-col structure
            for iValidStay = 1:length(stayWithPrecedingSwitchInds)
                % for each valid stay trial, load single-trial beta map of the
                % preceding switch trial
                betaNum = (r-1)*nRegsPerRun + switchWithSuccedingStayInds(iValidStay);
                betaNumStr = sprintf('%0*d', 4, betaNum);
                betaFile = fullfile(spmDirMni,['MNI_beta_' betaNumStr '.nii']);
                beta_nii = spm_vol(betaFile);
                beta_data = spm_read_vols(beta_nii);
                meanSignalInMask = mean(beta_data(mask_data ~= 0),'omitnan');
                % update in 3-col regressor
                stayWithPrecedingSwitch_switchMPfc(iValidStay,3) = meanSignalInMask;
            end
            % caculate the interaction regressor of the two regressors
            % get the 3-col structure and onsets
            stayWithPrecedingSwitch_switchMPfc_X_irrelDis = stayWithPrecedingSwitch_switchMPfc; 
            % demean the original regs and calculate interaction
            stayWithPrecedingSwitch_switchMPfc_X_irrelDis(:,3) = (stayWithPrecedingSwitch_switchMPfc(:,3)-mean(stayWithPrecedingSwitch_switchMPfc(:,3))) .* (stayWithPrecedingSwitch_irrelDist(:,3)-mean(stayWithPrecedingSwitch_irrelDist(:,3)));


            % save 3-col regressor txt files
            dlmwrite(fullfile(root,subStr,sessionStr,'behaviour','EVs',[subStr '_' sessionStr '_' runStr],'stayWithPrecedingSwitch_switchMPfc.txt'), stayWithPrecedingSwitch_switchMPfc, 'delimiter', ' ');
            dlmwrite(fullfile(root,subStr,sessionStr,'behaviour','EVs',[subStr '_' sessionStr '_' runStr],'stayWithPrecedingSwitch_irrelDist.txt'), stayWithPrecedingSwitch_irrelDist, 'delimiter', ' ');
            dlmwrite(fullfile(root,subStr,sessionStr,'behaviour','EVs',[subStr '_' sessionStr '_' runStr],'stayWithPrecedingSwitch_switchMPfc_X_irrelDis.txt'), stayWithPrecedingSwitch_switchMPfc_X_irrelDis, 'delimiter', ' ');
            dlmwrite(fullfile(root,subStr,sessionStr,'behaviour','EVs',[subStr '_' sessionStr '_' runStr],'stayWithoutPrecedingSwitch.txt'), onsetsStayWithPrecedingSwitch, 'delimiter', ' ');

        end
    end
end