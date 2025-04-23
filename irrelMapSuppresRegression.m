clear all
close all

subs = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14',...
    '16','19','20','21','22','23','24','25'};
v = nan(length(subs),2); % number of valid trials for Alon's regression analysis
b_bothDists_rt = nan(length(subs),2,3); % betas from regression RT, nSubs x nSess x nRegs
b_onlyIrrel_rt = nan(length(subs),2,2); % betas from regression RT, nSubs x nSess x nRegs
b_bothDists_cr = nan(length(subs),2,3); % betas from regression accuracy, nSubs x nSess x nRegs
b_onlyIrrel_cr = nan(length(subs),2,2); % betas from regression accuracy, nSubs x nSess x nRegs
b_numStayTrBeforeValidTr_cr = nan(length(subs),2,2); % betas from regression accuracy, nSubs x nSess x nRegs
b_numStayTrBeforeValidTr_rt = nan(length(subs),2,2); % betas from regression accuracy, nSubs x nSess x nRegs

for s = 1:length(subs)
    for session = 1:2;
        
        load(strtrim(ls(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/Subj_',subs{s},'/*_session_',num2str(session),'/data_',subs{s},'_',num2str(session),'.mat'])))
            
            
            %         load(['/vols/Scratch/mgarvert/ManyMaps/imagingData/Subj_',num2str(subjNo),'/session_',num2str(session),'/behaviour/data_Subj',num2str(subjNo),'_session',num2str(session),'.mat'])
            RT              = [];
            cr              = [];
            % regressors, we'll concatenate across blocks
            irrelDistInPrecedTr = [];
            relDistInPrecedTr   = [];

            numStayTrBeforeValidTr = [];

            
            for bl = 1:4
                
                d = data.scan{bl}.objDiff;
                % which object trials were followed by a choice trial (note
                % choice can be shorter than the number of obj trials
                % because it will end on the last obj trial followed by a
                % choice trial)
                choice = d.stimOn > 0 & d.RT(1:length(d.stimOn)) < 100;
                

                %         c2(session,bl,:) = regress(RT(session,(bl-1)*17+1:bl*17)',[ones(1,17); dRelCh(session,(bl-1)*17+1:bl*17);dRelUnch(session,(bl-1)*17+1:bl*17); dIrrelCh(session,(bl-1)*17+1:bl*17);dIrrelUnch(session,(bl-1)*17+1:bl*17); ]');
                
                % bhv data
                objTr = data.scan{bl}; % object trials
                chTr  = data.scan{bl}.objDiff; % choice trials

                % get distances on both maps

                % ix = find(objTr.map == m);

                % Exclude first trial: doesn't have a preceding object
                % ix(ix == 1) = [];
                ix = 2:length(objTr.map); % ignore first trial because it doesn't have a preceding object
                curr_obj = data.scan{bl}.seq(1,ix);
                prev_obj = data.scan{bl}.seq(1,ix-1);
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

                % get choice trials of interest - switch trials preceded by
                % stay trials

                % first take care of end cases - we want to look back 2
                % trials from the current trials, so we need to check what
                % happens if one of the first 2 trials are a choice trial
                if choice(1) & choice(2)
                    choice(1:2) = false;
                    switchChoiceTr = objTr.map(choice) ~= objTr.map(find(choice)-1);
                    switchChoiceTr = [false false switchChoiceTr];
                    noSwitchIn2trials = objTr.map(find(choice)-1) == objTr.map(find(choice)-2);
                    noSwitchIn2trials = [false false noSwitchIn2trials];
                    choice(1:2) = true;
                elseif choice(1) 
                    choice(1) = false;
                    switchChoiceTr = objTr.map(choice) ~= objTr.map(find(choice)-1);
                    switchChoiceTr = [false switchChoiceTr];
                    noSwitchIn2trials = objTr.map(find(choice)-1) == objTr.map(find(choice)-2);
                    noSwitchIn2trials = [false noSwitchIn2trials];
                    choice(1) = true;
                elseif choice(2)
                    switchChoiceTr = objTr.map(choice) ~= objTr.map(find(choice)-1);
                    choice(2)=false;
                    noSwitchIn2trials = objTr.map(find(choice)-1) == objTr.map(find(choice)-2);
                    noSwitchIn2trials = [false noSwitchIn2trials];      
                    choice(2)=true;
                else % most cases
                    switchChoiceTr = objTr.map(choice) ~= objTr.map(find(choice)-1);
                    noSwitchIn2trials = objTr.map(find(choice)-1) == objTr.map(find(choice)-2);
                end
                validChoiceTrialsForReg = switchChoiceTr & noSwitchIn2trials; % in indeces of choice trials in block
                choiceIndsinObjtr = find(choice);
                validChoiceIndsInObjTr = choiceIndsinObjtr(validChoiceTrialsForReg);
                numValidTrials(s,session,bl) = length(validChoiceIndsInObjTr); % to put in the title of fig 

                cr          = [cr d.cr(validChoiceIndsInObjTr)];
                RT          = [RT log(d.RT(validChoiceIndsInObjTr))];

                % now need to calculate distnaces of all objTr from their
                % preceding ObjTr for both relevant and irrelevant maps. 
                % The regressor of interest is the irrlevent dist on the trial preceding 
                % the switch trial. 
                % use relevant dist on that trial as confound?
                % use both distances on the switch trials themselves as
                % confounds?
                irrelDistInPrecedTr = [irrelDistInPrecedTr dIrrel(validChoiceIndsInObjTr-1)];
                relDistInPrecedTr = [relDistInPrecedTr dRel(validChoiceIndsInObjTr-1)];
                
                % get number of stay trials prior to switch trial
                d = [nan diff(objTr.map)];
                numStayTrBeforeValidTr_currentBlock = nan(length(validChoiceIndsInObjTr),1);
                for iValid = 1:length(validChoiceIndsInObjTr)
                    count = 0;
                    goBack = true;
                    iBack = 1;
                    while goBack                        
                        try % for when  validChoiceIndsInObjTr(iValid) <= iBack                         
                            if d(validChoiceIndsInObjTr(iValid)-iBack)==0
                                count = count + 1;
                                iBack = iBack+1;
                            else
                                goBack = false;
                            end
                        catch
                            goBack = false;
                        end                        
                    end
                    numStayTrBeforeValidTr = [numStayTrBeforeValidTr count];    
                end
            end

            [b_onlyIrrel_cr(s,session,:)] = regress(cr',[ones(length(irrelDistInPrecedTr),1),irrelDistInPrecedTr']);
            [b_bothDists_cr(s,session,:)] = regress(cr',[ones(length(irrelDistInPrecedTr),1),irrelDistInPrecedTr',relDistInPrecedTr']);
            [b_onlyIrrel_rt(s,session,:)] = regress(RT',[ones(length(irrelDistInPrecedTr),1),irrelDistInPrecedTr']);
            [b_bothDists_rt(s,session,:)] = regress(RT',[ones(length(irrelDistInPrecedTr),1),irrelDistInPrecedTr',relDistInPrecedTr']);
            [b_numStayTrBeforeValidTr_cr(s,session,:)] = regress(cr',[ones(length(numStayTrBeforeValidTr),1),numStayTrBeforeValidTr']);
            [b_numStayTrBeforeValidTr_rt(s,session,:)] = regress(RT',[ones(length(numStayTrBeforeValidTr),1),numStayTrBeforeValidTr']);
    end
end

%% 

figure() % accuracy
tcl = tiledlayout(2,1);
nexttile; hold on
y = b_onlyIrrel_cr(:,1,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
plot([0,2],[0,0],'k--')
title(sprintf('session 1, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,1,:),3)),P))
nexttile; hold on
y = b_onlyIrrel_cr(:,2,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
plot([0,2],[0,0],'k--')
title(sprintf('session 2, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,2,:),3)),P))
title(tcl, 'efct on correctness at choice+switch trl of irrel dist in preceding trl')

figure() % rt
tcl = tiledlayout(2,1);
nexttile; hold on
y = b_onlyIrrel_rt(:,1,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
plot([0,2],[0,0],'k--')
title(sprintf('session 1, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,1,:),3)),P))
nexttile; hold on
y = b_onlyIrrel_rt(:,2,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
title(sprintf('session 2, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,2,:),3)),P))
plot([0,2],[0,0],'k--')
title(tcl, 'efct on rt at choice+switch trl of irrel dist in preceding trl')

figure() % accuracy # preceding stay trials
tcl = tiledlayout(2,1);
nexttile; hold on
y = b_numStayTrBeforeValidTr_cr(:,1,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
plot([0,2],[0,0],'k--')
title(sprintf('session 1, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,1,:),3)),P))
nexttile; hold on
y = b_numStayTrBeforeValidTr_cr(:,2,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
plot([0,2],[0,0],'k--')
title(sprintf('session 2, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,2,:),3)),P))
title(tcl, 'efct on correctness at choice+switch trl of # preceding stay trials')

figure() % rt # preceding stay trials
tcl = tiledlayout(2,1);
nexttile; hold on
y = b_numStayTrBeforeValidTr_rt(:,1,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
plot([0,2],[0,0],'k--')
title(sprintf('session 1, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,1,:),3)),P))
nexttile; hold on
y = b_numStayTrBeforeValidTr_rt(:,2,2);
[~,P] = ttest(y);
bar(1,mean(y)); errorbar(1,mean(y),std(y)./sqrt(length(y)))
plot([0,2],[0,0],'k--')
title(sprintf('session 2, avg %0.2f trls p/sub, p=%0.2f',mean(sum(numValidTrials(:,2,:),3)),P))
title(tcl, 'efct on RT at choice+switch trl of # preceding stay trials')
