clear

rootScripts = '/vols/Scratch/mgarvert/ManyMaps/imagingData/scripts/alon'; % this includes a copy of rsatoolbox
addpath(genpath(rootScripts));
spmPath = '/vols/Scratch/abaram/MATLAB/spm12';
addpath(spmPath);

rootData    = '/vols/Scratch/mgarvert/ManyMaps/imagingData/';

ROIs = {'322_both_07_left_parahippoc_2p5_mask', '322_diff_04_mPFC_2p5'};

% only subjects with both sessions
subjects = {'Subj_1','Subj_2','Subj_3','Subj_4','Subj_5','Subj_6','Subj_7','Subj_8',...
    'Subj_9','Subj_10','Subj_11','Subj_12','Subj_13','Subj_14','Subj_16',...
    'Subj_19','Subj_20','Subj_21','Subj_22','Subj_23','Subj_24','Subj_25' };

sessions = {'session_1','session_2'};


%% Run RSA
for iROI = 1:length(ROIs)
    for iSub=1:length(subjects)
        for iSess=1:2
            disp([iSub iSess])
            rsa_roi_rdm_withcollapse(rootData, subjects{iSub}, ROIs{iROI}, 'correlation', sessions{iSess});
        end
    end
end

stats_roi(rootData, subjects, sessions, '322_both_07_left_parahippoc_2p5_mask', 'relIrrel_sameMap') % Fig S2G
stats_roi(rootData, subjects, sessions, '322_both_07_left_parahippoc_2p5_mask', 'relIrrel_sameMapNoSameCond') % Fig S2H
stats_roi(rootData, subjects, sessions, '322_diff_04_mPFC_2p5', 'distRel_diffMaps') % Fig 4H
