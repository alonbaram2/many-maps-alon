% these are from groupStatsWilcoxon
subjects = {'Subj_1','Subj_2','Subj_3','Subj_4','Subj_5','Subj_6','Subj_7','Subj_8',...
   'Subj_9','Subj_10','Subj_11','Subj_12','Subj_13','Subj_14','Subj_16',...
  'Subj_19','Subj_20','Subj_21','Subj_22','Subj_23','Subj_24','Subj_25' };
distType='correlation';
rootData    = '/vols/Scratch/mgarvert/ManyMaps/imagingData/';
example_subjDir = fullfile(rootData,'rsa_alon','Subj_1','statistics',distType,'session_1');
analyses_fnames = dir(fullfile(example_subjDir,'*.mgh'));
example_subjFile = MRIread(fullfile(example_subjDir,analyses_fnames(1).name));
mkdir(fullfile(rootData,'rsa_alon','groupStats',distType,'Wilcoxon'));
sessions = {'session_1','session_2','both','diff'};

%%

fname = ['distRel_bothMaps_within_smth5_rh.mgh']

iSess=1;
allSubjData = nan(length(subjects),length(example_subjFile.vol));
P = nan(1,length(example_subjFile.vol));
for iSubj= 1:length(subjects)
subjFile = MRIread(fullfile(rootData,'rsa_alon',subjects{iSubj},'statistics',distType,sessions{iSess},fname));
allSubjData(iSubj,:) = subjFile.vol;
end
figure; hold on
title(fname(1:end-4))
h = histogram(allSubjData(:)); hold on; 
plot([0,0],[0,max(h.Values)],'r--'); 
plot([median(allSubjData(:)),median(allSubjData(:))],[0,max(h.Values)],'g--')
legend({'kendalTauCoeffs','0','median'})
