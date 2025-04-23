function groupStatsWilcoxon (rootData,subjects,distType)

example_subjDir = fullfile(rootData,'rsa_alon','Subj_1','statistics',distType,'session_1');
analyses_fnames = dir(fullfile(example_subjDir,'*.mgh'));
example_subjFile = MRIread(fullfile(example_subjDir,analyses_fnames(1).name));
mkdir(fullfile(rootData,'rsa_alon','groupStats',distType,'Wilcoxon'));

sessions = {'session_1','session_2','both','diff'};
for iFiles = 1:length(analyses_fnames)
    for iSess = 1:length(sessions)
        allSubjData = nan(length(subjects),length(example_subjFile.vol));
        P = nan(1,length(example_subjFile.vol));
        for iSubj=1:length(subjects)
            subjFile = MRIread(fullfile(rootData,'rsa_alon',subjects{iSubj},'statistics',distType,sessions{iSess},analyses_fnames(iFiles).name));                
            allSubjData(iSubj,:) = subjFile.vol;
        end
        for iVer = 1:size(allSubjData,2)
            [P(iVer),~,~] = signrank(allSubjData(:,iVer),[],'tail','right');
        end
        fname_save = ['P_' sessions{iSess} '_' analyses_fnames(iFiles).name];
        fname_fullpath = fullfile(rootData,'rsa_alon','groupStats',distType,'Wilcoxon',fname_save);
        surfToSave = example_subjFile;
        surfToSave.vol = single(1-P);
        surfToSave.fspec = fname_fullpath;
        MRIwrite(surfToSave,fname_fullpath)
    end
end
