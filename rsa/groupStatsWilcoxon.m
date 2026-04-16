function groupStatsWilcoxon (rootData,subjects,sl)

example_subjDir = fullfile(rootData,'rsa_alon','Subj_1','statistics',sl,'session_1');
if strcmp(sl, 'WhBr_surf_r10_v100')
    analyses_fnames = dir(fullfile(example_subjDir,'*.mgh'));
else
    analyses_fnames = dir(fullfile(example_subjDir,'*.nii'));
end
example_subjFile = MRIread(fullfile(example_subjDir,analyses_fnames(1).name));
mkdir(fullfile(rootData,'rsa_alon','groupStats',sl,'Wilcoxon'));

sessions = {'session_1','session_2','both','diff'};
if strcmp(sl,'WhBr_surf_r10_v100')
    for iFiles = 1:length(analyses_fnames)
        for iSess = 1:length(sessions)
            allSubjData = nan(length(subjects),length(example_subjFile.vol));
            P = nan(1,length(example_subjFile.vol));
            for iSubj=1:length(subjects)
                subjFile = MRIread(fullfile(rootData,'rsa_alon',subjects{iSubj},'statistics',sl,sessions{iSess},analyses_fnames(iFiles).name));
                allSubjData(iSubj,:) = subjFile.vol;
            end
            for iVox = 1:size(allSubjData,2)
                [P(iVox),~,~] = signrank(allSubjData(:,iVox),[],'tail','right');
            end
            fname_save = ['P_' sessions{iSess} '_' analyses_fnames(iFiles).name];
            fname_fullpath = fullfile(rootData,'rsa_alon','groupStats',sl,'Wilcoxon',fname_save);
            surfToSave = example_subjFile;
            surfToSave.vol = single(1-P);
            surfToSave.fspec = fname_fullpath;
            MRIwrite(surfToSave,fname_fullpath)
        end
    end
else
    for iFiles = 1:length(analyses_fnames)
        for iSess = 1:length(sessions)
            allSubjData = nan(length(subjects),size(example_subjFile.vol,1),size(example_subjFile.vol,2),size(example_subjFile.vol,3));
            P = nan(size(example_subjFile.vol));
            for iSubj=1:1%length(subjects)
                subjFile = MRIread(fullfile(rootData,'rsa_alon',subjects{iSubj},'statistics',sl,sessions{iSess},analyses_fnames(iFiles).name));
                allSubjData(iSubj,:,:,:) = subjFile.vol;
            end
            voxVec = find(~isnan(example_subjFile.vol));
            for iVox = 1:length(voxVec)
                [P(iVox),~,~] = signrank(allSubjData(:, example_subjFile.vol, iVox),[],'tail','right');
            end
            fname_save = ['P_' sessions{iSess} '_' analyses_fnames(iFiles).name];
            fname_fullpath = fullfile(rootData,'rsa_alon','groupStats',sl,'Wilcoxon',fname_save);
            surfToSave = example_subjFile;
            surfToSave.vol = single(1-P);
            surfToSave.fspec = fname_fullpath;
            MRIwrite(surfToSave,fname_fullpath)
        end
    end
end