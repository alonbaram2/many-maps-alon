function stackSubjectsStatsSurfMaps(rootData,subjects)
system('module add freesurfer')
sessions = {'diff','both','session_1','session_2'};
% get analyses names from example subject
exampleInputDir = fullfile(rootData,'rsa_alon','Subj_1','statistics','correlation','session_1');
cd(exampleInputDir)
tmpFname = dir('*xRun1324_smth5*'); % find all files that are of the right cross-validation - this can be changed if needed.
fname = cell(size(tmpFname));
for iAnalysis = 1:length(fname)
    fname{iAnalysis} = tmpFname(iAnalysis).name(1:end-4); % get rid of .mgh suffix
end

for iSess=1:length(sessions)
    outputDir = fullfile(rootData,'rsa_alon','allSubjStacked','correlation',sessions{iSess});
    mkdir(outputDir);
    for iAn = 1:length(fname)
    % if strcmp(sl,'surf')
            str = ['mri_concat '];
            for iSub = 1:length(subjects)
                str = [str fullfile(rootData,'rsa_alon',subjects{iSub},'statistics','correlation',sessions{iSess},[fname{iAn} '.mgh']) ' '];
            end
            str = [str '--o ' fullfile(outputDir,[fname{iAn} '_allSubj.mgh'])];
            system (str);           
    % else            
    %         str = ['fslmerge -t ' fullfile(outputDir,fname{iAn}) ' '];
    %         for iSub = 1:length(subjects)
    %             str = [str fullfile(rootData,'RDMs',glm,subjects{iSub},'rsaContrasts',distType,sl,[fname{iAn} '.nii']) ' '];
    %         end
    %         system (str);
    %         gunzip(fullfile(outputDir,[fname{iAn} '.nii.gz']));
    %         delete(fullfile(outputDir,[fname{iAn} '.nii.gz']))
    % end 
    end
end