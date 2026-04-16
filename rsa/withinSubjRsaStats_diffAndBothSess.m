function withinSubjRsaStats_diffAndBothSess(rootData,subj, sl)

inputs_dir = fullfile(rootData,'rsa_alon',subj,'statistics',sl);
mkdir(fullfile(rootData,'rsa_alon',subj,'statistics',sl ,'both'));
mkdir(fullfile(rootData,'rsa_alon',subj,'statistics',sl,'diff'));

if strcmp(sl,'WhBr_surf_r10_v100') %surface
    fnames = dir(fullfile(inputs_dir,'session_1','*.mgh'));
else %ROI volume
    fnames = dir(fullfile(inputs_dir,'session_1','*.nii'));
end
for iFiles = 1:length(fnames)
    sess1 = MRIread(fullfile(inputs_dir,'session_1',fnames(iFiles).name));
    sess2 = MRIread(fullfile(inputs_dir,'session_2',fnames(iFiles).name));
    sess_both = sess1; % header info the same
    sess_both.vol = single(sess1.vol + sess2.vol / 2);
    sess_both.fspec = fullfile(rootData,'rsa_alon',subj,'statistics',sl,'both',fnames(iFiles).name);
    MRIwrite(sess_both,fullfile(rootData,'rsa_alon',subj,'statistics',sl,'both',fnames(iFiles).name));
    sess_diff = sess1; % header info the same
    sess_diff.vol = single(sess2.vol - sess1.vol);
    sess_diff.fspec = fullfile(rootData,'rsa_alon',subj,'statistics',sl,'diff',fnames(iFiles).name);
    MRIwrite(sess_diff,fullfile(rootData,'rsa_alon',subj,'statistics',sl,'diff',fnames(iFiles).name));
end

