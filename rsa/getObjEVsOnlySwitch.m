root = '/vols/Scratch/mgarvert/ManyMaps/imagingData/';

subjects = {'Subj_1','Subj_2','Subj_3','Subj_4','Subj_5','Subj_6','Subj_7','Subj_8',...
   'Subj_9','Subj_10','Subj_11','Subj_12','Subj_13','Subj_14','Subj_16',...
  'Subj_19','Subj_20','Subj_21','Subj_22','Subj_23','Subj_24','Subj_25' };

sessions = {'session_1','session_2'};

emptyRegs = false(length(subjects),length(sessions),4,34,2); % nSubj x nSess x nRuns x nObj x stay/switch

for switchTrial=0:1
    for iSubj = 1:length(subjects)
        for iSess = 1:2
            for iRun=1:4
                EVsDir = fullfile(root,subjects{iSubj},sessions{iSess},'behaviour','EVs',[subjects{iSubj} '_' sessions{iSess} '_run_' ,num2str(iRun)]);
                % Load reference files and extract first column
                ref1 = load(fullfile(EVsDir,['all_stimuli_map1switch_' num2str(switchTrial) '_distRel.txt']));
                ref2 = load(fullfile(EVsDir,['all_stimuli_map2switch_' num2str(switchTrial) '_distRel.txt']));
                ref1_vals = ref1(:,1);
                ref2_vals = ref2(:,1);
                for iObj=1:34
                    fname = fullfile(EVsDir,['obj_' num2str(iObj) '.txt']);
                    % Load the current obj_*.txt file
                    data = load(fname);

                    % Find values in first column that exist in either reference files
                    common_vals = data(ismember(data(:,1), ref1_vals) | ismember(data(:,1), ref2_vals), 1);

                    % Construct filtered output: keep only matching rows, fix cols 2 and 3
                    filtered_data = [common_vals, ...
                        repmat([2.000000, 1.000000], length(common_vals), 1)];

                    % Write new file "
                    output_fname = fullfile(EVsDir,['obj_' num2str(iObj) '_switch' num2str(switchTrial) '.txt']);
                    if isempty(filtered_data)
                        warning(['Empty reg: ' output_fname])
                        emptyRegs(iSubj,iSess,iRun,iObj,switchTrial+1) = true;
                    end
                    fid = fopen(output_fname, 'w');
                    fprintf(fid, '%.6f %.6f %.6f\n', filtered_data');
                    fclose(fid);
                end
            end
        end
    end
end

%%

% Collapse over the sessions dimension (size 2), using logical OR (any) - 
% if one of the sessions is empty that's enough to disqualify

B = any(emptyRegs, 2);  % Resulting size will be 22 x 1 x 4 x 34 x 2
% Assume A is 22 x 1 x 4 x 34 x 2

% Logical OR within each pair
pair1 = B(:,:,1,:,:) | B(:,:,3,:,:);  % size: 22 x 1 x 1 x 34 x 2
pair2 = B(:,:,2,:,:) | B(:,:,4,:,:);  % size: 22 x 1 x 1 x 34 x 2

% Logical AND across pairs (both must be true)
C = pair1 & pair2;  % size: 22 x 1 x 1 x 34 x 2

% Remove the singleton 3rd dimension
C = squeeze(C);  % Resulting size: 22 x 34 x 2
