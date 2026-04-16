root = '/vols/Scratch/mgarvert/ManyMaps/imagingData/';

subjects = {'Subj_1','Subj_2','Subj_3','Subj_4','Subj_5','Subj_6','Subj_7','Subj_8',...
   'Subj_9','Subj_10','Subj_11','Subj_12','Subj_13','Subj_14','Subj_16',...
  'Subj_19','Subj_20','Subj_21','Subj_22','Subj_23','Subj_24','Subj_25' };

sessions = {'session_1','session_2'};

emptyRegs = false(length(subjects),length(sessions),4,34,2); % nSubj x nSess x nRuns x nObj x stay/switch

for switchTrial=0:1
    clear bothMaps bothMaps_times
    for iSubj = 1:length(subjects)
        for iSess = 1:2
            for iRun=1:4
                EVsDir = fullfile(root,subjects{iSubj},sessions{iSess},'behaviour','EVs',[subjects{iSubj} '_' sessions{iSess} '_run_' ,num2str(iRun)]);
                % Load reference files and extract first column (we want to
                % collapse over maps, for a given switch/stay condition)        
                map1 = load(fullfile(EVsDir,['all_stimuli_map1switch_' num2str(switchTrial) '_distRel.txt']));
                map2 = load(fullfile(EVsDir,['all_stimuli_map2switch_' num2str(switchTrial) '_distRel.txt']));
                % first collapse over maps for all switch/stay trials an
                % save
                map1_times = map1(:,1);
                map2_times = map2(:,1);
                bothMaps_times = sortrows([map1_times;map2_times]);
                % fix cols 2 and 3. Note that col 2 (duration) is 2s
                % because for RSA analyses we take the entire duration
                % of the stim presentation (unlike for suppression
                % where duration i 1s)
                bothMaps = [bothMaps_times, ...
                    repmat([2.000000, 1.000000], length(bothMaps_times), 1)];
                output_fname = fullfile(EVsDir,['allObj_switch' num2str(switchTrial) '.txt']);
                fid = fopen(output_fname, 'w');
                fprintf(fid, '%.6f %.6f %.6f\n', bothMaps);
                fclose(fid);
                
                % now split to each obj
                for iObj=1:34
                    fname = fullfile(EVsDir,['obj_' num2str(iObj) '.txt']);
                    % Load the current obj_*.txt file
                    data = load(fname);

                    % Find values in first column that exist in either
                    % reference files: switch/stay trials of current obj
                    common_vals = data(ismember(data(:,1), map1_times) | ismember(data(:,1), map2_times), 1);

                    % Construct filtered output: keep only matching rows,
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



%% Find how many empty RDM elements per subj*obj*switch/stay

% Collapse over the sessions dimension (size 2), using logical OR (any) - 
% if one of the sessions is empty that's enough to disqualify

B = any(emptyRegs, 2);  % Resulting size will be 22 x 1 x 4 x 34 x 2
% Assume A is 22 x 1 x 4 x 34 x 2

% Logical OR within each pair of runs - e.g. if run 3 is empty, 
% tha'ts enough to count the runs pair 1-3 as empty
pair1 = B(:,:,1,:,:) | B(:,:,3,:,:);  % size: 22 x 1 x 1 x 34 x 2
pair2 = B(:,:,2,:,:) | B(:,:,4,:,:);  % size: 22 x 1 x 1 x 34 x 2

% Logical AND across runs pairs (only count as empty pairs objects where
% both pairs of runs 1-3 AND 2-4 are empty)
C = pair1 & pair2;  % size: 22 x 1 x 1 x 34 x 2

% Remove the singleton 3rd dimension
C = squeeze(C);  % Resulting size: 22 x 34 x 2
