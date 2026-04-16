clear; clc;

% ============================
% SETTINGS
% ============================

baseDir = "/vols/Scratch/mgarvert/ManyMaps/imagingData";
behDir  = "/vols/Scratch/mgarvert/ManyMaps/scan_1.1"; %#ok<NASGU>

subjPattern = "Subj_*";
sessPattern = "session_*";

firstLevelDirName = "1stLevel";
design = "design_322_fsl_";

nblocks = 4;          % runs per SPM

% =========================================================
% WHICH REGRESSORS TO CORRELATE?
% =========================================================
%
% In the GLM script, the all_objects pmods are named:
%   dist_rel_map<m>_switch<s>
%   dist_irrel_map<m>_switch<s>
%
% Define the 4 combinations
obj_pair_labels = ["map1_stay","map1_switch","map2_stay","map2_switch"];

% NOTE: "stay" corresponds to switch=0, "switch" corresponds to switch=1
rel_patterns   = ["dist_rel_map1_switch0",  "dist_rel_map1_switch1",  "dist_rel_map2_switch0",  "dist_rel_map2_switch1"];
irrel_patterns = ["dist_irrel_map1_switch0","dist_irrel_map1_switch1","dist_irrel_map2_switch0","dist_irrel_map2_switch1"];

% ============================
% FIND SUBJECTS
% ============================

subjDirs = dir(fullfile(baseDir, subjPattern));
subjDirs = subjDirs([subjDirs.isdir]);

corrMaps = nan(numel(subjDirs),1);

results = struct([]);
rawResults = struct([]);

row = 0;
rowRaw = 0;

fprintf("Found %d subject folders.\n", numel(subjDirs));

% ============================
% MAIN LOOP
% ============================

for s = 1:numel(subjDirs)

    subjName = string(subjDirs(s).name);

    % ============================
    % correlate distance maps
    % ============================
    try
        load(strtrim(ls(["/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/" + subjName + "/*_session_1/data_" + extractAfter(subjName,"Subj_") + "_1.mat"])))
        corrMaps(s,1) = corr(squareform(data.map{1,1})', squareform(data.map{1,2})');
    catch
        corrMaps(s,1) = nan;
    end

    % ============================
    % get data for correlating regressors
    % ============================

    subjPath = fullfile(baseDir, subjName);
    sessDirs = dir(fullfile(subjPath, sessPattern));
    sessDirs = sessDirs([sessDirs.isdir]);

    if isempty(sessDirs)
        fprintf("WARNING: No sessions found for %s\n", subjName);
        continue;
    end

    for ss = 1:numel(sessDirs)

        sessName = string(sessDirs(ss).name);
        sessPath = fullfile(subjPath, sessName);

        designPath = fullfile(sessPath, firstLevelDirName, design);
        spmMatPath = fullfile(designPath, "SPM.mat");

        if ~isfile(spmMatPath)
            fprintf("WARNING: Missing SPM.mat: %s\n", spmMatPath);
            continue;
        end

        S = load(spmMatPath);
        SPM = S.SPM;
        X = SPM.xX.X;

        if ~isfield(SPM, "xX") || ~isfield(SPM.xX, "name")
            fprintf("WARNING: SPM.xX.name missing in %s\n", spmMatPath);
            continue;
        end

        % Column names as string array for robust matching
        xnames = string(SPM.xX.name(:));

        for run = 1:nblocks

            runTag = "Sn(" + string(run) + ")";  % how SPM prefixes per-session columns

            for k = 1:numel(obj_pair_labels)

                % =========================================================
                % HRF-convolved parametric modulators
                % =========================================================
                g_rel = find(startsWith(xnames, runTag) & contains(lower(xnames), lower(rel_patterns(k))));
                g_irr = find(startsWith(xnames, runTag) & contains(lower(xnames), lower(irrel_patterns(k))));

                if isempty(g_rel) || isempty(g_irr)
                    fprintf("WARNING: Missing columns for %s | %s | run %d | %s (rel found=%d, irrel found=%d)\n", ...
                        subjName, sessName, run, obj_pair_labels(k), numel(g_rel), numel(g_irr));
                    r_obj = nan;
                else
                    if numel(g_rel) > 1
                        fprintf("WARNING: Multiple rel matches for %s | %s | run %d | %s. Using first.\n", ...
                            subjName, sessName, run, obj_pair_labels(k));
                        g_rel = g_rel(1);
                    end

                    if numel(g_irr) > 1
                        fprintf("WARNING: Multiple irrel matches for %s | %s | run %d | %s. Using first.\n", ...
                            subjName, sessName, run, obj_pair_labels(k));
                        g_irr = g_irr(1);
                    end

                    r_obj = corr(X(:,g_rel), X(:,g_irr), "type","Pearson", "rows","pairwise");
                end

                row = row + 1;
                results(row).subject = subjName;
                results(row).session = sessName;
                results(row).spmMat  = string(spmMatPath);
                results(row).run     = run;
                results(row).label   = obj_pair_labels(k);
                results(row).colRel  = g_rel;
                results(row).colIrr  = g_irr;
                results(row).r       = r_obj;

                % =========================================================
                % RAW parametric modulators (before convolution)
                % =========================================================
                r_raw = nan;
                raw_rel = [];
                raw_irr = [];

                if isfield(SPM, "Sess") && numel(SPM.Sess) >= run && isfield(SPM.Sess(run), "U") && ~isempty(SPM.Sess(run).U)
                    % Search over conditions within the run
                    for u = 1:numel(SPM.Sess(run).U)
                        U = SPM.Sess(run).U(u);

                        if ~isfield(U, "P") || isempty(U.P)
                            continue;
                        end

                        % P may be struct array with field "name"
                        if isfield(U.P, "name")
                            pNames = string({U.P.name});
                        else
                            continue;
                        end

                        relIdx = find(contains(lower(pNames), lower(rel_patterns(k))));
                        irrIdx = find(contains(lower(pNames), lower(irrel_patterns(k))));

                        if ~isempty(relIdx) && isempty(raw_rel)
                            if isfield(U.P(relIdx(1)), "P")
                                raw_rel = U.P(relIdx(1)).P(:);
                            end
                        end

                        if ~isempty(irrIdx) && isempty(raw_irr)
                            if isfield(U.P(irrIdx(1)), "P")
                                raw_irr = U.P(irrIdx(1)).P(:);
                            end
                        end
                    end

                    if ~isempty(raw_rel) && ~isempty(raw_irr)
                        if numel(raw_rel) == numel(raw_irr)
                            r_raw = corr(raw_rel, raw_irr, "type","Pearson", "rows","pairwise");
                        else
                            fprintf("WARNING: Raw pmod length mismatch for %s | %s | run %d | %s (rel=%d, irrel=%d)\n", ...
                                subjName, sessName, run, obj_pair_labels(k), numel(raw_rel), numel(raw_irr));
                        end
                    else
                        fprintf("WARNING: Missing raw pmods for %s | %s | run %d | %s\n", ...
                            subjName, sessName, run, obj_pair_labels(k));
                    end
                else
                    fprintf("WARNING: SPM.Sess(%d).U missing for %s | %s\n", run, subjName, sessName);
                end

                rowRaw = rowRaw + 1;
                rawResults(rowRaw).subject = subjName;
                rawResults(rowRaw).session = sessName;
                rawResults(rowRaw).spmMat  = string(spmMatPath);
                rawResults(rowRaw).run     = run;
                rawResults(rowRaw).label   = obj_pair_labels(k);
                rawResults(rowRaw).r       = r_raw;
            end
        end

        fprintf("OK: %s | %s (processed)\n", subjName, sessName);
    end
end

% ============================
% TABLES + SUMMARY
% ============================

T = struct2table(results);
Traw = struct2table(rawResults);

% ----------------------------
% HRF-convolved summary
% ----------------------------
if ~isempty(T)
    % Participant-level mean correlation (average across all available runs x labels)
    Gsub = groupsummary(T, "subject", "mean", "r");

    idx = strcmp(Gsub.Properties.VariableNames, 'mean_r');
    if any(idx)
        Gsub.Properties.VariableNames{idx} = 'r_mean_withinSubj';
    end

    % Overall summary across participants
    rsub = Gsub.r_mean_withinSubj;

    fprintf("\n============================\n");
    fprintf("PARTICIPANT-LEVEL SUMMARY (mean r per participant)\n");
    fprintf("============================\n");
    fprintf("N participants = %d\n", numel(rsub));
    fprintf("Mean = %.4f\n", mean(rsub, 'omitnan'));
    fprintf("Std  = %.4f\n", std(rsub,  'omitnan'));
    fprintf("Min  = %.4f\n", min(rsub,  [], 'omitnan'));
    fprintf("Max  = %.4f\n", max(rsub,  [], 'omitnan'));
    fprintf("Max |r| = %.4f\n", max(abs(rsub), [], 'omitnan'));

    % participant-level means separately per label
    GsubLabel = groupsummary(T, ["subject","label"], "mean", "r");
    GlabelAcrossSubj = groupsummary(GsubLabel, "label", ["mean","std","min","max"], "mean_r");

    fprintf("\n============================\n");
    fprintf("LABEL-WISE SUMMARY (participant means, then summarized)\n");
    fprintf("============================\n");
    for i = 1:height(GlabelAcrossSubj)
        fprintf("%s: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n", ...
            string(GlabelAcrossSubj.label(i)), ...
            GlabelAcrossSubj.mean_mean_r(i), ...
            GlabelAcrossSubj.std_mean_r(i), ...
            GlabelAcrossSubj.min_mean_r(i), ...
            GlabelAcrossSubj.max_mean_r(i));
    end

    % 2) Participant x session mean (average within each session, per participant)
    GsubSess = groupsummary(T, ["subject","session"], "mean", "r");  % -> mean_r

    % 3) Session-level summary across participants
    sessions = unique(GsubSess.session);

    fprintf("\n============================\n");
    fprintf("SESSION-LEVEL SUMMARY (participant means within session, then across participants)\n");
    fprintf("============================\n");

    for i = 1:numel(sessions)
        sess = sessions(i);
        ii = GsubSess.session == sess;
        vals = GsubSess.mean_r(ii);

        fprintf("Session %s | N participants = %d | Mean=%.4f | Std=%.4f | Min=%.4f | Max=%.4f | Max|r|=%.4f\n", ...
            sess, sum(ii), ...
            mean(vals,'omitnan'), std(vals,'omitnan'), ...
            min(vals,[],'omitnan'), max(vals,[],'omitnan'), ...
            max(abs(vals),[],'omitnan'));
    end
else
    Gsub = table();
    GsubSess = table();
    fprintf("No HRF-convolved correlations found.\n");
end

% ----------------------------
% RAW summary
% ----------------------------
if ~isempty(Traw)
    GsubRaw = groupsummary(Traw, "subject", "mean", "r");

    idx = strcmp(GsubRaw.Properties.VariableNames, 'mean_r');
    if any(idx)
        GsubRaw.Properties.VariableNames{idx} = 'r_mean_withinSubj';
    end

    GsubSessRaw = groupsummary(Traw, ["subject","session"], "mean", "r");

    rsubRaw = GsubRaw.r_mean_withinSubj;

    fprintf("\n============================\n");
    fprintf("RAW PARTICIPANT-LEVEL SUMMARY (mean r per participant)\n");
    fprintf("============================\n");
    fprintf("N participants = %d\n", numel(rsubRaw));
    fprintf("Mean = %.4f\n", mean(rsubRaw, 'omitnan'));
    fprintf("Std  = %.4f\n", std(rsubRaw,  'omitnan'));
    fprintf("Min  = %.4f\n", min(rsubRaw,  [], 'omitnan'));
    fprintf("Max  = %.4f\n", max(rsubRaw,  [], 'omitnan'));
    fprintf("Max |r| = %.4f\n", max(abs(rsubRaw), [], 'omitnan'));

    GsubLabelRaw = groupsummary(Traw, ["subject","label"], "mean", "r");
    GlabelAcrossSubjRaw = groupsummary(GsubLabelRaw, "label", ["mean","std","min","max"], "mean_r");

    fprintf("\n============================\n");
    fprintf("RAW LABEL-WISE SUMMARY (participant means, then summarized)\n");
    fprintf("============================\n");
    for i = 1:height(GlabelAcrossSubjRaw)
        fprintf("%s: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n", ...
            string(GlabelAcrossSubjRaw.label(i)), ...
            GlabelAcrossSubjRaw.mean_mean_r(i), ...
            GlabelAcrossSubjRaw.std_mean_r(i), ...
            GlabelAcrossSubjRaw.min_mean_r(i), ...
            GlabelAcrossSubjRaw.max_mean_r(i));
    end
else
    GsubRaw = table();
    GsubSessRaw = table();
    fprintf("No RAW correlations found.\n");
end

% ============================
% SAVE TABLES
% ============================

outCsv_all      = fullfile(baseDir, "pm_corr_all_objects_allRows.csv");
outCsv_sub      = fullfile(baseDir, "pm_corr_all_objects_participantMeans_allSessionsPooled.csv");
outCsv_subSess  = fullfile(baseDir, "pm_corr_all_objects_participantMeans_bySession.csv");

writetable(T, outCsv_all);
writetable(Gsub, outCsv_sub);
writetable(GsubSess, outCsv_subSess);

fprintf("\nSaved row-level table to:\n%s\n", outCsv_all);
fprintf("Saved participant-mean (pooled) table to:\n%s\n", outCsv_sub);
fprintf("Saved participant-mean (by session) table to:\n%s\n", outCsv_subSess);

if ~isempty(Traw)
    outCsv_raw_all     = fullfile(baseDir, "pm_corr_all_objects_RAW_allRows.csv");
    outCsv_raw_sub     = fullfile(baseDir, "pm_corr_all_objects_RAW_participantMeans_allSessionsPooled.csv");
    outCsv_raw_subSess = fullfile(baseDir, "pm_corr_all_objects_RAW_participantMeans_bySession.csv");

    writetable(Traw, outCsv_raw_all);
    writetable(GsubRaw, outCsv_raw_sub);
    writetable(GsubSessRaw, outCsv_raw_subSess);

    fprintf("Saved RAW row-level table to:\n%s\n", outCsv_raw_all);
    fprintf("Saved RAW participant-mean (pooled) table to:\n%s\n", outCsv_raw_sub);
    fprintf("Saved RAW participant-mean (by session) table to:\n%s\n", outCsv_raw_subSess);
end

%% 
% ============================
% PLOTTING
% ============================
try
    % --- Boxplot of corrMaps (distance-map correlations)
    figure;
    validIdx = ~isnan(corrMaps);
    boxplot(corrMaps(validIdx));
    hold on;
    x = ones(sum(validIdx),1);
    jitter = (rand(sum(validIdx),1)-0.5)*0.2;
    scatter(x + jitter, corrMaps(validIdx), 30, 'k','filled','MarkerFaceAlpha',0.6);
    hold off;
    ylabel('Pearson r between graph distances');
    ylim([-0.1, 0.1])
    xlim([0.7, 1.3])
    title('Graph distances were orthogonalised');
    grid off;

    % ---------------------------------------------------------
    % Histogram: HRF-convolved regressor correlations
    % (session-separated, using subject x session averages)
    % ---------------------------------------------------------
    if ~isempty(T) && ~isempty(GsubSess)

        figure;
        hold on;

        sessNames = unique(string(GsubSess.session));

        if isempty(sessNames)
            warning('No session means to plot.');
        else
            allr = GsubSess.mean_r;
            allr = allr(~isnan(allr));

            if isempty(allr)
                warning('No HRF-convolved regressor correlations to plot (all NaN).');
            else
                nbins = 20;
                edges = linspace(min(allr), max(allr), nbins+1);

                colors = lines(numel(sessNames));
                legendEntries = strings(0);

                for i = 1:numel(sessNames)
                    sess = sessNames(i);
                    vals = GsubSess.mean_r(string(GsubSess.session) == sess);
                    vals = vals(~isnan(vals));

                    if ~isempty(vals)
                        histogram(vals, edges, ...
                            'Normalization','probability', ...
                            'FaceAlpha',0.5, ...
                            'FaceColor', colors(i,:));
                        legendEntries(end+1) = sess; %#ok<AGROW>
                    end
                end

                xlabel('Mean Pearson r per subject \times session');
                ylabel('Proportion');
                title('HRF-convolved parametric modulator correlations');
                if ~isempty(legendEntries)
                    legend(cellstr(legendEntries), 'Location','best');
                end
                grid off;
            end
        end

        hold off;
    end

    % ---------------------------------------------------------
    % Histogram: RAW regressor correlations
    % (session-separated, using subject x session averages)
    % ---------------------------------------------------------
    if ~isempty(Traw) && ~isempty(GsubSessRaw)

        figure;
        hold on;

        sessNamesRaw = unique(string(GsubSessRaw.session));

        if isempty(sessNamesRaw)
            warning('No raw session means to plot.');
        else
            allr = GsubSessRaw.mean_r;
            allr = allr(~isnan(allr));

            if isempty(allr)
                warning('No RAW regressor correlations to plot (all NaN).');
            else
                nbins = 20;
                edges = linspace(min(allr), max(allr), nbins+1);

                colors = lines(numel(sessNamesRaw));
                legendEntries = strings(0);

                for i = 1:numel(sessNamesRaw)
                    sess = sessNamesRaw(i);
                    vals = GsubSessRaw.mean_r(string(GsubSessRaw.session) == sess);
                    vals = vals(~isnan(vals));

                    if ~isempty(vals)
                        histogram(vals, edges, ...
                            'Normalization','probability', ...
                            'FaceAlpha',0.5, ...
                            'FaceColor', colors(i,:));
                        legendEntries(end+1) = sess; %#ok<AGROW>
                    end
                end

                xlabel('Mean Pearson r per subject \times session');
                ylabel('Proportion');
                title('Raw parametric modulator correlations');
                if ~isempty(legendEntries)
                    legend(cellstr(legendEntries), 'Location','best');
                end
                grid off;
            end
        end

        hold off;
    end

catch ME
    warning("Plotting failed: %s", ME.message);
end