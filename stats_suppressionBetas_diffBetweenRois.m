%% =========================================================
%  ROI beta difference tests: mask1 - mask2
%  - 4 one-sided one-sample t-tests vs 0
%  - 6 paired t-tests between all pairs of cells
%  - Holm and Bonferroni corrections
%  - 2x2 repeated-measures ANOVA
%% =========================================================

clear; clc;
clear functions;

mask2 = 'harvardoxford_HPC_bilateral';
mask1 = 'juelich_bothERC_thr25';

root1 = ['/vols/Scratch/mgarvert/ManyMaps/imagingData/2ndLevel/design_322_fsl_/mask/' mask1];
root2 = ['/vols/Scratch/mgarvert/ManyMaps/imagingData/2ndLevel/design_322_fsl_/mask/' mask2];

rel_con = '02_rel_dist_stay';
irrel_con = '03_irrel_dist_stay';

nSub = 25;
nSess = 2;
nGraph = 2;

betas1 = nan(nSub, nSess, nGraph);
betas2 = nan(nSub, nSess, nGraph);

%% -----------------------------
%% Load beta values for both masks
%% -----------------------------
for iSub = 1:nSub
    for iSess = 1:nSess
        try
            rel_file1 = fullfile(root1, rel_con, ...
                [sprintf('%02d', iSub) '_' num2str(iSess) ...
                '_' mask1 '_' rel_con '.txt']);

            irrel_file1 = fullfile(root1, irrel_con, ...
                [sprintf('%02d', iSub) '_' num2str(iSess) ...
                '_' mask1 '_' irrel_con '.txt']);

            rel_file2 = fullfile(root2, rel_con, ...
                [sprintf('%02d', iSub) '_' num2str(iSess) ...
                '_' mask2 '_' rel_con '.txt']);

            irrel_file2 = fullfile(root2, irrel_con, ...
                [sprintf('%02d', iSub) '_' num2str(iSess) ...
                '_' mask2 '_' irrel_con '.txt']);

            betas1(iSub, iSess, 1) = load(rel_file1);
            betas1(iSub, iSess, 2) = load(irrel_file1);

            betas2(iSub, iSess, 1) = load(rel_file2);
            betas2(iSub, iSess, 2) = load(irrel_file2);

        catch
            % leave NaN if missing/unreadable
        end
    end
end

%% -----------------------------
%% Compute difference scores
%% -----------------------------
% Difference = mask1 - mask2
betas = betas1 - betas2;

%% -----------------------------
%% Reshape into 4 cells
%% Order:
%%   1 = S1-Rel
%%   2 = S1-Irrel
%%   3 = S2-Rel
%%   4 = S2-Irrel
%% -----------------------------
cell_data = {
    squeeze(betas(:,1,1)), ...
    squeeze(betas(:,1,2)), ...
    squeeze(betas(:,2,1)), ...
    squeeze(betas(:,2,2))
};

cell_names = {'S1-Rel','S1-Irrel','S2-Rel','S2-Irrel'};
nCells = numel(cell_names);

%% =========================================================
%% 1) One-sample one-sided t-tests vs 0
%% =========================================================
fprintf('\n=== One-sample t-tests on %s - %s vs 0 ===\n', mask1, mask2);

t_one = nan(nCells,1);
p_one = nan(nCells,1);
n_one = nan(nCells,1);

for i = 1:nCells
    x = cell_data{i};
    x = x(~isnan(x));
    n_one(i) = numel(x);

    [~, p, ~, stats] = ttest(x, 0, 'Tail', 'right'); % H1: mean > 0
    t_one(i) = stats.tstat;
    p_one(i) = p;

    fprintf('%s: t = %.4f, p = %.6f\n', cell_names{i}, t_one(i), p_one(i));
end

% Bonferroni
m1 = numel(p_one);
p_one_bonf = min(1, p_one * m1);

% Holm
[p_sorted, idx] = sort(p_one(:), 'ascend');
m = numel(p_sorted);
p_holm_sorted = min(1, p_sorted .* (m:-1:1)');

for k = 2:m
    p_holm_sorted(k) = max(p_holm_sorted(k), p_holm_sorted(k-1));
end

p_one_holm = nan(m,1);
p_one_holm(idx) = p_holm_sorted;
p_one_holm = reshape(p_one_holm, size(p_one));

fprintf('\nHolm-corrected p-values (vs 0):\n');
for i = 1:nCells
    fprintf('%s: p_corr = %.6f\n', cell_names{i}, p_one_holm(i));
end

fprintf('\nBonferroni-corrected p-values (vs 0):\n');
for i = 1:nCells
    fprintf('%s: p_corr = %.6f\n', cell_names{i}, p_one_bonf(i));
end

T_one = table(cell_names(:), n_one, t_one, p_one, p_one_holm, p_one_bonf, ...
    'VariableNames', {'Cell','N','tStat','pRaw','pHolm','pBonf'});

disp(' ');
disp('=== Summary: one-sample tests vs 0 ===');
disp(T_one);

%% =========================================================
%% 2) Pairwise paired t-tests between all 4 cells
%% =========================================================
fprintf('\n=== Pairwise paired t-tests on %s - %s ===\n', mask1, mask2);

pairs = nchoosek(1:nCells, 2);
nPairs = size(pairs,1);

t_pair = nan(nPairs,1);
p_pair = nan(nPairs,1);
n_pair = nan(nPairs,1);
pair_names = strings(nPairs,1);

for i = 1:nPairs
    a = pairs(i,1);
    b = pairs(i,2);

    x = cell_data{a};
    y = cell_data{b};

    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);
    n_pair(i) = numel(x);

    [~, p, ~, stats] = ttest(x, y); % paired, two-sided
    t_pair(i) = stats.tstat;
    p_pair(i) = p;

    pair_names(i) = cell_names{a} + " vs " + cell_names{b};
    fprintf('%s: t = %.4f, p = %.6f\n', pair_names(i), t_pair(i), p_pair(i));
end

% Bonferroni
m2 = numel(p_pair);
p_pair_bonf = min(1, p_pair * m2);

% Holm
[p_sorted, idx] = sort(p_pair(:), 'ascend');
m = numel(p_sorted);
p_holm_sorted = min(1, p_sorted .* (m:-1:1)');

for k = 2:m
    p_holm_sorted(k) = max(p_holm_sorted(k), p_holm_sorted(k-1));
end

p_pair_holm = nan(m,1);
p_pair_holm(idx) = p_holm_sorted;
p_pair_holm = reshape(p_pair_holm, size(p_pair));

fprintf('\nHolm-corrected p-values (pairwise):\n');
for i = 1:nPairs
    fprintf('%s: p_corr = %.6f\n', pair_names(i), p_pair_holm(i));
end

fprintf('\nBonferroni-corrected p-values (pairwise):\n');
for i = 1:nPairs
    fprintf('%s: p_corr = %.6f\n', pair_names(i), p_pair_bonf(i));
end

T_pair = table(pair_names, n_pair, t_pair, p_pair, p_pair_holm, p_pair_bonf, ...
    'VariableNames', {'Comparison','N','tStat','pRaw','pHolm','pBonf'});

disp(' ');
disp('=== Summary: pairwise tests ===');
disp(T_pair);

%% =========================================================
%% 3) 2x2 repeated-measures ANOVA
%%    Within-subject factors:
%%      Session = {S1, S2}
%%      Graph   = {Rel, Irrel}
%% =========================================================

Y = table( ...
    betas(:,1,1), betas(:,1,2), betas(:,2,1), betas(:,2,2), ...
    'VariableNames', {'S1_Rel','S1_Irrel','S2_Rel','S2_Irrel'} );

within = table( ...
    categorical({'S1'; 'S1'; 'S2'; 'S2'}), ...
    categorical({'Rel'; 'Irrel'; 'Rel'; 'Irrel'}), ...
    'VariableNames', {'Session','Graph'} );

rm = fitrm(Y, 'S1_Rel-S2_Irrel ~ 1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'Session*Graph');

disp(' ');
disp('=== 2x2 repeated-measures ANOVA ===');
disp(ranovatbl);

cell_means = squeeze(mean(betas, 'omitnan'));
disp(' ');
disp('=== Cell means (difference scores, subjects averaged) ===');
disp(array2table(cell_means(:)', ...
    'VariableNames', {'S1_Rel','S1_Irrel','S2_Rel','S2_Irrel'}));