clear
close all

% Remove subjects with no or only one session
ix = ones(25,1); ix([15,17,18]) = 0; ix = ix == 1;

whichPC = 'jalapeno';
if strcmp(whichPC,'jalapeno')
    basedir = '/home/fs0/mgarvert/scratch/ManyMaps/';
    basedir_beh = '/home/fs0/mgarvert/scratch/ManyMaps/imagingData';
    addpath(genpath(fullfile(basedir,'imagingData','scripts','alon','utilities')));
    dir = '/home/fs0/mgarvert/scratch/ManyMaps/imagingData/2ndLevel/design_322_fsl_/mask';    
elseif strcmp(whichPC,'oregano')
    basedir = '/home/raid1/garvert/tu_garvert_cloud/Projects/ManyMaps/';
    basedir_beh = '/home/raid1/garvert/tu_garvert_cloud/owncloud-gwdg/Projects/ManyMaps/';
    addpath (genpath('/data/tu_garvert_cloud/owncloud-gwdg/matlab_scripts'))
elseif strcmp(whichPC,'dell')
    basedir = 'C:/Users/mgarvert/ownCloudCBS/Projects/ManyMaps/';
    basedir_beh = 'C:/Users/mgarvert/ownCloudCBS/Projects/ManyMaps/';
    addpath C:/Users/mgarvert/ownCloudCBS/matlab_scripts
elseif strcmp(whichPC,'mac')
    basedir = '/data/tu_garvert_cloud/owncloud-gwdg/Projects/ManyMaps/';
    basedir_beh = '/data/tu_garvert_cloud//owncloud-gwdg/Projects/ManyMaps/';
    addpath (genpath('/data/tu_garvert_cloud/owncloud-gwdg/matlab_scripts'))
elseif strcmp(whichPC,'online')
    basedir = '/Users/garvert_UW/DropboxB/DateienMG/Privat/Projects/ManyMaps/';
    basedir_beh = "/Users/garvert_UW/DropboxB/DateienMG/Privat/Projects/ManyMaps/imagingData";
    addpath '/Users/garvert_UW/ownCloud/matlab_scripts'

    dir = '/Users/garvert_UW/DropboxB/DateienMG/Privat/Projects/ManyMaps/imagingData/design_322_fsl_/mask';
end

con{1} = '02_rel_dist_stay';
con{2} = '03_irrel_dist_stay';

roi{1} = '322_both_07_left_parahippoc_2p5_mask';
roi{2} = '322_1_07_left_parahippoc_2p5_mask';
roi{3} = '322_2_07_mPFC_2p5_mask';

% Number of ROIs
num_rois = length(roi);

% Placeholder for data
pe = nan(length(con), 25, 2, num_rois);

% Load data for each ROI
for r = 1:num_rois
    for subj = 1:25
        for session = 1:2
            for c = 1:size(con, 2)
                if exist(fullfile(dir, roi{r}, con{c}, sprintf('%d_%u_%s_%s.txt', subj, session, roi{r}, con{c})), 'file')
                    pe(c, subj, session, r) = load(fullfile(dir, roi{r}, con{c}, sprintf('%d_%u_%s_%s.txt', subj, session, roi{r}, con{c})));
                elseif exist(fullfile(dir, roi{r}, con{c}, sprintf('%02d_%u_%s_%s.txt', subj, session, roi{r}, con{c})), 'file')
                    pe(c, subj, session, r) = load(fullfile(dir, roi{r}, con{c}, sprintf('%02d_%u_%s_%s.txt', subj, session, roi{r}, con{c})));
                else
                    pe(c, subj, session, r) = nan;
                end
            end
        end
    end
end

% Plotting
colors = lbmap(7, 'RedBlue');
label = {'relevant', 'irrelevant'};


figure('Position', [100, 100, 1200, 400]);

for r = 1:num_rois
    subplot(1, num_rois, r);
    hold on;
    for session = 1:2
        for c = 1:size(con, 2)
            no_outlier = pe(c, :, session, r) < (nanmean(pe(c, :, session, r)) + 3*nanstd(pe(c, :, session, r))) & pe(c, :, session, r) > (nanmean(pe(c, :, session, r)) - 3*nanstd(pe(c, :, session, r)));
            bar(c + (session-1)*size(con, 2), nanmean(pe(c, no_outlier, session, r)), 'FaceColor', colors(c, :));
            v = c + (session-1)*size(con, 2) - 0.2 + 0.4*rand(sum(no_outlier), 1);
            scatter(v, pe(c, no_outlier, session, r), 35, colors(c, :), 'filled','MarkerEdgeColor', 'k');
            errorbar(c + (session-1)*size(con, 2), nanmean(pe(c, no_outlier, session, r)), nanstd(pe(c, no_outlier, session, r))/sqrt(22), 'k');
            [h, p] = ttest(pe(c, no_outlier, session, r));
        end
    end
    xticks(1:size(con, 2)*2);
    xtickangle(45);
    title(sprintf('ROI %s', roi{r}), 'Interpreter', 'none');
    ylabel('Effect size (a.u.)');
    ax = gca;
    axislength(r, :) = ax.YLim;

    ax.TickLabelInterpreter = 'none';
    ax.XTickLabel = {'relevant','irrelevant','relevant','irrelevant'};
end


% prepImg


%% Arena task

arena = nan(25,2);
overall_arena = nan(2,25,136);
position_arena = nan(2,25,17,2);
for subj = 1:25
    for map = 1:2
        if exist([basedir,'/scan_1.1/datafiles/arena/similarityJudgementData_Subj_',num2str(subj),'_map_',num2str(map),'_trial1.mat'],'file')
            a = load(([basedir,'/scan_1.1/datafiles/arena/similarityJudgementData_Subj_',num2str(subj),'_map_',num2str(map),'_trial1.mat']));
            load([basedir,'/scan_1.1/datafiles/alldata/data_',num2str(subj),'_',num2str(session),'.mat']);
            

            for i = 1:17
                order(i) = find(data.stimuli(2,:) == data.stimuli(1,i));
                order2(i) = find(data.stimuli(1,:) == data.stimuli(2,i));
            end
            
            if map == 1
                overall_arena(map,subj,:)   = a.distMat_ltv;
                pos_arena(map,subj,:,:)     = a.itemPositions;
                arena(subj,map) = corr(squareform(data.map{map})',a.distMat_ltv');
                
                sqfm = squareform(a.distMat_ltv);
                sqfm(order2,order2) = sqfm;                
%                 arena(subj,3) = corr(squareform(data.map{1})',squareform(sqfm)');
            
            else
                sqfm = squareform(a.distMat_ltv);
                sqfm(order,order) = sqfm;
                overall_arena(map,subj,:) = squareform(sqfm);
                pos_arena(map,subj,:,:)     = a.itemPositions(order,:);
                arena(subj,map) = corr(squareform(data.map{1})',squareform(sqfm)');
%                 arena(subj,4) = corr(squareform(data.map{1})',a.distMat_ltv');
            
            end
        end
    end
end

writematrix(arena, fullfile(basedir_beh, 'beh_measures','arena.csv'));

[h,p] = ttest(atanh(arena))


%% Overall performance as a function of the neural effect
counter = 0;
% type = {'rel','irrel'};
figure('Position', [100, 100, 1200, 200]);
session = 2;
roi = 3;
for c = 1:2
    for map = 1:2
        counter = counter+1;
        subplot(1,4,counter)
        scatter(arena(~isnan(arena(:,1)),map)', pe(c, ~isnan(arena(:,1)), session, roi),50, colors(c, :),'filled')
        [r,p] = corrcoef(arena(~isnan(arena(:,1)),map), pe(c, ~isnan(arena(:,1)), session, roi));
        
        lsline
        xlabel('Percent correct')
        ylabel('Effect size')
        title(sprintf('Arena %u, %s, %0.2f',map,con{c},p(1,2)), 'Interpreter','none')
        ylim([-0.8 1.2])
    end
end
set(gcf,'color','w')



