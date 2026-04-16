clear
close all

addpath(genpath(fullfile('/vols/Scratch/mgarvert/ManyMaps','imagingData','scripts','alon','utilities')));

arena = nan(25,2);
overall_arena = nan(2,25,136);
position_arena = nan(2,25,17,2);
for subj = 1:25
    for map = 1:2
        if exist(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/arena/similarityJudgementData_Subj_',num2str(subj),'_map_',num2str(map),'_trial1.mat'],'file')
            a = load((['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/arena/similarityJudgementData_Subj_',num2str(subj),'_map_',num2str(map),'_trial1.mat']));
            load(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/alldata/data_',num2str(subj),'_2.mat']);
            

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
            
            end
        end
    end
end


[h,p] = ttest(atanh(arena))

figure;
subplot(1,3,1)
Y=mdscale(data.map{1},2);
eudist = pdist(Y,'euclidean');
imagesc(squareform(eudist))
axis equal off
for map = 1:2
    subplot(1,3,map+1)
    imagesc(squareform(squeeze(nanmean(overall_arena(map,:,:),2))))
    axis equal off
end
colormap(lbmap(25,'RedBlue'))

%% Subpanel D
% Define colors
colors = lbmap(2, 'RedBlue');

% Data preparation
[h, p] = ttest(atanh(arena));
mean_arena = nanmean(arena);
std_arena = nanstd(arena);
n1 = sum(~isnan(arena(:,1))); % number of non-NaN values
n2 = sum(~isnan(arena(:,2))); % number of non-NaN values

% Plot
figure;
hold on;

% Bar plots
bar(1, mean_arena(1), 'FaceColor', colors(2, :));
bar(2, mean_arena(2), 'FaceColor', colors(1, :));

% Error bars
errorbar(1, mean_arena(1), std_arena(1) / sqrt(n1), 'k');
errorbar(2, mean_arena(2), std_arena(2) / sqrt(n2), 'k');

% Individual data points
scatter(ones(size(arena(:, 1)))- 0.2 + 0.4*rand(size(arena(:, 1),1), 1), arena(:, 1), 35, colors(2, :), 'filled','MarkerEdgeColor', 'k');
scatter(2 * ones(size(arena(:, 2)))- 0.2 + 0.4*rand(size(arena(:, 2),1), 1), arena(:, 2), 35, colors(1, :), 'filled','MarkerEdgeColor', 'k');

% Formatting
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Graph 1', 'Graph 2'});
ylabel('Correlation between arena distances and true distances');
box off;
set(gcf,'color','w');

%% Subpanel E
%% Analyse data from cross-map suppression

cr = nan(25,2,4);
RT_effect = nan(25,2,5);
cr_effect = nan(25,2,5);
cr_rel = nan(25,2);
cr_irrel = nan(25,2);
RT_all_effect = nan(25,5);
cr_all_effect = nan(25,5);
pool_across_sessions = 0;


h1 = figure;
for subj = 1:25
    sessiontype = []; choice = []; cr_cum = []; RT_cum = []; dist_rel = []; dist_irrel = []; dist_rel_ch = []; dist_irrel_ch = []; dist_rel_unch = []; dist_irrel_unch = [];
    trial = [];
        
    for session = 1:2
        if ~pool_across_sessions
            choice = []; cr_cum = []; RT_cum = []; dist_rel = []; dist_irrel = []; dist_rel_ch = []; dist_irrel_ch = []; dist_rel_unch = []; dist_irrel_unch = [];
        end
        try
                load(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/alldata/data_',num2str(subj),'_',num2str(session),'.mat']);
                disp(subj)
                %                 load(strtrim(ls([basedir_beh,'scan_1.1/datafiles/Subj_',num2str(subj),'/*_session_',num2str(session),'/data_',num2str(subj),'_',num2str(session),'.mat'])))
                %             load(strtrim(['C:/Users/mgarvert/Downloads/data/savedata/data_',num2str(subj),'_',num2str(session),'.mat']))
                for run = 1:4
                    % remove decoy choices
                    ix_ch = data.scan{run}.objDiff.distRel(:,1)~=data.scan{run}.objDiff.distRel(:,2);
                ix_ch = ix_ch(1:length(data.scan{run}.objDiff.distRelCh));
                
                cr(subj,session,run) =  sum(data.scan{run}.objDiff.cr(ix_ch)==1)/sum(data.scan{run}.objDiff.cr(ix_ch)~=0);
                
                % Regression analysis to determine the influence of
                % chosen and unchosen distances on behaviour
                choice = [choice; data.scan{run}.objDiff.choice(ix_ch)'];
                cr_cum = [cr_cum; (data.scan{run}.objDiff.cr(ix_ch)'+1)/2];
                RT_cum = [RT_cum; log(data.scan{run}.objDiff.RT(ix_ch)')];
                dist_rel_ch = [dist_rel_ch; data.scan{run}.objDiff.distRelCh(ix_ch)'];
                dist_irrel_ch = [dist_irrel_ch; data.scan{run}.objDiff.distRelUnch(ix_ch)'];
                dist_rel_unch = [dist_rel_unch; data.scan{run}.objDiff.distIrrelCh(ix_ch)'];
                dist_irrel_unch = [dist_irrel_unch; data.scan{run}.objDiff.distIrrelUnch(ix_ch)'];
                dist_rel = [dist_rel; data.scan{run}.objDiff.distRel(ix_ch,:)];
                dist_irrel = [dist_irrel; data.scan{run}.objDiff.distIrrel(ix_ch,:)];
                sessiontype = [sessiontype; repmat(session,sum(ix_ch),1)];
                trial = [trial; repmat(run,sum(ix_ch),1)];
            end
            if session == 2
                day(subj) = data.day;
            end
            
            
            %                 cr_cum = categorical(cr_cum);
            
            % How many choices are correct on the relevant / irrelevant
            % map?
            choice(choice == 1) = 2; choice(choice == -1) = 1;
            cr_rel(subj,session) = (sum(dist_rel(choice == 1,1) < dist_rel(choice == 1,2)) + sum(dist_rel(choice == 2,1) > dist_rel(choice == 2,2))) /sum(dist_rel(:,1) ~= dist_rel(:,2));
            cr_irrel(subj,session)  = (sum(dist_irrel(choice == 1,1) < dist_irrel(choice == 1,2)) + sum(dist_irrel(choice == 2,1) > dist_irrel(choice == 2,2))) /sum(dist_irrel(:,1) ~= dist_irrel(:,2));
            
            
            only_include_real_data = RT_cum>0;
            
            if ~pool_across_sessions
                RT_effect(subj,session,:) = regress(RT_cum(only_include_real_data),[ones(sum(only_include_real_data),1) dist_rel_ch(only_include_real_data),dist_irrel_ch(only_include_real_data) dist_rel_unch(only_include_real_data),dist_irrel_unch(only_include_real_data)]);
                cr_effect(subj,session,:) = glmfit([dist_rel_ch(only_include_real_data),dist_irrel_ch(only_include_real_data) dist_rel_unch(only_include_real_data),dist_irrel_unch(only_include_real_data)],cr_cum(only_include_real_data), 'binomial', 'link', 'logit');
            end
            
            RT_all_effect(subj,:) = regress(RT_cum(only_include_real_data),[ones(sum(only_include_real_data),1) dist_rel_ch(only_include_real_data),dist_irrel_ch(only_include_real_data) dist_rel_unch(only_include_real_data),dist_irrel_unch(only_include_real_data)]);
            cr_all_effect(subj,:) = glmfit([dist_rel_ch(only_include_real_data),dist_irrel_ch(only_include_real_data) dist_rel_unch(only_include_real_data),dist_irrel_unch(only_include_real_data)],cr_cum(only_include_real_data), 'binomial', 'link', 'logit');
        catch
        end
    end
end

cr_blocks = nanmean(cr, 3);
mean_cr = squeeze(nanmean(nanmean(cr, 3)));
std_cr = squeeze(nanstd(nanmean(cr, 3)));
n = sum(~isnan(squeeze(nanmean(cr, 3)))); % number of non-NaN values

figure, bar(mean_cr, 'FaceColor', [248/256, 206/256, 161/256]); hold on
errorbar(1:2, mean_cr, std_cr ./ sqrt(n), 'k','linestyle','none');
hold on;
for bl = 1:2
    v = bl - 0.2 + 0.4 * rand(25, 1);
    scatter(v, nanmean(cr(:, bl), 3), 35, [248/256, 206/256, 161/256], 'filled','MarkerEdgeColor', 'k');
end
[h, p] = ttest(cr_blocks(:, 1), cr_blocks(:, 2));
set(gca, 'XTickLabel', {'session 1', 'session 2'},'YLim',[0.4, 1]);
ylabel('Percent correct choices');
prepImg;


%% Subpanel F

figure;
scatter(nanmean(arena, 2), nanmean(cr_blocks,2)', 75, [248/256, 206/256, 161/256], 'filled','MarkerEdgeColor', 'k');
[r, p] = corrcoef(nanmean(arena, 2), nanmean(cr_blocks,2), 'rows', 'complete');
lsline;
title(sprintf('r = %0.2f, p = %0.3f', r(1,2), p(1, 2)));
xlabel('Correlation between arena distances and true distances');
ylabel('Percent correct choices');
prepImg

%% Subpanels G and H
colors = lbmap(7,'RedBlue');
colix = [1 2 6 7];

n = sum(~isnan(squeeze(RT_all_effect(:,2:end)))); % number of non-NaN values
h2 = figure;
subplot(1,2,1)
b = bar(squeeze(nanmean(RT_all_effect(:,2:end),1)),'FaceColor','flat','EdgeColor','none');
hold on, errorbar(1:4,squeeze(nanmean(RT_all_effect(:,2:end),1)),squeeze(nanstd(RT_all_effect(:,2:end),1))/sqrt(n),'k');
b.CData = colors([1 2 6 7],:);
hold on
for bl = 1:4
    v = bl-0.2+0.4*rand(25,1);
    scatter(v,RT_all_effect(:,bl+1),35,colors(colix(bl),:),'filled','MarkerEdgeColor', 'k');
end
hold on, errorbar(1:size(RT_all_effect,2)-1,squeeze(nanmean(RT_all_effect(:,2:end),1)),squeeze(nanstd(RT_all_effect(:,2:end))/sqrt(sum(~isnan(RT_all_effect(:,1))))),'k','linestyle','none');
title(['RT effect both sessions'])
xticklabels({'chosen rel','chosen irrel','unchosen rel', 'unchosen irrel'}), xtickangle(45)
xtickangle(45)
prepImg
disp(['cr both sessions'])
[h,p] = ttest(squeeze(RT_all_effect(:,2:end)))
ylim([-0.05 0.5])

colors = lbmap(7,'RedBlue');
subplot(1,2,2)
b = bar(squeeze(nanmean(cr_all_effect(:,2:end),1)),'FaceColor','flat','EdgeColor','none');
hold on, errorbar(1:4,squeeze(nanmean(cr_all_effect(:,2:end),1)),squeeze(nanstd(cr_all_effect(:,2:end),1))/sqrt(n),'k');
b.CData = colors([1 2 6 7],:);
hold on
for bl = 1:4
    v = bl-0.2+0.4*rand(25,1);
    scatter(v,cr_all_effect(:,bl+1),35,colors(colix(bl),:),'filled','MarkerEdgeColor', 'k');
end
hold on, errorbar(1:size(cr_all_effect,2)-1,squeeze(nanmean(cr_all_effect(:,2:end),1)),squeeze(nanstd(cr_all_effect(:,2:end))/sqrt(sum(~isnan(cr_all_effect(:,1))))),'k','linestyle','none');
title(['cr effect both sessions'])
prepImg
xtickangle(45)
xticklabels({'chosen rel','chosen irrel','unchosen rel', 'unchosen irrel'}), xtickangle(45)
disp(['cr both sessions'])
[h,p] = ttest(squeeze(cr_all_effect(:,2:end)))
ylim([-80 80])
set(gcf,'color','w')

%% Subpanel I-L

label = {'chosen rel','chosen irrel','unchosen rel', 'unchosen irrel'};
figure('Position', [100, 100, 1200, 400]);
for s = 1:4
    subplot(1,4,s)
    scatter(nanmean(arena, 2), RT_all_effect(:,s+1),50,colors(colix(s),:),'filled','MarkerEdgeColor', 'k');
    [r, p] = corrcoef(nanmean(arena, 2), RT_all_effect(:,s+1),'rows','complete');
    lsline;
    title(sprintf('%s, r = %0.2f, p = %0.2f', label{s}, r(1, 2), p(1, 2)), 'Interpreter', 'none');
    ylabel('RT effect')
    xlabel('Correlation between arena distances and true distances');
    prepImg
end


