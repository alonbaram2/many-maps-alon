function  ss=FscatJit2(identifiers, data, varargin)

% Plots data as both jittered points and an error bar, in the
% order the identifiers are passed to it. It preserves whatever order
% was in the fbmm structure. Also includes the option to define the circle
% size or use a default. Also allows someone to add bars to the plot.
%
% FscatJit(identifiers, data)
% Uses the default circle size, no bars
%
% FscatJit(identifiers, data, circleSize)
% Defines the circle size, no bars
%
% FscatJit(identifiers, data, circleSize, 'on')
% Defines the circle size, adds bars

% Adam Claridge-Chang 20120522 Takes numeric variables on the X-axis as
% well as nominal.
% Labels xaxis with the uidents instead of numbers Adam CC 20121113

% Sameer Aryal Jan 22, 2013. Can now compute mean difference and plot it
% next to the bar scatjits, along with a floating axis.

%% For testing_______
% Mamma mia
% % % % clear all
% % % % clc
% % % % close all
% % % % load('FscatJit_data.mat') %for two genotypes
% % load('FscatJit_data_five.mat') % for 5 genotypes
% % load ageingfbmm.mat % for debugging reverse axis
% % load NaLactaefbmm.mat
% load lightcurvefbmm.mat
% % for n-genotypes
% % clearvars -except fbmm
% % identifiers = fbmm.szGenotype;
% % data = fbmm.AngVbyEpoch(:,2);
%
% % % This option to test numeric identifiers
% % a=rand(length(identifiers),1);
% % b=a*4;
% % identifiers=ceil(b);
% % identifiers(1:45)=7;

% % Generate artificial data
% reps=20;
% id={'abrams', 'daddy'}';% , 'tank', 'blah', 'blech', 'cch', 'blargh', 't', 'y', 'fgh', 'sdgh', 'sryj', 'wetrh'}';
% % id={'abrams', 'daddy', 'tank', 'blah', 'blech', 'cch'}';%, 'blargh', 't', 'y', 'fgh', 'sdgh', 'sryj', 'wetrh'}';
% % id={'abrams', 'daddy', 'tank', 'blah', 'blech', 'cch', 'blargh', 't', 'y', 'fgh', 'sdgh', 'sryj', 'wetrh'}';
% identifiers=repmat(id, reps, 1);
% data=rand(length(identifiers), 1);

% Set circleSize


figure; 
%% Deal with the varargin options
% fprintf('Total number of inputs = %d\n',nargin);
nVarargs = length(varargin);
% fprintf('Inputs in varargin(%d):\n',nVarargs);
% If no circleSize is given, define a default
if nVarargs == 0
    lims = [];
    isPaired = 'N';
    circleSize=170;
    barstate='off';
elseif nVarargs == 1 && isfloat(varargin{1});
    lims = varargin{1};
    isPaired = 'N';
    circleSize= 170;
    barstate='off';
elseif nVarargs == 1 && ~isfloat(varargin{1}); % this is the case for paired
    lims = [];
    isPaired = varargin{1};
    circleSize=50;
    barstate='off';
elseif nVarargs == 2
    lims = varargin{1};
    isPaired = varargin{2};
    circleSize=170;
    barstate='off';
elseif nVarargs == 3
    lims = varargin{1};
    isPaired = varargin{2};
    circleSize=varargin{3};
    barstate='off';
elseif nVarargs == 4
    lims = varargin{1};
    isPaired = varargin{2};
    circleSize=varargin{3};
    barstate=varargin{4};
end

switch barstate
    case 'on'
        middle_bar='off';
    case 'off'
        middle_bar='on';
    otherwise
end



%% Repack into nan-padded columns if they are in vector form
if size(identifiers)==size(data)
    [celld, uidents] = repackData (identifiers, data);
else
    uidents = identifiers;
    celld = data;
end
% Define nex
nex=length(celld);


if length(uidents) > 2
    p = panel();
    p.pack([50,50], 1);
end

%% Define plotting
if iscell(identifiers)
    X=1:length(uidents);
else
    X=uidents;
end


if length(uidents) > 2
    p(1,1).select();
end
%% Plot bars as option
if strcmp(barstate, 'on')
    for idx=1:nex
        curDat=celld{idx};
        av(idx)=nanmean(curDat);
    end
    [b] = bar(av);
    mydata=(1:length(av));
    bar_child=get(b,'Children');
    set(bar_child,'CData',mydata);
    % colormap(summer)
    % beta=0.5;
    % brighten(beta)
end

%% Plot scatjits
jitFactor=0.2;
circleSize=circleSize./max(X);% So circleSize scales with n-data columns

if strcmp(barstate, 'off') && strcmp(isPaired, 'N')
    colors = lines(100);
    for idx=1:nex
        curDat=celld{idx};
        hold on
        [s1] = scatJit(curDat, jitFactor, X(idx) ,circleSize,colors(idx,:));
              
    end
end


set(gca,'XTick',X);
hold on;

barwidth=0.5;
linewidth = 1;
% Use 1.96X SEM (margin of error on z-distibution) as the error bars
for idx=1:nex
       
    curDat=celld{idx};
    [av(idx), moes, bci] = bootmoes(curDat);
    
    er(:, idx) = moes; 
    CI(idx, :) = bci;  
    Value(idx, :) = av(idx);
    N(idx, :) = length(curDat);
%     stdcd=nanstd(curDat);
%     er(idx)=1.96*(stdcd/sqrt(length(curDat)-1));
end
hold on;
stats = table(uidents, Value, CI, N, 'VariableNames',{'Group','Value','CIs','N'});

[e1] = tripleErrorBars(av, er, X, barwidth, linewidth, middle_bar);


%% Set ticks, contigent on whether it is 2 or some other number of datasets
if length(celld)==2;
    Xplus=horzcat(X, 3);
    disp(Xplus);
    set(gca,'XTick',Xplus)
%     mdidents=vertcat(uidents,' ');
%     mdidents=vertcat(uidents,'{\Delta}');
    mdidents=vertcat(uidents,[uidents{2} ' - ' uidents{1}]);
    set(gca, 'xtickLabel', mdidents);
    set(gca, 'XLim', [0 length(mdidents)+1]);
    set(gca,'FontSize',18,'FontName','Helvetica')
    ylabel('1-r','FontSize',18,'FontName','Helvetica');
%     ylabel('value','FontSize',18,'FontName','Helvetica');
else
    set(gca,'XTick',X);
    %     set(gca, 'xtickLabel', uidents);
    set(gca, 'xtickLabel', []);
    set(gca, 'XLim', [0 length(uidents)+1]);
    ylabel('value','FontSize',18,'FontName','Helvetica');
end



%% Insert mean difference and CIs on a different axis if there is a pair
if length(celld)==2;
    if strcmpi(isPaired, 'Y')
        colors = lines(100);
        hold off
        if length(celld{1}) == length (celld{2})
            curDat = [celld{1} celld{2}];
        else
            disp ('Number of flies do not match. Aborting');
            return
        end
        [rows, ~] = size(curDat);
        paired = plot(curDat(1, :));
        set(paired, 'Color', [0    0    0]);
        hold on;
        for idx = 2:rows
            paired = plot(curDat(idx, :));
            set(paired, 'Color', [0    0    0]);
        end
        set(gca,'XTick',Xplus)
        set(gca, 'xtickLabel', mdidents);
        set(gca, 'XLim', [0 length(mdidents)+1], 'box', 'off');
        ylabel('1-r','FontSize',18,'FontName','Helvetica');
%         ylabel('value','FontSize',18,'FontName','Arial');
        set(gca,'FontSize',18,'FontName','Helvetica')        
        tripleErrorBars(av, er, [.5 2.5], barwidth, linewidth, middle_bar);
    end


    
    % Get the mean difference and CIs
 
    esm='md';
    if strcmp(isPaired, 'Y')
        ss=mes(celld{2},celld{1},esm,'nBoot',10000, 'isDep', 1);
    else
        ss=mes(celld{2},celld{1},esm,'nBoot',10000);
    end
    avr=repmat(ss.md,2, 1);
    moes=abs(avr-ss.mdCi);
    delta_name = {[uidents{2}, ' minus ', uidents{1}]}; 
    stats_delta = table(delta_name, ss.md, ss.mdCi', NaN, 'VariableNames', {'Group','Value','CIs','N'});
    stats = [stats; stats_delta];

    % Position of the reference axes (the axes that hold the scatjit)
    refAxes = gca;
    if isempty(lims) == 0
        set(refAxes, 'YLim', lims);
        [num, num3, ax1Pos, x, y, x2, y2, yNew, ~] = setThirdAxis(refAxes, av, ss.md);
        if ss.md < 0
            while (num3*y2-(num3-1)*y) <= ax1Pos(2)
                [num, num3, ax1Pos, x, y, x2, y2, yNew, ~] = setThirdAxis(refAxes, av, ss.md);
            end
            %         elseif ss.md>0
            %             while (num3*y2-(num3-1)*y) >= ax1Pos(2) + ax1Pos(4);
            %                 [num, num3, ax1Pos, x, y, x2, y2, yNew, ~] = setThirdAxis(refAxes, av, ss.md);
            %             end
        end
    else
        [num, num3, ax1Pos, x, y, x2, y2, yNew, ~] = setThirdAxis(refAxes, av, ss.md);
    end
    
    
    if ss.md < 0
        
        while (num3*y2-(num3-1)*y) <= ax1Pos(2)
            refLims = get(refAxes, 'YLim');
            refDiff = 0.1*(refLims(2) - refLims(1));
            newRefLims = [refLims(1)-refDiff refLims(2)+refDiff];
            set(refAxes, 'YLim', newRefLims);
            [num, num3, ax1Pos, x, y, x2, y2, yNew, ~] = setThirdAxis(refAxes, av, ss.md);
        end
        
        
        
        % Errorbar axis
        errorBarAxis = axes('Position',[ax1Pos(1) num3*y2-(num3-1)*y ax1Pos(3) yNew-(num3*y2-(num3-1)*y)]);
        p3= errorbar(3 ,ss.md, moes(1),moes(2));
        
        axis([0 4 num3*ss.md (num-1)*abs(ss.md)]);
        
        
        
        % Dummy axis that emuluates errorbar axis, but is placed next
        % to the mean difference
        dummyAxis = axes('Position',[ax1Pos(1) + ax1Pos(3)-((x-ax1Pos(1))/2) num3*y2-(num3-1)*y .001 yNew-(num3*y2-(num3-1)*y)]);
        axis([0 4 num3*ss.md (num-1)*abs(ss.md)]);
        
        marker ='v';
    elseif ss.md > 0
        
        while (num3*y2-(num3-1)*y) >= ax1Pos(2) + ax1Pos(4);
            refLims = get(refAxes, 'YLim');
            refDiff = 0.1*(refLims(2) - refLims(1));
            newRefLims = [refLims(1)-refDiff refLims(2)+refDiff];
            set(refAxes, 'YLim', newRefLims);
            [num, num3, ax1Pos, x, y, x2, y2, ~, ~] = setThirdAxis(refAxes, av, ss.md);
        end
        
        
        % Errorbar axis
        errorBarAxis = axes('Position',[ax1Pos(1) num*y-(num-1)*y2 ax1Pos(3) ((num3*y2-(num3-1)*y) - (num*y-(num-1)*y2))]);
        p3= errorbar(3 ,ss.md, moes(1),moes(2));
        
        axis([0 4 (num-1)*(-1)*ss.md num3*(ss.md)]);
        
        % Dummy axis that emuluates the y values of hte errorbar axis, but is placed next
        % to the mean difference
        dummyAxis = axes('Position',[ax1Pos(1) + ax1Pos(3)-((x-ax1Pos(1))/2) num*y-(num-1)*y2 .001 ((num3*y2-(num3-1)*y) - (num*y-(num-1)*y2))]);
        axis([0 4 (num-1)*(-1)*ss.md num3*(ss.md)]);
        marker ='^';
        
    else
        
        errorBarAxis = axes('Position',[ax1Pos(1) y-0.3 ax1Pos(3) 0.6]);
        
        p3= errorbar(3 ,ss.md, moes(1),moes(2));
        axis([0 4 -.1 .1]);
        
        dummyAxis = axes('Position',[ax1Pos(1) + ax1Pos(3)-((x-ax1Pos(1))/2) y-0.3 .001 0.6]);
        axis([0 4 -0.1 0.1]);
        
        marker ='o';
    end
    
    if strcmp(isPaired, 'Y')
        colors = lines(100);
        [x,~] = dsxy2figxy(refAxes, .5, av(1));
        [x2, ~] = dsxy2figxy(refAxes, 2.5, av(2));
        line1= annotation('line', [x ax1Pos(1) + ax1Pos(3)-((x-ax1Pos(1)))], [y y]);
        line2=annotation('line', [x2 ax1Pos(1) + ax1Pos(3)-((x-ax1Pos(1)))], [y2 y2]);
        
    else
        % Plot the lines that join the means to the difference axis.
        line1= annotation('line', [x ax1Pos(1) + ax1Pos(3)-((x-ax1Pos(1))/2)], [y y]);
        line2=annotation('line', [x2 ax1Pos(1) + ax1Pos(3)-((x-ax1Pos(1))/2)], [y2 y2]);
    end
    
    set(line1, 'LineStyle', ':');
    set(line2, 'LineStyle', ':');
    marker = 'v';
    set(p3,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',5,...
        'Marker',marker,...
        'LineStyle','none');
    
    set(errorBarAxis, 'Visible', 'Off');
    set(dummyAxis, 'YAxisLocation', 'right');
    set(dummyAxis,'FontSize',18)
    %     set(dummyAxis, 'FontSize', 7);
    
else
    % To ensure consistency with the md plot so that each figure has 3
    % child axes
    refAxes = gca;
    errorBarAxis = axes('Position', get(refAxes, 'Position'));
    dummyAxis = axes('Position', get(refAxes, 'Position'));
    set(errorBarAxis,'Color','k', 'Visible', 'Off');
    set(dummyAxis,'Color','k', 'Visible', 'Off');
    marker = 'o';
    ylabel('value','FontSize',18,'FontName','Arial');
    if isempty(lims) == 0
        set(refAxes, 'YLim', lims);
    end
    
    % Get the mean difference and CIs
    esm='md';
    clear moes;
    
    avr = zeros(2, length(celld));
    moes = zeros(2, length(celld));
    ci = zeros(2, length(celld));
    N = NaN(length(celld)-1,1);
    
    for idx = 2:length(celld)
        ss=mes(celld{idx},celld{1},esm,'nBoot',10000);
        avr(:,idx)=repmat(ss.md,2, 1);
        moes(:,idx)=abs(avr(:,idx)-ss.mdCi);
        ci(:,idx) = ss.mdCi;
        delta_name(idx-1,:) = {[uidents{idx}, ' minus ', uidents{1}]};
    end
    
    stats_delta = table(delta_name, avr(1,2:end)', ci(:,2:end)', N, 'VariableNames', {'Group','Value','CIs','N'});
    stats = [stats; stats_delta];
    
    [ciMin, ~] = min(ci(1,:));
    [ciMax, ~] = max(ci(2,:));
    
    axLo = ciMin - .1*abs(ciMax-ciMin);
    axHi = ciMax + .1*abs(ciMax-ciMin);
    
    p(2,1).select();
    ylabel('delta value','FontSize',18,'FontName','Arial');
    p4= errorbar(1:length(avr(1,:)), avr(1,:), moes(1,:), moes(2,:));
    
    axis([0 length(avr(1,:))+1 axLo axHi]);
    
    marker ='o';
    
    set(p4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'Color','k',...
        'MarkerSize',5,...
        'Marker',marker,...
        'LineStyle','none');
    set(gca, 'Xtick', X, 'xtickLabel', uidents);
    
    line1 = refline(0,0);
    set(line1, 'lineStyle', ':')
    
    %% Multiple pairwise comparisons
    clearvars avr moes ci delta_name N;
    if mod(length(celld),2)==0
        figure;
        pwmd = panel();
        pwmd.pack([50,50], 1);
 
        marker = 'o';
        
        % Pairwise mean differences
        idx = 1;
        jdx = 1;
        
        while jdx < length(celld)
            ss=mes(celld{jdx+1},celld{jdx},esm,'nBoot',10000);
            avr(:,idx)=repmat(ss.md,2, 1);
            moes(:,idx)=abs(avr(:,idx)-ss.mdCi);
            ci(:,idx) = ss.mdCi;
            delta_name(idx,:) = {[uidents{jdx+1}, ' minus ', uidents{jdx}]};
            
            jdx = jdx+2;
            idx = idx + 1;
        end
        stats_delta = table(delta_name, avr(1,:)', ci', NaN(length(delta_name),1), 'VariableNames', {'Group','Value','CIs','N'});
        stats = [stats; stats_delta];
        
        count = 0;
        idx = 1;
        while idx <= length(X)
            for jdx = 1:2
                newX(idx) = X(idx)+(count);
                idx = idx + 1;
            end
            count = count + 1;
        end
        
        pwmd(1,1).select();
        refAxes = gca;
        marker = 'o';
        ylabel('value','FontSize',18,'FontName','Arial');
        if isempty(lims) == 0
            set(refAxes, 'YLim', lims);
        end
        set(refAxes, 'XLim', [0 max(newX)+1]);
        
        if strcmp(barstate, 'on')
            for idx=1:nex
                curDat=celld{idx};
                av(idx)=nanmean(curDat);
            end
            [b] = bar(newX, av);
            mydata=(1:length(av));
            bar_child=get(b,'Children');
            for idx = 1:nex
                if mod(idx, 2)==1
                    mydata(idx) = 0;
                else
                    mydata(idx) = 1;
                end
            end
            set(bar_child,'CData',mydata);
            cmap = [0     0     0; 0    0    0];
            set(k, 'EdgeColor', [0    0    0], 'LineWidth', 1.25);
            set(gca, 'box', 'off');
            colormap(cmap);

        end
        
        if strcmp(barstate, 'off')
            for idx = 1:nex
                curDat=celld{idx};
                hold on
                scatJit(curDat, jitFactor, newX(idx), circleSize,colors(idx,:));
            end
        end
        tripleErrorBars(av, er, newX,barwidth, linewidth, middle_bar);
        
        
        pwmd(2,1).select();
        ylabel('delta value','FontSize',18,'FontName','Arial');

        idx = 1;
        count = 1;
        esX = zeros(1, length(newX)/2);
        while count <= length(newX)/2
            esX(count) = (newX(idx) + newX(idx+1))/2;
            idx = idx + 2;
            count = count + 1;
        end
        
        p5= errorbar(esX, avr(1,:), moes(1,:), moes(2,:));
        
        marker ='o';
        set(p5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'Color','k',...
            'MarkerSize',5,...
            'Marker',marker,...
            'LineStyle','none'...
           );
        set(refAxes, 'Xtick', []);
        %         set(refAxes,'XTick',newX, 'xtickLabel', uidents);
        set(gca, 'Xtick', newX, 'xtickLabel', uidents);
        set(gca, 'XLim', get(refAxes, 'Xlim'));
        
        line1 = refline(0,0);
        set(line1, 'lineStyle', ':')
    end
    
end
ss = stats;
end
