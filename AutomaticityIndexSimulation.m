% close all; 
clear all; 
clc

plotIndVersions = true;
saveFigure = false;
shape = 'Linear'; %Sigmoid, Linear, InvU

simulationIterations = 100;

CaseScaledData = {};
impairLevels = 0:1/(simulationIterations-1):1; %fix impairment, inclusive on both ends, so inbetween is simulationIteration - 1
  % set up path, x range,
rng(2022);
if strcmp(shape,'Sigmoid')
    x = -10:0.2:10;
elseif strcmp(shape, 'Linear')
    x = 0:0.2:10;
end
if simulationIterations <= 5
    plotFigure = true;
else
    plotFigure = false;
end
scriptDir = fileparts(matlab.desktop.editor.getActiveFilename); 
%redefine version to be a list, so the current version 2 will become 1:2.
% -5: performance top only, more impaired, lower top; this
% makes impaired always better bc even though perf change is the same, the
% normalization factor changed and hence the scaling is always bigger for
% impaired
% -4: performance range only; more impaired, larger range = lower btm; 
% if both -5 and -4 present, need to use coupled case
% -3: PFC top only, more impaired, lower top
% -2: PFC btm only, more impaired, higher btm
% if both -3 and -2 present, need to use coupled case
% -1: vary inflection pt of the performance
% NOT USED 0: performance range + top; not useful anymore replaced by [-5,-4]
% NOT USED v1: top and bottom of PFC only, not used anymore, replaced by [-3,-2]
% v2: pfc inflection points; or pfc intercept for linear.
% v3: slope of impaired PFC < control slope;
% v4:slope of impaired performance > control performance
% v5: linear only, slope of impaired PFC > control slope.
versions = {[0:1],[0:2],[0:3],[0:4],[0,1,2,4],[0],[1],[2],[3],[4],[-1],[-1:4]}; %first 5 covered in previous notebook.
%Previous t v5 = [1,2,4]; current v4 = [1,2,3,4]

%% generate all possible combinations of version choices.
caseCoding = [5,3,2,4,-1]%for linear: [5,3,2,4,-1]; for sigmoid[-5:-1,2:4];
n = length(caseCoding); k =1;
versions = dec2bin(sum(nchoosek(2.^(0:n-1),k),2)) - '0'; 
for k = 2:n
    if strcmp()
    versions = [versions; dec2bin(sum(nchoosek(2.^(0:n-1),k),2)) - '0']; 
    %append binary array of options to do n chooses k, the binary can be
    %used for logical indexing to get correct case combinations for each n
    %chooses k option.
end
versions = logical(versions);
versionsFull = versions;


%% main simulation loop
close all; clc
% versions = {[-5],[-4],[-3],[-2],[-1],[2],[3],[4],[-5,-4],[-3,-2],[-5,-4,-1,4],[-3,-2,2,3],[-5,-4,-1,4,-3,-2,2,3]}; %all perf, all pfc, all
versions = {[3,4]}%[3],[5],[2],[4],[-1]}
FullScaledData = {};
for verIdx = 1:length(versions)%[1 2 3 5 6 11]%6:numel(versions)
    % clearvars -except versions verIdx
    verIdx
    version = versions{verIdx};
%     version = caseCoding(versions(verIdx,:))
%     close all;
    % Changing param/settings
    CaseScaledData = {};
    scaleFactors = [1,2,4,10];%[1, 2, 4, 10];
    for scaleFactor = scaleFactors
        perfInflectionRatioToPFCCtr = 0.85; %TODO: To test for sensitivity
        if strcmp(shape, 'Sigmoid')
            bPerfCtlVal = 2.5; % between 90 - 95%
        elseif strcmp(shape, 'Linear')
            bPerfCtlVal = 12;
        end
        btmPerfCtlVal = 0.75;

        saveDir = [scriptDir '/Data/SimulationResults/Version' num2str(verIdx) '/'];
        if strcmp(filesep,'\') %mac enviroment, replace \ with /
            saveDir = replace(saveDir,"/",filesep);
        end

        fullData = nan(simulationIterations, length(x),6); %y1,y2,y3,y4,y1scaled,y2scaled

        for i = 1:simulationIterations
            a1 = 1; %slope,  smaller a smaller slope
            b1 = 0; %shifts, larger b, curve shifted more towards left (negative end of x-axis)
            %b = x value where the response is halfway between max and min
            top1 = 1;
            bottom1 = 0;
            % plateauX1 = (log((top1-bottom1)/(0.9-bottom1) -1) -b1)/(-a1);

            %Control performance curves, same top, much smaller range, negative and sharper
            %slope, plateau (curve starts dropping before healthy PFC curve plauteau
            %(or before it reaches 70-80% of the max), but remains steady after pass
            %50% of PFC curve.
            topPerfCtl = 1; %arbitrary chosen set value
            btmPerfCtl = btmPerfCtlVal;%AutomaticityIndexHelper.generateRandValInRange(0,0.9); %0.75; %TODO: To test for sensitivity
            aPerfCtl = -1; %AutomaticityIndexHelper.generateRandValInRange(-5,-1);
            bPerfCtl = bPerfCtlVal;
        %     bPerfCtl = AutomaticityIndexHelper.inverseSigmoidValue(perfInflectionRatioToPFCCtr, a1,b1);
            
            %impairmentLevel: parameter to control the range, smaller value means more
            %impaired. in range (0,1), the top and bottom of impaired PFC and control
            %will be impairLevel * healthyRange
            impairLevel = impairLevels(i); %fixed level
        %     impairLevel = AutomaticityIndexHelper.generateRandValInRange(0,1);%random level

            %pfc for impaired, same slope, differ in range (centered around
            %0.5, shrink top and bottom, smaller range, scaled by impairLevel; differ
            %in inflection point: scaled by impairLevel (larger the impairment, smaller
            %then inflection point) b = (b1-minx) * (1-impairLevel) + minx =
            %b1*impairLevel + minx*impairLevel = minx*impairLevel
            %same slope.
            if any(ismember(version, -2)) && any(ismember(version, -3)) %coupled top and btm
                %this coupling is required bc top and btm needs to be bounded so when top changes, range needs to adjust
    %             i.e., if impair = 1, top = 0, range is unchanging (range = 1,
    %             this is an ill case)
                topPFCImpaired = top1 - impairLevel/2;
                btmPFCImpaired = bottom1+impairLevel/2;
            elseif any(ismember(version, -2)) %btm only, more impaired higher btm
                topPFCImpaired = top1;
                btmPFCImpaired = impairLevel; %no impair, btm = 0, full impair, btm = 1
            elseif any(ismember(version, -3)) %top only, more impaired, lower top
                topPFCImpaired = 1 - impairLevel;%no impair, top = 0, full impair, top = 0
                btmPFCImpaired = bottom1; 
            else
                topPFCImpaired = top1;
                btmPFCImpaired = bottom1;
            end
            bPFCImpaired = b1;
            if any(ismember(version,2)) %Includes 2, inflection points, more impaired, starts changing earlier
                if strcmp(shape, 'Sigmoid')
                    bPFCImpaired = b1 - (b1 - min(x)) * impairLevel; 
                elseif strcmp(shape, 'Linear')
                    bPFCImpaired = b1 +  impairLevel; 
                end
            end
            if any(ismember(version, 3)) %diff slope, more impaired, less steep slope
                aPFCImpaired = a1*(1-impairLevel);
            elseif any(ismember(version, 5)) && strcmp(shape, 'Linear') %linear, diff slope, more impaired, more steep slope (over-recruitment)
                aPFCImpaired = a1*(1+impairLevel);
            else
                aPFCImpaired = a1; 
            end

            %performance of impaired, same slope as the control performance, inflection
            %point is scaled by impairLevel (larger impairlevel, smaller ratio/starts changing earlier, 
            %increase PFC but won't be able to maintain the performance as long
            %top is scaled by impairLevel (smaller than control performance top, larger impair, smaller top), 
            %range is scaled by impairLevel (larger than control performance range, higher impair, larger range)
            %slope: more negative than control, more impaired, more negative
            if any(ismember(version, -1))
                if strcmp(shape,'Sigmoid')
            %         perfInflectionRatioToPFCImpaired = perfInflectionRatioToPFCCtr * (1 - impairLevel);
                    bPerfImpaired = bPerfCtl - (bPerfCtl - min(x)) * impairLevel; %unbounded [min, bCtrl]
                elseif strcmp(shape,'Linear')
                    bPerfImpaired = bPerfCtl - impairLevel; %unbounded [min, bCtrl]
                end
            else
        %         perfInflectionRatioToPFCImpaired = 0.85;
                bPerfImpaired = bPerfCtl;
            end
        %     bPerfImpaired = AutomaticityIndexHelper.inverseSigmoidValue(topPFCImpaired*perfInflectionRatioToPFCImpaired, aPFCImpaired,bPFCImpaired); %find inflection point

            if any(ismember(version, 4)) %steeper slope for more impaired ppl
                aPerfImpaired = aPerfCtl * (1+impairLevel); 
            else 
                aPerfImpaired = aPerfCtl; 
            end
            if any(ismember(version, -5)) && any(ismember(version, -4)) %coupled
        %             topPerfImpaired=(1-impairLevel/2) *topPerfCtl; %top is lower than control perf, higher impair, lower top
%                 topPerfImpaired=(1-impairLevel/10) *topPerfCtl; %top is
%                 lower than control perf, higher impair, lower top, this
%                 doesn't cover the full range
%                 btmPerfImpaired = topPerfImpaired - (topPerfCtl-btmPerfCtl)*(1+impairLevel); %higher impair, lower bottom
                topPerfImpaired = topPerfCtl - (topPerfCtl - btmPerfCtl)*impairLevel; 
                %i.e., when impairLevel = 1, top = btm of the control; this
                %is arbitrary, could also make top = 0.5 when impairLevel =
                %1 and adjust the range accordingly. To make things
                %bounded, have to choose where you want the top the be then
                %decide bottom (btm is a dependent var)
%                 rangeI = (topPerfCtl - btmPerfCtl) + 0.5*impairLevel;
%                 btmPerfImpaired = topPerfImpaired - rangeI;
                btmPerfImpaired = btmPerfCtl - btmPerfCtl*impairLevel; % this is equivalent as the 2 lines above
            elseif any(ismember(version, -5)) %top only, more impaired, lower top
                topPerfImpaired = topPerfCtl - btmPerfCtl * impairLevel; %at max impair, btm = 0
                btmPerfImpaired = topPerfImpaired - (topPerfCtl - btmPerfCtl); %range is always the same as control
            elseif any(ismember(version, -4)) %range only, more impaired, larger range (i.e., lower btm)
                topPerfImpaired = topPerfCtl;
                btmPerfImpaired = btmPerfCtl * (1-impairLevel);
            else
                topPerfImpaired = topPerfCtl;
                btmPerfImpaired = btmPerfCtl;
            end
            
            if strcmp(shape, 'Sigmoid')
                y1 = AutomaticityIndexHelper.sigmoidValue(x, a1, b1, top1, bottom1);
                yPerfCtl = AutomaticityIndexHelper.sigmoidValue(x, aPerfCtl, bPerfCtl, topPerfCtl, btmPerfCtl);
                yPFCImpaired = AutomaticityIndexHelper.sigmoidValue(x, aPFCImpaired, bPFCImpaired, topPFCImpaired, btmPFCImpaired);
                yPerfImpaired = AutomaticityIndexHelper.sigmoidValue(x, aPerfImpaired, bPerfImpaired, topPerfImpaired, btmPerfImpaired);
            elseif strcmp(shape, 'Linear')
                y1 = AutomaticityIndexHelper.linearValue(x, a1, b1);
                yPerfCtl = AutomaticityIndexHelper.linearValue(x, aPerfCtl, bPerfCtl);
                %pfc impaired, slope impaired < slope control (version 3), or slope > slope
                %control (version 5), or intercept > intercept control (version 2)
                %where bImpaired = b + p. 
                yPFCImpaired = AutomaticityIndexHelper.linearValue(x, aPFCImpaired, bPFCImpaired);
                %perf impaired. slope impaired (more negative) < slope
                %control (version 4) or intercept < intercept control
                %(version -1), where b=bCtr-p
                yPerfImpaired = AutomaticityIndexHelper.linearValue(x, aPerfImpaired, bPerfImpaired);
            end
            
            fullData = plotAndStoreSimulationResult(plotFigure, x, aPerfCtl,bPerfCtl,topPerfCtl,btmPerfCtl, yPerfCtl,...
                aPerfImpaired, bPerfImpaired, topPerfImpaired, btmPerfImpaired, yPerfImpaired, y1, yPFCImpaired, ...
                impairLevel, fullData, i, scaleFactor, shape);
            %TODO: this could be improved, no need to simulate for each scale
            %factor, just scale them.
        end
        CaseScaledData{end+1} = fullData;
    end

    if plotIndVersions%only plot individual case plots if version is a reasonable number
    %% plot simulated range
    %plot simualted range
    xmin = min(x);
    colorOrder=colororder;
    %plot simulated results
    avgYs = mean(fullData(:,:,:),1);
    std_dev = std(fullData(:,:,:),1);
    topBorder = max(fullData(:,:,:));%avgYs + std_dev;
    bottomBorder = min(fullData(:,:,:)); %avgYs - std_dev;
    x2 = [x, fliplr(x)];
    f1 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
    for i = 1:4
        plot(x,topBorder(:,:,i),'LineWidth',1,'Color',[colorOrder(i,:),0.3],'HandleVisibility','off');
        plot(x,bottomBorder(:,:,i),'LineWidth',1,'Color',[colorOrder(i,:),0.3],'HandleVisibility','off');
        inBetween = [topBorder(:,:,i), fliplr(bottomBorder(:,:,i))];
        if i == 1
            fill(x2, inBetween, colorOrder(i,:),'FaceAlpha',0.3);
        else
            fill(x2, inBetween, colorOrder(i,:),'FaceAlpha',0.3,'HandleVisibility','off');
        end
        plot(x, avgYs(:,:,i), 'Color',colorOrder(i,:), 'LineWidth', 5);
    end
    legend('Max/Min','Control','Impaired','ControlPerformance','ImpairedPerformance')
    xlim([min(x) max(x)])
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    xticks([])
    xlabel('Task Difficulty')
    yticks([])
    ylabel('PFC Activation or Performance')

    %plot the T stats before and after scaling
    % %hard to see bc of overlapping
    % figure(); hold on;
    % for i = [1:2 5:6]
    %     plot(x,topBorder(:,:,i),'LineWidth',3,'Color',[colorOrder(i+1,:),0.3]);
    %     plot(x,bottomBorder(:,:,i),'LineWidth',3,'Color',[colorOrder(i+1,:),0.3]);
    %     inBetween = [topBorder(:,:,i), fliplr(bottomBorder(:,:,i))];
    %     fill(x2, inBetween, colorOrder(i+1,:),'FaceAlpha',0.3);
    %     plot(x, avgYs(:,:,i), 'Color',colorOrder(i+1,:), 'LineWidth', 3);
    % end
    % legend('ControlPFC','ImpairedPFC','ScaledControlPFC','ScaledImpairedPFC')

    %% plot the difference between scaled Ts between healthy and impaired (should always be > 0)
    f2 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
    for scaleIdx = 1:length(CaseScaledData)
        subplot(2,length(CaseScaledData),scaleIdx); hold on;
        fullData = CaseScaledData{scaleIdx};
        scaledDiff = fullData(:,:,6) - fullData(:,:,5); %stimulationSample x task difficulty
        [row, col] = find(scaledDiff < 0);%should both be empty, row = iteration#,col=datapoint index of the x array that doesn't work
        percentTimeSuccess = 1-length(unique(row)) / simulationIterations;
        [sortedImpairLevel, sortIdex] = sort(fullData(:,1,7));
        c = parula(simulationIterations);
        for it = 1:simulationIterations
            plot(x, scaledDiff(sortIdex(it),:), 'LineWidth', 1.5,'Color',c(it,:));
        end
        plot([min(x) max(x)], [0 0], 'k--','LineWidth',3)
        xlim([min(x) max(x)])
        title(sprintf('a=x%d: %.3f Success', scaleFactors(scaleIdx),percentTimeSuccess))
        xticks([])
        if scaleIdx == 1
            ylabel('ScaledImpaired - ScaledControl')
        elseif scaleIdx == length(CaseScaledData)
            legend({'Iteration'})
            cb = colorbar('TickLabels',0:100:1000);
            cb.Label.String = 'Sorted Impairment (0=least impaired)';
        end

        % plot heatmap showing the slope from less impaired to more impaired 
        sortedY2Scaled = fullData(sortIdex,:,6);
        sortedY2ScaledDiff = diff(sortedY2Scaled); %by default, diff by row (from 1 impair to next)
        subplot(2,length(CaseScaledData),length(CaseScaledData)+scaleIdx); hold on;
        imagesc(sortedY2ScaledDiff >=0); 
        xlim([0,length(x)+1])
        [rowImp, ~] = find(sortedY2ScaledDiff < 0);%should both be empty, row = iteration#,col=datapoint index of the x array that doesn't work
        percentTimeSuccess = 1-length(unique(rowImp)) / length(sortedY2ScaledDiff);
        title(sprintf('x%d: %.3f Success', scaleFactors(scaleIdx),percentTimeSuccess))
        if scaleIdx==1
            ylabel('Sample Sorted (0 = 2nd - least impaired)')
        elseif scaleIdx == length(CaseScaledData)
            colorbar
        end
    end

    sgtitle(sprintf('ScaledImpaired - ScaledControl Per Simulation (Ifr=%.2f, btm=%.2f)',perfInflectionRatioToPFCCtr, btmPerfCtlVal))
    han=axes(f2,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Task Difficulty')

    set(findall(gcf,'-property','FontSize'),'FontSize',20)

    % save results
    if saveFigure
        if not(isfolder(saveDir))
            mkdir(saveDir)
        end
        saveSuffix = ['btm' num2str(btmPerfCtlVal*100) '_bCtr' num2str(bPerfCtlVal*100)];
%         saveas(f1,[saveDir 'SimulatedRange_' saveSuffix]);
        saveas(f1,[saveDir 'SimulatedRange_' saveSuffix '.png']);
%         saveas(f2,[saveDir 'ScaledDiff_' saveSuffix]);
        saveas(f2,[saveDir 'ScaledDiff_' saveSuffix '.png']);
        clear f1 f2
    end

    %% visualize slices of task difficulty
    % f3 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
    % subplot(2,2,1); hold on;
    % taskDifficultyToEval = 35;
    % scatter(fullData(:,1,7), fullData(:,taskDifficultyToEval,6));
    % title(sprintf('~35%% of Max Task Difficulty (x=%.2f)',x(taskDifficultyToEval)))
    % subplot(2,2,2);
    % taskDifficultyToEval = 51;
    % scatter(fullData(:,1,7), fullData(:,taskDifficultyToEval,6));
    % title(sprintf('~50%% of Max Task Difficulty (x=%.2f)',x(taskDifficultyToEval)))
    % subplot(2,2,3);
    % taskDifficultyToEval = 65;
    % scatter(fullData(:,1,7), fullData(:,taskDifficultyToEval,6));
    % title(sprintf('~65%% of Max Task Difficulty (x=%.2f)',x(taskDifficultyToEval)))
    % subplot(2,2,4);
    % taskDifficultyToEval = 91;
    % scatter(fullData(:,1,7), fullData(:,taskDifficultyToEval,6));
    % title(sprintf('~90%% of Max Task Difficulty (x=%.2f)',x(taskDifficultyToEval)))
    % legend({'Iteration'})
    % sgtitle(sprintf('(x%d),Ifr=%.2f, btm=%.2f', scaleFactor,perfInflectionRatioToPFCCtr, btmPerfCtlVal))
    % 
    % han=axes(f3,'visible','off'); 
    % han.Title.Visible='on';
    % han.XLabel.Visible='on';
    % han.YLabel.Visible='on';
    % xlabel(han,'ImpairmentLevel (Higher is More Impaired)')
    % ylabel(han,'Scaled PFC (Higher is Worse)')
    % set(findall(gcf,'-property','FontSize'),'FontSize',20)

    % %% plot heatmap showing the slope from less impaired to more impaired 
    % 
    % sortedY2Scaled = fullData(sortIdex,:,6);
    % sortedY2ScaledDiff = diff(sortedY2Scaled); %by default, diff by row (from 1 impair to next)
    % % f7 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
    % % subplot(1,2,1); hold on; imagesc(sortedY2ScaledDiff); xlabel('Task Difficulty'); colorbar
    % % ylabel('Simulation Sample Sorted by Impairment Level (0 = 2nd - least impaired)')
    % % xlim([0,length(x)+1])
    % % title('Scaled PFC Difference (ImpairmentHigh-Low)')
    % % subplot(1,2,2); 
    % hold on; imagesc(sortedY2ScaledDiff >=0); xlabel('Task Difficulty'); colorbar
    % xlim([0,length(x)+1])
    % [rowImp, ~] = find(sortedY2ScaledDiff < 0);%should both be empty, row = iteration#,col=datapoint index of the x array that doesn't work
    % percentTimeSuccess = 1-length(unique(rowImp)) / length(sortedY2ScaledDiff);
    % title(sprintf('(x%d) PFC High Impair > Low Impair (%.3f Success, Ifr=%.2f, btm=%.2f)', scaleFactor,percentTimeSuccess,perfInflectionRatioToPFCCtr, btmPerfCtlVal))
    % set(findall(gcf,'-property','FontSize'),'FontSize',20)
    
    %% plot working vs not working examples
    for scaleIdx = 1:length(CaseScaledData)
        fullData = CaseScaledData{scaleIdx};
        scaledDiff = fullData(:,:,6) - fullData(:,:,5); %stimulationSample x task difficulty
        [row, col] = find(scaledDiff < 0);
        % plot examples of the trials that failed to scale.
        % [row1, col1] = find(sortedY2ScaledDiff < 0);
        % row1 = unique(row1);
        % failedSample = union(unique(row),row1);
        failedSample = unique(row);
        if ~isempty(failedSample)
        failedSample = randsample(failedSample,min(8,length(failedSample)));
        totalSampleToShow = 16;
        for badIdx = 1:min(length(failedSample),totalSampleToShow) %show max the first 8
        % badIteration = randsample(row,1);
            if badIdx == 1
                f4 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
                sgtitle(sprintf('Failing Examples (x%d, Ifr=%.2f, btm=%.2f)',scaleFactors(scaleIdx),perfInflectionRatioToPFCCtr, btmPerfCtlVal))      
            elseif badIdx == totalSampleToShow/2 + 1
                f5 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
                sgtitle(sprintf('Failing Examples (x%d, Ifr=%.2f, btm=%.2f)',scaleFactors(scaleIdx),perfInflectionRatioToPFCCtr, btmPerfCtlVal))      
            end

            subplot(4,4,mod(badIdx*2-1,16)); hold on;
            badIteration = failedSample(badIdx);
            for i = 1:4
                plot(x, fullData(badIteration,:,i), 'Color',colorOrder(i,:), 'LineWidth', 3);
            end
            xlim([min(x) max(x)])
        %     xticks([])
        %     yticks([])

            if badIdx == 1 || badIdx == totalSampleToShow/2+1
                xlabel('Task Difficulty')
                ylabel('PFC Activation or Performance')
            end

            if badIdx == 8 || badIdx == 16
                legend('Control','Impaired','ControlPerformance','ImpairedPerformance','Location','southoutside')
                subplot(4,4,16);
            else
                subplot(4,4,mod(badIdx*2,16));
            end
            hold on;
            for i = 5:6
                plot(x, fullData(badIteration,:,i), 'Color',colorOrder(i-4,:), 'LineWidth', 3);
            end
            if badIdx == 8 || badIdx == 16
                legend('Control','Impaired','Location','southoutside')
            end

            if badIdx == 1 || badIdx == 9
                title('Scaled')
                xlim([min(x) max(x)])
        %         xticks([])
                xlabel('Task Difficulty')
        %         yticks([])
                ylabel('Scaled PFC Activation')
            end
            title(sprintf('Failing Example (Sample %d, Impair: %.2f)',badIteration,fullData(badIteration,1,7)))
            set(findall(gcf,'-property','FontSize'),'FontSize',15)
        end
        end
        % plot examples of trials that scaled properly
        workingSample = setdiff(1:simulationIterations, unique(row)); %set(A,B) return data in A but not in B
        if ~isempty(workingSample)
            goodIterations = randsample(workingSample,min(2,length(workingSample)));
            f6 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
            for goodIdx = 1:length(goodIterations) 
                goodIteration = goodIterations(goodIdx);
                subplot(2,4,(goodIdx-1)*4+1:(goodIdx-1)*4+2); hold on;
                for i = 1:4
                    plot(x, fullData(goodIteration,:,i), 'Color',colorOrder(i,:), 'LineWidth', 3);
                end
                xlim([min(x) max(x)])

                if goodIdx == 1
                    legend('Control','Impaired','ControlPerformance','ImpairedPerformance','location','bestoutside')
                    xlabel('Task Difficulty')
                    ylabel('PFC Activation or Performance')
                end

                subplot(2,4,(goodIdx-1)*4+3:(goodIdx-1)*4+4); hold on;
                for i = 5:6
                    plot(x, fullData(goodIteration,:,i), 'Color',colorOrder(i-4,:), 'LineWidth', 3);
                end
                if goodIdx == 1
                    legend('Control','Impaired','location','bestoutside')
                    title('Scaled')
                    xlabel('Task Difficulty')
                    ylabel('Scaled PFC Activation')
                end
                title(sprintf('Working Example (Sample %d, Impair: %.2f)',goodIteration,fullData(goodIteration,1,7)))
                set(findall(gcf,'-property','FontSize'),'FontSize',12)
            end
            sgtitle(sprintf('Working Examples (x%d, Ifr=%.2f, btm=%.2f)',scaleFactors(scaleIdx),perfInflectionRatioToPFCCtr, btmPerfCtlVal))      
        end
        if saveFigure
            saveSuffix = [num2str(scaleFactors(scaleIdx)) 'btm' num2str(btmPerfCtlVal*100) '_bCtr' num2str(bPerfCtlVal*100)];
            if exist('f4','var')
%                 saveas(f4,[saveDir 'FailExample_batch1_x' saveSuffix]);
                saveas(f4,[saveDir 'FailExample_batch1_x' saveSuffix '.png']);
            end
            if exist('f6','var')
%                 saveas(f6,[saveDir 'WrokingExample_x' saveSuffix]);
                saveas(f6,[saveDir 'WrokingExample_x' saveSuffix '.png']);
            end
            clear f4 f6
        end
    end
    end
    FullScaledData{end+1} = CaseScaledData;
end

%% plot all cases together
% versions = versions(versions(:,1),:); %select subset
if ~plotIndVersions
    scaleFactorIdx = 3; %to index from scaleFactors[1,2,4,10]
    %prepare data
    sortedY2ScaledDiffByImpairLevel = nan(length(FullScaledData),simulationIterations-1,length(x)+1); %cases x samples x task difficulty (later pad by -1 to separate cases)
    sortedY2ScaledDiffByDifficulty = nan(length(FullScaledData),simulationIterations,length(x)); %cases x sample x task difficulty; later sample dimension will be padded by -1 to separate cases
    for ii = 1:length(FullScaledData)
        fullData = FullScaledData{ii}{scaleFactorIdx};
        [sortedImpairLevel, sortIdex] = sort(fullData(:,1,7));
        sortedY2Scaled = fullData(sortIdex,:,6); %samples X difficulty X 7; y1=pfcControl,y2=pfcImpaired,y3=PerfCtr,y4=PerfImp,y1scaled,y2scaled,impairlevel
        sortedY2ScaledDiff = diff(sortedY2Scaled); %by default, diff by row (from 1 impair to next)
        sortedY2ScaledDiffByImpairLevel(ii,:,:) = [sortedY2ScaledDiff>=0,-1*ones(size(sortedY2ScaledDiff,1),1)];
        %pad in -1 to separate cases
        sortedY2ScaledDiffByDifficulty(ii,:,:) = [sortedY2ScaledDiff>=0;-1*ones(1,size(sortedY2ScaledDiff,2))];
    end
    
    
    %% plot the results
    for paramsToPlot = 1:8%-3:-2 %-1 doesn't plot well, show white spaces
%         -1: all, -2: pfc, -3: performance; 4: 1 param at a time
        f1 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
        f2 = figure('Units','Normalized','OuterPosition',[0 0 1 1]); hold on;
    %     performance: [-5,-4,-1,4], index [1,2,5,8]; PFC: [-3,-2,2,3], index
    %     [3,4,6,7]
        if paramsToPlot == -1
            titleSuffix = 'All';
            rows = 1; cols = 1;
            subCases = [0];
            rowtitles = {'All'};
        elseif paramsToPlot <= -2
            rows = 4; cols = 1;
            rowtitles = {'Top','Bottom','X_0','Slope'};
            if paramsToPlot == -2
                titleSuffix = 'PFC';
                subCases = [3 4 6 7]; %top, bottom, inflection, slope
            else
                titleSuffix = 'Performance';
                subCases = [1 2 5 8]; %top, bottom, inflection, slope
            end
        else %1:8
            rows = 8; cols = 1;
            rowtitles = cell(1,8);
            titleSuffixForEachParam = {'PerfTop','PerfBottom','PFCTop','PFCBottom','PerfX_0','PFCX_0','PFCSlope','PerfSlope'};
            titleSuffix = titleSuffixForEachParam{paramsToPlot};
%             rowsWithParam = versions(versions(:,paramsToPlot),:); %top, bottom, inflection, slope of PFC then performance
            numParamChanged = sum(versions,2);
            subCases = cell(1,8);
            for ii =1:8
                rowtitles{ii} = [num2str(ii) 'Param Changing'];
                subCases{ii} = find(versions(:,paramsToPlot) & numParamChanged == ii);
            end
        end
        titleSuffix = [titleSuffix 'ScaledX' num2str(scaleFactors(scaleFactorIdx))];
        for caseIdx = 1:length(subCases)%1:length(caseCoding)
            %select sub cases to plot
            if paramsToPlot == -1
                subCase = logical(ones(size(versions,1),1)); %select all
            elseif paramsToPlot <= -2
                subCase = versions(:,subCases(caseIdx)); %find the column of 
                %the corresponding params and use that column as a logical array to index 
                %cases that param is changing
            else
                subCase = subCases{caseIdx};
            end
            byImpairLevel = permute(sortedY2ScaledDiffByImpairLevel(subCase,:,:),[2,3,1]);
            byImpairLevel = reshape(byImpairLevel,simulationIterations-1,[],1);
            set(0,'CurrentFigure',f1); 
            subplot(rows,cols,caseIdx);
            imagesc(byImpairLevel); colorbar('YTick',[-1 0 1]);
            title(rowtitles{caseIdx});
%             xlim([0 size(byImpairLevel,2)+1]);
%             ylim([0 size(byImpairLevel,1)+1]);
            
            byDifficulty = permute(sortedY2ScaledDiffByDifficulty(subCase,:,:),[3,2,1]);
            byDifficulty = reshape(byDifficulty,length(x),[],1);
            set(0,'CurrentFigure',f2); 
            subplot(rows,cols,caseIdx);
            imagesc(byDifficulty); colorbar('YTick',[-1 0 1]);
            title(rowtitles{caseIdx});
%             xlim([0 size(byImpairLevel,2)+1]);
%             ylim([0 size(byImpairLevel,1)+1]);
        end
        for f = 1:2
            if f == 1
                han=axes(f1,'visible','off');
                set(0,'CurrentFigure',f1); 
                xlabelStr = 'Task Difficulty in Increasing Parameter Choice Complexity (1 Param Left)';
                ylabelStr = 'Sample Sorted (0 = 2nd - least impaired)';
            else
                han=axes(f2,'visible','off');
                set(0,'CurrentFigure',f2); 
                xlabelStr = 'Sorted Sample in Increasing Parameter Choice Complexity (1 Param Left)';
                ylabelStr = 'Task Difficulty';
            end
%             han.Title.Visible='on';
            han.XLabel.Visible='on';
            han.YLabel.Visible='on';
            xlabel(han,xlabelStr)
            ylabel(han, ylabelStr);
%             title(han,[titleSuffix ' Params By ' ylabelStr]);
            
            sgtitle([titleSuffix ' Params By ' ylabelStr]);
            set(findall(han,'-property','FontSize'),'FontSize',20)
        end

        if saveFigure
            saveDir = [scriptDir '/Data/SimulationResults/AllParams/'];
            saveas(f1,[saveDir titleSuffix 'ParamByImpairmentLevel.png']);
            saveas(f1,[saveDir titleSuffix 'ParamByImpairmentLevel.fig']);
            saveas(f2,[saveDir titleSuffix 'ParamByTaskDifficulty.png']);
            saveas(f2,[saveDir titleSuffix 'ParamByTaskDifficulty.fig']);
        end
    end
end
%% plot reference exp functions
% close all;
% f = figure(); hold on;
% xref = -1.5:0.01:1.5;
% subplot(2,2,1); hold on;
% plot(xref, exp(xref),'LineWidth',3);
% ylim([-5,25])
% subplot(2,2,2); hold on;
% plot(xref, exp(xref),'LineWidth',3);
% plot(xref, exp(2*xref),'LineWidth',3);
% ylim([-5,25])
% subplot(2,2,3); hold on;
% plot(xref, exp(xref),'LineWidth',3);
% plot(xref, exp(2*xref),'LineWidth',3);
% plot(xref, exp(4*xref),'LineWidth',3);
% ylim([-5,25])
% subplot(2,2,4); hold on;
% plot(xref, exp(xref),'LineWidth',3);
% plot(xref, exp(2*xref),'LineWidth',3);
% plot(xref, exp(4*xref),'LineWidth',3);
% plot(xref, exp(10*xref),'LineWidth',3);
% % plot(xref,xref,'.--','LineWidth',3);
% % plot(xref,4*xref,'.--','LineWidth',3);
% % plot(xref,10*xref,'.--','LineWidth',3);
% plot(xref,zeros(1,length(xref)),'k--','LineWidth',1.5);
% legend('a=1','a=2','a=4','a=10')
% % legend('a=1','a=2','a=4','a=10','y=x','y=4*x','y=10*x')
% ylim([-5,25])
% sgtitle('Change of y = exp(a*x) as a Changes')
% han=axes(f,'visible','off'); 
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'Scaled PFC');
% xlabel(han,'Performance Change Normalized by Baseline');
% set(findall(gcf,'-property','FontSize'),'FontSize',30)

%% plot working example demo
% figure()
% subplot(1,3,1); hold on;
% plot(x, fullData(goodIteration,:,1), 'b', 'LineWidth', 3);
% plot(x, fullData(goodIteration,:,2), 'r', 'LineWidth', 3);
% xlim([xmin -xmin])
% title('PFC')
% legend('S1','S2')
% xticks([])
% xlabel('Task Difficulty')
% 
% subplot(1,3,2); hold on;
% plot(x, fullData(goodIteration,:,3), 'b--', 'LineWidth', 3);
% plot(x, fullData(goodIteration,:,4), 'r--', 'LineWidth', 3);
% title('Performance')
% xticks([])
% xlabel('Task Difficulty')
% 
% subplot(1,3,3); hold on;
% plot(x, fullData(goodIteration,:,5), 'b', 'LineWidth', 3);
% plot(x, fullData(goodIteration,:,6), 'r', 'LineWidth', 3);
% xticks([])
% xlabel('Task Difficulty')
% title('Attentional Cost')
% 
% % title(sprintf('Working Example (Sample %d, Impair: %.2f)',goodIteration,fullData(goodIteration,1,7)))
% set(findall(gcf,'-property','FontSize'),'FontSize',25)

%% plotting and saving data helper funciton
function fullData = plotAndStoreSimulationResult(plotFigure, x, a3,b3,top3,bottom3, y3, a4, b4, top4, bottom4, y4, y1, y2, impairLevel, fullData, i, scaleFactor, shape)
    %y3 is performance of control; y4 = performance of impaired group
    %y1 = PFC of control; y2 = PFC of impaired.
    if plotFigure
        colorOrder=colororder;
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        subplot(1,2,1);
        plot(x,y1,'LineWidth',3);
        hold on;
        if strcmp(shape, 'Sigmoid')
            ylim([-0.1 1.1])
        end
    %     xlim([xmin,6])
        plot(x,y2,'LineWidth',3);
        plot(x,y3,'LineWidth',3);
        plot(x,y4,'LineWidth',3);
        legend('Control','Impaired','ControlPerformance','ImpairedPerformance')
        title('Raw Schematics')
    end
    
    %scale performance
    %get value at 5*min(x), usually -50, to make sure get the true
    %top/baseline performance (if at -10 at small inflection value, the
    %data might be cutoff)
    if strcmp(shape, 'Sigmoid')
        perfBase3 = AutomaticityIndexHelper.sigmoidValue(20*min(x), a3, b3, top3, bottom3);
        perfBase4 = AutomaticityIndexHelper.sigmoidValue(20*min(x), a4, b4, top4, bottom4);
    elseif strcmp(shape, 'Linear')
        perfBase3 = AutomaticityIndexHelper.linearValue(min(x), a3, b3);
        perfBase4 = AutomaticityIndexHelper.linearValue(min(x), a4, b4);
    end
    deltaPerf3 = (perfBase3 - y3) / perfBase3;
    scale3 = (1/2 + (exp(scaleFactor*deltaPerf3))/2);
    y1Scaled = scale3 .* y1;

    deltaPerf4 = (perfBase4 - y4) / perfBase4; 
    scale4 = (1/2 + (exp(scaleFactor*deltaPerf4))/2);
    y2Scaled = scale4 .* y2;

    if plotFigure
        subplot(1,2,2);
        plot(x,y1,'LineWidth',3);
    %     xlim([xmin, 6])
        hold on;
        plot(x,y2,'LineWidth',3);
        plot(x,y1Scaled,'--','LineWidth',3, 'Color',colorOrder(1,:));
        plot(x,y2Scaled,'--','LineWidth',3, 'Color',colorOrder(2,:));
        legend('ControlPFC','ImpairedPFC','ScaledControlPFC','ScaledImpairedPFC')
        title('Scaled PFC(T-value)')
    end

    fullData(i,:,1) = y1;
    fullData(i,:,2) = y2;
    fullData(i,:,3) = y3;
    fullData(i,:,4) = y4;
    fullData(i,:,5) = y1Scaled;
    fullData(i,:,6) = y2Scaled;
    fullData(i,:,7) = impairLevel;
end


%% older version + version notes
%% version 5: impairment is defined by top and bottom + b (inflection point) + slope of performance only
%Assumptions: 1. the performance decays to mid point (inflection point) when PFC reaches 85% of max, true for
%control and impaired
% 2. Control performance have same top as control PFC (1), and bottom in
% range (0.75,0.9)
% 3. Impaired and control performance have different slope (control: between (-5,-1), 
% faster slope than PFC curves) and ****** impaired decays faster than control******
% Different in range. For impaired, larger the impairLevel, larger then range 
% and smaller the top performance level, but the
% top performance is kept above 0.5 for all impairment level. 
% 4. Impaired and control PFC have different
% slope, range and inflection point. ****** Slope: impaired has less steep slope
% than healthy. ******
% Range: Impaired is centered around 0.5 (mid y value of control) and range
%varies between 1-(0.2,0.9). Inflection point: impaired start
%to increase early (inflection point b < 0, b = (b1-minx) * (1-impairLevel) + minx =
%b1*impairLevel + minx*impairLevel = minx*impairLevel, smaller impairLevel, 
% less negative b, starts increase later)
%% version 4: impairment is defined by top and bottom + b (inflection point) + slope of PFC AND performance
%Assumptions: 1. the performance decays to mid point (inflection point) when PFC reaches 85% of max, true for
%control and impaired
% 2. Control performance have same top as control PFC (1), and bottom in
% range (0.75,0.9)
% 3. Impaired and control performance have different slope (control: between (-5,-1), 
% faster slope than PFC curves) and ****** impaired decays faster than control******
% Different in range. For impaired, larger the impairLevel, larger then range 
% and smaller the top performance level, but the
% top performance is kept above 0.5 for all impairment level. 
% 4. Impaired and control PFC have different
% slope, range and inflection point. ****** Slope: impaired has less steep slope
% than healthy. ******
% Range: Impaired is centered around 0.5 (mid y value of control) and range
%varies between 1-(0.2,0.9). Inflection point: impaired start
%to increase early (inflection point b < 0, b = (b1-minx) * (1-impairLevel) + minx =
%b1*impairLevel + minx*impairLevel = minx*impairLevel, smaller impairLevel, 
% less negative b, starts increase later)


%% version 3: impairment is defined by top and bottom + b (inflection point) + slope of PFC
%Assumptions: 1. the performance decays to mid point (inflection point) when PFC reaches 85% of max, true for
%control and impaired
% 2. Control performance have same top as control PFC (1), and bottom in
% range (0.75,0.9)
% *** 3. Impaired and control performance have different slope (control: between (-5,-1), 
% faster slope than PFC curves). Different in range. For impaired, larger the
% impairLevel, larger then range and smaller the top performance level, but the
% top performance is kept above 0.5 for all impairment level. 
% 4. Impaired and control PFC have different
% slope, range and inflection point. *** Slope: impaired has less steep slope
% than healthy. ***
% Range: Impaired is centered around 0.5 (mid y value of control) and range
%varies between 1-(0.2,0.9). Inflection point: impaired start
%to increase early (inflection point b < 0, b = (b1-minx) * (1-impairLevel) + minx =
%b1*impairLevel + minx*impairLevel = minx*impairLevel, smaller impairLevel, 
% less negative b, starts increase later)

%% version 2: impairment is defined by top and bottom + b (inflection point)
%Assumptions: 1. the performance decays to mid point (inflection point) when PFC reaches 85% of max, true for
%control and impaired
% 2. Control performance have same top as control PFC (1), and bottom in
% range (0.75,0.9)
% 3. Impaired and control performance have the same slope (between (-5,-1), 
% faster slope than PFC curves). Different in range. For impaired, larger the
% impairLevel, larger then range and smaller the top performance level, but the
% top performance is kept above 0.5 for all impairment level.
% Changing Assumption: 4. Impaired and control PFC have the same slope (a), differ in
%range and inflection point. Range: Impaired is centered around 0.5 (mid y value of control) and range
%varies between 1-(0.2,0.9). Inflection point: impaired start
%to increase early (inflection point b < 0, b = (b1-minx) * (1-impairLevel) + minx =
%b1*impairLevel + minx*impairLevel = minx*impairLevel, smaller impairLevel, 
% less negative b, starts increase later)

%% version 1: impairment is defined by top and bottom only
%Assumptions: 1. the performance decays to mid point (inflection point) when PFC reaches 85% of max, true for
%control and impaired
% 2. Control performance have same top as control PFC (1), and bottom is
% fixed to 0.75. (reference control performance)
% 3. Impaired and control performance have the same slope (between (-5,-1), 
% faster slope than PFC curves). Different in range. For impaired, larger the
% impairLevel, larger then range and smaller the top performance level, but the
% top performance is kept above 0.5 for all impairment level.
% **Changing Assumption 3. Impaired and control PFC have the same slope and b, only differ in
%range. Impaired is centered around 0.5 (mid y value of control) and range
%varies between 1-(0.2,0.9).

%% version 0 - first attempt
% if version == -10
%     for i = 1:simulationIterations
%     a1 = 1; %slope,  smaller a smaller slope
%     b1 = 0; %shifts, larger b, curve shifted more towards left (negative end of x-axis)
%     %b = x value where the response is halfway between max and min
%     top1 = 1;
%     bottom1 = 0;
%     y1 = AutomaticityIndexHelper.sigmoidValue(x, a1, b1, top1, bottom1);
%     plateauX1 = (log((top1-bottom1)/(0.9-bottom1) -1) -b1)/(-a1);
% 
%     %smaller top, smaller y range, shallower slope (55 - 95% of healthy),
%     %plateau before healthy reaches 70-80% (assume it happens at x = 1)
%     topPerfCtl = AutomaticityIndexHelper.generateRandValInRange(0.5,0.95); %(0.5,0.9)
%     btmPerfCtl = topPerfCtl - 0.5; %always -0.5 of top
%     aPerfCtl = AutomaticityIndexHelper.generateRandValInRange(0.5,1);%at least half of the healthy curve slope(0.55,0.95)
%     plateauX2 = AutomaticityIndexHelper.generateRandValInRange(-3,1); 
%     bPerfCtl = AutomaticityIndexHelper.findMatchingb(topPerfCtl, btmPerfCtl, aPerfCtl, plateauX2);
%     yPerfCtl = AutomaticityIndexHelper.sigmoidValue(x, aPerfCtl, bPerfCtl, topPerfCtl, btmPerfCtl);
% 
%     %performance curves, same top, much smaller range, negative and sharper
%     %slope, plateau (curve starts dropping before healthy PFC curve plauteau
%     %(or before it reaches 70-80% of the max), but remains steady after pass
%     %50% of PFC curve.
%     topPFCImpaired = 1; %arbitrary chosen set value
%     btmPFCImpaired = AutomaticityIndexHelper.generateRandValInRange(0.75,0.9);
%     aPFCImpaired = AutomaticityIndexHelper.generateRandValInRange(-5,-1);%(-10,-1)
%     % plateauX3 = AutomaticityIndexHelper.generateRandValInRange(0,1);
%     plateauX3 = plateauX1 - AutomaticityIndexHelper.generateRandValInRange(0,0.25); %starts dropping 0.5-1 before the pfc plateau
%     bPFCImpaired = AutomaticityIndexHelper.findMatchingb(topPFCImpaired, btmPFCImpaired, aPFCImpaired, plateauX3);
%     yPFCImpaired = AutomaticityIndexHelper.sigmoidValue(x, aPFCImpaired, bPFCImpaired, topPFCImpaired, btmPFCImpaired);
% 
%     %performance curve for impaired, same top as PFC of impaired, larger y
%     %range (0.6 instead of 0.5), shallower slope than healthy perf but not too
%     %slow (0.5-1) * slopeperfhealthy, plateau 0.5-1 before the impaired PFC
%     topPerfImpaired = AutomaticityIndexHelper.generateRandValInRange(0.6,0.9); 
%     btmPerfImpaired = topPerfImpaired - 0.6; %should be all above 0
%     aPerfImpaired = AutomaticityIndexHelper.generateRandValInRange(-1,-0.5);%(-10,-1)
%     % a4 = a3*AutomaticityIndexHelper.generateRandValInRange0.2,0.6);%(-10,-1)
%     plateauX4 = plateauX2 - AutomaticityIndexHelper.generateRandValInRange(0.5,1);
%     bPerfImpaired = AutomaticityIndexHelper.findMatchingb(topPerfImpaired, btmPerfImpaired, aPerfImpaired, plateauX4);
%     yPerfImpaired = AutomaticityIndexHelper.sigmoidValue(x, aPerfImpaired, bPerfImpaired, topPerfImpaired, btmPerfImpaired);
%     xmin = -20;
% 
%     if plotFigure
%         fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
%         subplot(1,2,1);
%         plot(x,y1,'LineWidth',3);
%         hold on;
%         ylim([-0.1 1.1])
%         xlim([xmin,6])
%         plot(x,yPerfCtl,'LineWidth',3);
%         plot(x,yPFCImpaired,'LineWidth',3);
%         plot(x,yPerfImpaired,'LineWidth',3);
%         legend('Control','Impaired','ControlPerformance','ImpairedPerformance')
%         title('Raw Schematics')
%     end
% 
%     %healthy perf
%     perfBase3 = AutomaticityIndexHelper.sigmoidValue(xmin, aPFCImpaired, bPFCImpaired, topPFCImpaired, btmPFCImpaired);
%     deltaPerf3 = (perfBase3 - yPFCImpaired) / perfBase3;
%     scale3 = (1/2 + (exp(100*deltaPerf3))/2);
%     y1Scaled = scale3 .* y1;
% 
%     perfBase4 = AutomaticityIndexHelper.sigmoidValue(xmin, aPerfImpaired, bPerfImpaired, topPerfImpaired, btmPerfImpaired);
%     deltaPerf4 = (perfBase4 - yPerfImpaired) / perfBase4; 
%     scale4 = (1/2 + (exp(100*deltaPerf4))/2);
%     y2Scaled = scale4 .* yPerfCtl;
% 
%     if plotFigure
%         subplot(1,2,2);
%         plot(x,y1,'LineWidth',3);
%         xlim([xmin, 6])
%         hold on;
%         plot(x,yPerfCtl,'LineWidth',3);
%         plot(x,y1Scaled,'--','LineWidth',3);
%         plot(x,y2Scaled,'--','LineWidth',3);
%         legend('ControlPFC','ImpairedPFC','ScaledControlPFC','ScaledImpairedPFC')
%         title('Scaled PFC(T-value)')
%     end
% 
%     fullData(i,:,1) = y1;
%     fullData(i,:,2) = yPerfCtl;
%     fullData(i,:,3) = yPFCImpaired;
%     fullData(i,:,4) = yPerfImpaired;
%     fullData(i,:,5) = y1Scaled;
%     fullData(i,:,6) = y2Scaled;
% 
%     % subplot(1,3,3); hold on;
%     % plot(x,scale3,'LineWidth',3)
%     % plot(x,scale4,'LineWidth',3)
%     % legend('ScaleHealthy','ScaleImpaired')
%     % title('Scaling Factor')
% 
%     end
% end


%% older version
% % the impaired group should have shallower slope (smaller a), shifted to
% % the left (platue sooner, b>b1), smaller range (scale y) & larger
% % intercept (shift up)
% aScale = 0.5; %between (0,1)
% bScale = 5;% >1
% rangeScale = 0.5; %between (0,1)
% interceptShift = 0.15; % between (0,1-rangeScale)

% y1 = 1./(1+exp(-a1*x+b1));
% aMatch1 = findMatchingSlope(1, 0, 1, 3.19)
% aMatch2 = findMatchingSlope(top, bottom, b2, -1)
% 
% a2 = a1*aScale;
% % b2 = b1*bScale;
% top = 0.8; %(bottom,1)
% bottom = 0.2; %(0,top)
% b2 = -7; %halfway y reached before 9, [-5,0)
% % y2 = interceptShift + rangeScale * 1./(1+exp(-(a2*x+b2)));
% y2 = bottom + (top - bottom)./(1+exp(-a2*x+b2));
% plot(x,y2,'LineWidth',3);

% 
% %performance should start to drop slower than PFC increase and platau later
% %(shift to right)
% rangeScalePerf = 0.2; %much smaller scale, way << 1
% bScalePerf = -2; % <1
% perfy1 = 1 -rangeScalePerf * 1./(1+exp(-(a1*x+b1*bScalePerf)));
% plot(x,perfy1,'LineWidth',3);
% 
% % performance of impaired group, should plateau later, drop similar or
% % smaller/flatter, similar or larger range than PFC of impaired
% range2ScalePerf = 1; %similar or larger, [1, 1.5)
% b2ScalePerf = -0.1; % plateau later, but not too late, smaller b, (-0.1,0.7)
% a2ScalePerf = 0.7; %(0,1), more likely [0.5,1) %this change is linked to b2 scale
% perfy2 = 1 - (interceptShift + rangeScale *range2ScalePerf * 1./(1+exp(-(a2*a2ScalePerf*x+b2*b2ScalePerf))));
% plot(x,perfy2,'LineWidth',3);
% 
% % range2ScalePerf = 2; %large (>1)
% % b2ScalePerf = 0.5; % plateau later, but not too late, shifts to the right b>1
% % a2ScalePerf = 0.3; %(0,1) slower decrease, probably (0,0.5)
% % interceptShiftPerf = 0.3; %(0, 1-rangeScalePerf * range2ScalePerf-interceptShift
% % perfy2 = 1 - rangeScalePerf * range2ScalePerf * 1./(1+exp(-(a1*a2ScalePerf*x+b1*bScalePerf*b2ScalePerf))) - interceptShiftPerf;
% % plot(x,perfy2,'LineWidth',3);
% 