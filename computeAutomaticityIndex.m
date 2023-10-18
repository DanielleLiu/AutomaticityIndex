function [subjectID, hboData, hbrData, scaledData, responders, perfDataArray, deltaPerf, deltaPerfGaitOrCog, scaleFactors, missingDataMask]=computeAutomaticityIndex(perfDataPath, pfcDataPath, studyID, visitNum, saveResAndFigure, saveDir, scaleMethod, logData, plotNormalizeDiff)
% (perfDataPath, pfcDataPath, studyID, visitNum, saveResAndFigure, saveDir, scaleMethod, logData, plotNormalizeDiff)
if nargin < 3 || isempty(visitNum) || isnan(visitNum)
    visitNum = 1; %Default visit 1.
end
% default scaleMethod = 'Exp', this arg has to use keywords exactly as Exp
% or Linear
if nargin <= 6 || isempty(scaleMethod) || any(isnan(scaleMethod))
    scaleMethod = 'Exp'; %Linear or Exp, default Exp
end
%logData default false.
%true if shift unscaled to positive and log y-axis for all; false o.w.
if nargin <= 7 || isempty(logData) || any(isnan(logData))
    logData = false; %Linear or Exp, default Exp
end

if nargin <= 8 || isempty(plotNormalizeDiff) || any(isnan(plotNormalizeDiff))
    plotNormalizeDiff = false; %default false. this is for plotting to visualize data only.
    %calculation of index always used normalized data perf
end

%% load perf data
% perfData = readcell('/Users/mac/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofPittsburgh/BEAM lab - fNIRS Workgroup/Rosso_fNIRS_code/Matlab/Emma_SampleDataNMCM_forShuqi/NMCM_PRIMA_MOVE_combined_wide_19Nov2021.csv');
% /Users/mac/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SML/Projects/fNIR Project/Code_NIRS_Automaticity/Data/EmmaDataset/PRIMA_051722
% load 19Nov2021 Data
% perfData = readcell([scriptDir 'EmmaDataset' filesep 'NMCM_PRIMA_MOVE_combined_wide_19Nov2021.csv']);
perfData = readcell(perfDataPath);
% readcell([scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'Performance_GaitSpeed_Alphabet' filesep 'NMCM_PRIMA_MMH_Visit1_performance_subjavg.csv']);

perfDataHeader = perfData(1,:);
studyIdCol = find(strcmp(perfDataHeader,'Cohort'));
perfDataFull = perfData(contains(perfData(:,studyIdCol),studyID),:);
visitCol = strcmp(perfDataHeader,'Visit');
perfDataFull = perfDataFull(strcmp(perfDataFull(:,visitCol),['Visit ' num2str(visitNum)]) |...
    any(cell2mat(perfDataFull(:,visitCol)) == visitNum,2),:);

%in order: ID (numerical version, 10xxx, 2xxx, 21xxx, study, task, surface, gait, alphabet
%the provided array are in order of the data header, so the indexed sub
%data will have headers in the same order.
relevantCols = ismember(perfDataHeader,{'MergeID','Cohort','cond','GaitSpeed_mps','RCRT'});
perfData = perfDataFull(:,relevantCols); %the perfData now will contain the 5 columns above in the same order listed above.
% load old data 19Nov2021
% perfData = perfDataFull(:,[1:4,gaitSpeedIdx , alphabetIdx]);
%the index here is the sub column (following order in line 44.
gaitSpeedCol = 4; alphabetCol = 5; condCol = 3; idCol = 1; 
walkData = perfData(strcmp(perfData(:,condCol),'Even'),:);
walk2Data = perfData(strcmp(perfData(:,condCol),'Even_ABC'),:);
s2Data = perfData(strcmp(perfData(:,condCol),'Standing_ABC'),:);
walkUneven2Data = perfData(strcmp(perfData(:,condCol),'Uneven_ABC'),:);
checkID = isequal(walkData(:,idCol),walk2Data(:,idCol),s2Data(:,idCol),walkUneven2Data(:,idCol));
if ~checkID
    warning('Gait and cognitive performance data ID did NOT match.')
end

%Remove missing data (PRIMA 2 subject gait speed missing, remove this person for all)
missingDataMask = cellfun(@ismissing, [walkData(:,gaitSpeedCol),walk2Data(:,[alphabetCol,gaitSpeedCol]),s2Data(:,alphabetCol),walkUneven2Data(:,[alphabetCol,gaitSpeedCol])]);
missingDataMask = any(missingDataMask,2);
missingID = walkData(missingDataMask,idCol) %number version: 10XXX for NMCM, 21xxx for PRIMA and 30xxx for MHT
%previously missing 21180, 21282 PRIMA. but value exists in 01272023Data
%If using string version of iD, converts to the string id version: xxxxNIRS
% for mIdx = 1:numel(missingID) 
%     mid = missingID{mIdx};
%     missingID{mIdx} = [num2str(mid-20000) 'NIRS'];
% end
walkData(missingDataMask,:) = [];
walk2Data(missingDataMask,:) = [];
walkUneven2Data(missingDataMask,:) = [];
s2Data(missingDataMask,:) = [];

subjectID = walkData(:,1);

subjectID = cellstr(num2str(cell2mat(subjectID))); 
perfDataArray = cell2mat([walkData(:,gaitSpeedCol),s2Data(:,alphabetCol),walk2Data(:,[gaitSpeedCol,alphabetCol]),walkUneven2Data(:,[gaitSpeedCol,alphabetCol])]);


%% load PFC data
% pfcData = readcell('/Users/mac/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofPittsburgh/BEAM lab - fNIRS Workgroup/Rosso_fNIRS_code/Matlab/Emma_SampleDataNMCM_forShuqi/SubjStats_ROIw_PFC.csv');
% %load previously processed data. Version 051722
% if iscell(studyID)
%     pfcDataFull = readcell([scriptDir 'EmmaDataset' filesep 'SubjStats_ROIw_PFC_NMCA_PRIMA_051722.csv']);
% else
%     pfcDataFull = readcell([scriptDir 'EmmaDataset' filesep studyID '_051722' filesep 'SubjStats_ROIw_PFC_' studyID '_051722.csv']);
% end
%always load all studies. Version 12/01/2022
% pfcData = readcell([scriptDir 'EmmaDataset' filesep 'SubjStats_ROIw_PFC_NMCM_112222.csv']);
% pfcDataFull = [pfcData, repmat("NMCM", size(pfcData,1),1)];
% pfcDataFull(1,end) = 'Cohort';
% pfcData = readcell([scriptDir 'EmmaDataset' filesep 'SubjStats_ROIw_PFC_PRIMA_111822.csv']);
% pfcData = [pfcData, repmat("PRIMA", size(pfcData,1),1)];
% pfcDataFull = [pfcDataFull; pfcData];
% pfcData = readcell([scriptDir 'EmmaDataset' filesep 'SubjStats_ROIw_PFC_MMH_111822.csv']);
% pfcData = [pfcData, repmat("MMH", size(pfcData,1),1)];
% pfcDataFull = [pfcDataFull; pfcData];
%load new data, 1/27/2022
pfcDataFull = readcell(pfcDataPath);
%readcell([scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'fNIRS' filesep 'NMCM_PRIMA_MMH_Visit1_fNIRS_subjavg.csv']);
pfcDataHeader = pfcDataFull(1,:);

% if loaded all data, (version 12/01/2022)
%subselect the current studies of interest.
cohortCol = find(strcmp(pfcDataHeader,'Cohort'));

if iscell(studyID) %multiple studies
    pfcData = pfcDataFull(contains(pfcDataFull(:,cohortCol),studyID),:);
else %single study, options: Move MYHAT, NMCM, PRIMA.
    pfcData = pfcDataFull(strcmp(pfcDataFull(:,cohortCol),studyID),:);
end
% % if loading specific studies only (version 051722)
% pfcData = pfcDataFull;

%sub select session 1 only and whole PFC only.
visitCol = strcmp(pfcDataHeader,'Visit');
pfcData = pfcData(strcmp(pfcData(:,visitCol),['Visit ' num2str(visitNum)]) |...
    any(cell2mat(pfcData(:,visitCol)) == visitNum,2),:);

roiCol = find(strcmp(pfcDataHeader,'ROI'));
pfcData = pfcData(strcmp(pfcData(:,roiCol),'PFC'),:); %select whole brain PFC result only.
% tasks = {'Even','Standing_ABC','Even_ABC','Uneven_ABC'};
% pfcData = pfcData(ismember(pfcData(:,taskCol),tasks),:);

%need to get walk, s2, w2, walkUneven2; for hbo and hbr
typeCol = find(strcmp(pfcDataHeader,'type'));
taskCol = find(strcmp(pfcDataHeader,'cond')); %find(strcmp(pfcDataHeader,'Contrast'));
tCol = find(strcmp(pfcDataHeader,'tstat'));%find(strcmp(pfcDataHeader,'T'));

%remove missing data
idCol = find(strcmp(pfcDataHeader,'MergeID'));
pfcData(ismember(cell2mat(pfcData(:,idCol)),cell2mat(missingID)),:)=[]; %remove rows for the person with missing data.

%check if pfc and perf data match.
walkData = pfcData(strcmp(pfcData(:,taskCol),'Even'),:);
s2Data = pfcData(strcmp(pfcData(:,taskCol),'Standing_ABC'),:);
walk2Data = pfcData(strcmp(pfcData(:,taskCol),'Even_ABC'),:);
walkUneven2Data = pfcData(strcmp(pfcData(:,taskCol),'Uneven_ABC'),:);
checkID = isequal(walkData(:,idCol),walk2Data(:,idCol),s2Data(:,idCol),walkUneven2Data(:,idCol));
if ~checkID
    warning('ID of PFC data did NOT match');
end
checkType = isequal(walkData(:,typeCol),walk2Data(:,typeCol),s2Data(:,typeCol),walkUneven2Data(:,typeCol));
if ~checkType
    warning('ID of data type hbo/hbr did NOT match');
end

checkID = [];
% mergeIdCol = find(strcmp(pfcDataHeader,'MergeID'));

%remove subjects with misssing PFC or performance data.
%find missing PFC data, subjectID comes from performance
[diffId, idxPerf] = setdiff(str2num(cell2mat(subjectID)),cell2mat(walkData(1:2:end,idCol)))
perfDataArray(idxPerf,:)=[]; %remove entry that has performance but not PFC.
subjectID(idxPerf)=[];
% %TODO: do this for things in PFC but not in performance; and this is
% probably better than the for loop to check for IDs&if there is also
% gaurantee of same ordering.

for id = 1:length(subjectID)
%     %check ID when ID is in XXXNIRS format
%     if str2num(subjectID{id}) < 20000 %naming convension for NMCM: 10xxx in combined, or xxxNIRS where xxx is subject id
%         checkID(end+1) = str2num(subjectID{id}) == 10000+str2num(walk2Data{id*2,idCol}(1:3));
%     elseif str2num(subjectID{id}) > 21000 && str2num(subjectID{id}) < 30000  %PRIMA naming convention: 2xxx in combined, or xxxxNIRS where xxxx is subject ID
%         checkID(end+1) = str2num(subjectID{id}) == 20000+str2num(walk2Data{id*2,idCol}(1:4));
%     elseif str2num(subjectID{id}) >= 30000 %MMH naming convention: 30xxx in combined, and # 1-117 in combined perf data
%         checkID(end+1) = str2num(subjectID{id}) - 30000 == str2num(walk2Data{id*2,idCol});
%     end
% check ID when both IDs have the numerical format
    checkID(end+1) = str2num(subjectID{id}) == walk2Data{id*2,idCol};
end
if ~all(checkID)
    warning('PFC data and performance data subject ID did NOT match')
end

pfcDataTstats = [walkData(:,[typeCol,tCol]),s2Data(:,tCol),walk2Data(:,tCol),walkUneven2Data(:,tCol)];
%in subject x task format, task order: walkeven, s2, w2,
%walkuneven2, 1st column list the type (Hbo or Hbr)

%% setup variabile of performance and PFC in various format that will be used later.
% %if load full data (version 12/01/2022)
% hboData = str2double(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbo'),2:end));
% hbrData = -str2double(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbr'),2:end)); %hbr should be oppositive of hbo (more negative means more activation)
% %if load specific study 05/16/2022 version or load full new data 01/27/23
% version
hboData = cell2mat(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbo'),2:end));
hbrData = -cell2mat(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbr'),2:end)); %hbr should be oppositive of hbo (more negative means more activation)

deltaPerf = (-perfDataArray(:,1)+perfDataArray(:,[1,3,5])) ./ perfDataArray(:,1) + (-perfDataArray(:,2)+perfDataArray(:,[2,4,6])) ./ perfDataArray(:,2) ;
%deltaWalk for W2, W2Uneven, then deltaCog for W2, W2Unevent
deltaPerfGaitOrCog = [(-perfDataArray(:,1)+perfDataArray(:,[3,5])) ./ perfDataArray(:,1), (-perfDataArray(:,2)+perfDataArray(:,[4,6])) ./ perfDataArray(:,2)];
%deltaPer has all 0s for 1st column (baseline), then w2, w2u; pad another
%column of 0 for 2 baselines (walk only, alphabet only)
deltaPerf = [zeros(size(deltaPerf,1),1),deltaPerf]; 
% subject x task, tasks: w2 and walkUneven2, 
%only keep data for S2, Walkeven, Walk2, Walkuneven2. Remove WalkUneven since no composite
%performance exists.
%(walk- walk2) / walk + (s2 - walk2)/s2
%% scale and plot other dataset
% clc;

if plotNormalizeDiff
    normSaveSuffix = '_NormDiff';
else
    normSaveSuffix = '';
end
colorOrder = colororder;
scaledData = cell(2,4); %rows: title, then data; columns: hbo scaled both, hbr scaled by both
responders={[]}; %record non responders, scaleFactors x 2, columns = hbo, hbr
%rows = 1 col per scaleFactor in order of scaleFactors, each cell 
%contains 2 column vectors (could be different length) of nonresponder indexes, The 2 column vectors are
%in order: w2, w2u
exampleFactorAroundOptimalToPlot = 0; %how many examples +- optimal to plot. default = 0, only plot the optimal one. 
% if 1, will plot 2 examples (one 0.5 below optimal, one 0.5 above
% optimal); if 2 will plot 4 examples (2 below and 2 above).
totalSubPltCol = 2 * (exampleFactorAroundOptimalToPlot + 2); %2 x (examples + optimal + perf or pfc unscaled)
totalRow = 2;
blueColor = [0 0.4470 0.7410]; %blue
greyColor = [0.7 0.7 0.7]; %grey.
for i = 1:2 %1=hbo, 2=hbr
    fPerf = figure('units','normalized','outerposition',[0 0 1 1]);
    if i == 1
        pfcdata = hboData;
        titleStr = 'Hbo';
    else
        pfcdata = hbrData;
        titleStr = 'Hbr';
    end
    % shift to be in all positive range (not necessarily right, losing info about increase/decrease compared to rest)
    shiftAmount = abs(min(pfcdata(:,3:4),[],'all')) + 1;
    dataPfcHboVsBaseInOrderShifted = pfcdata + shiftAmount;
    
    %plot performance change normalized by baseline perf.
    s1=subplot(totalRow,totalSubPltCol,1);
    axis square;
    hold on;
%     plot(1:3,deltaPerf(:,2:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', 'k');%plot all with 1 color

    increaseLogicalPerf = deltaPerf(:,4) - deltaPerf(:,3) > 0; %W2U - W2
    increaseRatioPerf = sum(increaseLogicalPerf) / size(deltaPerf,1);
%     plot(1:2,deltaPerf(:,2:3),'o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%decrease group
    plot(2:3,deltaPerf(increaseLogicalPerf,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%gray (unexpected)
%     ylim([-0.65 0.67]); 
    ylabel('\DeltaPerformance');
    xticks([2, 3]); xticklabels({'evenABC','unevenABC'})
    title(sprintf('evenABC -> unevenABC\nIncrease (unexpected): %.2f', increaseRatioPerf));
    
    s1=subplot(totalRow,totalSubPltCol,totalSubPltCol+1); %plot on 2nd row.
    axis square;
    hold on;
    plot(2:3,deltaPerf(~increaseLogicalPerf,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', blueColor);%blue (expected, non-change or decrease)
%     ylim([-0.65 0.67])
%     lineItems = findobj(s1,'Type','Line');
%     legend(lineItems([1,end]),'Performance Decrease (Expected)','Performance Increase (Unexpected)')
        
%     %plot each subject as its own color
%     for dIdx = 1:size(deltaPerf,1)
% %         plot([1:size(deltaPerf,2)],deltaPerf(dIdx,3:4),'o-','LineWidth',2.5,'MarkerSize',5,'DisplayName',subjectID{dIdx,1});
%         plot(1:3,deltaPerf(dIdx,2:4),'o-','LineWidth',1.5,'MarkerSize',5,'DisplayName',subjectID{dIdx,1});
%     end
    xticks([2, 3])
%     xticklabels({'Wk','S2','W2','W2Uneven'})
    xticklabels({'evenABC','unevenABC'})
    ylabel('\DeltaPerformance')
    % ylabel('(SpeedTask-SpeedWalk)/SpeedWalk + (AlphaTask - AlphaS2) / AlphaS2')

    %plot bar graph of the %subjects
%     figure(fCounts); subplot(1,totalSubPltCol,1); hold on; 
%     axis square;
    title(sprintf('evenABC -> unevenABC\nDecrease (expected): %.2f', 1-increaseRatioPerf));
%     bar(1,increaseRatio,1,'FaceColor',greyColor,'DisplayName','Unexpected (Perf Increase)');
%     bar(2,-100+increaseRatio,1,'FaceColor',blueColor,'DisplayName','Expected (Perf Decrease)');
%     b = bar([increaseRatio -100+increaseRatio; 55 45; 77 23],'Stacked')
%     xlim([0.5 2.5])
%    

    %plot perf of 2 groups in CI.
%     fCI = figure('units','normalized','outerposition',[0 0 1 1]);
%     disp('Perf W2U-W2')
    %plot on the 3rd row
    %%Plot the CI of the whole group.
%     subplot(totalRow,totalSubPltCol*2,4*totalSubPltCol+1); axis square; hold on;
%     PlotHelper.plotCI(1, deltaPerf(:,4) - deltaPerf(:,3),'k','Full', true)
%     xlim([0.8,1.2]); xticks([]); 
%     title(sprintf('Perf W2U-W2 Dec: %.2f', 1-increaseRatioPerf)); legend();
    %plot on the 2nd column (unexpected first row, expected 2nd row)
    subplot(totalRow,totalSubPltCol,2); axis square;grid on; hold on;
%     %plot CI per group
%     PlotHelper.plotCI(1, deltaPerf(increaseLogicalPerf,4) - deltaPerf(increaseLogicalPerf,3),greyColor,'Inc(Unexp)', false) %gray, unexpected increase
%     xlim([0.8 1.2]); xticks([]); legend(); ylim([0 0.5]); yticks([0 0.2 0.4]); grid on; grid minor
    %%plot the difference between W2U-WU as scatter.   
    if plotNormalizeDiff%normalized version
        dataToplot = (deltaPerf(increaseLogicalPerf,4) - deltaPerf(increaseLogicalPerf,3))./abs(deltaPerf(increaseLogicalPerf,3));
        ylabel('(unevenABC - uevenABC) / uevenABC');
    else %%non-normalized version of the difference 
        dataToplot = deltaPerf(increaseLogicalPerf,4) - deltaPerf(increaseLogicalPerf,3);
        ylabel('uevenABC - evenABC');
    end
    plot(ones(size(dataToplot)),dataToplot,'o','Color',greyColor, 'LineWidth',2,'MarkerSize',10,'DisplayName','Inc(unexp)');
    xticks([]); grid on; xlabel('Inc(unexp)'); ylim([0 1.1*max(dataToplot)]); %ylim([0 0.5]); 
    title('Amount of increase')
%     legend(); %legend('location','northeastoutside'); 
    %plot on 2nd row 2nd column (expected perf)
    subplot(totalRow,totalSubPltCol,totalSubPltCol+2); axis square;
%     PlotHelper.plotCI(1, deltaPerf(~increaseLogicalPerf,4) - deltaPerf(~increaseLogicalPerf,3),blueColor,'Dec(Exp)', false) %blue, expected decrease
%     xlim([0.8 1.2]); xticks([]); legend(); ylim([-0.5 0]); grid on;
    if plotNormalizeDiff%normalized version
        dataToplot = (deltaPerf(~increaseLogicalPerf,4) - deltaPerf(~increaseLogicalPerf,3))./abs(deltaPerf(~increaseLogicalPerf,3));
        ylabel('(unevenABC - uevenABC) / uevenABC');
    else %Non-normalized version of the difference 
        dataToplot = deltaPerf(~increaseLogicalPerf,4) - deltaPerf(~increaseLogicalPerf,3);
        ylabel('uevenABC - evenABC'); %FIXME: this ylabel will disappear after plotting.
    end
    plot(ones(size(dataToplot)),dataToplot,'o','Color',blueColor, 'LineWidth',2,'MarkerSize',10,'DisplayName','Dec(unexp)')
    xticks([]); grid on; xlabel('Dec(exp)'); 
    ylim([1.1*min(dataToplot) 0]); %ylim([-0.5 0]); 
    title('Amount of decrease');
%     legend(); %legend('location','northeastoutside'); 
    
    %plot PFC increaser vs decreaser
    fPFC = figure('units','normalized','outerposition',[0 0 1 1]);
    %plot unscaled PFC.
    % {{'StandAndAlphabet2RC';'StandAndAlphabet3RC';'WalkRestCorrected'; 'WalkAndAlphabet2RC';'WalkAndAlphabet3RC'},...
    subplot(totalRow,totalSubPltCol,1); axis square; hold on;
    increaseLogicalPFC = dataPfcHboVsBaseInOrderShifted(:,4) - dataPfcHboVsBaseInOrderShifted(:,3) >= 0;
    increaseRatioRawPFC = sum(increaseLogicalPFC) / size(dataPfcHboVsBaseInOrderShifted,1);

    plot(dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%decrease group
%     ylim([0 15]); 
    ylabel(['PFCActivation_{' titleStr '}']);
    title(sprintf('evenABC -> unevenABC\nDecrease (unexpected): %.2f', 1-increaseRatioRawPFC));
    xticks([1:2]); xticklabels({'evenABC','unevenABC'});

    s1 = subplot(totalRow,totalSubPltCol,totalSubPltCol+1); %plot on row 2
    axis square; hold on;
    plot(dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', blueColor);%increase or non change group
%     lineItems = findobj(s1,'Type','Line');
%     legend(lineItems([1]),'W2<=W2U (Expected)','location','northeastoutside')
    
    %plot each subject with its own color
%     for dIdx = 1:size(deltaPerf,1)
%         plot(1:2,dataPfcHboVsBaseInOrderShifted(dIdx,3:4),'o-','LineWidth',2.5,'MarkerSize',5,'DisplayName',subjectID{dIdx,1});
%     end
    %plot the line indicating the amount of shift.
%     plot(1:2,repmat(shiftAmount,1,2),'k-','LineWidth',1.5,'DisplayName',['Shifted ' num2str(shiftAmount)])
%         xlim([0.75,2.25])
    xticks([1:2]); xticklabels({'evenABC','unevenABC'});
%     ylim([0 15])
    ylabel(['PFCActivation_{' titleStr '}']);

    title(sprintf('evenABC -> unevenABC\nIncrease (expected): %.2f', increaseRatioRawPFC));
    outlierCounts = sum(isoutlier(dataPfcHboVsBaseInOrderShifted(:,3:4)));
%     legend('NumColumns',1,'Location','northeastoutside')
    
    %plot CI and unevenABC - evenABC scatter for unscaledPFC
    disp([titleStr 'Raw W2U-W2'])
%     subplot(totalRow,totalSubPltCol*2,totalSubPltCol*4+1); %plot on row 3
%     axis square; hold on;
%plot full CI of the whole group.
%     PlotHelper.plotCI(1, dataPfcHboVsBaseInOrderShifted(:,4) - dataPfcHboVsBaseInOrderShifted(:,3),'k','Full', true);
%     xlim([0.8,1.2]); xticks([]); %ylim([-30 20])
%     title(sprintf('Raw W2U-W2 Inc: %.2f', increaseRatioRawPFC)); %legend(); 
    %2nd column, row 1 unexpected
    subplot(totalRow,totalSubPltCol,2); axis square; hold on;%plot on row 3
%plot ci 
%     PlotHelper.plotCI(1, dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,4) - dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,3),greyColor,'Dec(unexp)', false); %gray (unexpected group)
%     xlim([0.8 1.2]); xticks([]); legend(); grid on; ylim([-4.5 0]); 
    %%plot the difference between W2U-WU.
    if plotNormalizeDiff%normalized version
        dataToplot = (dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,4) - dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,3))./abs(dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,3)); %normalized
        ylabel('(unevenABC - uevenABC) / uevenABC');
    else
        dataToplot = dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,4) - dataPfcHboVsBaseInOrderShifted(~increaseLogicalPFC,3);
        ylabel('uevenABC - evenABC');
    end
    plot(ones(size(dataToplot)),dataToplot,'o','Color',greyColor, 'LineWidth',2,'MarkerSize',10,'DisplayName','Dec(unexp)');
    xticks([]); grid on; xlabel('Dec(unexp)'); ylim([1.1*min(dataToplot) 0]); %ylim([-4.5 0]);
    title('Amount of decrease');
%     legend(); %legend('location','northeastoutside'); 
    %2nd column, 2nd row expected direction.
    subplot(totalRow,totalSubPltCol,totalSubPltCol+2); axis square; hold on;%plot on row 3
%     PlotHelper.plotCI(1, dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,4) - dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,3),blueColor,'Inc(exp)', false); %blue
%     xlim([0.8 1.2]); xticks([]); legend(); grid on; ylim([0 4.5]); 
    %%plot the difference between W2U-WU.
    if plotNormalizeDiff%normalized version
        dataToplot = (dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,4) - dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,3))./abs(dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,3));
        ylabel('(unevenABC - uevenABC) / uevenABC');
    else
        dataToplot = dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,4) - dataPfcHboVsBaseInOrderShifted(increaseLogicalPFC,3);
        ylabel('uevenABC - evenABC');
    end
    plot(ones(size(dataToplot)),dataToplot,'o','Color',blueColor, 'LineWidth',2,'MarkerSize',10,'DisplayName','Inc(exp)')
    xticks([]); grid on; xlabel('Inc(exp)'); ylim([0 1.1*max(dataToplot) ]); %ylim([0 4.5]); 
%     legend(); %legend('location','northeastoutside'); 
    title('Amount of increase');
    
    % determine optimal scaleFactor
    %objective is to maximize the %increase from W2 to W2Uneven and
    %minimize the % of outliers in both tasks (%outlier calculated by
    %outlier#/totalSubNum pertask and add up the % for both tasks;
    %outlier is defined as value outside of 3MAD of mean)
    %MAD = median absolute deviation, a robust statistic that gives the average distance of the data points from the median.
    objectiveFcnVal = [];%nan(1,length(1:0.1:5)+1); %size matching scale factors to try + 1 (no scale)
    scaleRange = [0:0.1:5];%[0,1:0.1:5]; %1st entry is no scale
    for scaleFactor = scaleRange
        if strcmp(scaleMethod, 'Exp')
%                 scale = 1/2 + exp(-scaleFactor*deltaPerf)/2; %exponential scaling
            scale = exp(-scaleFactor*deltaPerf); %exponential without division
        elseif strcmp(scaleMethod, 'Linear')
            scale = 1-scaleFactor*deltaPerf;
        end
        pfcHboVsBaseScaled = scale .* dataPfcHboVsBaseInOrderShifted;
        increaseRatio = sum((pfcHboVsBaseScaled(:,4) - pfcHboVsBaseScaled(:,3)) >=0) / size(pfcHboVsBaseScaled,1);
%             outlierCounts = sum(isoutlier(pfcHboVsBaseScaled(:,3:4)));
        objectiveFcnVal(end+1) = increaseRatio-sum(sum(isoutlier(pfcHboVsBaseScaled(:,3:4),'mean')))/(2*size(pfcHboVsBaseScaled,1));
    end
    fObjectiveFcn = figure();hold on;
    plot(scaleRange, objectiveFcnVal,'LineWidth',3);
    [bestObjVal, bestIdx] = max(objectiveFcnVal)
    xline(scaleRange(bestIdx),'k--','LineWidth',3);
    legend('ObjectiveFcn',['BestScale: ' num2str(scaleRange(bestIdx))],'Location','southoutside')
    xlabel('\alpha');
    ylabel('ObjectiveFcnValue')
    title([titleStr 'Objective Fcn'])
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    if saveResAndFigure
        saveas(fObjectiveFcn,[saveDir titleStr 'VsRest_ScaleByGaitAndCognitive_ScaleFactorObjFcn'])
        set(gcf,'renderer','painters')
        saveas(fObjectiveFcn,[saveDir titleStr 'VsRest_ScaleByGaitAndCognitive_ScaleFactorObjFcn.png'])
    end

    % plot the optimal scaleFactor
    figure(fPFC)
    scaleIdx = 1;
    scaleFactors = scaleRange(bestIdx); %could change this to show multiple example scenarios.
%         pfcVsBaseScaledPerFactor = nan([size(scaleFactors,2),size(dataPfcHboVsBaseInOrderShifted)]);
    if exampleFactorAroundOptimalToPlot >=1
        scaleFactors = [scaleFactors - (exampleFactorAroundOptimalToPlot/2:-0.5:0.5) scaleFactors scaleFactors + (0.5:0.5:exampleFactorAroundOptimalToPlot/2)]
    end
    %to plot examples of our choise.
%     figure(); hold on;; totalSubPltCol = 5; scaleFactors = 1:5; scaleIdx = -1;
    for scaleFactor = scaleFactors
        if strcmp(scaleMethod, 'Exp')
%             scale = 1/2 + exp(-scaleFactor*deltaPerf)/2; %exponential scaling
            scale = exp(-scaleFactor*deltaPerf); %exponential without division
        elseif strcmp(scaleMethod, 'Linear')
            scale = 1-scaleFactor*deltaPerf;
        end
        pfcHboVsBaseScaled = scale .* dataPfcHboVsBaseInOrderShifted;
        increaseLogical = pfcHboVsBaseScaled(:,4) - pfcHboVsBaseScaled(:,3) >= 0;
        increaseRatioScaled = sum(increaseLogical) / size(pfcHboVsBaseScaled,1);
%         
%         %% plot increase vs decrease group that tracks across perf, PFC, and scaled data.
%         %8 Colors are hard to track. 
%         fAll = figure('units','normalized','outerposition',[0 0 1 1]);
%         colorOrder = colororder('default');
%         colorOrder(end+1,:) = greyColor;
%         colorOrder(6,:) = [0 0 0];
%         logicalExpected = [~increaseLogicalPerf, increaseLogicalPFC]; %true is expected direction
%         logicalGroupOptions = truth_table(2);
%         dataAll = {deltaPerf, dataPfcHboVsBaseInOrderShifted,pfcHboVsBaseScaled};
%         for row = 1:2 %exp, then unexp
%             for col = 1:2
%                 subplot(totalRow,totalSubPltCol*2,totalSubPltCol*2*(row-1)+col*2-1:totalSubPltCol*2*(row-1)+col*2); hold on; %perf, (expected)
%                 perfRow = find(logicalGroupOptions(:,col));
%                 if row == 2
%                     perfRow = find(~logicalGroupOptions(:,col));
%                 end
%                 for r = perfRow'
%                     subIndex = all(logicalExpected == logicalGroupOptions(r,:),2);
%                     if any(subIndex)
%                         plot(2:3,dataAll{col}(subIndex,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', colorOrder(r,:));%gray (unexpected)
%                     end
%                 end
%             end
%         end
%         subplot(3,6,1:2); title('\DeltaPerformance Decreaser (Expected)')
%         subplot(3,6,3:4); title('Hbo Increaser (Expected)')
%         subplot(3,6,5:6); title('ScaledHbo Increaser (Expected)')
%         subplot(3,6,7:8); title('\DeltaPerformance Increaser (Unexpected)')
%         subplot(3,6,9:10); title('Hbo Decreaser (Unexpected)')
%         subplot(3,6,11:12); title('ScaledHbo Decreaser (Unexpected)')
        
%         %% plot PFC vs Performance and track differences
%         fPFCPerf = figure('units','normalized','outerposition',[0 0 1 1]);
%         s1=subplot(totalRow,totalSubPltCol*2,1:2);
%         axis square;
%         hold on;
% 
%         plot(2:3,deltaPerf(increaseLogicalPerf,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%gray (unexpected)
%         ylim([-0.65 0.67])
%         s1=subplot(totalRow,totalSubPltCol*2,totalSubPltCol*2+[1:2]); %plot on 2nd row.
%         axis square;
%         hold on;
%         plot(2:3,deltaPerf(~increaseLogicalPerf,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', blueColor);%blue (expected, non-change or decrease)
%         ylim([-0.65 0.67])
%         xticks([1:size(deltaPerf,2)])
%         xticklabels({'Baseline','W2','W2Uneven'})
%         ylabel('(Baseline - Task)/Baseline (Combined)')
%         increaseRatioPerf = sum(increaseLogicalPerf) / size(deltaPerf,1);
%         title(sprintf('Perf W2toW2U Dec: %.4f', 1-increaseRatioPerf));
%         
%         subplot(totalRow,totalSubPltCol*2,3:4); axis square; hold on;
%         plot(dataPfcHboVsBaseInOrderShifted((~increaseLogicalPFC)&(increaseLogicalPerf),3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%people whoc come from grey performance group (increased performance, unexpected)
%         plot(dataPfcHboVsBaseInOrderShifted((~increaseLogicalPFC)&(~increaseLogicalPerf),3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', blueColor);%people who come from blue perf group (decrease perf) %decrease group
%         ylim([0 15])
%         s1 = subplot(totalRow,totalSubPltCol*2,totalSubPltCol*2+[3:4]); %plot on row 2
%         axis square; hold on;
%         plot(dataPfcHboVsBaseInOrderShifted((increaseLogicalPFC)&(~increaseLogicalPerf),3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', blueColor);%people who come from blue perf group (decrease perf, expected)
%         plot(dataPfcHboVsBaseInOrderShifted((increaseLogicalPFC)&(increaseLogicalPerf),3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%people whoc come from grey performance group (increased performance, unexpected)
%         lineItems = findobj(s1,'Type','Line');
%         legend(lineItems([1]),'W2<=W2U (Expected)');%,'location','northeastoutside')
%         ylim([0 15])
%         ylabel(titleStr); title('Raw W2toW2U Dec')
% 
%         plot(1:2,repmat(shiftAmount,1,2),'k-','LineWidth',1.5,'DisplayName',['Shifted ' num2str(shiftAmount)])
%         ylabel([titleStr ' vs Rest (Shifted to be above 0)'])
%         xticks([1:2])
%         xticklabels({'W2','W2Uneven'})
%         title(sprintf('Raw W2toW2U Inc: %.2f', increaseRatioRawPFC));
% %         legend('NumColumns',1,'Location','northeastoutside')

        %plot increase vs decrease group scaled values with changes from
        %unscaled PFC, then with changes from delta perf
        for j = 1:2 %1=plot on fPFC, 2 = plot on fPerf.
            if j == 1
                figure(fPFC);
                refLogical = increaseLogicalPFC;
                displaySuffix = 'PFC';
            else
                figure(fPerf);
                refLogical = ~increaseLogicalPerf; %ref logic is expected performance (perf should dec)
                displaySuffix = 'Perf';
            end
            %3rd column or 5,7, ... if more examples are plotted, 1st row,
            %decrease index (unexpected)
            subplot(totalRow,totalSubPltCol,2+scaleIdx*2-1); axis square; hold on; 
            plot(pfcHboVsBaseScaled((~increaseLogical)&(~refLogical),3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%decrease in pfc then at index (grey)
            plot(pfcHboVsBaseScaled((~increaseLogical)&(refLogical),3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', blueColor);%increase in pfc (blue), now decrease at index
            ylim([0 max(pfcHboVsBaseScaled,[],'all')+1]);
            ylabel('Automaticity Index'); xticks([1 2]);xticklabels({'evenABC','unevenABC'});          
            title(sprintf('evenABC -> unevenABC ScaledX%.1f\nDecrease (unexpected): %.2f', scaleFactor, 1-increaseRatioScaled));
            
            %3rd (or other odd columns) 2nd row (increase index, expected)
            subplot(totalRow,totalSubPltCol,totalSubPltCol+2+scaleIdx*2-1); axis square; hold on;      
            plot(pfcHboVsBaseScaled(increaseLogical & refLogical,3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', blueColor);%>= in PFC (blue), then increase or non change group in index (blue)
            plot(pfcHboVsBaseScaled(increaseLogical & (~refLogical),3:4)','o-','LineWidth',1.5,'MarkerSize',5,'Color', greyColor);%< in PFC (grey), now increase or non change group (blue)
            ylim([0 max(pfcHboVsBaseScaled,[],'all')+1])

            %Plot individual colors
            %         for dIdx = 1:size(deltaPerf,1)
    %             plot(1:2,pfcHboVsBaseScaled(dIdx,3:4),'o-','LineWidth',2.5,'MarkerSize',5,'DisplayName',subjectID{dIdx,1});
    %         end
            ylabel('Automaticity Index')
            xticks([1:2]);
    %         xlim([1 2])
            xticklabels({'evenABC','unevenABC'});
            title(sprintf('evenABC -> unevenABC ScaledX%.1f\nIncrease (expected): %.2f', scaleFactor, increaseRatioScaled));

            sgtitle([titleStr 'VsRest ComparedTo ' displaySuffix])
              
            %plot CI of differences
    %         subplot(2,totalSubPltCol,2+scaleIdx); axis square; hold on;
    %         subplotTitleStr = sprintf('ScaledX%.1f W2U-WU Inc: %.2f', scaleFactor, increaseRatioScaled);
    %         disp([titleStr subplotTitleStr])
    %         PlotHelper.plotCI(1, pfcHboVsBaseScaled(:,4) - pfcHboVsBaseScaled(:,3),'k','Full', true);
    %         xlim([0.8,1.2]); xticks([]); 
    %         title(subplotTitleStr); %legend();
            %odd # columns starting from 4th, 1st row, increaser diff.
            subplot(totalRow,totalSubPltCol,2+scaleIdx*2); %row 3
            axis square; hold on;
            %plot CI of increaes and decrease;
%             PlotHelper.plotCI(1, pfcHboVsBaseScaled(~increaseLogical,4) - pfcHboVsBaseScaled(~increaseLogical,3),greyColor,'Dec(unexp)', false); %gray group
%             xlim([0.8 1.2]); xticks([]); grid on; ylim([-5.5 0]); yticks([-4 -2 0]); grid minor;
            %plot the difference (W2U - WU) as scatter points only
            if plotNormalizeDiff%normalized version
                indexDiffVal = (pfcHboVsBaseScaled(:,4) - pfcHboVsBaseScaled(:,3))./pfcHboVsBaseScaled(:,3);
                ylabel('(unevenABC - uevenABC) / uevenABC');
            else
                indexDiffVal = (pfcHboVsBaseScaled(:,4) - pfcHboVsBaseScaled(:,3));
                ylabel('uevenABC - evenABC');       
            end
            dataToplot = indexDiffVal((~increaseLogical)&(~refLogical)); %decreaser with unexpected perf (inc) or pfc (decrease) (grey)
            plot(ones(size(dataToplot)),dataToplot,'o','Color',greyColor, 'LineWidth',2,'MarkerSize',10,'DisplayName',['Unexpected ' displaySuffix]); 
            dataToplot = indexDiffVal((~increaseLogical)&refLogical); %decreaser with expected perf (decrease) or pfc (increase) (blue)
            plot(ones(size(dataToplot))+0.5,dataToplot,'o','Color',blueColor, 'LineWidth',2,'MarkerSize',10,'DisplayName',['Expected '  displaySuffix]); 
            xlim([0.5,2]);grid on; xticks({}); %xticks([1 1.5]); xticklabels({['Unexpected' displaySuffix],['Expected' displaySuffix]});
            ylim([min(indexDiffVal,[],'all')*1.1 0]);legend(); %legend('location','northeastoutside'); 
            title('Amount of decrease');
            
            %odd #columns (start from 4th) 2nd row, decreaser diff
            subplot(totalRow,totalSubPltCol,totalSubPltCol+2+scaleIdx*2); axis square; hold on;
            %plot CI of increase and decrease
%             PlotHelper.plotCI(1, pfcHboVsBaseScaled(increaseLogical,4) - pfcHboVsBaseScaled(increaseLogical,3),blueColor,'Inc(exp)', false);
%             xlim([0.8 1.2]); xticks([]); grid on; ylim([0 5.5]);yticks([0 2 4]);grid minor;%legend(); 
            %plot the data points only
            dataToplot = indexDiffVal(increaseLogical&(~refLogical)); %increaser with unexpected perf (increase) or pfc (decrease) (grey)
            plot(ones(size(dataToplot)),dataToplot,'o','Color',greyColor, 'LineWidth',2,'MarkerSize',10,'DisplayName',['Unexpected '  displaySuffix]); 
            dataToplot = indexDiffVal(increaseLogical&refLogical); %increaser with expected perf (dec) or pfc (increase) (blue)
            plot(ones(size(dataToplot))+0.5,dataToplot,'o','Color',blueColor, 'LineWidth',2,'MarkerSize',10,'DisplayName',['Expected '  displaySuffix]); 
            xlim([0.5,2]); grid on; xticks({});%xticks([1 1.5]); xticklabels({['Unexpected' displaySuffix],['Expected' displaySuffix]});
            ylim([0 max(indexDiffVal,[],'all')*1.1]); legend(); %legend('location','northeastoutside'); 
            title('Amount of increase');
        end
        
        outlierCounts = sum(isoutlier(pfcHboVsBaseScaled(:,3:4)))
        objectiveFcn = increaseRatioScaled-max(outlierCounts)/size(dataPfcHboVsBaseInOrderShifted,1)
%             pfcVsBaseScaledPerFactor(scaleIdx, :, :) = pfcHboVsBaseScaled;           
        scaledData{1,(scaleIdx - 1)*2 + i} = [titleStr ' vs Rest Scaled by Gait and Cog x' num2str(scaleFactor)];
        scaledData{2,(scaleIdx - 1)*2 + i} = pfcHboVsBaseScaled;
        %save data  column order (i.e., order of the task)
        scaledData{3,(scaleIdx - 1)*2 + i} = {'Wk','S2','W2','W2Uneven'};
        scaleIdx = scaleIdx + 1;
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    
    fCounts = figure('units','normalized','outerposition',[0 0 1 1]);
    %subplot(1,totalSubPltCol,1); 
    hold on; 
%     axis square;
%     increaseRatioPerf = sum(increaseLogical) / size(deltaPerf,1) * 100;
% %     bar(1,increaseRatio,1,'FaceColor',greyColor,'DisplayName','Unexpected (Perf Increase)');
%     bar(2,-100+increaseRatio,1,'FaceColor',blueColor,'DisplayName','Expected (Perf Decrease)');
    b = bar([increaseRatioPerf*100 -100+increaseRatioPerf*100; increaseRatioRawPFC*100 -100+increaseRatioRawPFC*100; increaseRatioScaled*100 -100+increaseRatioScaled*100])
    b(1).FaceColor = blueColor;
    b(2).FaceColor = [0.7 0.7 0.7];
    legend({'Increase','Decrease'})
    ylim([-100,100]);
    xticks([1,2,3]);
    xticklabels({'\DeltaPerformance','RawPFC',['Scaledx' num2str(scaleFactor,'%.1f')]});
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    if saveResAndFigure
        figure(fCounts);
        set(gcf,'renderer','painters')
        saveas(fCounts,[saveDir titleStr 'VsRest_IncDecSubjCount.fig'])
%         s = findobj('type','legend'); delete(s)
        saveas(fCounts,[saveDir titleStr 'VsRest_IncDecSubjCount.png'])
        
        figure(fPerf);
        set(findall(gcf,'-property','FontSize'),'FontSize',20); set(gcf,'renderer','painters'); 
        saveas(fPerf,[saveDir titleStr 'VsRest_Scaled' num2str(scaleFactor,'%.1f') 'ByGaitAndCognitive_VsPerf' normSaveSuffix '.fig'])
%         s = findobj('type','legend'); delete(s)
        saveas(fPerf,[saveDir titleStr 'VsRest_Scaled' num2str(scaleFactor,'%.1f') 'ByGaitAndCognitive_VsPerf' normSaveSuffix '.png'])
        
        figure(fPFC);
        set(findall(gcf,'-property','FontSize'),'FontSize',20); set(gcf,'renderer','painters');
        saveas(fPFC,[saveDir titleStr 'VsRest_Scaled' num2str(scaleFactor,'%.1f') 'ByGaitAndCognitive_VsPFC' normSaveSuffix '.fig'])
%         s = findobj('type','legend'); delete(s)
        saveas(fPFC,[saveDir titleStr 'VsRest_Scaled' num2str(scaleFactor,'%.1f') 'ByGaitAndCognitive_VsPFC' normSaveSuffix '.png'])
    end
        
    % plot ppl who didn't scale incrementally
    scaleIdx = 1;
    for scaleFactor = scaleFactors
        responders{i,scaleIdx} = []; %FIXME: this might not be working.
        increasingColor = [0.7 0.7 0.7]; %grey
        increaseLogical = pfcHboVsBaseScaled(:,4) - pfcHboVsBaseScaled(:,3) >= 0;
        cols = 5;
        f2 = figure('units','normalized','outerposition',[0 0 1 1]);
        for subGroup = 1:2
            if subGroup == 1
                subGroupIdx = find(increaseLogical);
                titleId = 'IncreasingDTOnly';
            else %non responder
                subGroupIdx = find(~increaseLogical);
                titleId = 'NonIncreasingDTOnly';
                %previous version where hbo and hbr are stored in 1 cell
                %array.
%                 removeNonResponder{i,scaleIdx}{(i-1)*2+1} = subGroupIdx;
%                 removeNonResponder{i,scaleIdx}{(i-1)*2+2} = subGroupIdx; %repeat for W2U
                responders{i,scaleIdx} = increaseLogical;
%                 responders{i,scaleIdx}{2} = increaseLogical; %repeat for W2U
                strToPrint = [titleStr ' Non Responder x' num2str(scaleFactor)];
                for j = 1:numel(subGroupIdx)
                    strToPrint = [strToPrint '  ' subjectID{subGroupIdx(j)} ' '];
                end
            end
            titleId = 'IncreaseVsNonIncrease';
            
            s1 = subplot(3,4,1);
            hold on;
            if subGroup == 1
                plot(deltaPerf(subGroupIdx,3:4)','o-','LineWidth',2.5,'MarkerSize',5,'Color', increasingColor);%,'DisplayName',subjectID{dIdx,1});
            else
                plot(deltaPerf(subGroupIdx,3:4)','o-','LineWidth',2.5,'MarkerSize',5);%,'DisplayName',subjectID{dIdx,1});
                lineItems = findobj(s1,'Type','Line');
                legend(lineItems([1,end]),'NonIncreasing','Increasing')
            end
            xticks([1:size(deltaPerf,2)])
%             xticklabels({'Wk','S2','W2','W2U'})
            xticklabels({'W2','W2Uneven'});
            ylabel('Combined \delta Performance')
            title('Combined \delta Performance')

            subplot(3,4,2)
            hold on;
            if subGroup == 1
                plot(dataPfcHboVsBaseInOrderShifted(subGroupIdx,3:4)','o-','LineWidth',2.5,'MarkerSize',5,'Color', increasingColor);%,'DisplayName',subjectID{dIdx,1});
            else
                plot(dataPfcHboVsBaseInOrderShifted(subGroupIdx,3:4)','o-','LineWidth',2.5,'MarkerSize',5);%,'DisplayName',subjectID{dIdx,1});
            end
            plot([1:2],repmat(shiftAmount,1,2),'k-','LineWidth',1.5,'DisplayName',['Shifted ' num2str(shiftAmount)])
%             xlim([0,size(deltaPerf,2)+1])
            ylabel([titleStr ' vs Rest (Shifted>0=)'])
            xticks([1:size(deltaPerf,2)])
%             xticklabels({'Walk','S2','W2','W2U'})
            xticklabels({'W2','W2Uneven'});
            title('Raw PFC');

            subplot(3,4,3:4); hold on;
            if subGroup == 1
                plot(pfcHboVsBaseScaled(subGroupIdx,3:4)','o-','LineWidth',2.5,'MarkerSize',5,'Color',increasingColor);
            else
                plot(pfcHboVsBaseScaled(subGroupIdx,3:4)','o-','LineWidth',2.5,'MarkerSize',5);
            end
            
            ylabel(['Scaled ',titleStr,' vs Rest'])
            xticks([1:size(deltaPerf,2)])
%             xticklabels({'Walk','S2','W2','W2Uneven'})
            xticklabels({'W2','W2Uneven'});
        
            subplot(3,4,5:6); hold on; %gait only
            if subGroup == 1
                plot(perfDataArray(subGroupIdx,[1,3,5])','o-','LineWidth',2.5,'MarkerSize',5,'Color',increasingColor);
                decreaseRatio_gait1 = sum(perfDataArray(subGroupIdx,3) >= perfDataArray(subGroupIdx,[5])) / length(subGroupIdx);
            else
                plot(perfDataArray(subGroupIdx,[1,3,5])','o-','LineWidth',2.5,'MarkerSize',5);
                decreaseRatio_gait2 = sum(perfDataArray(subGroupIdx,3) >= perfDataArray(subGroupIdx,[5])) / length(subGroupIdx);
            end
            title(sprintf('Gait Performance.'));
            xticks([1:3]);
            xticklabels({'Walk','W2','W2Uneven'});
            %should decrease or unchange
%             xticklabels({'W2','W2Uneven'});
            ylabel('Gait Speed (m/s)')

            subplot(3,4,7:8); hold on; %cognitive only
            if subGroup == 1
                plot(perfDataArray(subGroupIdx,[2,4,6])','o-','LineWidth',2.5,'MarkerSize',5,'Color',increasingColor);
                decreaseRatio_alpha1 = sum(perfDataArray(subGroupIdx,4) >= perfDataArray(subGroupIdx,[6])) / length(subGroupIdx);
            else
                plot(perfDataArray(subGroupIdx,[2,4,6])','o-','LineWidth',2.5,'MarkerSize',5);
                decreaseRatio_alpha2 = sum(perfDataArray(subGroupIdx,4) >= perfDataArray(subGroupIdx,[6])) / length(subGroupIdx);
            end
            title(sprintf('Alphabet Performance'));
            xticks([1:3]);
            xticklabels({'S2','W2','W2Uneven'});
%             xticklabels({'W2','W2Uneven'});
            ylabel('Correct Letter / s)');   
                        
%             subplot(3,4,10); %gait speed only
%             plot((perfDataArray(subGroupIdx,[1,3,5])./perfDataArray(subGroupIdx,1))','o-','LineWidth',2.5,'MarkerSize',5);
%             title('Gait Speed as % of Base');
%             xticks([1:3]);
%             xticklabels({'Wk','W2','W2Uneven'});
% %             xticklabels({'W2','W2Uneven'});
%             ylabel('Gait%');
            subplot(3,4,9); 
            hold on;
            diffData = perfDataArray(subGroupIdx,5) - perfDataArray(subGroupIdx,3);
            bar(subGroup,mean(diffData));
%             text(subGroup-0.5,mean(diffData)+3*std(diffData),sprintf('<0:%.2f',decreaseRatio_gait))
            if subGroup == 1 %skip legend
                errorbar(subGroup,mean(diffData),std(diffData),'k.','LineWidth',2,'HandleVisibility','off')
            else
                errorbar(subGroup,mean(diffData),std(diffData),'k.','LineWidth',2)
            end
            scatter(subGroup*ones(size(diffData))-0.2,diffData,65,'k','filled','HandleVisibility','off')
            xticks([1:2]);
            xticklabels({'Inc','NonInc'});
            title('W2U-W2');
            if subGroup == 2
                legend({sprintf('Inc SlowDown: %.2f',decreaseRatio_gait1),sprintf('NonInc SlowDown: %.2f',decreaseRatio_gait2),'SE'},'NumColumns',2)
            end
            
            subplot(3,4,10); hold on;%gait speed difference in %
            if subGroup == 1
                plot(((-perfDataArray(subGroupIdx,1) + perfDataArray(subGroupIdx,[1,3,5]))./perfDataArray(subGroupIdx,1))','o-','LineWidth',2.5,'MarkerSize',5,'Color',increasingColor);
            else
                plot(((-perfDataArray(subGroupIdx,1) + perfDataArray(subGroupIdx,[1,3,5]))./perfDataArray(subGroupIdx,1))','o-','LineWidth',2.5,'MarkerSize',5);
            end
            title('(DT - GaitBase) / Base');
            xticks([1:3]);
            xticklabels({'Wk','W2','W2Uneven'});
%             xticklabels({'W2','W2Uneven'});
            ylabel('\delta Gait%');
                            
            %bar graph showing difference between W2Uneven-W2 (most should
            %be negative - slow down)
%             subplot(3,4,12); %cognitive only
%             plot((perfDataArray(subGroupIdx,[2,4,6])./perfDataArray(subGroupIdx,2))','o-','LineWidth',2.5,'MarkerSize',5);
%             title('Alphabet as % of Base');
%             xticks([1:3]);
%             xticklabels({'S2','W2','W2Uneven'});
% %             xticklabels({'W2','W2Uneven'});
%             ylabel('Alpha%');
            subplot(3,4,11); 
            hold on;
            diffData = perfDataArray(subGroupIdx,6) - perfDataArray(subGroupIdx,4);
            bar(subGroup,mean(diffData));
%             text(subGroup-0.5,mean(diffData)+3*std(diffData),sprintf('<0:%.2f',decreaseRatio_alpha))
            if subGroup == 1 %skip legend
                errorbar(subGroup,mean(diffData),std(diffData),'k','LineWidth',2,'HandleVisibility','off')
            else
                errorbar(subGroup,mean(diffData),std(diffData),'k','LineWidth',2)
            end
            scatter(subGroup*ones(size(diffData))-0.2,diffData,65,'k','filled','HandleVisibility','off')
            xticks([1:2]);
            xticklabels({'Inc','NonInc'});
            title('W2U-W2');
            if subGroup == 2
                legend({sprintf('Inc LessAlpha: %.2f',decreaseRatio_alpha1),sprintf('NonInc LessAlpha: %.2f',decreaseRatio_alpha2),'SE'},'NumColumns',2)
            end
            
            subplot(3,4,12); hold on; %cognitive difference in %
            if subGroup == 1
                plot(((-perfDataArray(subGroupIdx,2) + perfDataArray(subGroupIdx,[2,4,6]))./perfDataArray(subGroupIdx,2))','o-','LineWidth',2.5,'MarkerSize',5,'Color',increasingColor,'HandleVisibility','Off');
            else
                plot(((-perfDataArray(subGroupIdx,2) + perfDataArray(subGroupIdx,[2,4,6]))./perfDataArray(subGroupIdx,2))','o-','LineWidth',2.5,'MarkerSize',5);
            end
            title('(DT - AlphaBase) / Base');
            xticks([1:3]);
            xticklabels({'S2','W2','W2Uneven'});
%             xticklabels({'W2','W2Uneven'});
            ylabel('\delta Alpha%');
            
            if subGroup == 2
                legend(subjectID(subGroupIdx,:))
            end
            sgtitle([titleId ' ' titleStr ' Subjects x' num2str(scaleFactor)])
            set(findall(gcf,'-property','FontSize'),'FontSize',20)
        end
        if saveResAndFigure
            saveas(f2,[saveDir titleStr 'VsRest_Scaled' num2str(scaleFactor,'%.1f') 'ByGaitAndCognitive_' titleId '.fig'])
            set(gcf,'renderer','painters')
            saveas(f2,[saveDir titleStr 'VsRest_Scaled' num2str(scaleFactor,'%.1f') 'ByGaitAndCognitive_' titleId '.png'])
        end
    end
end



end