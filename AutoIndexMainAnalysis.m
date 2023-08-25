%%
close all; clc; clear all;
scaleMethod = 'Exp'; %Linear or Exp
logData = false; %true if shift unscaled to positive and log y-axis for all; false o.w.
studyID = {'NMCM','PRIMA','Move MYHAT'};%{'NMCM','PRIMA','Move MYHAT'}; %NMCM, PRIMA, MMH
scriptDir = fileparts(matlab.desktop.editor.getActiveFilename); 
scriptDir = strrep(scriptDir,'AutoIndexAnalysis','');
scriptDir = [scriptDir 'Data' filesep];

saveDir = [scriptDir 'IndexAppliedResults' filesep '012723Data' filesep];
if iscell(studyID)
    saveDir = [saveDir strjoin(studyID, '_') scaleMethod filesep]
else
    saveDir = [saveDir studyID scaleMethod filesep]
end
% if logData
%     saveCorrDir = [saveDir 'LogY_']
% else
%     saveCorrDir = [saveDir 'RegularY_']
% end
if not(isfolder(saveDir))
    mkdir(saveDir)
end
saveResAndFigure = false;
%% load perf data
% perfData = readcell('/Users/mac/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofPittsburgh/BEAM lab - fNIRS Workgroup/Rosso_fNIRS_code/Matlab/Emma_SampleDataNMCM_forShuqi/NMCM_PRIMA_MOVE_combined_wide_19Nov2021.csv');
% /Users/mac/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SML/Projects/fNIR Project/Code_NIRS_Automaticity/Data/EmmaDataset/PRIMA_051722
% load 19Nov2021 Data
% perfData = readcell([scriptDir 'EmmaDataset' filesep 'NMCM_PRIMA_MOVE_combined_wide_19Nov2021.csv']);
perfData = readcell([scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'Performance_GaitSpeed_Alphabet' filesep 'NMCM_PRIMA_MMH_Visit1_performance_subjavg.csv']);

perfDataHeader = perfData(1,:);
studyIdCol = find(strcmp(perfDataHeader,'Cohort'));
perfDataFull = perfData(contains(perfData(:,studyIdCol),studyID),:);

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
% for i = 1:length(subjectID)
%     subjectIDStr{end+1} = num2str(subjectID);
% end
% subjectID = subjectIDStr;
subjectID = cellstr(num2str(cell2mat(subjectID))); 
perfDataArray = cell2mat([walkData(:,gaitSpeedCol),s2Data(:,alphabetCol),walk2Data(:,[gaitSpeedCol,alphabetCol]),walkUneven2Data(:,[gaitSpeedCol,alphabetCol])]);
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
pfcDataFull = readcell([scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'fNIRS' filesep 'NMCM_PRIMA_MMH_Visit1_fNIRS_subjavg.csv']);
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
visitCol = find(strcmp(pfcDataHeader,'Visit'));
pfcData = pfcData(strcmp(pfcData(:,visitCol),'Visit 1'),:);
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
%walkuneven2, additional column for hbo and hbr
% %if load full data (version 12/01/2022)
% hboData = str2double(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbo'),2:end));
% hbrData = -str2double(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbr'),2:end)); %hbr should be oppositive of hbo (more negative means more activation)
% %if load specific study 05/16/2022 version or load full new data 01/27/23
% version
hboData = cell2mat(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbo'),2:end));
hbrData = -cell2mat(pfcDataTstats(strcmp(pfcDataTstats(:,1),'hbr'),2:end)); %hbr should be oppositive of hbo (more negative means more activation)

%% scale and plot other dataset
% clc;
normalizeDiff = false;
if normalizeDiff
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
    if normalizeDiff%normalized version
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
    if normalizeDiff%normalized version
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
    if normalizeDiff%normalized version
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
    if normalizeDiff%normalized version
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
            if normalizeDiff%normalized version
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

%% find and prepare relevant cognitive and age variables to correlate with =
%load cognitive data, version 01272023
cogData = readcell([scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'Demographics_Other' filesep 'NMCM_PRIMA_MMH_Visit1_demographic_other_080323.csv']);
cogDataHeader = cogData(1,:);
%include physical fatigubility score and life space assessment score (maybe
%more like physical measures)
% cogCol = find(ismember(cogDataHeader,{'MergeID','Cohort','age','sexF','eduyr','PFS','traila','trailb','MMSE_total','Overall_LSA'}));
cogCol = find(ismember(cogDataHeader,{'MergeID','Cohort','age','sexF','eduyr','traila','trailb','MMSE_total','race'}));
%replace my hat participant to have same name as fNIRS and performance
%data.
cogData(strcmp(cogData(:,2),'Move MyHAT'),2) = {'Move MYHAT'};
cogData = cogData(contains(cogData(:,2),studyID),:);
cogData = cogData(:,cogCol);
missingCogDataMask = strcmp(cogData(:,3),'NA'); %have age missing
cogData = cogData(~missingCogDataMask,:); %this should make the cognitive data and fnirs performance full data match.
cogData = cogData(~missingDataMask,:); %remove rows with missing performance of fnirs data.
checkID = isequal(str2double(subjectID), cell2mat(cogData(:,1))); %make sure every subject is present
% %data in previous version 051722
% cogCol = find(ismember(perfDataHeader,{'age','sex','eduyr','traila','trailb','MMSE_total'}));
% cogData = perfDataFull(:,[1:4,cogCol]); %should have 5 repeats per subject
% cogData = cogData(strcmp(cogData(:,3),'Even'),:); %take a random one, should have 29 entries, 1 persubject
% cogData = cogData(~missingDataMask,:);
% checkID = isequal(char(subjectID), num2str(cell2mat(cogData(:,1)))); %make sure every subject is present
% cogDataFull = cogData;
% cogData = cogData(:,5:end);
if ~checkID
    warning('Cognitive data have missing/mismatch subjects.')
end
%missing 21275, 21276 PRIMA trialB, and edu, trailAB, MMSE for 56 MMH
%subjects.
otherMissingVal = ~cellfun(@isnumeric,cogData(:,3:end-1)); %check missing entry of demographics. starting from column 3
CogDataMissingID = cogData(any(otherMissingVal,2),:)
cogData = cogData(:,3:end); %leave only the data columns.
cogData(otherMissingVal) = {NaN}; %missing data was 'NA' make them nan to support number operation

%handle race
race = categorical(cogData(:,end));
[raceCat,~,raceNumerical] = unique(race);
% cogData(:,end) = raceNumerical;
cogData = [cell2mat(cogData(:,1:end-1)), raceNumerical]; 
%in format: subject x ['age','sex','eduyr','traila','trailb','MMSE_total','race']


%% summarize the data (demographics and performance).
clc
% old data loading version 051722
% cogLabels = perfDataHeader(cogCol);
% New data label is in a separate file for cog/demographic data. 01272023
cogLabels = cogDataHeader(cogCol(3:end));

% merge education and race categories to avoid single observation in a
% category
eduColIdx = find(strcmp('eduyr',cogLabels)); %education is the 3rd column, value 2 means < high school
dataToMerge = cogData(:,eduColIdx)==2;
disp('Less than high school cognitive data row index:')
find(dataToMerge)
cogData(dataToMerge,eduColIdx) = 3; %reassign value 2 to 3, merge into: high school or less

raceColIdx = find(strcmp('race',cogLabels)); %race is the last column, 
dataToMerge = ismember(cogData(:,raceColIdx),[1,2,4]); %value 1,2,4 AI, AS(American indian or alaska native), 4=other/mix
disp('Unique race row index:')
find(dataToMerge)
cogData(dataToMerge,raceColIdx) = 3; %reassign 1,2,4 to 3 (black), merge into: white vs nonwhite 

groups = floor(cellfun(@str2num,subjectID)/10000); 

summaryStats = cell(27,2*4+1);
summaryStats(:,1)=[{'Study','n'},cogLabels(1:end-1),{'W2Hbo','W2UHbo','W2Hbr','W2UHbr',...
    'GaitSpeedWalk','GaitSpeedW2','GaitSpeedW2U','AlphaS2','AlphaW2','AlphaW2U',...
    'Education2','Education3','Education4','Education5','Race1','Race2','Race3','Race4','Race5'}]';
for study=1:4
    if study == 4
        studyMask = groups<=3;
        fprintf('\n\nAll\n')
        summaryStats{1,study*2}='All';
    else
        studyMask = groups==study;
        fprintf(['\n\n',studyID{study}, '\n'])
        summaryStats{1,study*2}=studyID{study};
    end
    fprintf('n=%d',sum(studyMask))
    summaryStats{2,study*2}=sum(studyMask);
    disp(cogLabels)
    fprintf('Mean+SD    ')
    %\x00B1 is +- sign.
%     fprintf('%.1f \x00B1 %.1f
%     ',nanmean(cogData(studyMask,:)),nanstd(cogData(studyMask,:)))
%     %incorrect formatting, will do mean1+-mean2, go through 2 list by
%     finish 1st list then 2nd
    fprintf('\nSex (female) %d (%.1f%%)',sum(cogData(studyMask,2)),100*sum(cogData(studyMask,2))/length(cogData(studyMask,2)))
    fprintf('\nEducation (by categories) \n')
    summary(categorical(cogData(studyMask,3)))
    summaryStats(3:8,study*2:study*2+1)=num2cell([nanmean(cogData(studyMask,1:end-1));nanstd(cogData(studyMask,1:end-1))])';
    summaryStats(4,study*2:study*2+1)=num2cell([sum(cogData(studyMask,2)),100*sum(cogData(studyMask,2))/length(cogData(studyMask,2))]);
    %append to end 19-22 rows for different levels of education counts and
    %percentage (3-high school/equivalent; 4-college; 5-post graduate; 2-less than high school or others)
    eduCats = cellfun(@str2num,categories(categorical(cogData(studyMask,3)))); %education categories counts
    eduByCat = countcats(categorical(cogData(studyMask,3))); %education categories counts
    summaryStats(17+eduCats,study*2:study*2+1)=num2cell([eduByCat, eduByCat./sum(eduByCat)*100]);
    %append to end 23-27 rows for different race counts and
    %percentage (1-AI: American Indian or Alaskan Native, 2-AS: Asian; 3-B: Black; 4- others, pt then specified 'Mixed', 5-W:White
    raceCats = cellfun(@str2num,categories(categorical(cogData(studyMask,end)))); %race categories counts
    raceByCat = countcats(categorical(cogData(studyMask,end))); %race categories counts
    summaryStats(22+raceCats,study*2:study*2+1)=num2cell([raceByCat, raceByCat./sum(raceByCat)*100]);
    fprintf('\nRace (by categories) \n')
    summary(categorical(cogData(studyMask,end)))
    
    fprintf('\n\nHbo   walkEvenHbo,   walkUnevenHbo\n   ')
    fprintf('%.2f + %.2f    ',nanmean(hboData(studyMask,3:4)),nanstd(hboData(studyMask,3:4)))
    summaryStats(9:10,study*2:study*2+1)=num2cell([nanmean(hboData(studyMask,3:4));nanstd(hboData(studyMask,3:4))])';
    
    fprintf('\n\nHbr  walkEvenHbr,   walkUnevenHbr\n    ')
    fprintf('%.2f + %.2f    ',nanmean(hbrData(studyMask,3:4)),nanstd(hbrData(studyMask,3:4)))
    summaryStats(11:12,study*2:study*2+1)=num2cell([nanmean(hbrData(studyMask,3:4));nanstd(hbrData(studyMask,3:4))])';

    fprintf('\n\nGaitSpeed    Walk,   WalkEven,   WalkUneven\n    ')
    fprintf('%.2f + %.2f    ',nanmean(perfDataArray(studyMask,[1,3,5])),nanstd(perfDataArray(studyMask,[1,3,5])))
    summaryStats(13:15,study*2:study*2+1)=num2cell([nanmean(perfDataArray(studyMask,[1,3,5]));nanstd(perfDataArray(studyMask,[1,3,5]))])';
    
    fprintf('\n\nAlphabetPerf    S2,   W2,   W2Uneven\n    ')
    fprintf('%.2f + %.2f    ',nanmean(perfDataArray(studyMask,[2,4,6])),nanstd(perfDataArray(studyMask,[2,4,6])))
    summaryStats(16:18,study*2:study*2+1)=num2cell([nanmean(perfDataArray(studyMask,[2,4,6]));nanstd(perfDataArray(studyMask,[2,4,6]))])';
end

summaryStats = [{'Variables','Mean','SD','Mean','SD','Mean','SD','Mean','SD'};summaryStats];
if saveResAndFigure
    summaryStats = cell2table(summaryStats(1:end,:));
    % Write the table to a CSV file
    writetable(summaryStats,[saveDir 'DataSummaryStats' strjoin(studyID) '.csv']);
end
%% set up data for correlation
%load only 2 DT conditions, prepare pfcDatAll as a cell with 1st entry
%unscaled, 2nd entry scaled in order of scaleFactors. 
condIdx = 3:4; %w2 and w2uneven
pfcDataAll = {hboData(:,condIdx), hbrData(:,condIdx)};%{deltaPerf(:,3:4),[hboData(:,condIdx), hbrData(:,condIdx)]};

if logData
    for i = 1:numel(pfcDataAll)
        %shift data up by the min +1
        pfcDataAll{i} = pfcDataAll{i} + abs(min(pfcDataAll{i}))+1;
        pfcDataAll{i} = log(pfcDataAll{i});
    end
end

saveStr = {'HboUnscaled','HbrUnscaled'};%{'Performance','Unscaled'};
pfcTaskLabels = {'W2','W2Uneven'};%{'W2 Hbo','W2Uneven Hbo','W2 Hbr','W2Uneven Hbr'};

% scaledDataToAdd = [1,2];
for toAdd = 1:length(scaleFactors)
    for hboType = 1:2 %hbo, then hbr
        pfcDataAll{end+1} = scaledData{2,(toAdd-1)*2+hboType}(:,condIdx);
    %     scaledData{2,(toAdd-1)*2+2}(:,condIdx)
        if logData
            pfcDataAll{end} = log(pfcDataAll{end});
        end
        %save string should include hbo/hbr, scaled by both or 1perf, scale
        %factor
        saveStr{end+1} = [scaledData{1,(toAdd-1)*2+hboType}(1:3) 'ScaledByBoth' scaledData{1,(toAdd-1)*2+hboType}(end-3:end)]
%         scaledData{1,(toAdd-1)*2+1}
    end
end

% find the max/visual outliers for scaled W2 and scaled W2Uneven
% [sortData,sortIdx] = sort(pfcDataAll{end}); %hbo, W2: top3: 10022,10020,10019; W2Uneven top: 10020
%from outlier test maybe remove:  {[10019]}    {[10020]}    {[10022]}    {[10027]}
% their index is 18,19,21,25

% append pfcdataall with the last cell as the responder subjects.
pfcDataResponder = {};
dataInCell = pfcDataAll;
for dataIdx = 3:numel(pfcDataAll) %skip the unscaled one.
    %separate the matrix of pfcdataall into a cell of 4 column vectors bc
    %different amount of responders will be removed per condition and
    %hbo/hbr condition.
    responderIdx = responders{dataIdx-2}; %1x2 cell arrays of 2 column logical vectors.
    if ~all(responderIdx) %nonresponders exist.
        saveStr{end+1} = [saveStr{dataIdx} 'Responder'];
        currData = pfcDataAll{dataIdx};
        pfcDataResponder{end+1} = currData(responderIdx,:);
    end
%     toRemove = responders{dataIdx-2};%index shift by 2 bc responders start from scaled data.
%     if ~isempty(toRemove)
%         pfcDataResponder{end+1} = {};
%         saveStr{end+1} = [saveStr{dataIdx} 'Responder'];
%         for j = 1:numel(toRemove)
%             currData = pfcDataAll{dataIdx}(:,j);
%             currData(toRemove{j})= [];
%             pfcDataResponder{end}{j} = currData;
%         end
%     end
end

% %convert the unscaled and scaled data from a matrix to a cell of 4 column
% vectors to be consistent with the responder data.
% for dataIdx = 1:numel(pfcDataAll) 
%     dataInCell{dataIdx} = mat2cell(pfcDataAll{dataIdx},[size(pfcDataAll{dataIdx},1)],ones([1,size(pfcDataAll{dataIdx},2)]));
% end
%append responders to the end of data, 1 responder group per scaledData
pfcDataAll = [dataInCell, pfcDataResponder];
%% test for normality
% clc;
% [18,19,21];
%not normal: for scaled with W2,
removeMeanOutlier={};
removeToNormal = {};
for dataIdx = 1:numel(pfcDataAll)
    removeMeanOutlier{end+1} = {[],[],[],[]};
    removeToNormal{end+1} = {[],[],[],[]};
    pfcDataCurr = pfcDataAll{dataIdx};
    for pfcIdx = 1:size(pfcDataCurr,2) %w2 and w3       
        [h,pnorm] = kstest(normalize(pfcDataCurr(:,pfcIdx)));
        fprintf('\nPFCData: %s %s, normal(0 = normal): %d, p = %.2f',saveStr{dataIdx}, pfcTaskLabels{pfcIdx}, h, pnorm)
        [h,pnorm] = adtest(pfcDataCurr(:,pfcIdx));
        fprintf('\nPFCData: %s %s, adtest normal(0 = normal): %d, p = %.2f',saveStr{dataIdx}, pfcTaskLabels{pfcIdx}, h, pnorm)

        tf = find(isoutlier(pfcDataCurr(:,pfcIdx),'mean'));
        if ~isempty(tf)
            removeMeanOutlier{dataIdx}{pfcIdx} = tf;
            removeToNormal{dataIdx}{pfcIdx} = tf;
            if ~logData && strcmp(saveStr{dataIdx},'ScaledByBothX4')
                if pfcIdx == 1 
                    tf = [19;21;18;67;57]; %remove to normal
                elseif pfcIdx == 2 
                    tf=[19;45;61;42;67];
                end
                removeToNormal{dataIdx}{pfcIdx} = tf;
            end
            nonOutlierData = pfcDataCurr(:,pfcIdx);
            nonOutlierData(tf) = [];
            [h,pnorm] = kstest(normalize(nonOutlierData));
            fprintf('\nAfterRemovingOutlier PFCData: %s %s, normal(0 = normal): %d, p = %.2f\n',saveStr{dataIdx}, pfcTaskLabels{pfcIdx}, h, pnorm)
            disp(['Outilers: ' num2cell(tf)' subjectID(tf)']);
        end
        
        %test normality within responders
        
    end
end
% outlier: x2 w2 hbo and hbr, index 21-10022; always normal
% x4 w2 hbo: 21 - 10022 (not normal even afte removal; normal after removal 21,19 - 10020,18 - 10019); 
% x4 w2uneven hbo: normal always, outlier: 19 - 10020
% x4 w2 hbr: 21 - 10022, normal after removal
% x4 w2uneven hbr: 19 - 10020, always normal;

%% normality test for cog idx
%all normal
removeMeanOutlierX = {};

for cogIdx = 1:size(cogData,2)
%     normalitytest(cogData(:,cogIdx)')
    [h,pnorm] = kstest(normalize(cogData(:,cogIdx)));
    fprintf('\nCogData: %s, normal(0 = normal): %d, p = %.2f',cogLabels{cogIdx}, h, pnorm)
    tf = find(isoutlier(cogData(:,cogIdx),'mean'));
    if ~isempty(tf)
        removeMeanOutlierX{cogIdx} = tf; %add in outliers
        nonOutlierData = cogData(:,cogIdx);
        nonOutlierData(tf) = [];
        [h,pnorm] = kstest(normalize(nonOutlierData));
        fprintf('\nAfterRemovingOutlier CogData: %s , normal(0 = normal): %d, p = %.2f',cogLabels{cogIdx}, h, pnorm)
        disp(['Outilers: ' subjectID(tf)']);
    end
end

%outliers using mean
% PFCData: HboScabledByBoth W2, normal(0 = normal): 1, p = 0.00
% AfterRemovingOutlier PFCData: HboScabledByBoth W2, normal(0 = normal): 1, p = 0.02
%     {'Outilers: '}    {[10022]}
% PFCData: HboScabledByBoth W2Uneven, normal(0 = normal): 0, p = 0.18
% AfterRemovingOutlier PFCData: HboScabledByBoth W2Uneven, normal(0 = normal): 0, p = 0.24
%     {'Outilers: '}    {[10020]}
% PFCData: HbrScaledByBoth W2, normal(0 = normal): 1, p = 0.01
% AfterRemovingOutlier PFCData: HbrScaledByBoth W2, normal(0 = normal): 0, p = 0.16
%     {'Outilers: '}    {[10022]}
% PFCData: HbrScaledByBoth W2Uneven, normal(0 = normal): 0, p = 0.41
% AfterRemovingOutlier PFCData: HbrScaledByBoth W2Uneven, normal(0 = normal): 0, p = 0.94
%     {'Outilers: '}    {[10020]}

% outlier found using median
% PFCData: HboUnscaled W2, normal(0 = normal): 0, p = 0.79
% AfterRemovingOutlier PFCData: HboUnscaled W2, normal(0 = normal): 0, p = 0.74
%     {'Outilers: '}    {[10002]}    {[10022]}    {[10029]}
% PFCData: HboUnscaled W2Uneven, normal(0 = normal): 0, p = 0.81
% PFCData: HbrUnscaled W2, normal(0 = normal): 0, p = 0.72
% PFCData: HbrUnscaled W2Uneven, normal(0 = normal): 0, p = 0.65
% AfterRemovingOutlier PFCData: HbrUnscaled W2Uneven, normal(0 = normal): 0, p = 0.93
%     {'Outilers: '}    {[10002]}    {[10011]}    {[10021]}
% PFCData: HboScabledByBoth W2, normal(0 = normal): 1, p = 0.00
% AfterRemovingOutlier PFCData: HboScabledByBoth W2, normal(0 = normal): 0, p = 0.31
%     {'Outilers: '}    {[10019]}    {[10020]}    {[10022]}    {[10027]}
% PFCData: HboScabledByBoth W2Uneven, normal(0 = normal): 0, p = 0.18
% AfterRemovingOutlier PFCData: HboScabledByBoth W2Uneven, normal(0 = normal): 0, p = 0.24
%     {'Outilers: '}    {[10020]}
% PFCData: HbrScaledByBoth W2, normal(0 = normal): 1, p = 0.01
% AfterRemovingOutlier PFCData: HbrScaledByBoth W2, normal(0 = normal): 0, p = 0.56
%     {'Outilers: '}    {[10019]}    {[10020]}    {[10022]}
% PFCData: HbrScaledByBoth W2Uneven, normal(0 = normal): 0, p = 0.41
% AfterRemovingOutlier PFCData: HbrScaledByBoth W2Uneven, normal(0 = normal): 0, p = 0.94
%     {'Outilers: '}    {[10020]}

%% cog var correlations
f = figure('units','normalized','outerposition',[0 0 1 1]);
cogData = [cogData,cogData(:,5) - cogData(:,4)];
cogLabels = [cogLabels,'TMTBMinusTMTA'];
[cogR, cogP] = corrplot(cogData, 'Type','Spearman','TestR',"on",'VarNames',cogLabels)
if saveResAndFigure
    saveas(f, [saveDir 'CogVar_CorrelationMatrix.fig'])
    saveas(f, [saveDir 'CogVar_CorrelationMatrix.png'])
end

%% anova to compare if the 3 studies data are significantly different
close all; clc
groups = floor(cellfun(@str2num,subjectID)/10000); %create 1,2,3,arrays to identify subject groups
fprintf('\nW2HboScaled\n')
[p,tbl,stats] = anova1(pfcDataAll{3}(:,1),groups)
fprintf('\nW2UHboScaled\n')
[p,tbl,stats] = anova1(pfcDataAll{3}(:,2),groups)
fprintf('\nW2HbrScaled\n')
[p,tbl,stats] = anova1(pfcDataAll{4}(:,1),groups)
fprintf('\nW2UHbrScaled\n')
[p,tbl,stats] = anova1(pfcDataAll{4}(:,2),groups)

%% check correlation among the demographic vars
f = figure();
cogDataTable = array2table(cogData);
cogDataTable.Properties.VariableNames = strrep(strrep(cogLabels,'_',''),'il',''); %the figure can only display up to 5 chars for labels, so make traila/b traa/b
[r, p] = corrplot(cogDataTable(:,[1,2,3,6]),'Type','Spearman','TestR','on')
if saveResAndFigure
    saveas(f, [saveDir 'Demographic_CorrelationMatrix.png'])
end

%% initialize a big cell array for model summaries.
mdlSummary = {'ResponseName','Predictors','SigRegressor','R2Ordinary','R2Adjusted','MMSEBeta','pValue','AIC','BIC','ResidualAbnormal(0Normal)','LackOfFit(1lack)','RobustOpt','NumSub','Model'};
%TODO: 
%4. add task as a regressor and task x mmse to see if predicts 1 task
%better than the other.

%% regression model of the index vs demographic variables.
clc; close all;
robustOption = 'off'; %default on use bisquare
%%run regression model of performance using PFC + demographic variables
lmTable = array2table([pfcDataAll{1},pfcDataAll{2},pfcDataAll{3},pfcDataAll{4},cogData,deltaPerf(:,3:4),responders{1},responders{2}]);
lmTable.Properties.VariableNames = [{'W2Hbo','W2UHbo','W2Hbr','W2UHbr','W2HboScaled','W2UHboScaled','W2HbrScaled','W2UHbrScaled'},cogLabels,'DeltaPerfW2','DeltaPerfW2U','ResponderHbo','ResponderHbr'];
missingTableDataMask = any(isnan(lmTable{:,:}),2); %find rows that contain nan.
% if using trail b, has 2 missing data.
% lmTable(missingTableDataMask,:) = []; %remove the nan rows? 53,54 has nan for trail b value when NMCM and PRIMA together
%lowest education (outlier: 78); max performance error: 104; race outlier:
%18,63,134 %[104,78,18,63,135]
%104,130 poor performer
% lmTable([104 130],:) = []; %remove poor performers.
% warning('Removing sample')
%set up MMSE into categories
%Option1. 18-23 MCI, 24-30 normal
% lmTable.MMSE_total = (lmTable.MMSE_total >= 24) +1;
%option2. above 29: normal, above 25 separates AD and MCI
% lmTable.MMSE_total = 1*(lmTable.MMSE_total >= 29) + 2*(lmTable.MMSE_total >= 25 & lmTable.MMSE_total <= 28);

%set up a large table with both tasks 
lmTableMerged =  [[lmTable.W2Hbo; lmTable.W2UHbo],[lmTable.W2Hbr; lmTable.W2UHbr],[lmTable.W2HboScaled; lmTable.W2UHboScaled],...
    [lmTable.W2HbrScaled;lmTable.W2UHbrScaled],[lmTable.DeltaPerfW2;lmTable.DeltaPerfW2U],...
    [ones(size(lmTable.DeltaPerfW2));2*ones(size(lmTable.DeltaPerfW2))]];
lmTableMerged = array2table(lmTableMerged);
lmTableMerged.Properties.VariableNames = [{'Hbo'},{'Hbr'},{'HboScaled'},{'HbrScaled'},{'DeltaPerf'},{'TaskDifficulty'}];
lmTableMerged = [lmTableMerged,repmat(lmTable(:,[9:15]),2,1)];

%create customized zscore ignore nan: https://www.mathworks.com/matlabcentral/answers/249566-zscore-a-matrix-with-nan
%standard zscore function from matlab, if any column contains nan this will make the whole column nan
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan')); 

lmTableMerged{:,:} = zscor_xnan(lmTableMerged{:,:}); 
lmTable{:,:} = zscor_xnan(lmTable{:,:});
%fixed model or upper model if using stepwise.
% cogTerms = 'MMSE_total + age + sexF + eduyr + TMTBMinusTMTA'; 
% cogTerms = 'MMSE_total + age + sexF + eduyr + race';%keep demographic only
cogTerms = 'MMSE_total';
%terms to choose in stepwise regression
predictorTerms = split(cogTerms,' + ');
% cogTerms = 'age + MMSE_total + sexF + eduyr + traila + trailb + TMTBMinusTMTA';
% https://www.mathworks.com/matlabcentral/fileexchange/60551-vif-x?tab=discussions;
% vif is the diagonal elements of the inverse of the correlation
% coefficient matrix.
% https://stats.stackexchange.com/questions/558013/is-there-any-way-to-determine-vif-of-some-variable-that-is-included-in-the-datas
cogColForVIF = [9:11 14 15] %exclude traila, b
lmTable.Properties.VariableNames(cogColForVIF)
vif(lmTable{:, cogColForVIF}) %Cog only; this by default uses pearson's correlation, results are not the same but within similar range if using spearman 
lmTable.Properties.VariableNames([16, cogColForVIF])
vif(lmTable{:,[16, cogColForVIF]}) %cog & trailB-A
lmTable.Properties.VariableNames([17, cogColForVIF])
vif(lmTable{:,[17, cogColForVIF]}) %cog & perf W2
lmTable.Properties.VariableNames([18, cogColForVIF])
vif(lmTable{:,[18, cogColForVIF]}) %cog & perf W2U
lmTable.Properties.VariableNames([1,3, cogColForVIF])
vif(lmTable{:,[1,3, cogColForVIF]}) %cog & PFC W2 Hbo and Hbr
lmTable.Properties.VariableNames([2,4, cogColForVIF])
vif(lmTable{:,[2,4, cogColForVIF]}) %cog & PFC W2U Hbo and Hbr
lmTable.Properties.VariableNames([1,cogColForVIF])
vif(lmTable{:,[1,cogColForVIF]}) %cog & PFC W2 Hbo only 
lmTable.Properties.VariableNames([2,cogColForVIF])
vif(lmTable{:,[2,cogColForVIF]}) %cog & PFC W2U Hbo only
lmTable.Properties.VariableNames([3, cogColForVIF])
vif(lmTable{:,[3, cogColForVIF]}) %cog & PFC W2 Hbr only
lmTable.Properties.VariableNames([4, cogColForVIF])
vif(lmTable{:,[4, cogColForVIF]}) %cog & PFC W2U Hbr only
lmTable.Properties.VariableNames([9 14])
vif(lmTable{:,[9 14]}) %MMSE and age only
% lmTable.Properties.VariableNames([9 13])
% vif(lmTable{:,[9 13]}) %MMSE and trailB only
% lmTable.Properties.VariableNames([9:11 14 15])
% vif(lmTable{:,[9:11 14 15]}) %MMSE, age, edu, sex
% lmTable.Properties.VariableNames([12 13])
% vif(lmTable{:,[12 13]}) %trail A,B
% lmTable.Properties.VariableNames([9:13])
% vif(lmTable{:,[9:13]}) %trailAB, age, edu sex

lmTable.sexF = categorical(lmTable.sexF);
lmTable.eduyr = categorical(lmTable.eduyr);
lmTable.race = categorical(lmTable.race);
% lmTable.MMSE_total = categorical(lmTable.MMSE_total);

%by default matlab throws away nan rows for fitlm.
yVars = {'DeltaPerfW2','DeltaPerfW2U'};
for y = 1:2
    if y == 1
%         xTerm = 'W2Hbo + W2Hbr';
        xTerm = 'W2Hbo';
%        xTerm = [];
    else
%         xTerm = 'W2UHbo + W2UHbr';
        xTerm = 'W2UHbo';
%         xTerm = []; %remove Hbo from perf regressor
    end
    for intercept = 1%0:1
        if ~intercept %0, no intercept
            interceptTerm = '-1';
        else
            interceptTerm = '';
        end
        mdl = fitlm(lmTable, [yVars{y} ' ~ ' xTerm ' + ' cogTerms interceptTerm],'RobustOpts',robustOption)
%         mdl = fitlm(lmTable, [yVars{y} ' ~ ' cogTerms interceptTerm])
%         mdl = stepwiselm(lmTable, 'ResponseVar',yVars{y}, 'upper',[yVars{y} ' ~ ' xTerm ' + ' cogTerms interceptTerm],...
%             'PredictorVars',[xTerm, predictorTerms],'Verbose',2,'intercept',logical(intercept))
        Rsqured = mdl.Rsquared
        aicbic = [mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC]
        anova(mdl, 'summary')
        PlotHelper.plotRegressionDiagnostics(mdl, eval(['lmTable.' yVars{y}]), saveResAndFigure, [saveDir 'RegDiagnosis_' yVars{y} '_Intercep_' num2str(intercept)])
        mdlSummary = PlotHelper.compileModelSummaries(mdl,mdlSummary);
        %perform the model without PFC as a regressor
        mdl = fitlm(lmTable, [yVars{y} ' ~ ' cogTerms interceptTerm],'RobustOpts',robustOption)
        aicbic = [mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC]
        anova(mdl, 'summary')
        mdlSummary = PlotHelper.compileModelSummaries(mdl,mdlSummary);
    end
    fprintf('\n\n\n')
end

%%run regression model of PFC using demographic + performance variables.
yVars = {'W2Hbo','W2UHbo','W2Hbr','W2UHbr'};
% yVars = {'W2Hbo','W2UHbo'};
for y = 1:numel(yVars)
    if ismember(y,[1,3])
        xTerm = 'DeltaPerfW2';
%         xTerm = []; %remove performance from Hbo regressors
    else
        xTerm = 'DeltaPerfW2U';
%         xTerm = [];
    end
    for intercept = 1%0:1
        if ~intercept %0, no intercept
            interceptTerm = '-1';
        else
            interceptTerm = '';
        end
        mdl = fitlm(lmTable, [yVars{y} ' ~ ' xTerm ' + ' cogTerms interceptTerm],'RobustOpts',robustOption)
%         mdl = fitlm(lmTable, [yVars{y} ' ~ ' cogTerms interceptTerm])
%         mdl = stepwiselm(lmTable, 'ResponseVar',yVars{y}, 'upper',[yVars{y} ' ~ ' cogTerms interceptTerm],...
%             'PredictorVars',[xTerm, predictorTerms],'Verbose',2,'intercept',logical(intercept))
%         mdl = stepwiselm(lmTable, 'ResponseVar',yVars{y}, 'upper',[yVars{y} ' ~ ' xTerm ' + ' cogTerms interceptTerm],...
%             'PredictorVars',[xTerm, predictorTerms],'Verbose',2,'intercept',logical(intercept))
        Rsqured = mdl.Rsquared
        aicbic = [mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC]
        anova(mdl, 'summary')
%         PlotHelper.plotRegressionDiagnostics(mdl, eval(['lmTable.' yVars{y}]), saveResAndFigure, [saveDir 'RegDiagnosis_' yVars{y} '_Intercep_' num2str(intercept)])
%         mdlSummary = PlotHelper.compileModelSummaries(mdl,mdlSummary);
        %run the model without performance as a regressor
        mdl = fitlm(lmTable, [yVars{y} ' ~ ' cogTerms interceptTerm],'RobustOpts',robustOption)
        aicbic = [mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC]
        anova(mdl, 'summary')
        mdlSummary = PlotHelper.compileModelSummaries(mdl,mdlSummary);
    end
    fprintf('\n\n\n')
end

%%run regression model of scaled index vs demographic
yVars = {'W2HboScaled','W2UHboScaled','W2HbrScaled','W2UHbrScaled'};
for y = 1:numel(yVars)
    for intercept = 1%0:1
        if ~intercept %0, no intercept
            interceptTerm = '-1';
        else
            interceptTerm = '';
        end
        mdl = fitlm(lmTable, [yVars{y} ' ~ ' cogTerms interceptTerm],'RobustOpts',robustOption)
%         mdl = stepwiselm(lmTable, 'ResponseVar',yVars{y}, 'upper',[yVars{y} ' ~ ' cogTerms interceptTerm],...
%             'PredictorVars',predictorTerms,'Verbose',2,'intercept',logical(intercept))
        Rsqured = mdl.Rsquared
        aicbic = [mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC]
        anova(mdl, 'summary')
%         PlotHelper.plotRegressionDiagnostics(mdl, eval(['lmTable.' yVars{y}]), saveResAndFigure, [saveDir 'RegDiagnosis_' yVars{y} '_Intercep_' num2str(intercept)])
        mdlSummary = PlotHelper.compileModelSummaries(mdl,mdlSummary);
    end
    fprintf('\n\n\n')
end

%% fit all metrics together with dummy variable encoding
close all; clc
lmTableMerged =  [[lmTable.W2HboScaled; lmTable.W2Hbo; -lmTable.DeltaPerfW2],...
    [lmTable.W2UHboScaled; lmTable.W2UHbo; -lmTable.DeltaPerfW2U],...
    [ones(size(lmTable.DeltaPerfW2));2*ones(size(lmTable.DeltaPerfW2));3*ones(size(lmTable.DeltaPerfW2))]];
lmTableMerged = array2table(lmTableMerged);
lmTableMerged.Properties.VariableNames = [{'W2'},{'W2U'},{'Metrics'}];
lmTableMerged = [lmTableMerged,repmat(lmTable(:,[9:16]),3,1)];
lmTableMerged.Metrics = categorical(lmTableMerged.Metrics); %1=index, 2=Hbo, 3=Perf
adjustedTerm = ' + age + sexF + eduyr + race';
adjustedWtInteraction = [adjustedTerm '  + age:Metrics + sexF:Metrics + eduyr:Metrics + race:Metrics'];
% run regression with y = 1+MMSE+MMSExCat + Cat, unadjusted
yVars = {'W2','W2U'};
for adjusted = 0:2%1%0:1
    for y = 1:2
        if adjusted == 0 %not adjusted
            mdl = fitlm(lmTableMerged, [yVars{y} '~MMSE_total+Metrics+MMSE_total:Metrics'],'RobustOpts',robustOption)
        elseif adjusted == 1 %adjusted no interaction
            mdl = fitlm(lmTableMerged, [yVars{y} '~MMSE_total+Metrics+MMSE_total:Metrics' adjustedTerm],'RobustOpts',robustOption)
        else %adjusted and interaction
            mdl = fitlm(lmTableMerged, [yVars{y} '~MMSE_total+Metrics+MMSE_total:Metrics' adjustedWtInteraction],'RobustOpts',robustOption)
        end
        Rsqured = mdl.Rsquared
        aicbic = [mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC]
        anova(mdl, 'summary')
        PlotHelper.plotRegressionDiagnostics(mdl, eval(['lmTableMerged.' yVars{y}]), saveResAndFigure, ...
            [saveDir 'RegDiagnosis_Merged_adjusted_' num2str(adjusted)])
        mdlSummary = PlotHelper.compileModelSummaries(mdl,mdlSummary);
    end
end

%% save the mdl summary after running all combinations
if saveResAndFigure
    save([saveDir 'Reg_ModelZscoreSummary' strjoin(studyID) '_Merged_173_WtSensAnalsys_WtInteraction'],'mdlSummary');
    % Convert cell to a table and use first row as variable names
    mdlSumTable = cell2table(mdlSummary(2:end,:),'VariableNames',mdlSummary(1,:));
%     mdlSumZScoreTable = cell2table(mdlSummaryZscore(2:end,:),'VariableNames',mdlSummaryZscore(1,:));
    % Write the table to a CSV file
    writetable(mdlSumTable,[saveDir 'Reg_ModelZscoreSummary' strjoin(studyID) '_Merged_173_WtSensAnalsys_WtInteraction.csv']);
%     writetable(mdlSumZScoreTable,[saveDir 'Reg_ModelZscoreSummary' strjoin(studyID) 'MMSEHboHbr_Demog_Remove78_104.csv']);
end

%% run the merged model with task difficulty as a regressor
lmTableMerged.TaskDifficulty = categorical(lmTableMerged.TaskDifficulty);
cogTerms = 'MMSE_total';
yvars = {'DeltaPerf','Hbo','HboScaled','Hbr','HbrScaled'};
for y = yvars
    mdl = fitlm(lmTableMerged, [y{1} ' ~ TaskDifficulty * ' cogTerms '-1']) %no intercept
    Rsqured = mdl.Rsquared
    anova(mdl, 'summary')
%     mdlSummary = PlotHelper.compileModelSummaries(mdl,mdlSummary);
end

%% plot index vs performance, Hbo fitted vs actual comparison
%specify index for W2.
% indexMdlRow = 22;  hboMdlRow = 18; PerfMdlRow = 15; %172
% indexMdlRow = 10;  hboMdlRow = 6; PerfMdlRow = 3; %173 huber
% indexMdlRow = 46;  hboMdlRow = 42; PerfMdlRow = 39; %173 default
% indexMdlRow = 34;  hboMdlRow = 30; PerfMdlRow = 27; %168 remove all
% indexMdlRow = 58;  hboMdlRow = 54; PerfMdlRow = 51; %173 no z score
% indexMdlRow = 58;  hboMdlRow = 54; PerfMdlRow = 51; %173 MMSE only, huber
% indexMdlRow = 70;  hboMdlRow = 66; PerfMdlRow = 63; %173 MMSE only, nonrobust
% indexMdlRow = 82;  hboMdlRow = 78; PerfMdlRow = 75; %168 adjusted, nonrobust
% indexMdlRow = 94;  hboMdlRow = 90; PerfMdlRow = 87; %169 adjusted, nonrobust, kept low performance
indexMdlRow = 58;  hboMdlRow = 54; PerfMdlRow = 51; 
saveSuffix = 'non_171_MMSE_Only_merged';
removedOutlierIdx = [];%[3 4];

% set up labels and line locations for problematic points.
%lowest education (outlier: 78); max performance error: 104; race outlier:
%18,63,134 
outlierIdxs = [];%[18,63,78,104,130,135];%[104-3];%[104,78,18,63,135] %[78,18,63,135-1]%
outlierLegends = {' '};%{' ','Race:AS','Race:Other','LowEd','BadPerf','BadPerfAll','Race:AI'};
outlierLegends(removedOutlierIdx+1)=[] %remove label for sample removed.
for rmIdx = removedOutlierIdx
    outlierIdxs(rmIdx+1:end) = outlierIdxs(rmIdx+1:end) - 1;
end
outlierIdxs(removedOutlierIdx)=[]

for taskIdx = 1:2 %W2 then W2u
    close all; clc
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    if taskIdx == 2 
        %always assum indx and hbo is next row below, performance is 2 rows
        %below
        indexMdlRow = indexMdlRow+1;  hboMdlRow = hboMdlRow+1; PerfMdlRow = PerfMdlRow+2;
        varNames = {'W2UHboScaled','W2UHbo','DeltaPerfW2U'};
        taskName = 'W2U';
    else
        taskName = 'W2';
    end
    
    mdlRows = [indexMdlRow,hboMdlRow,PerfMdlRow];%index, perf, hbo in order.
    varNames = {'HboScaled','Hbo','DeltaPerf'};
    varNames = [strcat(taskName,varNames(1:2)), strcat(varNames(3),taskName)]; %create var naming format matching the regression var name in the lmTable
    for plotIdx = 1:3
        figure(1);
        subplot(2,3,plotIdx); hold on;
        mdlSummary{mdlRows(plotIdx),end}.Formula
        scatter(mdlSummary{mdlRows(plotIdx),end}.Fitted, eval (['mdlSummary{mdlRows(plotIdx),end}.Variables.' varNames{plotIdx}]),'k');
    %     [temp_rho, temp_p] = corr(mdlSummary{mdlRows(plotIdx),end}.Fitted, eval (['mdlSummary{mdlRows(plotIdx),end}.Variables.' varNames{plotIdx}]))
    %     [temp_rho, temp_p] = corr(eval (['mdlSummary{mdlRows(plotIdx),end}.Variables.' varNames{plotIdx}]),mdlSummary{mdlRows(plotIdx),end}.Variables.trailb)
        for outlierLoc = outlierIdxs
            scatter(mdlSummary{mdlRows(plotIdx),end}.Fitted(outlierLoc), eval (['mdlSummary{mdlRows(plotIdx),end}.Variables.' varNames{plotIdx} '(outlierLoc)']),'*','LineWidth',3);
        end
    
        xlabel(['Fitted ' mdlSummary{mdlRows(plotIdx),1}]); ylabel(['Actual ' mdlSummary{mdlRows(plotIdx),1}]);
        hold on; plot(xlim, xlim,'DisplayName','y=x'); 
        if plotIdx == 1
            legend([outlierLegends,'y=x'])
        end
        axis square;
        title(sprintf('R^2=%.2f, p=%.3f, BIC=%.2f,n=%d',mdlSummary{mdlRows(plotIdx),4},mdlSummary{mdlRows(plotIdx),7},mdlSummary{mdlRows(plotIdx),9},mdlSummary{mdlRows(plotIdx),13}))
        hold off;

        subplot(2,3,3+plotIdx);
        plot(mdlSummary{mdlRows(plotIdx),end},'Marker','.')
        mdlSummary{mdlRows(plotIdx),end}.plotPartialDependence('MMSE_total'); hold on;
        scatter(mdlSummary{mdlRows(plotIdx),end}.Variables.MMSE_total,eval(['mdlSummary{mdlRows(plotIdx),end}.Variables.' varNames{plotIdx}]))
        axis square;
        
        figure(2);
        if ~isempty(mdlSummary{mdlRows(plotIdx),end}.Robust)
            spCol = 3;
            subplot(3,spCol,plotIdx*spCol-2); hold on;
            plot(mdlSummary{mdlRows(plotIdx),end}.Robust.Weights,'x')
            for outlierLoc =outlierIdxs  
                xline(outlierLoc)
            end
            ylabel('SampleWeights')
            title(['Weights ' mdlSummary{mdlRows(plotIdx),1}]);
            xlabel('SampleIdx')
        else
            spCol = 2;
        end
        
        subplot(3,spCol,plotIdx*spCol-1) %diagnostics
        plotDiagnostics(mdlSummary{mdlRows(plotIdx),end},'cookd')
        hold on;
        for outlierLoc = outlierIdxs%[104,78,18,63,135] %[78,18,63,135-1]%
            xline(outlierLoc)
        end
        
        subplot(3,spCol,plotIdx*spCol) %diagnostics
        plotDiagnostics(mdlSummary{mdlRows(plotIdx),end})
        hold on;
        for outlierLoc = outlierIdxs %[78,18,63,135-1]%
            xline(outlierLoc)
        end
        if plotIdx == 1
            legend(outlierLegends)
        end
    end
    
    %display key model info.
    mdlSummary([1, indexMdlRow,hboMdlRow,PerfMdlRow],[1,2,4,6,7,9,13])
    
    if saveResAndFigure
        figure(1)
        set(findall(gcf,'-property','FontSize'),'FontSize',19)
        saveas(f,[saveDir 'zScore_MdlFitComparison_Robust_' saveSuffix '_' taskName ''])
        set(gcf,'renderer','painters')
        saveas(f,[saveDir 'zScore_MdlFitComparison_Robust_' saveSuffix '_' taskName '.png'])
        figure(2)
        set(findall(gcf,'-property','FontSize'),'FontSize',19)
%         saveas(f2,[saveDir 'zScore_MdlFitComparison_Robust_' saveSuffix '_' taskName '_diagnostics.fig'])
        set(gcf,'renderer','painters')
        saveas(f2,[saveDir 'zScore_MdlFitComparison_Robust_' saveSuffix '_' taskName '_diagnostics.png'])
    end
end
%% Exploratory: find unique combinations of age and MMSE and rows with duplicate value (the duplicate will cause model under fit)
tempData = [lmTable.age, lmTable.MMSE_total];
[tempDataUnique, Iunique,~] = unique(tempData, 'rows','first');
size(tempDataUnique,1) < size(tempData,1)
ixDupRows = setdiff(1:size(tempData,1),Iunique);
dupRowVals = tempData(ixDupRows,:)

%% Exploratory: logistic regression to see if the responder and non-responder are different
% responderY = responders{1};
% responderY(missingTableDataMask) = []; %remove missing data.
% responderY = categorical(responderY);
% responseCategories = categories(responderY) %last one will be ref category
% fprintf('Ref category: %s', responseCategories{end})
% predictorVar = lmTable{:,[9,12:14]};
% % predictorVar =
% % [predictorVar,dummyvar(lmTable.sexF),dummyvar(lmTable.eduyr)] %ill
% % conditioned here
% lmTable.Properties.VariableNames([9,12:14,10,11])
% predictorVar = [predictorVar, double(lmTable.sexF),double(lmTable.eduyr)];
% [B,dev,stats] = mnrfit(predictorVar,responderY);
% B
% dev %deviance, best is 1, smaller the better.
% stats.p
lmTable.ResponderHbo = categorical(lmTable.ResponderHbo);
lmTable.ResponderHbr = categorical(lmTable.ResponderHbr);
lgFit=fitglm(lmTable,['ResponderHbo ~ ' cogTerms],'Distribution','binomial','Link','logit')
lgFit.Rsquared
lgFit=fitglm(lmTable,['ResponderHbo ~ ' cogTerms '-1'],'Distribution','binomial','Link','logit')
lgFit.Rsquared
lgFit.devianceTest

%% Exploratory: plot significant ones.
% w2hbo = pfcDataAll{3}; %W2 and W2U
% w2hbo = w2hbo(:,1);
% PlotHelper.computeAndPlotCorrelations(cogData(:,1),w2hbo,subjectID,['Correlation Between Age vs ' saveStr{3} ' W2'],...,
%     cogLabels(1), {'W2Hbo'}, saveResAndFigure, [saveDir 'Corr_Age_vs_' saveStr{3} '_W2'],true);
% PlotHelper.computeAndPlotCorrelations(cogData(:,end),w2hbo,subjectID,['Correlation Between MMSE vs ' saveStr{3} ' W2'],...,
%     cogLabels(end), {'W2Hbo'}, saveResAndFigure, [saveDir 'Corr_MMSE_vs_' saveStr{3} '_W2'],true);
mdl = fitlm(lmTable, 'W2HboScaled ~ 1 + age + sexF + eduyr + traila + trailb + MMSE_total')
figure(); mdl.plotPartialDependence('MMSE_total')
hold on; scatter(lmTable.MMSE_total,lmTable.W2HboScaled)

mdl = fitlm(lmTable, 'W2HboScaled ~ age + sexF + eduyr + traila + trailb + MMSE_total -1 ')
figure(); mdl.plotPartialDependence('age')
hold on; scatter(lmTable.age,lmTable.W2HboScaled)
figure(); mdl.plotPartialDependence('MMSE_total')
hold on; scatter(lmTable.MMSE_total,lmTable.W2HboScaled)

% plotFitX = xlim;
% plotFitY = mdl.Coefficients.Estimate(end) * plotFitX + mdl.Coefficients.Estimate(1);
% plot(xlim,plotFitY,'k','LineWidth',2.5,'handleVisibility','off');
% PlotHelper.computeAndPlotCorrelations(cogData,[pfcDataAll{4}],subjectID,['Correlation Between Demographic vs ' saveStr{4}],...,
%     cogLabels, {'W2Hbr','W2UHbr'}, saveResAndFigure, [saveDir 'Corr_Cog_vs_' saveStr{4}],true);

%% MCI Analysis: Compare between subjects without or without CDR for move participants
if strcmp(studyID, 'Move MYHAT') 
    cdrData = readcell([scriptDir filesep 'Data' filesep 'EmmaDataset' filesep '01272023Data' filesep 'MYHAT_CDR_matchToHMBLv1' filesep 'MMH_HMBLv1_CDRfromMYHATstudy.csv']);
    cdrDataHeader = cdrData(1,:);
    cdrData = cdrData(2:end,:);
    missingCdrMask = cellfun(@(x) any(isa(x,'missing')), cdrData);
    cdrData(missingCdrMask) = {nan};
    cdrData = cell2mat(cdrData);
    cdrCol = find(ismember(cdrDataHeader,{'MOVE_ID','CDR_MYHATcycle_matchHMBL1','maxCDR_allMYHATcycles_beforeHMBL1'}));
    cdrIdCol = cdrCol(1); cdrCurrCol = cdrCol(2); cdrPastCol = cdrCol(3); 
    currMCIMask = cdrData(:,cdrCurrCol) > 0;
    currMCI = cdrData(currMCIMask,:);
    pastMCIMask = cdrData(:,cdrPastCol) > 0 & cdrData(:,cdrCurrCol) == 0;
    pastMCI = cdrData(pastMCIMask,:);
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    varsToPlot = {'W2HboScaled','W2Hbo','DeltaPerfW2','W2UHboScaled','W2UHbo','DeltaPerfW2U'};
    varsToPlot = {'age','sexF','eduyr'}
    for varIdx = 1:numel(varsToPlot);
        subplot(2,3,varIdx); hold on;
        currVar = varsToPlot{varIdx};
        PlotHelper.plotCI(1,eval(['lmTable.' currVar '(currMCIMask)']),'k','CurrentMCI',true)
        PlotHelper.plotCI(2,eval(['lmTable.' currVar '(pastMCIMask)']),'k','PastMCI',false)
        PlotHelper.plotCI(3,eval(['lmTable.' currVar '((~currMCIMask) & (~pastMCIMask))']),'k','NonMCI',false)
        xticks([1 2 3]); xticklabels({'CurrentMCI','PastMCI','NonMCI'});
        title(currVar)
    end
    legend();
    if saveResAndFigure
        set(gcf,'renderer','painters')
        saveas(f,[saveDir 'MCIvsControlIndexComp' strjoin(studyID) '.fig'])
        saveas(f,[saveDir 'MCIvsControlIndexComp' strjoin(studyID) '.png'])
    end
end

%% Schematics of the scaling procedure
sampleSubjIndex = [7,23,21];%7: deltaPerf = 0 at W2 and decrease at W2U, scaled up; 23: deltaper always positive, but decrease from W2 to W2U and scaled up
%21: deltaPerf improved. 
sampleColor = {'#0072BD','#D95319','#EDB120'}; %colors to use for sample subject, should match the # of sample subejcts index chosen
otherSubColor = 'k'; %main color of all the subjects
fScalingProcess = figure('units','normalized','outerposition',[0 0 1 1]); subplot(2,3,1); hold on; %plot the original delta perf
plot(1:2,deltaPerf(:,3:4),'o-','LineWidth',1,'MarkerSize',5,'Color', otherSubColor);
plot(2,deltaPerf(:,4),'.','MarkerSize',15,'Color', otherSubColor,'MarkerFaceColor',otherSubColor);
ylabel('\DeltaPerformance'); xlim([0.5,2.5]);xticks([1,2]); xticklabels({'W2','W2U'})
for s=1:numel(sampleSubjIndex)
    plot(1:2,deltaPerf(sampleSubjIndex(s),3:4),'o-','LineWidth',2,'MarkerSize',5,'Color',sampleColor{s}); %sample subject
    plot(2,deltaPerf(sampleSubjIndex(s),4),'.','MarkerSize',15,'Color',sampleColor{s},'MarkerFaceColor',sampleColor{s}); %sample subject
end
title('\DeltaPerformance')
scaleFactor = str2num(scaledData{1}(end-2:end))
scale = exp(-scaleFactor*deltaPerf);
subplot(2,3,2); hold on; %Plot the exp transformation
plot(deltaPerf(:,3),scale(:,3),'o','MarkerSize',8,'Color', otherSubColor,'HandleVisibility','off')
plot(deltaPerf(:,4),scale(:,4),'.','MarkerSize',24,'Color', otherSubColor,'HandleVisibility','off')
xlabel('\DeltaPerformance'); ylabel('Gain');%ylabel('exp(-\alpha*\DeltaPerf)');
for s=1:numel(sampleSubjIndex)
    plot(deltaPerf(sampleSubjIndex(s),3),scale(sampleSubjIndex(s),3),'o','MarkerSize',8,'LineWidth',2,'Color',sampleColor{s}); %sample subject
    plot(deltaPerf(sampleSubjIndex(s),4),scale(sampleSubjIndex(s),4),'.','MarkerSize',24,'LineWidth',2,'Color',sampleColor{s}); %sample subject
end
title('Gaint = f(\DeltaPerformance)');
legend({'SampleSubject1(evenABC)','SampleSubject1(UnevenABC)'})
% subplot(2,4,3); hold on; %plot transformed result
% plot(1:2,scale(:,3:4),'o-','LineWidth',1,'MarkerSize',5,'Color', otherSubColor);
% ylabel('exp(-\alpha*\DeltaPerf)'); xlim([0.5,2.5]);xticks([1,2]); xticklabels({'W2','W2U'})
% for s=1:numel(sampleSubjIndex)
%     plot(1:2,scale(sampleSubjIndex(s),3:4),'o-','LineWidth',2,'MarkerSize',5,'Color',sampleColor{s});%sample subject
% end
% title('exp(-\alpha*\DeltaPerf)');

subplot(2,3,4); hold on; %raw hbo
grid on; yline(1,'--','LineWidth',2,'Color',"k",'DisplayName','ShiftAmount'); 
plot(1:2,hboData(:,3:4),'o-','LineWidth',1,'MarkerSize',5,'Color', otherSubColor);
plot(2,hboData(:,4),'.','LineWidth',1,'MarkerSize',15,'Color', otherSubColor);
ylabel('Hbo'); xlim([0.5,2.5]);xticks([1,2]); xticklabels({'W2','W2U'});
yticks([-5 0 1 5 10 15]);
for s=1:numel(sampleSubjIndex)
    plot(1:2,hboData(sampleSubjIndex(s),3:4),'o-','LineWidth',2,'MarkerSize',5,'Color',sampleColor{s});%sample subject
    plot(2,hboData(sampleSubjIndex(s),4),'.','LineWidth',2,'MarkerSize',15,'Color',sampleColor{s});%sample subject
end
yRangeOriginal = ylim; 
title('Hbo');
subplot(2,3,5); hold on; %plot the shift
% shift to be in all positive range (not necessarily right, losing info about increase/decrease compared to rest)
shiftAmount = abs(min(hboData(:,3:4),[],'all')) + 1;
hboShifted = hboData(:,3:4) + shiftAmount;
grid on;yline(1,'--','LineWidth',2,'Color',"k",'DisplayName','\epsilon=1'); 
plot(1:2,hboShifted,'o-','LineWidth',1,'MarkerSize',5,'Color', otherSubColor,'HandleVisibility','off');
plot(2,hboShifted(:,2),'.','LineWidth',1,'MarkerSize',15,'Color', otherSubColor,'HandleVisibility','off');
ylabel('PFCActivation_{Hbo}'); xlim([0.5,2.5]);xticks([1,2]); xticklabels({'W2','W2U'});
yticks([-5 0 1 5 10 15]);
for s=1:numel(sampleSubjIndex)
    plot(1:2,hboShifted(sampleSubjIndex(s),:),'o-','LineWidth',2,'MarkerSize',5,'Color',sampleColor{s},'HandleVisibility','off');%sample subject
    plot(2,hboShifted(sampleSubjIndex(s),2),'.','LineWidth',2,'MarkerSize',15,'Color',sampleColor{s},'HandleVisibility','off');%sample subject
end
title('PFCActivation_{Hbo}'); %legend();
%make the 2 sub graph hae the same y-axis range (take the min from the original and max from the shifted).
yRangeShifted = ylim; ylim([yRangeOriginal(1) yRangeShifted(2)]);
subplot(2,3,4);ylim([yRangeOriginal(1) yRangeShifted(2)]); %shift axis of raw hbo

subplot(2,3,[3,6]); hold on; %raw hbo
plot(1:2,pfcDataAll{3},'o-','LineWidth',1,'MarkerSize',5,'Color', otherSubColor,'handlevisibility','off');
plot(2,pfcDataAll{3}(:,2),'.','LineWidth',1,'MarkerSize',15,'Color', otherSubColor,'handlevisibility','off');
ylabel('AutomaticityIndex'); xlim([0.5,2.5]);xticks([1,2]); xticklabels({'W2','W2U'});
for s=1:numel(sampleSubjIndex)
    plot(1:2,pfcDataAll{3}(sampleSubjIndex(s),:),'o-','LineWidth',2,'MarkerSize',5,'DisplayName','SampleSubject','Color',sampleColor{s});%sample subject
    plot(2,pfcDataAll{3}(sampleSubjIndex(s),2),'.','LineWidth',2,'MarkerSize',15,'DisplayName','SampleSubject','Color',sampleColor{s});%sample subject
end
title('AutomaticityIndex');
xlabel('Index = Gaint * PFCActivation')
sgtitle('Computation Process of Automaticity Index')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf,'renderer','painters')
if saveResAndFigure
    saveas(fScalingProcess,[saveDir studyID 'HboScalingSchematicsWithSampleSub3'])
    saveas(fScalingProcess,[saveDir studyID 'HboScalingSchematicsWithSampleSub3.png'])
end

%% descriptive figures of the variables that contribute to the indexes (raw gait speed of wk, w2,w2u, raw alphabet, PFC)
%bargraph or box plot.

%% Exploratory: check correlation of PFC and performance change - uncorrelated with deltaPerf or raw perf - replaced by regression analysis.
% %pfcDataAll first 2 cells are unscaled hbo and hbr, deltaPerf matrix with 4
% %columns for wk, s2, then w2, w2uneven
% PlotHelper.computeAndPlotCorrelations([deltaPerf(:,3:4),deltaPerfGaitOrCog],pfcDataAll{1},subjectID,['Correlation Between Hbo Unscaled And PerformanceChange in NIRS Task'],...,
%     {'DeltaPerfW2','DeltaPerfW2U','(W2-W)/W','(W2U-W)/W','(W2-S2)/S2','(W2U-S2)/S2'}, {'HboW2','HboW2U'}, saveResAndFigure, [saveDir 'Corr_UnscaledHbo_PerfChange'],true);
% PlotHelper.computeAndPlotCorrelations([deltaPerf(:,3:4),deltaPerfGaitOrCog],pfcDataAll{2},subjectID,['Correlation Between Hbr Unscaled And PerformanceChange in NIRS Task'],...,
%     {'DeltaPerfW2','DeltaPerfW2U','(W2-W)/W','(W2U-W)/W','(W2-S2)/S2','(W2U-S2)/S2'}, {'HbrW2','HbrW2U'}, saveResAndFigure, [saveDir 'Corr_UnscaledHbr_PerfChange'],true);
% % PFC vs raw performance
% PlotHelper.computeAndPlotCorrelations(perfDataArray,pfcDataAll{1},subjectID,['Correlation Between Hbo Unscaled And Performance in NIRS Task'],...,
%     {'Wk','S2','W2Speed','W2Alph','W2USpeed','W2UAlpha'}, {'HboW2','HboW2U'}, saveResAndFigure, [saveDir 'Corr_UnscaledHbo_Perf'],true);
% PlotHelper.computeAndPlotCorrelations(perfDataArray,pfcDataAll{2},subjectID,['Correlation Between Hbr Unscaled And Performance in NIRS Task'],...,
%     {'Wk','S2','W2Speed','W2Alph','W2USpeed','W2UAlpha'}, {'HbrW2','HbrW2U'}, saveResAndFigure, [saveDir 'Corr_UnscaledHbr_Perf'],true);
% 
% %% check if change in gait and change in alphabet are related?
% PlotHelper.computeAndPlotCorrelations(deltaPerfGaitOrCog(:,[1,2]),deltaPerfGaitOrCog(:,[3,4]),subjectID,['Correlation Between Change in Gait and Change in Cognitive Performance'],...,
%     {'(W2-W)/W','(W2U-W)/W'}, {'(W2-S2)/S2','(W2U-S2)/S2'}, saveResAndFigure, [saveDir 'Corr_Gait_vs_Cog'],true);
% 
% %% cognitive & demographics variable vs performance change, partial control for PFC.
% PlotHelper.computeAndPlotCorrelations(cogData,[deltaPerf(:,3),deltaPerfGaitOrCog(:,[1,3])],subjectID,['PartialControlHbo Correlation Between Demographic vs PerformanceChange W2'],...,
%     cogLabels, {'DeltaPerfW2','(W2-W)/W','(W2-S2)/S2'}, saveResAndFigure, [saveDir 'PartialCorrCtrHbo_Cog_vs_PerfChange_W2'],true, pfcDataAll{1}(:,1));
% PlotHelper.computeAndPlotCorrelations(cogData,[deltaPerf(:,4),deltaPerfGaitOrCog(:,[2,4])],subjectID,['PartialControlHbo Correlation Between Demographic vs PerformanceChange W2Uneven'],...,
%     cogLabels, {'DeltaPerfW2U','(W2U-W)/W','(W2U-S2)/S2'}, saveResAndFigure, [saveDir 'PartialCorrCtrHbo_Cog_vs_PerfChange_W2Uneven'],true, pfcDataAll{1}(:,2));
PlotHelper.computeAndPlotCorrelations(cogData,[deltaPerf(:,3:4),deltaPerfGaitOrCog(:,[1,3,2,4])],subjectID,['NoPartial Correlation Between Demographic vs PerformanceChange'],...,
    cogLabels, {'DeltaPerfW2','DeltaPerfW2U','(W2-W)/W','(W2-S2)/S2','(W2U-W)/W','(W2U-S2)/S2'}, saveResAndFigure, [saveDir 'NoPartial_Cog_vs_PerfChange'],true);
% 
% % demographics vs raw perf control for pfc
% PlotHelper.computeAndPlotCorrelations(cogData,perfDataArray(:,[3,4]),subjectID,['PartialControlHbo Correlation Between Demographic vs Performance W2'],...,
%     cogLabels, {'W2Speed','W2Alpha'}, saveResAndFigure, [saveDir 'PartialCorrCtrHbo_Cog_vs_Perf_W2'],true, pfcDataAll{1}(:,1));
% PlotHelper.computeAndPlotCorrelations(cogData,perfDataArray(:,[5,6]),subjectID,['PartialControlHbo Correlation Between Demographic vs Performance W2Uneven'],...,
%     cogLabels, {'W2USpeed','W2UAlpha'}, saveResAndFigure, [saveDir 'PartialCorrCtrHbo_Cog_vs_Perf_W2Uneven'],true, pfcDataAll{1}(:,2));
% 
% %% cognitive & demographics vs PFC, partial control for performance change
% %pfcdataall 1, 2 is unscaled
% PlotHelper.computeAndPlotCorrelations(cogData,[pfcDataAll{1}(:,1),pfcDataAll{2}(:,1)],subjectID,['PartialControlDeltaPerf Correlation Between Demographic vs Hbo and Hbr W2'],...,
%     cogLabels, {'W2Hbo','W2Hbr'}, saveResAndFigure, [saveDir 'PartialCorrCtrDeltaPerf_Cog_vs_HboHbr_W2'],true, deltaPerf(:,3));
% PlotHelper.computeAndPlotCorrelations(cogData,[pfcDataAll{1}(:,2),pfcDataAll{2}(:,2)],subjectID,['PartialControlDeltaPerf Correlation Between Demographic vs Hbo and Hbr W2Uneven'],...,
%     cogLabels, {'W2UHbo','W2UHbr'}, saveResAndFigure, [saveDir 'PartialCorrCtrDeltaPerf_Cog_vs_HboHbr_W2Uneven'],true, deltaPerf(:,4));
% %% cognitive & demographics vs PFC, without partial control
PlotHelper.computeAndPlotCorrelations(cogData,[pfcDataAll{1},pfcDataAll{2}],subjectID,['NoPartial Correlation Between Demographic vs Hbo and Hbr'],...,
    cogLabels, {'W2Hbo','W2UHbo','W2Hbr','W2UHbr'}, saveResAndFigure, [saveDir 'NoPartial_Cog_vs_HboHbr'],true);
% 
% %% Correlation btw cognitive & demographics vs PFC scaled
% %pfcdataAll 3&4 is scaled (hbo then hbr)
PlotHelper.computeAndPlotCorrelations(cogData,[pfcDataAll{3}],subjectID,['Correlation Between Demographic vs ' saveStr{3}],...,
    cogLabels, {'W2HboScaled','W2UHboScaled'}, saveResAndFigure, [saveDir 'Corr_Cog_vs_' saveStr{3}],true);
% PlotHelper.computeAndPlotCorrelations(cogData,[pfcDataAll{4}],subjectID,['Correlation Between Demographic vs ' saveStr{4}],...,
%     cogLabels, {'W2Hbr','W2UHbr'}, saveResAndFigure, [saveDir 'Corr_Cog_vs_' saveStr{4}],true);
% 
% %% parial control for age and MMSE
% %scaled
% PlotHelper.computeAndPlotCorrelations(cogData(:,1),[pfcDataAll{3},pfcDataAll{4}],subjectID,['PartialMMSE Correlation Between Age vs Scaled'],...,
%    cogLabels(1), {'W2Hbo','W2UHbo','W2Hbr','W2UHbr'},  saveResAndFigure, [saveDir 'PartialMMSE_Age_Scaled'],true,cogData(:,6));
% PlotHelper.computeAndPlotCorrelations(cogData(:,6),[pfcDataAll{3},pfcDataAll{4}],subjectID,['PartialAge Correlation Between MMSE vs Scaled'],...,
%     cogLabels(6),{'W2Hbo','W2UHbo','W2Hbr','W2UHbr'},  saveResAndFigure, [saveDir 'PartialAge_MMSE_Scaled'],true,cogData(:,1));
% 
% %perf
% PlotHelper.computeAndPlotCorrelations(cogData(:,1),[deltaPerf(:,3),deltaPerfGaitOrCog(:,[1,3])],subjectID,['PartialMMSEHbo Correlation Between Age vs PerformanceChange W2'],...,
%     cogLabels(1), {'DeltaPerfW2','(W2-W)/W','(W2-S2)/S2'}, saveResAndFigure, [saveDir 'PartialMMSEHbo_Age_vs_PerfChange_W2'],true, [pfcDataAll{1}(:,1),cogData(:,6)]);
% PlotHelper.computeAndPlotCorrelations(cogData(:,1),[deltaPerf(:,4),deltaPerfGaitOrCog(:,[2,4])],subjectID,['PartialMMSEHbo Correlation Between Age vs PerformanceChange W2Uneven'],...,
%     cogLabels(1), {'DeltaPerfW2U','(W2U-W)/W','(W2U-S2)/S2'}, saveResAndFigure, [saveDir 'PartialMMSEHbo_Age_vs_PerfChange_W2Uneven'],true, [pfcDataAll{1}(:,2),cogData(:,6)]);
% %vs MMSE control for Age
% PlotHelper.computeAndPlotCorrelations(cogData(:,6),[deltaPerf(:,3),deltaPerfGaitOrCog(:,[1,3])],subjectID,['PartialAgeHbo Correlation Between MMSE vs PerformanceChange W2'],...,
%     cogLabels(6), {'DeltaPerfW2','(W2-W)/W','(W2-S2)/S2'}, saveResAndFigure, [saveDir 'PartialAgeHbo_MMSE_vs_PerfChange_W2'],true, [pfcDataAll{1}(:,1),cogData(:,1)]);
% PlotHelper.computeAndPlotCorrelations(cogData(:,6),[deltaPerf(:,4),deltaPerfGaitOrCog(:,[2,4])],subjectID,['PartialAgeHbo Correlation Between MMSE vs PerformanceChange W2Uneven'],...,
%     cogLabels(6), {'DeltaPerfW2U','(W2U-W)/W','(W2U-S2)/S2'}, saveResAndFigure, [saveDir 'PartialAgeHbo_MMSE_vs_PerfChange_W2Uneven'],true, [pfcDataAll{1}(:,2),cogData(:,1)]);
% 
% %PFC
% %Age vs, control for MMSE, 
% PlotHelper.computeAndPlotCorrelations(cogData(:,1),[pfcDataAll{1}(:,1),pfcDataAll{2}(:,1)],subjectID,['PartialMMSEDeltaPerf Correlation Between Age vs Hbo and Hbr W2'],...,
%     cogLabels(1),{'W2Hbo','W2Hbr'},  saveResAndFigure, [saveDir 'PartialMMSEDeltaPerf_Age_vs_HboHbr_W2'],true, [deltaPerf(:,3),cogData(:,6)]);
% PlotHelper.computeAndPlotCorrelations(cogData(:,1),[pfcDataAll{1}(:,2),pfcDataAll{2}(:,2)],subjectID,['PartialMMSEDeltaPerf Correlation Between Age vs Hbo and Hbr W2Uneven'],...,
%     cogLabels(1), {'W2UHbo','W2UHbr'}, saveResAndFigure, [saveDir 'PartialMMSEDeltaPerf_Age_vs_HboHbr_W2Uneven'],true, [deltaPerf(:,4),cogData(:,6)]);
% %MMSE vs, control for Age, 
% PlotHelper.computeAndPlotCorrelations(cogData(:,6),[pfcDataAll{1}(:,1),pfcDataAll{2}(:,1)],subjectID,['PartialAgeDeltaPerf Correlation Between MMSE vs Hbo and Hbr W2'],...,
%     cogLabels(6),{'W2Hbo','W2Hbr'},   saveResAndFigure, [saveDir 'PartialAgeDeltaPerf_MMSE_vs_HboHbr_W2'],true, [deltaPerf(:,3),cogData(:,1)]);
% PlotHelper.computeAndPlotCorrelations(cogData(:,6),[pfcDataAll{1}(:,2),pfcDataAll{2}(:,2)],subjectID,['PartialAgeDeltaPerf Correlation Between MMSE vs Hbo and Hbr W2Uneven'],...,
%     cogLabels(6),{'W2UHbo','W2UHbr'}, saveResAndFigure, [saveDir 'PartialAgeDeltaPerf_MMSE_vs_HboHbr_W2Uneven'],true, [deltaPerf(:,4),cogData(:,1)]);
% 
% %% correlation btw cognitive and demographics vs PFC scaled responders only
%  %find responder index of cogdata to regress with
% pfcDataResponder = {};
% for dataIdx = 5:6 %responder for Hbo, then Hbr
%     responderIdx = responders{dataIdx-4}; %1x2 cell arrays of 2 column logical vectors.
%     if ~all(responderIdx) %nonresponders exist.
%         PlotHelper.computeAndPlotCorrelations(cogData(responderIdx,:),[pfcDataAll{dataIdx}],subjectID(responderIdx),['Correlation Between Demographic vs ' saveStr{dataIdx}],...,
%             cogLabels, {'W2Hbo','W2UHbo'}, saveResAndFigure, [saveDir 'Corr_Cog_vs_' saveStr{dataIdx}],true,nan,nan,nan,true);
%     end
% end



%% Correlate scaled PFC vs Cognitive and Demographic Variables.
% % removeToNormal = {[],[],{[21,19,18],[],[21],[]}}; %unscaled, scale x2, scale x4, in each cell: W2hbo, W2Uhbo, W2Hbr, W2UHbr
% % removeMeanOutlier = {[],{[21],[],[21],[]},{[21],[19],[21],[19]}};
% % removeNonResponder = {[],{[2,15,16,18:22,25,29],[2,15,16,18:22,25,29],[18,20,21,22,25,27,28],[18,20,21,22,25,27,28]},{[18,20,21,22,25,29],[18,20,21,22,25,29],[18,20,21,22,25,27,28],[18,20,21,22,25,27,28]}};
% % non responder x2: hbo: 3, 16,17,19-23, 27, 31; hbr: 19, 21-23, 27, 29, 30
% % non responder x4: hbo: 19,21,22,23,27,31; hbr: 19,21,22,23,27,29,30
% outlineRemovalMethodStr = {'NoRemoval','ToNormal','Outlier3Mean','RmNonResponders'};
% 
% %FIXME: PRIMA total MMSE has non integer input?
% % close all;
% % clc;
% % sigThreshold = 0.05;
% 
% for removalMethodIdx = [1]%[1,3,4] %1 = no removal
%     if removalMethodIdx == 2
%         outlierToRemove = removeToNormal;%([1,scaledDataToAdd+1]);% [18,19,21]; %25 maybe
%         outlierToRemoveX = removeMeanOutlierX;
%     elseif removalMethodIdx == 3
%         outlierToRemove = removeMeanOutlier;
%         outlierToRemoveX = removeMeanOutlierX;
% %     elseif removalMethodIdx == 4 %not relevant any more. taken care of
% %     above.
% %         outlierToRemove = removeNonResponder; %([1,scaledDataToAdd+1]);
% %         outlierToRemoveX = cell(1,size(cogData,2));
%     end
% 
%     for dataIdx = 1:numel(pfcDataAll) %unscaled, scaled by factors
%         % Manual calculation
%         pfcDataCurr = pfcDataAll{dataIdx}; %4 columns: w2hbo, w2u hbo, w2hr, w2u hbr
% 
%         % Self-implementation of the correlation plots
%         f = figure('units','normalized','outerposition',[0 0 1 1]);%('Position', get(0, 'Screensize'));
%         markerOrder = {'o','+','*','^','x','_','|','.','d','s'};
% 
%         for pfcIdx = 1:size(pfcDataCurr,2) %w2 and w3
%             for cogIdx = 1:size(cogData,2)%col
%                 subplot(size(pfcDataCurr,2),size(cogData,2),(pfcIdx-1)*size(cogData,2)+cogIdx);
%                 hold on;
%                 yToPlot = pfcDataCurr{:,pfcIdx};
%                 xToPlot = cogData(:,cogIdx);
%                 subjectIDCurr = subjectID;
%                 if contains(saveStr{dataIdx},'Responder')
%                     nonResponderToRemove = [responders{dataIdx-2}{pfcIdx}]; %shifted by 2 unscale doesn't have responder only group
%                     xToPlot(nonResponderToRemove,:) = [];
%                     subjectIDCurr(nonResponderToRemove,:) = [];
%                 end
%                 if removalMethodIdx~=1 && ((~isempty(outlierToRemove{dataIdx})) || (~isempty(outlierToRemoveX{cogIdx})))
%                     idxToRemove = [outlierToRemove{dataIdx}{pfcIdx};outlierToRemoveX{cogIdx}];
%                     yToPlot(idxToRemove,:) = [];
%                     xToPlot(idxToRemove,:) = [];
%                     subjectIDCurr(idxToRemove,:) = [];
%                 end
% %                 if dataIdx>=2
% %                     yToPlot = log(yToPlot);
% %                 end
% %                 if cogIdx == 1
% %                     [h,pnorm] = kstest(normalize(yToPlot));
% %                     tf = isoutlier(yToPlot,'mean');
% %                     fprintf('%s, %s, normal(0=normal), %d, p=%.2f outlider: ',saveStr{dataIdx},pfcTaskLabels{pfcIdx},h,pnorm)
% %                     disp(['Outilers: ' subjectID(tf)'])
% %                 end
%                 [rho,p] = corr(yToPlot,xToPlot,'Rows','complete'); %DT only vs cog, coefficient or correlation: rho; rows=complete: ignore nan values
%                 [rhoSpearman,pSpearman] = corr(yToPlot,xToPlot,'Type','Spearman','Rows','complete');
%     %                     xSig = QValueTable{:,nirsIdx};
%     %                     xSig = xSig <= sigThreshold;
%     %                     sigMarkerLabelled = false;
%                 %plot ind subject's dots
%                 for subIdx = 1:length(subjectIDCurr)
%     %                 plot(xToPlot(subIdx),yToPlot(subIdx),'o','Color',colorOrder(mod(subIdx,7),:),'LineWidth',2,'MarkerSize',7,'DisplayName',(subjectID{subIdx}));
%                     plot(xToPlot(subIdx),yToPlot(subIdx),'o','Color',colorOrder(1,:),'LineWidth',2,'MarkerSize',7,'DisplayName',(subjectIDCurr{subIdx}));
%                 end
%                 % Do the regression with an intercept of 0 and plot the line
%                 linFit = fitlm(xToPlot,yToPlot);
%                 plotFitX = xlim;
%                 plotFitY = linFit.Coefficients.Estimate(2) * plotFitX + linFit.Coefficients.Estimate(1);
%                 plot(xlim,plotFitY,'k','LineWidth',5,'handleVisibility','off');
%                 
%                 txtY = ylim;
%                 if p < 0.05 %show the text in red
% %                     text(plotFitX(1),txtY(2)-range(txtY)*0.25,sprintf('P:\\rho=%.2f,p=%.2f',rho(pfcIdx,cogIdx),p(pfcIdx,cogIdx)),'FontSize',40,'Color','r')
%                     text(plotFitX(1),txtY(2)-range(txtY)*0.95,sprintf('P:\\rho=%.4f,p=%.4f',rho,p),'FontSize',20,'Color','r')
%                 else
%                     text(plotFitX(1),txtY(2)-range(txtY)*0.95,sprintf('P:\\rho=%.4f,p=%.4f',rho,p),'FontSize',20,'Color','k')
%                 end              
% %                  text(plotFitX(1),txtY(2)-range(txtY)*0.35,sprintf('R^2=%.2f',linFit.Rsquared.Ordinary),'FontSize',20,'Color','k')
%                 if pSpearman < 0.05 %show the text in red
%                     text(plotFitX(1),txtY(2)-range(txtY)*0.85,sprintf('S:\\rho=%.4f,p=%.4f',rhoSpearman,pSpearman),'FontSize',40,'Color','r')
%                 else
%                     text(plotFitX(1),txtY(2)-range(txtY)*0.85,sprintf('S:\\rho=%.4f,p=%.4f',rhoSpearman,pSpearman),'FontSize',40,'Color','k')
%                 end
%                 if cogIdx == size(cogData,2)
%                     %handle categorical var
% %                     xToPlotCat = categorical(xToPlot);
% %                     disp([saveStr{dataIdx} pfcTaskLabels{pfcIdx} ' Removed Outliers ' outlineRemovalMethodStr{removalMethodIdx}])
% %                     linFitCat = fitlm(xToPlotCat,yToPlot)
%                     %kruskalwallis(x) returns the p-value for the null hypothesis that the data in each column of the matrix x comes from the same distributi
%                     kwp = kruskalwallis(yToPlot,xToPlot,'off');
%                     if kwp < 0.05
%                         textClr = 'r';
%                     else
%                         textClr = 'k';
%                     end
%                     text(plotFitX(1),txtY(2)-range(txtY)*0.1,sprintf('K-W p=%.4f',kwp),'FontSize',20,'Color',textClr)
%                 end
%             
%                 if cogIdx == 1 %col1, y label
%                     if logData
%                         ylabel(['Log ' pfcTaskLabels{pfcIdx}]);
%                     else
%                         ylabel(pfcTaskLabels{pfcIdx});
%                     end
%                 end
%                 if pfcIdx == size(pfcDataCurr,2) %last row, x label
%                     xlabel(cogLabels{cogIdx})
%                 end
%     %                 set(gca,'FontSize',20)
%                 xlim(plotFitX)
%                 axis square
%             end
%         end
% %         legend()
%         sgtitle(['Correlation Between PFC and Cog Measure ' saveStr{dataIdx} ' Removed: ' outlineRemovalMethodStr{removalMethodIdx}])
%         set(findall(gcf,'-property','FontSize'),'FontSize',15)
% 
%         %other ways to plot regression results
%         % refline(linFit.Coefficients.Estimate(2),linFit.Coefficients.Estimate(1),'k'); %coeff first then intercept
%         % plot(linFit)
%         if saveResAndFigure
% %             save([saveDir saveStr{dataIdx} '_' outlineRemovalMethodStr{removalMethodIdx} '_CorrDataAndResults.mat'],'rho','p','rhoSpearman','pSpearman','pfcDataCurr','cogData','-v7.3')
%             saveas(f, [saveCorrDir saveStr{dataIdx} '_' outlineRemovalMethodStr{removalMethodIdx} '_CorrelationRef.fig'])
%     %                 s = findobj('type','legend'); delete(s)
%             saveas(f, [saveCorrDir saveStr{dataIdx} '_' outlineRemovalMethodStr{removalMethodIdx} '_CorrelationRef.png'])
%         end
%     end
% end

% %% run regression models of the PFC data vs all cognitive and demographic variables.
% % clc
% cogCol = find(ismember(perfDataHeader,{'age','sex','eduyr','traila','trailb','MMSE_total'}));
% cogData = perfDataFull(:,[1:4,cogCol]); %should have 5 repeats per subject
% cogData = cogData(strcmp(cogData(:,3),'Even'),:); %take a random one, should have 29 entries, 1 persubject
% cogData = cogData(~missingDataMask,:);
% checkID = isequal(char(subjectID), num2str(cell2mat(cogData(:,1)))) %make sure every subject is present
% cogDataFull = cogData;
% cogData = cogData(:,5:end);
% otherMissingVal = ~cellfun(@isnumeric,cogData);
% cogData(otherMissingVal) = {NaN};
% cogData = cell2mat(cogData); %in format: subject x [triala, trialb, mmse, age]
% 
% for dataIdx = 1:numel(pfcDataAll)
%     for pfcIdx = 1:4
%         fprintf('%s, %s\n',saveStr{dataIdx}, pfcTaskLabels{pfcIdx})
%         ys = pfcDataAll{dataIdx}{:,pfcIdx}; %scaled
%         ys(removeMeanOutlier{dataIdx}{pfcIdx})=[];
%         Xs = cogData;
%         if contains(saveStr{dataIdx},'Responder')
%             nonResponderToRemove = [responders{dataIdx-2}{pfcIdx}]; %shifted by 2 unscale doesn't have responder only group
%             Xs(nonResponderToRemove,:) = [];
%         end
%         Xs(removeMeanOutlier{dataIdx}{pfcIdx},:)=[];
%         mdlFull = fitlm(Xs,ys)
%         mdlMMSEOnly = fitlm(Xs(:,6),ys)
%         mdlAgeOnly= fitlm(Xs(:,1),ys)
%     end
% end
% 
% %unscaled Hbo vs the 4 outcome variables.
% 
