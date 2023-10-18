%% This script is the main script that reproduces all the results in Liu et al., 2023
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10418087/
%1. computes auomaticity index on the given study, 
%2. load and summarizes demographic and cognitive data of the given study
%3. performs normality checks of the automaticity, pfc, performance, and
%cognitive data.
%4. run regressions to test relationship btw automaticity index, pfc,
%performance and demographic variables and cognitive variables.
%% Load perf, pfc data and compute automaticity index.
close all; clc; clear all;
studyID = {'Move MYHAT'};%{'Move MYHAT'};%{'NMCM','PRIMA','Move MYHAT'}; %NMCM, PRIMA, MMH
scriptDir = fileparts(matlab.desktop.editor.getActiveFilename); 
scriptDir = strrep(scriptDir,'AutoIndexAnalysis','');
scriptDir = [scriptDir 'Data' filesep];
saveDir = [scriptDir 'IndexAppliedResults' filesep '012723Data' filesep];
scaleMethod = 'Exp';
if iscell(studyID)
    saveDir = [saveDir strjoin(studyID, '_') scaleMethod filesep]
else
    saveDir = [saveDir studyID scaleMethod filesep]
end
if not(isfolder(saveDir))
    mkdir(saveDir)
end
saveResAndFigure = false;

perfDataPath=[scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'Performance_GaitSpeed_Alphabet' filesep 'NMCM_PRIMA_MMH_Visit1_performance_subjavg.csv'];
pfcDataPath = [scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'fNIRS' filesep 'NMCM_PRIMA_MMH_Visit1_fNIRS_subjavg.csv'];

visitNum = 1;
logData = false;

%load pfc, perf data and compute index.
[subjectID, hboData, hbrData, scaledData, responders, perfDataArray, deltaPerf,...
    deltaPerfGaitOrCog, scaleFactors,missingDataMask]=computeAutomaticityIndex(perfDataPath, ...
    pfcDataPath, studyID, visitNum, saveResAndFigure, saveDir, scaleMethod, logData);

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
[p,tbl,stats] = anova1(pfcDataAll{3}(:,1),groups);
p
fprintf('\nW2UHboScaled\n')
[p,tbl,stats] = anova1(pfcDataAll{3}(:,2),groups);
p
fprintf('\nW2HbrScaled\n')
[p,tbl,stats] = anova1(pfcDataAll{4}(:,1),groups);
p
fprintf('\nW2UHbrScaled\n')
[p,tbl,stats] = anova1(pfcDataAll{4}(:,2),groups);
p

%compare Hbo/Hbr
fprintf('\nW2Hbo\n')
[p,tbl,stats] = anova1(pfcDataAll{1}(:,1),groups); p %print p only
fprintf('\nW2UHbo\n')
[p,tbl,stats] = anova1(pfcDataAll{1}(:,2),groups); p %print p only
fprintf('\nW2Hbr\n')
[p,tbl,stats] = anova1(pfcDataAll{2}(:,1),groups); p %print p only
fprintf('\nW2UHbr\n')
[p,tbl,stats] = anova1(pfcDataAll{2}(:,2),groups); p %print p only

%compare combined change performance in gait and cognitive
fprintf('\n\DeltaPerf W2\n')
[p,tbl,stats] = anova1(deltaPerf(:,3),groups); p %print p only
fprintf('\n\DeltaPerf W2U\n')
[p,tbl,stats] = anova1(deltaPerf(:,4),groups); p %print p only

%compare age, sex, race, education?

%% check correlation among the demographic vars
f = figure();
cogDataTable = array2table(cogData);
cogDataTable.Properties.VariableNames = strrep(strrep(cogLabels,'_',''),'il',''); %the figure can only display up to 5 chars for labels, so make traila/b traa/b
[r, p] = corrplot(cogDataTable(:,:),'Type','Spearman','TestR','on')
if saveResAndFigure
    saveas(f, [saveDir 'Demographic_CorrelationMatrix.png'])
end

%% initialize a big cell array for model summaries.
mdlSummary = {'ResponseName','Predictors','SigRegressor','R2Ordinary','R2Adjusted','MMSEBeta','pValue','AIC','BIC','ResidualAbnormal(0Normal)','LackOfFit(1lack)','RobustOpt','NumSub','Model'};
%4. add task as a regressor and task x mmse to see if predicts 1 task
%better than the other. --> no significant interactions

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
cogTerms = 'trailb + age + sexF + eduyr';
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
%         mdlopen = fitlm(lmTable, [yVars{y} ' ~ ' cogTerms interceptTerm])
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
indexMdlRow = 10;  hboMdlRow = 6; PerfMdlRow = 3; 
saveSuffix = 'non_171_tmtb_demoAdjusted';
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
        title(sprintf('R^2=%.2f, p=%.3f, Beta=%.2f,n=%d',mdlSummary{mdlRows(plotIdx),4},mdlSummary{mdlRows(plotIdx),7},mdlSummary{mdlRows(plotIdx),6},mdlSummary{mdlRows(plotIdx),13}))
        hold off;

        subplot(2,3,3+plotIdx);
        plot(mdlSummary{mdlRows(plotIdx),end},'Marker','.')
%         mdlSummary{mdlRows(plotIdx),end}.plotPartialDependence('MMSE_total'); hold on;
% scatter(mdlSummary{mdlRows(plotIdx),end}.Variables.MMSE_total,eval(['mdlSummary{mdlRows(plotIdx),end}.Variables.' varNames{plotIdx}]))
        mdlSummary{mdlRows(plotIdx),end}.plotPartialDependence('trailb'); hold on;
        scatter(mdlSummary{mdlRows(plotIdx),end}.Variables.trailb,eval(['mdlSummary{mdlRows(plotIdx),end}.Variables.' varNames{plotIdx}]))
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
%         saveas(f,[saveDir 'zScore_MdlFitComp_Robust_' saveSuffix '_' taskName ''])
        set(gcf,'renderer','painters')
        saveas(f,[saveDir 'zScore_MdlFitComp_Robust_' saveSuffix '_' taskName '.png'])
        figure(2)
        set(findall(gcf,'-property','FontSize'),'FontSize',19)
%         saveas(f2,[saveDir 'zScore_MdlFitComp_Robust_' saveSuffix '_' taskName '_diagnostics.fig'])
        set(gcf,'renderer','painters')
%         saveas(f2,[saveDir 'zScore_MdlFitCompn_Robust_' saveSuffix '_' taskName '_diagnostics.png'])
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
    cdrData = readcell([scriptDir filesep 'EmmaDataset' filesep '01272023Data' filesep 'MYHAT_CDR_matchToHMBLv1' filesep 'MMH_HMBLv1_CDRfromMYHATstudy.csv']);
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
    
    % find age and sex matched controls for the currentMCI
    mciDemographics = lmTable(currMCIMask,:);
    nonMCI = lmTable((~currMCIMask) & (~pastMCIMask),:);
    matchSubjIdx = [];
    for i = 1:height(mciDemographics)
        matchSexIdx = find(nonMCI.sexF==mciDemographics.sexF(i));
        [minAgeDiff(i),matchAgeIdx(i)] = min(abs(nonMCI.age(matchSexIdx) - mciDemographics.age(i)));
        matchSubjIdx(i) = matchSexIdx(matchAgeIdx(i));
    end
    minAgeDiff    
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    varsToPlot = {'W2HboScaled','W2Hbo','DeltaPerfW2','W2UHboScaled','W2UHbo','DeltaPerfW2U'};
%     varsToPlot = {'age','sexF','eduyr'}
    for varIdx = [1]%1:numel(varsToPlot);
%         subplot(2,3,varIdx); hold on;
        hold on;
        currVar = varsToPlot{varIdx};
        %TODO: plotCI needs to have an additional argument for what kind of
        %CI are we plotting (inverse t or bootstrp?)
%         PlotHelper.plotCI(1,eval(['lmTable.' currVar '((~currMCIMask) & (~pastMCIMask))']),'k','NonMCI',false)
%         PlotHelper.plotCI(2,eval(['lmTable.' currVar '(pastMCIMask)']),'k','PastMCI',false)
        PlotHelper.plotCI(3,eval(['lmTable.' currVar '(currMCIMask)']),PlotHelper.colorOrder(1,:),'CurrentMCI',true)
        PlotHelper.plotCI(4,eval(['lmTable.' currVar '(matchSubjIdx)']),PlotHelper.colorOrder(2,:),'MatchedCtr',true)
%         xticks([1 2 3 4]); xticklabels({'NonMCI','PastMCI','CurrentMCI','AgeSexMatched NonMCI'});
        xticks([1 3 4]); xticklabels({'NonMCI','CurrentMCI','AgeSexMatched NonMCI'});
        title(currVar)
        dataToConnect = [eval(['lmTable.' currVar '(currMCIMask)']),eval(['lmTable.' currVar '(matchSubjIdx)'])];
        plot([3.4 4.4],dataToConnect','k','MarkerSize',0.5,'HandleVisibility','off')
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
PlotHelper.computeAndPlotCorrelations(cogData(),[pfcDataAll{3}],subjectID,['Correlation Between Demographic vs ' saveStr{3}],...,
    cogLabels, {'W2HboScaled','W2UHboScaled'}, saveResAndFigure, [saveDir 'Corr_Cog_vs_' saveStr{3}],true);
% PlotHelper.computeAndPlotCorrelations(cogData,[pfcDataAll{4}],subjectID,['Correlation Between Demographic vs ' saveStr{4}],...,
%     cogLabels, {'W2Hbr','W2UHbr'}, saveResAndFigure, [saveDir 'Corr_Cog_vs_' saveStr{4}],true);
%MMSE only (for grant power analysis)
% PlotHelper.computeAndPlotCorrelations(cogData(:,[6]),[pfcDataAll{3}],subjectID,['Correlation Between Demographic vs ' saveStr{3}],...,
%     cogLabels([6]), {'W2HboScaled','W2UHboScaled'}, saveResAndFigure, [saveDir 'Corr_MMSE_vs_' saveStr{3}],true);
%trials for presentations/exploratory analysis
PlotHelper.computeAndPlotCorrelations([pfcDataAll{3},pfcDataAll{1},deltaPerf(:,3:4)],cogData(:,[6]),subjectID,['Correlation Between Demographic vs ' saveStr{3}],...,
     {'W2HboScaled','W2UHboScaled','W2Hbo','W2UHbo','DeltaPerfW2','DeltaPerfW2U'},cogLabels([6]), saveResAndFigure, [saveDir 'Corr_trails_vs_' saveStr{3}],true);
PlotHelper.computeAndPlotCorrelations([pfcDataAll{3},pfcDataAll{1},deltaPerf(:,3:4)],cogData(:,[6]),subjectID,['Partialed Correlation Between Demographic vs ' saveStr{3}],...,
     {'W2HboScaled','W2UHboScaled','W2Hbo','W2UHbo','DeltaPerfW2','DeltaPerfW2U'},cogLabels([6]), saveResAndFigure, [saveDir 'PartDemCorr_trails_vs_' saveStr{3}],true,cogData(:,[1 2 3 7]));

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
