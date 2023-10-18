close all; clc; clear all;
studyID = {'PRIMA'};%{'Move MYHAT'};%{'NMCM','PRIMA','Move MYHAT'}; %NMCM, PRIMA, MMH
scriptDir = fileparts(matlab.desktop.editor.getActiveFilename); 
scriptDir = strrep(scriptDir,'AutoIndexAnalysis','');
scriptDir = [scriptDir 'Data' filesep];

saveDir = [scriptDir 'IndexAppliedResults' filesep '022223PRIMAAllSession' filesep];
if not(isfolder(saveDir))
    mkdir(saveDir)
end
saveResAndFigure = false;
perfDataPath=[scriptDir 'EmmaDataset' filesep 'PRIMAFullData02222023' filesep 'Performance_GaitSpeed_Alphabet' filesep 'PRIMA_AllVisits_performance_subjavg.csv'];
pfcDataPath = [scriptDir 'EmmaDataset' filesep 'PRIMAFullData02222023' filesep 'fNIRS' filesep 'PRIMA_AllVisits_fNIRS_subjavg.csv'];
subjectIdByVisit = {};
scaledHboByVisit = {};

%% get data over sessions 
DataByVisit = struct();
for visitNum = 1:4
    [subjectID, hboData, hbrData, scaledData, responders, perfDataArray, deltaPerf,...
        deltaPerfGaitOrCog, scaleFactors,missingDataMask]=computeAutomaticityIndex(perfDataPath, pfcDataPath, studyID, visitNum, saveResAndFigure, saveDir);
    DataByVisit.ID{visitNum} = subjectID;
    DataByVisit.scaledHbo{visitNum} = [scaledData{2,1}(:,3:4)]; %w2, w2uneven
    DataByVisit.scaledHbr{visitNum} = [scaledData{2,2}(:,3:4)]; %w2, w2uneven
    DataByVisit.deltaPerf{visitNum} = deltaPerf(:,3:4); %w2, w2uneven
    DataByVisit.Hbo{visitNum} = [hboData(:,3:4)]; %w2, w2uneven
    DataByVisit.Hbr{visitNum} = [hbrData(:,3:4)]; %w2, w2uneven
    DataByVisit.gaitSpeed{visitNum} = perfDataArray(:,[1,3,5]); %wk, s2, w2speed, w2alpha, w2uspeed, w2ualpha
    DataByVisit.alphaRate{visitNum} = perfDataArray(:,[2,4,6]);
%perfDataArray = wk, s2, w2speed, w2alpha, w2uspeed, w2ualpha
end

%% plot PFC evolution over time, 
% allScaledHbo = scaledHboByVisit{1}; %w2even: col1,3,5,7 for visit 1-4; w2uneven: column2,4,6,8 for visit1-4
% allScaledHbr = scaledHbrByVisit{1}; %like Hbo
% allDeltaPerf = scaledHboByVisit{1}; %like Hbo
numSubjV1 = length(DataByVisit.ID{1});
variablesToFill = {'scaledHbo','scaledHbr','deltaPerf','Hbo','Hbr','gaitSpeed','alphaRate'};
allData = struct();

%initialize with visit 1.
for vars = 1:length(variablesToFill)
    eval(['allData.' variablesToFill{vars} '= DataByVisit.' variablesToFill{vars} '{1}']);
end

%load data for visit2-4 with nan padding for missing subject.
for visitNum = 2:4
    commonSubj = ismember(DataByVisit.ID{1},DataByVisit.ID{visitNum});
    for vars = 1:length(variablesToFill)
        curVar = eval(['allData.' variablesToFill{vars}]);
        dataColNum = eval(['size(DataByVisit.' variablesToFill{vars} '{visitNum},2)']);
        curVar(:,(visitNum-1)*dataColNum+1:visitNum*dataColNum) = nan(numSubj,dataColNum);
        curVar(commonSubj,(visitNum-1)*dataColNum+1:visitNum*dataColNum) = eval(['DataByVisit.' variablesToFill{vars} '{visitNum}']);
        eval(['allData.' variablesToFill{vars} ' = curVar;'])
    end
%     allScaledHbr(:,visitNum*2-1:visitNum*2) = nan(numSubj,2);
%     allScaledHbr(commonSubj,visitNum*2-1:visitNum*2) = scaledHbrByVisit{visitNum};
end

%find subjects that consistently come back
commonSubj = all(~isnan(allData.scaledHbo)');
commonSubjId = DataByVisit.ID{1}(commonSubj);

%% plot data over session
for vars = [1,4,3]
%Index, Hbo, DeltaPerf
%     PlotHelper.barPlotWithIndiv(eval(['allData.' variablesToFill{vars} '(commonSubj,[1,3,5,7])'])',commonSubjId,{'V1','V2','V3','V4'},variablesToFill{vars},'PRIMA Over Session W2 (12wk apart)',true,[variablesToFill{vars} '_W2'])
%     PlotHelper.barPlotWithIndiv(eval(['allData.' variablesToFill{vars} '(commonSubj,[2,4,6,8])'])',commonSubjId,{'V1','V2','V3','V4'},variablesToFill{vars},'PRIMA Over Session W2Uneven (12wk apart)',true,[variablesToFill{vars} '_W2U'])
    PlotHelper.barPlotWithIndiv(eval(['allData.' variablesToFill{vars} '(commonSubj,[3,5,7]) - allData.' variablesToFill{vars} '(commonSubj,1)'])',commonSubjId,{'V2','V3','V4'},['Diff from V1 ' variablesToFill{vars} ],'PRIMA Diff from V1, W2 (12wk apart)',true,[variablesToFill{vars} '_W2_Diff'])
    PlotHelper.barPlotWithIndiv(eval(['allData.' variablesToFill{vars} '(commonSubj,[4,6,8]) - allData.' variablesToFill{vars} '(commonSubj,2)'])',commonSubjId,{'V2','V3','V4'},['Diff from V1 ' variablesToFill{vars} ],'PRIMA Diff from V1, W2Uneven (12wk apart)',true,[variablesToFill{vars} '_W2U_Diff'])
end
%change from v1 to v2,v3,v4



