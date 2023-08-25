%%
close all; clc; clear all;
scaleMethod = 'Exp'; %Linear or Exp
logData = false; %true if shift unscaled to positive and log y-axis for all; false o.w.
studyID = {'NMCM','PRIMA','Move MYHAT'};%{'NMCM','PRIMA','Move MYHAT'}; %NMCM, PRIMA, MMH
scriptDir = fileparts(matlab.desktop.editor.getActiveFilename); 
scriptDir = strrep(scriptDir,'AutoIndexAnalysis','');
scriptDir = [scriptDir filesep 'Data' filesep];

saveDir = [scriptDir 'IndexAppliedResults' filesep '012723Data' filesep];
if iscell(studyID)
    saveDir = [saveDir strjoin(studyID, '_') scaleMethod filesep]
else
    saveDir = [saveDir studyID scaleMethod filesep]
end
saveDir = [saveDir 'CogPerfCheck' filesep]
if not(isfolder(saveDir))
    mkdir(saveDir)
end
saveResAndFigure = true;
%% load perf data
% perfData = readcell('/Users/mac/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofPittsburgh/BEAM lab - fNIRS Workgroup/Rosso_fNIRS_code/Matlab/Emma_SampleDataNMCM_forShuqi/NMCM_PRIMA_MOVE_combined_wide_19Nov2021.csv');
% /Users/mac/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SML/Projects/fNIR Project/Code_NIRS_Automaticity/Data/EmmaDataset/PRIMA_051722
% load 19Nov2021 Data
% perfData = readcell([scriptDir 'EmmaDataset' filesep 'NMCM_PRIMA_MOVE_combined_wide_19Nov2021.csv']);
perfData = readcell([scriptDir 'EmmaDataset' filesep '01272023Data' filesep 'Performance_GaitSpeed_Alphabet' filesep 'NMCM_PRIMA_MMH_Visit1_performance_eaTrial.csv']);

perfDataHeaderFull = perfData(1,:);
studyIdCol = find(strcmp(perfDataHeaderFull,'Cohort'));
perfDataFull = perfData(contains(perfData(:,studyIdCol),studyID),:); %this will also get rid of the header

%in order: ID (numerical version, 10xxx, 2xxx, 21xxx, study, task, surface, gait, alphabet
%the provided array are in order of the data header, so the indexed sub
%data will have headers in the same order.
relevantCols = ismember(perfDataHeaderFull,{'MergeID','Cohort','Test','cond','RCRT','RGEN','NumGen','NumCorrect'});
perfData = perfDataFull(:,relevantCols); %the perfData now will contain the 5 columns above in the same order listed above.
perfDataHeader = perfDataHeaderFull(:,relevantCols);
for col= 5:8 %numerical columns, change NA to NaN
    perfData(strcmp(perfData(:,col),'NA'),col) = {NaN};
end

taskCol = 4; 
tasks = unique(perfData(:,taskCol));
tasks = tasks(contains(tasks,'ABC'))

%% plot NumCorrect vs NumTotal per trial per task
clc; close all;
badPerformers = [];subjWtLarge234FromR1 = []; subjWtLarge34FromR12=[]
for tk = 1:numel(tasks)
    task = tasks{tk}; 
    currData = cell2mat(perfData(strcmp(perfData(:,taskCol),task),[1,3,5:8]));
    repCol = 2; numGenCol = 5; numCorrectCol = 6; 

    groupedID = accumarray((currData(:,repCol)), (currData(:,1)), [], @(v){v} ) ; %cell array of size 4 (4 reps) for each subject
    groupedNumGen = accumarray((currData(:,repCol)), (currData(:,numGenCol)), [], @(v){v} ) ; %cell array of size 4 (4 reps) for each subject
    groupedNumCorct = accumarray((currData(:,repCol)), (currData(:,numCorrectCol)), [], @(v){v} ) ; %cell array of size 4 (4 reps) for each subject

    % visualization, numcorrect vs total per trial
    f = figure('units','normalized','outerposition',[0 0 1 1]); 
    for rep = 1:4
        subplot(3,4,[rep,rep+4]); hold on;
    %     dataToPlot = [groupedNumCorct{rep}(eval(['fullDataIdx' num2str(rep)])),groupedNumGen{rep}(eval(['fullDataIdx' num2str(rep)]))];
        dataToPlot = [groupedNumCorct{rep},groupedNumGen{rep}];
        plot([-0.25,0.25],dataToPlot)
        xticks([-0.2 0.2])
        xticklabels({'NumCorct','NumGen'})
        title(['Rep' num2str(rep)])

        subplot(3,4,rep+8); hold on; %3rd row
        changeVal = dataToPlot(:,2) - dataToPlot(:,1);
        plot(changeVal)
        title('NumGen - NumCorrect')
        xlabel('dataIdx')
        
        [sorted,sortIdx] = sort(changeVal,'descend');
        sortIdxOutlier = sortIdx(sorted >=10);
        badPerformer = groupedID{rep}(sortIdxOutlier);
        legend(num2str(badPerformer'),'Location','southoutside');
        
        %put the info into a separate table
        badPerformer = ismember(cell2mat(perfDataFull(:,1)),badPerformer) & (cell2mat(perfDataFull(:,5)) == rep) & strcmp(perfDataFull(:,6),task);
        badPerformers = [badPerformers;perfDataFull(badPerformer,:)];
    end
    sgtitle(task)
    set(findall(gcf,'-property','FontSize'),'FontSize',18)

    if saveResAndFigure
        saveas(f, [saveDir 'numGenVsNumCorrect_eaTrial' task '.fig'])
        saveas(f, [saveDir 'numGenVsNumCorrect_eaTrial' task '.png'])
    end
    
    % visualize per trial changes
    f = figure('units','normalized','outerposition',[0 0 1 1]); 
%     subplot(1,3,1); hold on;
    [ComId, Idx1,Idx2]=intersect(groupedID{1},groupedID{2});
    ComId=intersect(ComId,groupedID{4}); %intersect the shorter list
    ComId=intersect(ComId,groupedID{3}); %intersect the shorter list
    [ComId, ~,Idx1]=intersect(ComId,groupedID{1}); %intersect the full list again to get index from group1
    [ComId, ~,Idx2]=intersect(ComId,groupedID{2}); %intersect the full list again to get index from group1
    [ComId, ~,Idx3]=intersect(ComId,groupedID{3}); %intersect the full list again to get index from group1
    [ComId, ~,Idx4]=intersect(ComId,groupedID{4}); %intersect the full list again to get index from group1
    
    dataToPlot = [groupedNumCorct{2}(Idx2)-groupedNumCorct{1}(Idx1),...
        groupedNumCorct{3}(Idx3)-groupedNumCorct{1}(Idx1),...
        groupedNumCorct{4}(Idx4)-groupedNumCorct{1}(Idx1)];
    plot([1,2,3],dataToPlot,'o-','MarkerSize',10)
    xticks([1 2 3])
    xticklabels({'Rep2-1','Rep3-1','Rep4-1'})
    ylabel('Rep-Rep1 Num Correct (>0: improved)')
    sgtitle([task ' NumCorrect'])
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    if saveResAndFigure
        saveas(f, [saveDir 'numCorrect_vsTrail1' task '.fig'])
        saveas(f, [saveDir 'numCorrect_vsTrail1' task '.png'])
    end
    
    %find data that had |change| > 5 conssitently
    largeChange = all(abs(dataToPlot) > 5,2);
    largeChange = ComId(largeChange);
    largeChange = ismember(cell2mat(perfDataFull(:,1)),largeChange) & strcmp(perfDataFull(:,6),task);
    subjWtLarge234FromR1 = [subjWtLarge234FromR1;perfDataFull(largeChange,:)];
    
    largeChange = all(abs(dataToPlot(:,2:3)) > 5,2) & (~all(abs(dataToPlot) > 5,2)); %grab rep3 and 4 only
    largeChange = ComId(largeChange);
    largeChange = ismember(cell2mat(perfDataFull(:,1)),largeChange) & strcmp(perfDataFull(:,6),task);
    subjWtLarge34FromR12 = [subjWtLarge34FromR12;perfDataFull(largeChange,:)];
%     numCorrect over trials
%     fullNumCorct = [groupedNumCorct{rep}]
%     figure(); hold on;
%     ids = unique(currData(:,1));
%     for i = 1:numel(ids)
%         currSubMask = currData(:,1) == ids(i);
%         currXs = currData(currSubMask,repCol);
%         currYs = currData(currSubMask,numCorrectCol);
%         plot(currXs(2:end),currYs(2:end) - currYs(1));
%     end
%     plot(cell2mat(currData(:,repCol)),cell2mat(currData(:,numCorrectCol)))

end

%% visualize RCRT and RGEN per subject per trial to identify very poor performance (thresholding)
f = figure('units','normalized','outerposition',[0 0 1 1]); 
rcrtCol = 5;
rgenCol = 6;
rateCols = [5,6]; %rcrt, then rgen
titles = {'RCRT','RGEN'};
outlierIdx = {};%1st RCRT, then RGEN
for subPh = 1:2
    subplot(2,1,subPh); hold on;
    dataToPlot = cell2mat(perfData(:,rateCols(subPh)));
    plot(dataToPlot,'x') %rcrt
    yline(nanmean(dataToPlot));
    outliers = isoutlier(dataToPlot);
    plot(find(outliers), dataToPlot(outliers),'ro') %rcrt
    legend({'data','mean','outliers'})
    title(titles{subPh})
    outlierIdx{subPh} = outliers;  
end
set(findall(gcf,'-property','FontSize'),'FontSize',18)

if saveResAndFigure
    saveas(f, [saveDir 'RCRTVsRGEN.fig'])
    saveas(f, [saveDir 'RCRTVsRGEN.png'])
end
    
%% find out who and what trials are these outliers. 1) ppl who are bad at
%both, 2) who did high RGEN but low RCRT, or 3) low RGEN but high RCRT (?
%these might be ok, super performers.)
ComRow=outlierIdx{1} & outlierIdx{2};
rcrtOnlyOutlier = outlierIdx{1}~=ComRow; %find values in RCRT only
rgenOnlyOutlier = outlierIdx{2}~=ComRow; %find values in RCRT only

badBoth = perfDataFull([ComRow],:);
badRcrt = perfDataFull([rcrtOnlyOutlier],:);
badRgen = perfDataFull([rgenOnlyOutlier],:);

%%

if saveResAndFigure
    % Write the table to a CSV file  
    writetable(cell2table(badPerformers,'VariableNames',perfDataHeaderFull),[saveDir 'lowCorrectHighGen.csv']);
    writetable(cell2table(subjWtLarge234FromR1,'VariableNames',perfDataHeaderFull),[saveDir 'largeChangeFromRep1.csv']);
    writetable(cell2table(subjWtLarge34FromR12,'VariableNames',perfDataHeaderFull),[saveDir 'largeChangeFromRep1And2.csv']);
%     writetable(cell2table(badPerformers,'VariableNames',perfDataHeaderFull),[saveDir 'lowCorrectHighGen.csv']);
    writetable(cell2table(badBoth,'VariableNames',perfDataHeaderFull),[saveDir 'outlierRCRT_outlierRGEN.csv']);
    writetable(cell2table(badRcrt,'VariableNames',perfDataHeaderFull),[saveDir 'outlierRCRTOnly.csv']);
    writetable(cell2table(badRgen,'VariableNames',perfDataHeaderFull),[saveDir 'outlierRGENOnly.csv']);
end

%%
clc
temp = table(cell2mat(subjWtLarge234FromR1(:,[1])),subjWtLarge234FromR1(:,[6]))
unique(temp,'rows')

clc
temp = table(cell2mat(subjWtLarge34FromR12(:,[1])),subjWtLarge34FromR12(:,[6]))
unique(temp,'rows')

clc
temp = table(cell2mat(badPerformers(:,[1])),badPerformers(:,[6]))
unique(temp,'rows')
unique(temp.Var1)
length(find(strcmp(temp.Var2,'Even_ABC')))
length(find(strcmp(temp.Var2,'Uneven_ABC')))
length(find(contains(temp.Var2,'Standing_ABC')))

%%
clc
badPIdAndTask = table(cell2mat(badPerformers(:,[1])),badPerformers(:,[6]))
badBothIdAndTask = table(cell2mat(badRcrt(:,[1])),badRcrt(:,[6]));
badRCRTIdAndTask = table(cell2mat(badBoth(:,[1])),badBoth(:,[6]));
badP1 = intersect(badPIdAndTask,badBothIdAndTask)
badP2 = intersect(badPIdAndTask,badRCRTIdAndTask)

imp234 = table(cell2mat(subjWtLarge234FromR1(:,[1])),subjWtLarge234FromR1(:,[6]));
imp34 = table(cell2mat(subjWtLarge34FromR12(:,[1])),subjWtLarge34FromR12(:,[6]));
intersect(badP1,imp234)
intersect(badP2,imp234)

intersect(badP1,imp34)
intersect(badP2,imp34)

intersect(imp234,imp34)

