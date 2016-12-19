[x,xi] = sort(plotTable.returnData{'loom_stim'});
vidList = plotTable.videoList{'loom_stim'}(xi);

analysisDir = 'Y:\Data_pez3000_analyzed';

for lineRef = 6:numel(vidList)
% 2 - male, good tho
videoID = vidList{lineRef};
disp(videoID)
strParts = strsplit(videoID,'_');
exptID = strParts{4}(5:end);
dateID = strParts{3};
runID = [strParts{1} '_' strParts{2} '_' strParts{3}];
exptResultsRefDir = fullfile(analysisDir,exptID);
assessmentName = [exptID '_rawDataAssessment.mat'];
assessmentPath = fullfile(exptResultsRefDir,assessmentName);
if exist(assessmentPath,'file') == 2
    assessTable_import = load(assessmentPath);
    dataname = fieldnames(assessTable_import);
    assessTable = assessTable_import.(dataname{1});
else
    %         continue
end
vidPath = assessTable.Video_Path{videoID};
implay(vidPath)
set(findall(0,'tag','spcui_scope_framework'),'position',[2200 50 600 950]);
figure
uiwait(gcf)
end