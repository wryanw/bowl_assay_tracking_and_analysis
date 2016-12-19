
% funDir = '/Volumes/card/Ming_RubinLab/tracking_code';
funDir = 'Y:\Ming_RubinLab\tracking_code';
% this is where the function 'flyTracker_mingBowl_multiFly.m' can be found
cd(funDir)

% dataDir = '/track';
% this is where the video files are

runMode = 4;
% 1 - initialize fly positions
% 2 - generate indicator vector and background frames
% 3 - track flies
% 4 - run after all flies have been tracked
% 5 - diagnostic mode

showMovie = 0;
% relevant to modes 3 to 5

skipPrev = 0; % set '1' if you wish to skip previously saved data

dataDirOps = {'Y:\Ming_RubinLab\CsChrimson_circular_arena_assay'
    'Y:\Ming_RubinLab\CsChrimson_circular_arena_assay_new\Figure4'
    'Y:\Ming_RubinLab\CsChrimson_circular_arena_assay_new\Figure5'};
% dataDirOps = {'Y:\Ming_RubinLab\LC16vCtrl_CsChrimson_redux'};
for iterD = 1:numel(dataDirOps)
    dataDir = dataDirOps{iterD};
    flyTracker_mingBowl_multiFly(dataDir,showMovie,runMode,skipPrev)
end

%% reviewing the indicator traces

for iterD = 1:numel(dataDirOps)
    dataDir = dataDirOps{iterD};
    vidList = dir(fullfile(dataDir,'*.avi'));
    vidList = {vidList(:).name};
    [~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
    fileCt = numel(fileNames);
    for iterF = 1:fileCt
        fileRef = iterF;
        videoID = fileNames{fileRef};
        vidPath = fullfile(dataDir,[videoID '.avi']);
        indicatorPath = fullfile(dataDir,[videoID 'indicator.mat']);
        load(indicatorPath)
        ndc8rVec = indicatorVector;
        ndc8rVec = (ndc8rVec-min(ndc8rVec))/range(ndc8rVec);
        figure
        plot(ndc8rVec)
        hold all
        plot(deltaLights)
        uiwait(gcf)
    end
end