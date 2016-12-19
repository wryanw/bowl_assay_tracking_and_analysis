function generateStimVector(dataDir,skipPrev)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
if nargin == 0 || isempty(mfilename)
    dataDir = 'D:\rai_tracking\Rai2Ryan_all';
    skipPrev = 0;
end
vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileCt = numel(fileNames);
parfor iterF = 1:fileCt
    fileRef = iterF;
    videoID = fileNames{fileRef};
    vidPath = fullfile(dataDir,[videoID '.avi']);
    indicatorPath = fullfile(dataDir,[videoID 'indicator.mat']);
    backgroundPath = fullfile(dataDir,[videoID 'backgroundFrame.tif']);
    maxPath = fullfile(dataDir,[videoID 'maxFrame.tif']);
    indicatorVars = load(indicatorPath);
    if ~isempty(indicatorVars.deltaLights) && skipPrev
        continue
    end
    frect = indicatorVars.indicatorPos;
    try
        vidObjP = VideoReader(vidPath);
    catch ME
        getReport(ME)
        disp(vidPath)
    end
    frmCt = round(vidObjP.Duration*vidObjP.FrameRate);
    indicatorVector = zeros(1,frmCt);
    for iterFrm = 1:frmCt
        frmRaw = read(vidObjP,iterFrm);
        frmRaw = squeeze(frmRaw(:,:,1,:));
        if iterFrm == 1
            backgroundFrm = frmRaw;
            maxFrm = frmRaw;
        else
            backgroundFrm = min(backgroundFrm,frmRaw);
            maxFrm = max(maxFrm,frmRaw);
        end
        ndc8rFrm = frmRaw(frect(2):frect(2)+frect(4),frect(1):frect(1)+frect(3));
        indicatorVector(iterFrm) = sum(ndc8rFrm(:));
    end
    ndc8rVec = indicatorVector;
    ndc8rVec = (ndc8rVec-min(ndc8rVec))/range(ndc8rVec);
    deltaLights = zeros(1,frmCt)+0.5;
    for iterFrm = 1:frmCt-1
        if ndc8rVec(iterFrm) < 0.5 && ndc8rVec(iterFrm+1) > 0.5
            deltaLights(iterFrm) = 1;
        elseif ndc8rVec(iterFrm) > 0.5 && ndc8rVec(iterFrm+1) < 0.5
            deltaLights(iterFrm) = 0;
        end
    end
%     plot(ndc8rVec)
%     hold all
%     plot(deltaLights)
    parsave_indicatorVars(indicatorPath,indicatorVector,deltaLights,frect);
    parsave_backgroundFrame(backgroundPath,backgroundFrm);
    parsave_maxFrame(maxPath,maxFrm);
end

end

