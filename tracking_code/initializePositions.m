function initializePositions(dataDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
if nargin == 0 || isempty(mfilename)
    dataDir = 'D:\rai_tracking\Rai2Ryan_all';
end
vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileCt = numel(fileNames);
flySize = 50;
for iterF = 1:fileCt

    fileRef = iterF;
    videoID = fileNames{fileRef};
    vidPath = fullfile(dataDir,[videoID '.avi']);
    metaTablePath = fullfile(dataDir,[videoID 'metaInformation.mat']);
    indicatorPath = fullfile(dataDir,[videoID 'indicator.mat']);
    try
        vidObj = VideoReader(vidPath);
    catch ME
        getReport(ME)
        disp(vidPath)
    end
    frmRaw = read(vidObj,1);
    frmRaw = squeeze(frmRaw(:,:,1,:));
    runningFlag = true;
    while runningFlag
        imshow(frmRaw)
        hold on
        figPos = get(gcf,'Position');
        set(gcf,'Position',[figPos(1) round(figPos(2)*0.1) figPos(3:4)])
        if exist(metaTablePath,'file')
            metaTable = load(metaTablePath);
            metaTable = metaTable.metaTable;
            rowNames = metaTable.Properties.RowNames;
            for iterB = 1:numel(rowNames)
                oldID = rowNames{iterB};
                oldRect = str2double(strsplit(oldID,'_'));
                oldRectX = [oldRect(1) oldRect(1) oldRect(3) oldRect(3) oldRect(1)];
                oldRectY = [oldRect(4) oldRect(2) oldRect(2) oldRect(4) oldRect(4)];
                plot(oldRectX,oldRectY,'color','g')
            end
        else
            disp('new table')
        end
        if ~exist(indicatorPath,'file')
            choice = 'Set Indicator';
        else
            indicatorVars = load(indicatorPath);
            oldRect = indicatorVars.indicatorPos;
            oldRect(3) = oldRect(3)+oldRect(1);
            oldRect(4) = oldRect(4)+oldRect(2);
            oldRectX = [oldRect(1) oldRect(1) oldRect(3) oldRect(3) oldRect(1)];
            oldRectY = [oldRect(4) oldRect(2) oldRect(2) oldRect(4) oldRect(4)];
            plot(oldRectX,oldRectY,'color','r')
            drawnow
            choice = questdlg('Select from the following:','What to do?','Set Indicator','New Fly','Remove Selection','New Fly');
        end
        switch choice
            case 'Set Indicator'
                text(0,-20,'Drag a box from top,left to bottom,right over the indicator','color',[0 0 1])
                hrect = imrect;
                indicatorPos = round(getPosition(hrect));
                indicatorVector = [];
                deltaLights = [];
                parsave_indicatorVars(indicatorPath,indicatorVector,deltaLights,indicatorPos);
            case 'New Fly'
                text(0,-20,'Click once on a new fly and press enter/return','color',[0 0 1])
                try
                    [x, y] = getpts(gcf);
                    frect = round([x-round(flySize/2) y-round(flySize/2) x+round(flySize/2) y+round(flySize/2)]);
                    flyIDa = cellfun(@(x) num2str(x),num2cell(frect'),'uniformoutput',false);
                    flyIDa = [flyIDa repmat({'_'},4,1)]';
                    flyIDa = cat(2,flyIDa{1:end-1});
                    if ~exist('metaTable','var')
                        varNames = {'headPos','tailPos','bot_centroid','bot_theta','mvmnt_ref',...
                            'range_trk','change_trk','filtered_pos','filtered_vel','filtered_accel'};
                        metaTable = table(zeros(1,2),zeros(1,2),cell(1,1),cell(1,1),cell(1,1),...
                            cell(1,1),cell(1,1),cell(1,1),cell(1,1),cell(1,1),'VariableNames',varNames,...
                            'RowNames',{flyIDa});
                    else
                        rowNames = metaTable.Properties.RowNames;
                        for iterB = 1:numel(rowNames)
                            oldID = rowNames{iterB};
                            oldRect = str2double(strsplit(oldID,'_'));
                            testA = oldRect(1) >= frect(1) && oldRect(1) <= frect(3);
                            testB = oldRect(3) >= frect(1) && oldRect(3) <= frect(3);
                            testC = oldRect(2) >= frect(2) && oldRect(2) <= frect(4);
                            testD = oldRect(4) >= frect(2) && oldRect(4) <= frect(4);
                            cumTest = testA && (testC || testD) || testB && (testC || testD);
                            if cumTest
                                metaTable(oldID,:) = [];
                            end
                        end
                    end
                    frm_crop = frmRaw(frect(2):frect(4),frect(1):frect(3));
                    close all
                    imshow(frm_crop,'initialmagnification',500)
                    text(0,-20,'Click once on the fly head then the tail and press enter/return','color',[0 0 1])
                    [x, y] = getpts(gcf);
                    x = x+frect(1);
                    y = y+frect(2);
                    metaTable.headPos(flyIDa,:) = [x(1) y(1)];
                    metaTable.tailPos(flyIDa,:) = [x(2) y(2)];
                catch ME
                    getReport(ME)
                    runningFlag = false;
                end
            case 'Remove Selection'
                if exist('metaTable','var')
                    text(0,-20,'Click once on a square to remove and press enter/return','color',[0 0 1])
                    try
                        [x, y] = getpts(gcf);
                        frect = round([x-round(flySize/2) y-round(flySize/2) x+round(flySize/2) y+round(flySize/2)]);
                        rowNames = metaTable.Properties.RowNames;
                        for iterB = 1:numel(rowNames)
                            oldID = rowNames{iterB};
                            oldRect = str2double(strsplit(oldID,'_'));
                            testA = oldRect(1) >= frect(1) && oldRect(1) <= frect(3);
                            testB = oldRect(3) >= frect(1) && oldRect(3) <= frect(3);
                            testC = oldRect(2) >= frect(2) && oldRect(2) <= frect(4);
                            testD = oldRect(4) >= frect(2) && oldRect(4) <= frect(4);
                            cumTest = testA && (testC || testD) || testB && (testC || testD);
                            if cumTest
                                metaTable(oldID,:) = [];
                            end
                        end
                    catch
                        runningFlag = false;
                    end
                end
        end
        if exist('metaTable','var')
            save(metaTablePath,'metaTable')
            clearvars metaTable
        end
        close all
        disp(['progress: ' num2str(round(fileRef/fileCt*100)) ' %'])
    end
end
close all
end

