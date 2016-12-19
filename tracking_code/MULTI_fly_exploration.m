clear all
% close all
clc
rng('default');
rng(19);

runMode = 2;
if runMode == 1
    dataDir = 'Z:\Ming_RubinLab\CsChrimson_circular_arena_assay';
else
    dataDir = 'Z:\Ming_RubinLab\CsChrimson_circular_arena_assay_new\Figure5';
end
dataList = dir(fullfile(dataDir,'*.mat'));
dataList = {dataList(:).name}';

vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileNames = fileNames(1:end)';
fileCt = numel(fileNames);
%%
graphMode = 3;
graphOps = {'turning','progression','regression'};
graphStr = graphOps{graphMode};
showIndicatorPlots = 0;
showThetaScroll = 0;
meanPlot = 0;
dataCell = cell(fileCt,1);
dataLabels = cell(fileCt,1);
ydataA = [];
ydataB = [];
longMotFullCell = cell(fileCt,1);
turnMotFullCell = cell(fileCt,1);
for iterD = 1:fileCt
    videoID = fileNames{iterD};
    metaTablePath = fullfile(dataDir,[videoID 'metaInformation.mat']);
    metaTable = load(metaTablePath);
    metaTable = metaTable.metaTable;
    indicatorPath = fullfile(dataDir,[videoID 'indicator.mat']);
    indicatorVars = load(indicatorPath);
    indicatorVector = indicatorVars.indicatorVector;
    deltaLights = indicatorVars.deltaLights;
    indicatorPos = indicatorVars.indicatorPos;
    lightsOn = find(deltaLights == 1);
    lightsOff = find(deltaLights == 0);
    if lightsOn(1) > lightsOff(1)
        lightsOff(1) = [];
    end
    if numel(lightsOn) > numel(lightsOff)
        lightsOn = lightsOn(1:numel(lightsOff));
    elseif numel(lightsOn) < numel(lightsOff)
        lightsOff = lightsOff(1:numel(lightsOn));
    end
    if runMode == 2
        lightsOn = lightsOn(1:5);%restricts experiment count to 5
        lightsOff = lightsOff(1:5);
    end
    onDur = abs(min(lightsOn-lightsOff));
%     onDur = 150;
    
    offDur = abs(min(lightsOff-lightsOn));
    if showIndicatorPlots
        disp(videoID)
        plot(deltaLights)
        hold all
        indicatorVector = (indicatorVector-min(indicatorVector(:)))/range(indicatorVector);
        plot(indicatorVector)
        uiwait(gcf)
        continue
    end
    flyCt = size(metaTable,1);
    latMotionCell = cell(flyCt,numel(lightsOn));
    longMotionCell = cell(flyCt,numel(lightsOn));
    turningCell = cell(flyCt,numel(lightsOn));
    finalPosition = cell(flyCt,numel(lightsOn));
    rowNames = metaTable.Properties.RowNames;
    for iterF = 1:flyCt
        flyID = rowNames{iterF};
        flyRect = str2double(strsplit(flyID,'_'));
        flyRectY = [flyRect(4) flyRect(2)];
        
        posArray = metaTable.filtered_pos{flyID};
        xPos = posArray(:,1);
        yPos = posArray(:,2);
        thetaFilt = posArray(:,3);
        thetaRaw = metaTable.bot_theta{flyID};
        posRaw = metaTable.bot_centroid{flyID};
        
        hPosInit = metaTable.headPos(iterF,:);
        tPosInit = metaTable.tailPos(iterF,:);
        dists = sqrt(sum([(hPosInit(:,1)-tPosInit(:,1)).^2,(hPosInit(:,2)-tPosInit(:,2)).^2],2));
        avgL = mean(dists);
        pixL = 0.11;
        %assuming 2.5mm avg fly length
        %     pix2mm = @(x) x*(2.5/avgL);
        pix2mm = @(x) x*pixL;
        deg2rad = @(x) x.*(pi/180);
        rad2deg = @(x) x.*(180/pi);
        %50 fps recording rate
        frm2sec = @(x) x*50;
        origFilterWin = 9;
        
                    part2cut = (origFilterWin-1)/2+1;
%             dataVel = metaTable.filtered_vel{flyID};
%             dataPos = metaTable.filtered_pos{flyID};
            dataPos = [posRaw thetaRaw];
            dataVel = cat(1,[0 0 0],diff(dataPos));
%             dataVel = metaTable.filtered_vel{flyID};
            
            dataVel(1:part2cut-1,:) = NaN;
            dataVel(end-part2cut:end,:) = NaN;
            dataPos(1:part2cut-1,:) = NaN;
            dataPos(end-part2cut:end,:) = NaN;
            [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
            [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
            relXvel = (pix2mm(relXvel));
            relYvel = (pix2mm(relYvel));
            forwardBack = (relXvel);
%             plot(forwardBack)
%             hold all
            goodndcs = find(abs(forwardBack) < 1);
            forwardBack = interp1(goodndcs,forwardBack(goodndcs),(1:numel(forwardBack)),'linear','extrap');
%             forwardBack(abs(forwardBack) > 2) = 0;
            forwardBack = frm2sec(forwardBack);
            sideways = (relYvel);
%             sideways(abs(sideways) > 2) = 0;
            turning = rad2deg((dataVel(:,3)));
%             sum(abs(turning) >= 30)
            goodndcs = find(abs(turning) < 30);
            turning = interp1(goodndcs,turning(goodndcs),(1:numel(turning)),'linear','extrap');
            turning = frm2sec(turning);
            
%         part2cut = (origFilterWin-1)/2+1;
%                     dataVel = metaTable.filtered_vel{flyID};
%             dataPos = [posRaw thetaRaw];
%         dataVel(1:part2cut-1,:) = NaN;
%         dataVel(end-part2cut:end,:) = NaN;
%         dataPos(1:part2cut-1,:) = NaN;
%         dataPos(end-part2cut:end,:) = NaN;
%         [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
%         [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
%         relXvel = (pix2mm(relXvel));
%         relYvel = (pix2mm(relYvel));
%         forwardBack = (relXvel);
%         forwardBack = frm2sec(forwardBack);
%         sideways = (relYvel);
% %%%%%%%%%%%%%%% distance threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         turning = rad2deg(abs(dataVel(:,3)));
        turning(forwardBack > 60) = 0;
        turning(forwardBack < -30) = 0;
        forwardBack(forwardBack > 60) = 0;
        forwardBack(forwardBack < -30) = 0;
% %%%%%%%%%%%%%%% distance threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         turning(abs(turning) > 15) = 0;
%         turning = frm2sec(turning);
%         ydataA = [ydataA;turning];
        for iterL = 1:numel(lightsOn)
            if iterL > numel(lightsOff)
                continue
            end
            dataBegin = lightsOn(iterL);
            dataEnd = dataBegin+onDur;
%             dataBegin = lightsOn(iterL)+onDur;
%             dataEnd = dataBegin+onDur;
            if dataEnd > numel(forwardBack)
                continue
            end
            
            if ~isempty(strfind(videoID,'testflea_OL0077B'))
                
            elseif ~isempty(strfind(videoID,'pBDPGAL4U'))
                ydataB = [ydataB;relXvel(:)];
            end
            
            if meanPlot == 1
                longMotionCell{iterF,iterL} = mean(forwardBack(dataBegin:dataEnd));
%                 longMotionCell{iterF,iterL} = mean(diff(forwardBack(dataBegin:dataEnd)));
                latMotionCell{iterF,iterL} = mean(sideways(dataBegin:dataEnd));
                turningCell{iterF,iterL} = abs(mean(turning(dataBegin:dataEnd)));
            else
                dataBegin = dataBegin-150;
                dataEnd = dataEnd+150;
%                 dataEnd = round((dataEnd-dataBegin)/2)+dataEnd;
                if dataEnd > numel(forwardBack)
                    continue
                end
                if dataBegin < 1
                    continue
                end
                longMotionCell{iterF,iterL} = forwardBack(dataBegin:dataEnd)';
                latMotionCell{iterF,iterL} = sideways(dataBegin:dataEnd)';
                turningCell{iterF,iterL} = turning(dataBegin:dataEnd)';
            end
            finalPosition{iterF,iterL} = dataPos(dataEnd,1:2);
            %             subNames = metaTable.Properties.RowNames;
            %             exceltable.distance(subNames) = round(forwardBack*10)/10;
            %             exceltable.turning(subNames) = round(turning);
        end
        
    end
    expts2keep = true(flyCt,numel(lightsOn));
    for iterL = 1:numel(lightsOn)
        allPos = cell2mat(finalPosition(:,iterL));
        for iterF = 1:flyCt
            if expts2keep(iterF,iterL) == false
                continue
            end
            thisPos = finalPosition{iterF,iterL};
            distTest = find(pdist2(thisPos,allPos) < 10);
            if numel(distTest) == 1
                continue
            end
            for iterT = distTest(2:end)
                expts2keep(iterT,iterL:end) = false;
            end
        end
    end
    if strcmp(graphStr,'turning')
        dataCell{iterD} = cell2mat(turningCell(expts2keep));
    else
        dataCell{iterD} = cell2mat(longMotionCell(expts2keep));
    end
    strBrks = strfind(videoID,'_');
    if runMode == 1
        label = videoID(strBrks(1)+1:strfind(videoID,'_att')-1);
    else
        label = videoID(strBrks(1)+1:strfind(videoID,'_continuous')-1);
    end
    label = regexprep(label,'_UAS-CsChrimson','');
    dataLabels{iterD} = label;
    if meanPlot == 1
        fillData = NaN;
    else
%         fillData = NaN(dataEnd-dataBegin+1,1);
        fillData = [];
    end
    longMotionCell(~expts2keep) = {fillData};
    turningCell(~expts2keep) = {fillData};
    longMotFullCell{iterD} = longMotionCell;
    turnMotFullCell{iterD} = turningCell;
end

% xdata = rand(numel(ydataA),1);
% plot(xdata,ydataA,'.')
% hold on
% xdata = rand(numel(ydataB),1)+1;
% plot(xdata,ydataB,'.')
%% set labels and data order
excelDir = 'Z:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2_swapped.xlsx';
if runMode == 1
    excelTable = readtable(fullfile(excelDir,excelName),'Sheet','heatmap_bowl','readrownames',true);
    ctrlLine = {'pBDPGAL4U'};
else
    excelTable = readtable(fullfile(excelDir,excelName),'Sheet','new_bowl_data','readrownames',true);
    ctrlLine = {'pBDPGAL4U_w+_norpA_attP18_male_ret+'};
%         '26A03-AD_attP40_attP18_ret+'
%         '54A05-DBD_attP2_attP18_ret+'};
end
labelA = strtrim(excelTable.Properties.RowNames);
labelA = labelA(logical(excelTable.Keepers));
sortdata = excelTable.Order(labelA);
if runMode == 2
    labelColors = excelTable.color(labelA);
    labelB = strtrim(excelTable.table_reference(labelA));
else
    labelB = strtrim(excelTable.CellTypes(labelA));
end
dataCellB = cell(numel(labelA),1);
longMotFullCellB = cell(numel(labelA),1);
turnMotFullCellB = cell(numel(labelA),1);
data2keep = true(numel(labelA),1);
for iterD = 1:numel(labelA)
    dataRef = find(strcmp(dataLabels,labelA{iterD}));
    if isempty(dataRef)
        data2keep(iterD) = false;
        continue
    end
    dataCellB(iterD) = dataCell(dataRef);
    longMotFullCellB(iterD) = longMotFullCell(dataRef);
    turnMotFullCellB(iterD) = turnMotFullCell(dataRef);
end
dataCellB = dataCellB(data2keep);
labelA = labelA(data2keep);
labelB = labelB(data2keep);
if runMode == 2
    labelColors = labelColors(data2keep);
end
longMotFullCellB = longMotFullCellB(data2keep);
turnMotFullCellB = turnMotFullCellB(data2keep);
dataorder = excelTable.Order(labelA);
%% revisions A
if ~exist('hFigR','var')
    hFigR = figure('pos',[2086 600 1650 500],'color','w');
end
clf
saveFig = 1;

hax = zeros(5,1);
for iterA = 1:5
    hax(iterA) = subplot(1,5,iterA);
    set(hax(iterA),'nextplot','add');%,'xcolor','none','ycolor','none','color','none');
end
dataPreFull = [];
dataFirstFull = [];
negCtrlStrings = cell(0,1);
exptStrings = cell(0,1);
h4legA = [];
h4legB = [];
pooledCompare = [];
poolLabel = zeros(4,1);
for iterD = 1:numel(dataorder)
    dataref = dataorder(iterD);
    exptDataCell = longMotFullCellB{iterD};
    yplot = [];
    xplot = [];
    dataPre = [];
    offset = NaN;
    dataMeanPre = NaN(1,5);
    for iterL = 1:5
        plot([0 0]+(iterL-1)*500+150,[dataref-1 dataref+1],'color','k','parent',hax(1))
        plot([0 0]+(iterL-1)*500+201,[dataref-1 dataref+1],'color','r','parent',hax(1))
        data2add = cat(2,exptDataCell{:,iterL});
        if isempty(data2add), continue, end
        if isnan(offset)
            offset = median(mean(data2add(1:150,:),1));
            if ~isnan(offset)
                dataFirstFull = cat(1,dataFirstFull,[mean(data2add(1:150),1)' zeros(150,1)+dataref]);
            end
        end
        yplot = cat(1,yplot,smooth(median(data2add,2))*0.03,NaN(50*10-351,1));
        xplot = cat(1,xplot,(iterL-1)*50*10+(1:50*10)');
        data2add = mean(data2add(1:150,:),1);
        if max(strcmp({'gray'},labelColors{iterD}))
            if iterL == 1
                pooledCompare = cat(1,pooledCompare,[data2add(:) zeros(numel(data2add),1)+1]);
                if poolLabel(1) == 0
                    poolLabel(1) = text(1,-2,'first pre-trial experimental',...
                        'horizontalalignment','right','rotation',90,'parent',hax(5));
                end
            elseif iterL == 5
                pooledCompare = cat(1,pooledCompare,[data2add(:) zeros(numel(data2add),1)+3]);
                if poolLabel(3) == 0
                    poolLabel(3) = text(3,-2,'final pre-trial experimental',...
                        'horizontalalignment','right','rotation',90,'parent',hax(5));
                end
            end
        elseif max(strcmp({'red','green'},labelColors{iterD}))
            if iterL == 1
                pooledCompare = cat(1,pooledCompare,[data2add(:) zeros(numel(data2add),1)+2]);
                if poolLabel(2) == 0
                    poolLabel(2) = text(2,-2,'first pre-trial control',...
                        'horizontalalignment','right','rotation',90,'parent',hax(5));
                end
            elseif iterL == 5
                pooledCompare = cat(1,pooledCompare,[data2add(:) zeros(numel(data2add),1)+4]);
                if poolLabel(4) == 0
                    poolLabel(4) = text(4,-2,'final pre-trial control',...
                        'horizontalalignment','right','rotation',90,'parent',hax(5));
                end
            end
        end
        data2add = data2add-offset;
        dataMeanPre(iterL) = median(data2add);
        dataPreFull = cat(1,dataPreFull,[data2add' zeros(size(data2add,2),1)+dataref]);
        dataPre = cat(1,dataPre,[data2add'*0.03+dataref zeros(size(data2add,2),1)+iterL]);
    end
    if max(strcmp({'gray'},labelColors{iterD}))
        hp = plot((1:5),dataMeanPre,'parent',hax(3));
%         plot((1:5),dataMeanPre,'.','color','k','markersize',18,'parent',hax(3))
        h4legA = cat(1,h4legA,hp);
        exptStrings = cat(1,exptStrings,labelB(iterD));
    elseif max(strcmp({'red','green'},labelColors{iterD}))
        hp = plot((1:5),dataMeanPre,'parent',hax(4));
%         plot((1:5),dataMeanPre,'.','color','k','markersize',18,'parent',hax(4))
        negCtrlStrings = cat(1,negCtrlStrings,labelB(iterD));
        h4legB = cat(1,h4legB,hp);
    end
    if isempty(yplot), continue, end
    yplot(end,:) = NaN;
    
    
%     boxplot(hax(2),dataPre(:,1),dataPre(:,2))
%     plot([0 6],[dataref dataref],'color','k','parent',hax(2))
    
    plot(xplot(:),yplot(:)+dataref,'parent',hax(1),'color','k')
    plot([0 size(yplot,1)],[dataref dataref],'color','k','parent',hax(1))
    text(-50,dataref+0.5,labelB{iterD},'horizontalalignment','right','parent',hax(1))
    
    text(dataref,-20,labelB{iterD},'horizontalalignment','right','rotation',90,'parent',hax(2))
end
timeline = (0:10:10*5);
for iterX = 1:numel(timeline)
    text(timeline(iterX)*50+150,0,num2str(timeline(iterX)),'parent',hax(1),...
        'horizontalalignment','right','color','k')
end
text(median(timeline)*50,-0.5,'time of lights on - seconds','parent',hax(1),...
    'horizontalalignment','center')
hleg = legend(h4legA,exptStrings);
legPos = get(hleg,'pos');
legPos(2) =  0.2;
set(hleg,'pos',legPos)
legend(h4legB,negCtrlStrings)
for iterA = 3:4
    set(hax(iterA),'ylim',[-10 10])
    plot([1 5],[0 0],'color','k','parent',hax(iterA))
    plot([1.15 1 1 1.15],[-5 -5 5 5],'color','k','parent',hax(iterA))
    ylabels = [-5 0 5];
    for iterYL = 1:numel(ylabels)
        text(0.75,ylabels(iterYL),num2str(ylabels(iterYL)),'parent',hax(iterA),...
            'horizontalalignment','center','color','k')
    end
end


xlim = [0.5 8.5];
ylim = [-40 40];
ylabel = [0 20 40];
for iterY = 1:numel(ylabel)
    plot(xlim,ylabel(iterY)+[0 0],'parent',hax(2),'color','k')
    text(0.2,ylabel(iterY),num2str(ylabel(iterY)),'parent',hax(2),...
        'horizontalalignment','center','color','k')
end
text(-0.75,20,'mm/s','parent',hax(2),...
    'horizontalalignment','center','rotation',90)
hbox = boxplot(hax(2),dataFirstFull(:,1),dataFirstFull(:,2),'boxstyle','outline',...
    'colors','kkkkrrrr','notch','on','widths',0.75);
set(hbox(7,:),'marker','.','markeredgecolor','none','markersize',12)
set(hax(2),'ycolor','none','xcolor','none','color','none',...
    'ylim',ylim,'xlim',xlim)
for iterL = 1:4
    pstring = get(poolLabel(iterL),'string');
    dataCt = sum(pooledCompare(:,2) == iterL);
    pstring = cat(2,pstring,' - n = ',num2str(dataCt));
    set(poolLabel(iterL),'string',pstring)
end
pvalFirst = ranksum(pooledCompare(pooledCompare(:,2) == 1,1),pooledCompare(pooledCompare(:,2) == 2,1));
pvalLast = ranksum(pooledCompare(pooledCompare(:,2) == 4,1),pooledCompare(pooledCompare(:,2) == 3,1));
text(1.5,25,['p = ' num2str(pvalFirst,2)],'parent',hax(5),...
    'horizontalalignment','center','rotation',0)
text(3.5,25,['p = ' num2str(pvalLast,2)],'parent',hax(5),...
    'horizontalalignment','center','rotation',0)

hbox = boxplot(hax(5),pooledCompare(:,1),pooledCompare(:,2),'boxstyle','outline',...
    'colors','kkkkrrrr','notch','on','widths',0.75);
set(hbox(7,:),'marker','.','markeredgecolor','none','markersize',12)
ylim = [-30 30];
xlim = [0.5 4.5];
set(hax(5),'ycolor','none','xcolor','none','color','none',...
    'ylim',ylim,'xlim',xlim)
ylabel = [0 15 30];
for iterY = 1:numel(ylabel)
    plot(xlim,ylabel(iterY)+[0 0],'parent',hax(5),'color','k')
    text(0.2,ylabel(iterY),num2str(ylabel(iterY)),'parent',hax(5),...
        'horizontalalignment','center','color','k')
end
text(-0.75,20,'mm/s','parent',hax(5),...
    'horizontalalignment','center','rotation',90)

if saveFig == 1
    pdfName = 'revisions_increased_locomotion_before_lightsOn';
    writePath = fullfile('Z:\Ming_RubinLab',[pdfName '.eps']);
    export_fig(writePath,'-nocrop')
end
%% revisions B
if ~exist('hFigR','var')
    hFigR = figure('pos',[2086 200 1650 850],'color','w');
end
clf
saveFig = 0;

hax = zeros(6,1);
for iterA = 1:6
    hax(iterA) = subplot(2,3,iterA);
%     set(hax(iterA),'nextplot','add','xcolor','none','ycolor','none','color','none');
    set(hax(iterA),'nextplot','add');
end
for iterE = 1:2
    expData = [];
    triData = [];
    flyData = [];
    if iterE == 1
        ylim = [-15 30];
        ylab = [-15 0 15];
    else
        ylim = [0 360];
        ylab = [0 180 360];
    end
    
    for iterD = 1:numel(dataorder)
        dataref = dataorder(iterD);
        if iterE == 1
            exptDataCell = longMotFullCellB{iterD};
        else
            exptDataCell = turnMotFullCellB{iterD};
        end
        for iterL = 1:size(exptDataCell,2)
            data2add = cat(2,exptDataCell{:,iterL});
            if isempty(data2add), continue, end
            if strcmp(labelB{iterD},'LC16')
                flyData = cat(1,flyData,[data2add;(1:numel(data2add))]');
                triData = cat(1,triData,[data2add;zeros(size(data2add))+iterL]');
                axRef = 2+(iterE-1)*3;
                dataCt = sum(~isnan(data2add));
                text(iterL,ylim(1)-range(ylim)*0.03,['trial ' num2str(iterL) ' - n = ' num2str(dataCt)],...
                    'horizontalalignment','right','rotation',90,'parent',hax(axRef))
            end
        end
        data2add = cell2mat(exptDataCell);
        data2add = data2add(:)';
        expData = cat(1,expData,[data2add;zeros(size(data2add))+iterD]');
        axRef = 1+(iterE-1)*3;
        dataCt = sum(~isnan(data2add));
        text(dataref,ylim(1)-range(ylim)*0.03,[labelB{iterD} ' - n = ' num2str(dataCt)],...
            'horizontalalignment','right','rotation',90,'parent',hax(axRef))
    end
    flyData(isnan(flyData(:,1)),:) = [];
    for i = 1:max(flyData(:,2))
        if sum(flyData(:,2) == i) < 4
            flyData(flyData(:,2) == i,1) = NaN;
        end
    end
    flyData(isnan(flyData(:,1)),:) = [];
    axRef = 3+(iterE-1)*3;
    uniqdata = unique(flyData(:,2));
    for i = 1:numel(uniqdata)
        dataCt = sum(flyData(:,2) == uniqdata(i));
        text(i,ylim(1)-range(ylim)*0.03,['fly' num2str(i) ' - n = ' num2str(dataCt)],...
            'horizontalalignment','right','rotation',90,'parent',hax(axRef))
    end
    boxData = {expData,triData,flyData};
    for iA = 1:3
        axRef = iA+(iterE-1)*3;
        hbox = boxplot(hax(axRef),boxData{iA}(:,1),boxData{iA}(:,2),'boxstyle','outline',...
            'colors','k','notch','off','widths',0.75);
        set(hbox(7,:),'marker','.','markeredgecolor','none','markersize',12)
        haxPos = get(hax(axRef),'pos');
        haxPos(2) = haxPos(2)+0.1;
        haxPos(4) = haxPos(4)-0.1;
        set(hax(axRef),'ycolor','none','xcolor','none','color','none',...
            'ylim',ylim,'pos',haxPos)
        for iterY = 1:numel(ylab)
            plot([0.5 max(boxData{iA}(:,2))+0.5],ylab(iterY)+[0 0],'color','k',...
                'parent',hax(axRef))
            text(0,ylab(iterY),num2str(ylab(iterY)),'parent',hax(axRef),...
                'horizontalalignment','right','color','k')
        end
    end
end
if saveFig == 1
    pdfName = 'revisions_individual_versus_population';
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.eps']);
    export_fig(writePath,'-nocrop')
end
%% find control and visualize scatter
colorCode = {'orange','green','blue','red','gray'};
colorOps = [255,127,0
    77,175,74
    55,126,184
    228,26,28
    120 120 120];
saveFig = 0;
figure
% [~,dataorder] = sort(sortdata);
% dataorder = (1:numel(labelA));
if strcmp(graphStr,'turning')
    ylines = [0 180 360];
else
    ylines = [-15 0 15 30];
end
for iterY = 1:numel(ylines)
    plot([0 numel(dataorder)]+0.5,[ylines(iterY) ylines(iterY)],'linewidth',1,'color',[1 1 1]*0.7)
    text(0,ylines(iterY),num2str(ylines(iterY)),'fontsize',14,'horizontalalignment','right')
    hold on
end
yoff = -45;
if ~strcmp(graphStr,'turning')
    yoff = -20;
end
ctrldata = [];
txtColr = 'k';
for iterD = 1:numel(dataorder)
    dataref = dataorder(iterD);
    ydata = dataCellB{iterD};
    ydata(ydata > max(ylines)) = NaN;
    ydata(ydata < min(ylines)) = NaN;
    if max(strcmp(ctrlLine,labelA{iterD}))
        ctrldata = cat(1,ctrldata(:),ydata(:));
    end
end
posThresh = prctile(ctrldata,97.7);
negThresh = prctile(ctrldata,2.3);
for iterD = 1:numel(dataorder)
    dataref = dataorder(iterD);
    ydata = dataCellB{iterD};
    xdata = (rand(numel(ydata),1)-0.5)*0.5+dataref;
    ydata(ydata > max(ylines)) = NaN;
    ydata(ydata < min(ylines)) = NaN;
    
    if ~isempty(strfind(graphStr,'turning'))
        ydata = abs(ydata);
        midDataX = xdata(ydata < posThresh);
        midDataY = ydata(ydata < posThresh);
    else
        midDataX = xdata(ydata < posThresh & ydata > negThresh);
        midDataY = ydata(ydata < posThresh & ydata > negThresh);
    end
    plot([dataref dataref],[yoff*0.95 max(ylines)],'linewidth',1,'color',[1 1 1]*0.7)
    plot(xdata(:),ydata(:),'.','color','k')
    plot(midDataX(:),midDataY(:),'.','color',[0 0 0]+0.5)
    
    
%     labelA{dataref} = regexprep(labelA{dataref},'_male','');
%     labelA{dataref} = regexprep(labelA{dataref},'_attP18','');
%     labelA{dataref} = regexprep(labelA{dataref},'_attP40','');
%     labelA{dataref} = regexprep(labelA{dataref},'_attP2','');
%     text(dataref,xoff,labelA{dataref},'rotation',45,'horizontalalignment','left','interpreter','none')
    if strcmp(graphStr,'turning')
        labely = yoff*1.5;
    else
        labely = yoff;
    end
    if runMode == 2
        txtColr = colorOps(strcmp(colorCode,labelColors{iterD}),:)/255;
    end
    text(dataref,labely,labelB{iterD},'rotation',45,'horizontalalignment','right',...
        'color',txtColr,'interpreter','none','fontsize',14)
    if max(strcmp(ctrlLine,labelA{iterD}))
        ctrldata = cat(1,ctrldata(:),ydata(:));
    end
end
plot([0 numel(dataorder)]+0.5,repmat(posThresh,1,2),'color',[0.8 0 0],'linewidth',1,'linestyle','--')
if ~strcmp(graphStr,'turning')
    plot([0 numel(dataorder)]+0.5,repmat(negThresh,1,2),'color',[0.8 0 0],'linewidth',1,'linestyle','--')
end

ylim = get(gca,'ylim');
ylim(1) = yoff*1.5;
if strcmp(graphStr,'turning')
    text(0,ylim(2)+diff(ylim)*0.1,'Arena','fontsize',18)
    text(-1,median(ylines),'Turning  (deg/sec)','rotation',90,'fontsize',14,'horizontalalignment','center')
else
    text(0,ylim(2)+diff(ylim)*0.1,'Arena','fontsize',18)
    text(-1,median(ylines),'Walking  (mm/sec)','rotation',90,'fontsize',14,'horizontalalignment','center')
end
haxPos = get(gca,'position');
if strcmp(graphStr,'turning')
    if runMode == 2
        growShiftAx = [-0.5,-0.5,0,0.0];%grow width, grow height, shift X, shift Y
    else
        growShiftAx = [-0,-0.4,0,0.0];%grow width, grow height, shift X, shift Y
    end
else
    if runMode == 2
        growShiftAx = [-0.5,-0.35,0,0.0];%grow width, grow height, shift X, shift Y
    else
        growShiftAx = [-0,-0.1,0,0.0];%grow width, grow height, shift X, shift Y
    end
end
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(gca,'ylim',ylim,'box','off','xcolor','w','ycolor','w','pos',haxPos)
set(gcf,'color','w')
if runMode == 2
    set(gcf,'position',[101 278 1750 640])
else
    set(gcf,'position',[101 278 1750 440])
end
if saveFig
    if strcmp(graphStr,'turning')
        pdfName = 'Turning_bowl_assay_jitterPlot.pdf';
    else
        pdfName = 'Longitudinal_motion_bowl_assay_jitterPlot.pdf';
    end
    if runMode == 2
        pdfName = [pdfName(1:end-4) '_LC16controls.pdf'];
    end
    writePath = fullfile('Y:\Ming_RubinLab',pdfName);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName(1:end-4) '.eps']);
    export_fig(writePath,'-nocrop')
end

%% summarize binary results
for iterD = 1:numel(dataorder)
    dataref = dataorder(iterD);
    ydata = dataCellB{dataref};
    testInf = sum(isinf(ydata));
    if testInf > 0, testInf, end
    testNan = sum(isnan(ydata));
    if testNan > 0, testNan, end
    dataCount = numel(ydata);
    if strcmp(graphStr,'turning')
        ydata = sum(ydata > posThresh)/numel(ydata)*100;
    elseif strcmp(graphStr,'progression')
        ydata = sum(ydata > posThresh)/numel(ydata)*100;
    elseif strcmp(graphStr,'regression')
        ydata = sum(ydata < negThresh)/numel(ydata)*100;
    end
    
    
    xdata = iterD;
    plot(xdata(:),ydata(:),'.','markersize',22)
    hold on
    %     text(dataref,0,labelA{dataref},'rotation',45,'horizontalalignment','left')
    text(iterD,60,labelB{dataref},'rotation',90,'horizontalalignment','left')
    
    dataVal = ydata/100;
    if strcmp(graphStr,'turning')
        excelTable.BowlAssay_turning_prct(labelA{dataref}) = dataVal;
    elseif strcmp(graphStr,'progression')
        excelTable.BowlAssay_progressing_prct(labelA{dataref}) = dataVal;
    elseif strcmp(graphStr,'regression')
        excelTable.BowlAssay_regressing_prct(labelA{dataref}) = dataVal;
    end
    excelTable.TotalFlies(labelA{dataref}) = dataCount;
end
set(gca,'ylim',[0 100],'xlim',[0 numel(dataorder)+1])
set(gcf,'position',[101 678 1503 420])
if runMode == 1
%     writetable(excelTable,fullfile(excelDir,excelName),'Sheet','heatmap_bowl','writerownames',true);
else
%     writetable(excelTable,fullfile(excelDir,excelName),'Sheet','new_bowl_data','writerownames',true);
end
