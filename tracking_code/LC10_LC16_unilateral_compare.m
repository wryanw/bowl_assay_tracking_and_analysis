clear all
close all
clc

dataDir = 'Y:\Ming_RubinLab\LC10_quantification_newAvi';
dataList = dir(fullfile(dataDir,'*.mat'));
dataList = {dataList(:).name}';

vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileNames = fileNames(1:end-3);
fileCt = numel(fileNames);
% exceldir = 'D:\rai_tracking\Rai2Ryan_all';
% excelname = 'Stochastic_MDN.xlsx';
exceldir = 'Y:\Ming_RubinLab\LC10_quantification_newAvi';
excelname = 'LC10_unilateral_labeling_050715_Ming.xlsx';
exceltable = readtable(fullfile(exceldir,excelname),'Sheet','Sheet2','ReadRowNames',true);
patternOps = {'/'
    'Both'
    'Left'
    'Left > Right'
    'Left >> Right'
    'Right'
    'Right > Left'
    'Right >>Left'
    'no LC10'};
patternCode = {'NA'
    'B'
    'L'
    'L'
    'L'
    'R'
    'R'
    'R'
    'E'};
exceltable.pattern_Nern = strtrim(exceltable.pattern_Nern);
for iterC = 1:numel(patternOps)
    exceltable.pattern_Nern(strcmp(exceltable.pattern_Nern,patternOps{iterC})) = patternCode(iterC);
end
pattern2flip = ones(numel(exceltable.pattern_Nern),1);
pattern2flip(strcmp(exceltable.pattern_Nern,'L')) = -1;

showIndicatorPlots = 0;
showThetaScroll = 0;
sumPlot = 1;
saveFigs = 0;

plotMode = 2;
if plotMode == 1
    grpList = strtrim(exceltable.pattern_8categories);
    nametag = '_8categories';
elseif plotMode == 2
    grpList = strtrim(exceltable.pattern_4categories);
    nametag = '_4categories';
else
    grpList = strtrim(exceltable.pattern_3categories);
    nametag = '_3categories';
end
grpOps = unique(grpList);
grpOps(strcmp(grpOps,'')) = [];
grpOps(strcmp(grpOps,'Others')) = [];
% grpOps = {'E','B','U'};
exptList = exceltable.Properties.RowNames;
grpCt = numel(grpOps);
dataCell = cell(grpCt,3);
%
for iterD = 1:fileCt
    videoID = fileNames{iterD};
    vidRefPart = videoID(strfind(videoID,'_v')+2:end);
    vidRef = vidRefPart(1:strfind(vidRefPart,'_')-1);
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
    
    onDur = abs(min(lightsOn-lightsOff));
    onDur = round(onDur/2);
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
    rowNames = metaTable.Properties.RowNames;
    for iterF = 1:flyCt
        flyID = rowNames{iterF};
        flyRect = str2double(strsplit(flyID,'_'));
        flyRectY = [flyRect(4) flyRect(2)];
        otherID = rowNames{abs(iterF-3)};
        flyRect = str2double(strsplit(otherID,'_'));
        otherRectY = [flyRect(4) flyRect(2)];
        if max(flyRectY) < min(otherRectY)
            exptName = [vidRef 'a'];
        else
            exptName = [vidRef 'b'];
        end
        exptRef = find(strcmp(exptList,exptName));
        if isempty(exptRef)
            continue
        end
        grpRef = find(strcmp(grpOps,grpList{exptRef}));
        if isempty(grpRef)
            continue
        end
        
        posArray = metaTable.filtered_pos{flyID};
        xPos = posArray(:,1);
        yPos = posArray(:,2);
        thetaFilt = posArray(:,3);
        thetaRaw = metaTable.bot_theta{flyID};
        
        hPosInit = metaTable.headPos(iterF,:);
        tPosInit = metaTable.tailPos(iterF,:);
        dists = sqrt(sum([(hPosInit(:,1)-tPosInit(:,1)).^2,(hPosInit(:,2)-tPosInit(:,2)).^2],2));
        avgL = mean(dists);
        pixL = 0.11;
        %assuming 2.5mm avg fly length
        %     pix2mm = @(x) x*(2.5/avgL);
        pix2mm = @(x) x*pixL;
        %50 fps recording rate
        frm2sec = @(x) x*50;
        origFilterWin = 9;
        part2cut = (origFilterWin-1)/2+1;
        dataVel = metaTable.filtered_vel{flyID};
        dataPos = metaTable.filtered_pos{flyID};
        dataVel(1:part2cut-1,:) = NaN;
        dataVel(end-part2cut:end,:) = NaN;
        dataPos(1:part2cut-1,:) = NaN;
        dataPos(end-part2cut:end,:) = NaN;
        [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
        [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
        relXvel = (pix2mm(relXvel));
        relYvel = (pix2mm(relYvel));
        forwardBack = (relXvel);
        sideways = (relYvel)*pattern2flip(exptRef);
        turning = (dataVel(:,3)/(pi/180))*pattern2flip(exptRef);
        forwardBack(abs(forwardBack) > 2) = 0;
        sideways(abs(sideways) > 2) = 0;
        turning(turning > 15) = 0;
        instanceList = table2cell(exceltable(exptRef,:));
        
        for iterL = 1:numel(lightsOn)
            if iterL > numel(lightsOff)
                continue
            end
            if isnan(instanceList{iterL})
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'NA')
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'j')
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'c')
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'/')
                continue
            end
            dataBegin = lightsOn(iterL);
            dataEnd = dataBegin+onDur;
            if dataEnd > numel(forwardBack)
                continue
            end
            dataCell{grpRef,1} = [dataCell{grpRef,1};sum(forwardBack(dataBegin:dataEnd))];
            dataCell{grpRef,2} = [dataCell{grpRef,2};sum(sideways(dataBegin:dataEnd))];
            dataCell{grpRef,3} = [dataCell{grpRef,3};sum(turning(dataBegin:dataEnd))];
        end
    end
end
dataLabels = {'LC10'};
dataPlot = dataCell(3,3);
%%

dataDir = 'Y:\Ming_RubinLab\LC16_unilateral_quantification\all_examples_crop';
dataFiles = dir(fullfile(dataDir,'*filteredData.mat'));
dataFiles = {dataFiles(:).name};
dataLoading = cellfun(@(x) load(fullfile(dataDir,x)),dataFiles,'uniformoutput',false);
dataLoading = cellfun(@(x) x.analyzer_record,dataLoading,'uniformoutput',false);
dataSet = cat(1,dataLoading{:});
dataCt = size(dataSet,1);
initPosPath = fullfile(dataDir,'initial_positions.mat');
load(initPosPath)
exceldir = 'Y:\Ming_RubinLab\LC16_unilateral_quantification';
excelname = 'LC16_unilateral_videos_notes.xlsx';
exceltable = readtable(fullfile(exceldir,excelname),'Sheet','display_video','ReadRowNames',true);
rowNames = exceltable.Properties.RowNames;
% rowNames = cellfun(@(x) x(1:end-4),rowNames,'uniformoutput',false);
exceltable.Properties.RowNames = rowNames;
dataSet = [initPos exceltable dataSet];

patternTable = table2cell(readtable(fullfile(exceldir,excelname),'Sheet','pattern'));
keyTable = table2cell(readtable(fullfile(exceldir,excelname),'Sheet','key'));
patternCt = size(patternTable,1);

dataCell = cell(3,3);
colrs = {'r','g','b'};
for iterP = 1:patternCt
    grpRef = patternTable{iterP,2};
    rowNames = dataSet.Properties.RowNames;
    pTest = cellfun(@(x) ~isempty(strfind(x,['crop_' patternTable{iterP,1}])),rowNames);
    subSet = dataSet(pTest,:);
    subKeep = strcmp(subSet.tracking_notes,'good');
    subSet = subSet(subKeep,:);
    hPosInit = subSet.headPos;
    tPosInit = subSet.tailPos;
    dists = sqrt(sum([(hPosInit(:,1)-tPosInit(:,1)).^2,(hPosInit(:,2)-tPosInit(:,2)).^2],2));
    avgL = mean(dists);
    pixL = 0.11;
    %assuming 2.5mm avg fly length
%     pix2mm = @(x) x*(2.5/avgL);
    pix2mm = @(x) x*pixL;
    %50 fps recording rate
    frm2sec = @(x) x*50;
    
    dataCt = size(subSet,1);
%     return
    sideways = zeros(dataCt,1);
    forwardBack = zeros(dataCt,1);
    turning = zeros(dataCt,1);
    origFilterWin = 9;
    part2cut = (origFilterWin-1)/2+1;
    for iterV = 1:dataCt
        dataVel = subSet.filtered_vel{iterV};
        dataPos = subSet.filtered_pos{iterV};
        dataVel(1:part2cut-1,:) = [];
        dataVel(end-part2cut:end,:) = [];
        dataPos(1:part2cut-1,:) = [];
        dataPos(end-part2cut:end,:) = [];
        [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
        [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
        relXvel = (pix2mm(relXvel));
        relYvel = (pix2mm(relYvel));
        forwardBack(iterV) = sum(relXvel);
        sideways(iterV) = sum(relYvel)*patternTable{iterP,3};
        turning(iterV) = sum((dataVel(:,3)/(pi/180)))*patternTable{iterP,3};
        
%         forwardBack(iterV) = sum(smooth(relXvel,25));
%         sideways(iterV) = sum(smooth(abs(relYvel),45));
%         turning(iterV) = (sum((smooth((dataVel(:,3)/(pi/180)),25))));
    end
%     ydata = turning;
%     xdata = (rand(numel(ydata),1)-0.5)/2+iterP;
%     plot(xdata,ydata,'.','markersize',20,'color',colrs{grpRef})
%     hold on
    
    dataCell{grpRef,1} = [dataCell{grpRef,1};forwardBack];
    dataCell{grpRef,2} = [dataCell{grpRef,2};sideways];
    dataCell{grpRef,3} = [dataCell{grpRef,3};turning];
    subNames = subSet.Properties.RowNames;
    exceltable.distance(subNames) = round(forwardBack*10)/10;
    exceltable.turning(subNames) = round(turning);
end
dataLabels(2) = {'LC16'};
dataPlot(2) = dataCell(1,3);

%%
saveFig = 0;
datarefs = [1 2];
figure
backC = [1 1 1];
set(gcf,'units','normalized')
set(gcf,'position',[  0.1620    0.4150    0.4974    0.4750])
figPos = get(gcf,'position');
set(gcf,'color',backC,'position',figPos)
% hax = zeros(1,dataCt);
hax(1) = axes;
haxPos = get(hax(1),'position');
growShiftAx = [-0.5,-0.6,0.1,-0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on
xtickpos = [-360 0 360];
ydrop = -0.1;
ytickpos = [1 2 3]-0.5;
fontC = [0 0 0];
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[-0.2 2],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','w','xcolor','w')

plot([0 0],[0 3],'color',[0 0 0],'linewidth',2)
xplot = cat(1,xtickpos,xtickpos,NaN(1,numel(xtickpos)));
yplot = cat(1,zeros(1,numel(xtickpos))+ydrop/2,zeros(1,numel(xtickpos))+ydrop,zeros(1,numel(xtickpos)));
plot(xplot,yplot,'k','linewidth',2,'parent',hax(1))
plot([xtickpos(1) xtickpos(end)],[ydrop ydrop]/2,'k','linewidth',2,'parent',hax(1))

for iterL = 1:2
    fitdata = dataPlot{datarefs(iterL)}*(-1);
    ydata = fitdata;
    xdata = (rand(numel(ydata),1)-0.5)/2+iterL-0.5;
    hbox = boxplot(hax(1),ydata,'positions',iterL-0.5,'boxstyle','outline',...
        'colors','k','notch','on','widths',0.75,'orientation','horizontal');
    set(hbox(7,:),'marker','.','markeredgecolor','k','markersize',12)
    text(xtickpos(1)-(0.5)*range(xtickpos),ytickpos(iterL),['n = ' num2str(numel(ydata))],'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',16);
%     plot(ydata,xdata,'.','color',[0 0 0],'markersize',16)
%     medData = median(fitdata);
%     plot([medData medData],[iterL-1 iterL-.2],'color',[0 0 0],'linewidth',5)
%     err = 1.57*(prctile(fitdata,75)-prctile(fitdata,25))/sqrt(numel(fitdata));
%     plot([medData-err medData+err],[iterL-1 iterL-1],'color',[0 0 0],'linewidth',5)
end
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[-0.2 2],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','w','xcolor','w')
obj = findobj(gcf,'Tag','Upper Adjacent Value');
set(obj,'visible','off')
obj = findobj(gcf,'Tag','Lower Adjacent Value');
set(obj,'visible','off')
obj = findobj(gcf,'Tag','Upper Whisker');
set(obj,'LineStyle','-','LineWidth',2)%,'color','k');
obj = findobj(gcf,'Tag','Lower Whisker');
set(obj,'LineStyle','-','LineWidth',2)%,'color','k');
obj = findobj(gcf,'Tag','Box');
set(obj,'LineWidth',2)
obj = findobj(gcf,'Tag','Median');
set(obj,'LineWidth',2)

xtickval = cellstr(num2str(xtickpos'));
set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) -0.5];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),strtrim(xtickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20,'parent',hax);
end
for iterL = 1:numel(dataLabels)
    text(xtickpos(1)-(0.1)*range(xtickpos),ytickpos(iterL),(dataLabels{datarefs(iterL)}),'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',20);
end
pdfName = 'LC16_LC10_unilateral_compare.pdf';
text(xtickpos(1),2.5,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
if saveFig
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab',pdfName);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName(1:end-4) '.eps']);
    export_fig(writePath,'-nocrop')
end