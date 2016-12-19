clear all
close all
clc
dataOp = 3;
if dataOp == 1
    dataDir = 'Y:\Ming_RubinLab\LC16vCtrl_CsChrimson_redux';
else
    dataDir = 'Y:\Ming_RubinLab\LC16vCtrl_revisions';
end
vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileNames = fileNames(1:end)';
fileCt = numel(fileNames);
%%
meanPlot = 0;

for iterW = 1:2
    if iterW == 1
        graphMode = 1;
    else
        graphMode = 3;
    end
    graphOps = {'turning','progression','regression'};
    graphStr = graphOps{graphMode};
    showIndicatorPlots = 0;
    showThetaScroll = 0;
    dataCell = cell(fileCt,1);
    dataLabels = cell(fileCt,1);
    flyCounts = zeros(fileCt,1);
    ydataA = [];
    ydataB = [];
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
        lightsOn = lightsOn(1:5);%restricts experiment count to 5
        lightsOff = lightsOff(1:5);
        onDur = abs(min(lightsOn-lightsOff));
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
        if dataOp == 2
            lightsOn = lightsOn(1);
            lightsOff = lightsOff(1);
        elseif dataOp == 3
            lightsOn = lightsOn(5);
            lightsOff = lightsOff(5);
        end
        flyCt = size(metaTable,1);
        latMotionCell = cell(flyCt,numel(lightsOn));
        longMotionCell = cell(flyCt,numel(lightsOn));
        turningCell = cell(flyCt,numel(lightsOn));
        finalPosition = cell(flyCt,numel(lightsOn));
        rowNames = metaTable.Properties.RowNames;
        longMotFullCell = cell(flyCt,1);
        turnMotFullCell = cell(flyCt,1);
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
            if strcmp(graphStr,'turning')
                dataVel = metaTable.filtered_vel{flyID};
                dataPos = metaTable.filtered_pos{flyID};
            else
                dataPos = [posRaw thetaRaw];
                dataVel = cat(1,[0 0 0],diff(dataPos));
            end
            dataVel(1:part2cut-1,:) = NaN;
            dataVel(end-part2cut:end,:) = NaN;
            dataPos(1:part2cut-1,:) = NaN;
            dataPos(end-part2cut:end,:) = NaN;
            [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
            [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
            relXvel = (pix2mm(relXvel));
            relYvel = (pix2mm(relYvel));
            forwardBack = (relXvel);
            goodndcs = find(abs(forwardBack) < 0.5);
%             round(numel(goodndcs)/numel(turning)*100)
            forwardBack = interp1(goodndcs,forwardBack(goodndcs),(1:numel(forwardBack)),'linear','extrap');
            forwardBack = frm2sec(forwardBack);
            
            sideways = (relYvel);
            sideways(abs(sideways) > 2) = 0;
            
            turning = rad2deg(dataVel(:,3));
            goodndcs = find(abs(turning) < 30);
%             round(numel(goodndcs)/numel(turning)*100)
            turning = interp1(goodndcs,turning(goodndcs),(1:numel(turning)),'linear','extrap');
            turning = frm2sec(turning);
            
%             round(sum(forwardBack > 60)/numel(forwardBack)*100)
%             round(sum(forwardBack < -30)/numel(forwardBack)*100)
            turning(forwardBack > 60) = 0;
            turning(forwardBack < -30) = 0;
            forwardBack(forwardBack > 60) = 0;
            forwardBack(forwardBack < -30) = 0;
            
            for iterL = 1:numel(lightsOn)
                dataBegin = lightsOn(iterL)-round(onDur);
                dataEnd = lightsOff(iterL)+onDur;
                if dataEnd > numel(forwardBack)
                    continue
                end
                forwardBackS = forwardBack(dataBegin:dataEnd);
                sidewaysS = sideways(dataBegin:dataEnd);
                turningS = turning(dataBegin:dataEnd);
                if meanPlot ~= 1
                    smoothwin = 3;
                    forwardBackS = (forwardBackS');
                    sidewaysS = smooth(sidewaysS,smoothwin);
                    smoothwin = 9;
                    turningS = abs(smooth(turningS,smoothwin));
                end
                turningS = abs(turningS);
                if meanPlot == 1
                    meanFun = @(x) median(x);
                    longMotionCell{iterF,iterL} = [meanFun(forwardBackS(1:onDur))
                        meanFun(forwardBackS(onDur+1:onDur*2))
                        meanFun(forwardBackS(onDur*2+1:end))]';
                    latMotionCell{iterF,iterL} = [meanFun(sidewaysS(1:onDur))
                        meanFun(sidewaysS(onDur+1:onDur*2))
                        meanFun(sidewaysS(onDur*2+1:end))]';
                    turningCell{iterF,iterL} = ([meanFun(turningS(1:onDur))
                        meanFun(turningS(onDur+1:onDur*2))
                        meanFun(turningS(onDur*2+1:end))])';
                else
                    longMotionCell{iterF,iterL} = forwardBackS';
                    latMotionCell{iterF,iterL} = sidewaysS';
                    turningCell{iterF,iterL} = turningS';
                end
                finalPosition{iterF,iterL} = dataPos(dataEnd,1:2);
            end
            longMotFullCell{iterF} = forwardBack;
            turnMotFullCell{iterF} = turning;
        end
        
        expts2keep = true(flyCt,numel(lightsOn));
        for iterL = numel(lightsOn)
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
            data2keep = turningCell(expts2keep);
        else
            data2keep = longMotionCell(expts2keep);
        end
        data2keep = data2keep(:);
        data2keep = cat(1,data2keep{:});
        dataCell{iterD} = (data2keep);
        strBrks = strfind(videoID,'_');
        dataLabels{iterD} = videoID(strBrks(1)+1:strBrks(2)-1);
        flyCounts(iterD) = flyCt;
    end
    %
    % xdata = rand(numel(ydataA),1);
    % plot(xdata,ydataA,'.')
    % hold on
    % xdata = rand(numel(ydataB),1)+1;
    % plot(xdata,ydataB,'.')
    % set labels and data order
    if ~exist('excelTable','var')
        excelDir = 'Y:\Ming_RubinLab';
        excelName = 'Table1_v5_WRW_MW_v2.xlsx';
        excelTable = readtable(fullfile(excelDir,excelName),'Sheet','heatmap_bowl');
        lineNames = strtrim(excelTable.Driver);
        lineLabels = strtrim(excelTable.Driver);
        cellTypes = strtrim(excelTable.CellTypes);
    end
    ctrlLine = 'pBDPGAL4U';
    sortdata = zeros(fileCt,1);
    data2keep = true(fileCt,1);
    labelA = cell(fileCt,1);
    labelB = cell(fileCt,1);
    for iterD = 1:fileCt
        labelA{iterD} = lineLabels{cellfun(@(x) ~isempty(strfind(dataLabels{iterD},x)),lineNames)};
        labelB{iterD} = cellTypes{cellfun(@(x) ~isempty(strfind(dataLabels{iterD},x)),lineNames)};
        sortdata(iterD) = excelTable.Order(strcmp(excelTable.Driver,dataLabels{iterD}));
        data2keep(iterD) = logical(excelTable.Keepers(strcmp(excelTable.Driver,dataLabels{iterD})));
    end
%     expt46logical = find(cellfun(@(x) ~isempty(strfind(x,'testflea_OL0046B')),fileNames));
%     dataCell(expt46logical(1)) = {cat(1,dataCell{expt46logical})};
%     dataCell(expt46logical(2)) = [];
%     labelA(expt46logical(2)) = [];
%     labelB(expt46logical(2)) = [];
%     sortdata(expt46logical(2)) = [];
%     data2keep(expt46logical(2)) = [];
%     flyCounts(expt46logical(2)) = [];
    if iterW == 1
        turningDataCell = dataCell;
    else
        backingDataCell = dataCell;
    end
    plot(smooth(mean(cat(1,longMotFullCell{:}))))
    hold on
end
%%
close all
ydrop = -0.1;
figure
set(gcf,'pos',[2062 264 1566 633],'color','w')
hax = axes;
haxPosA = get(hax(1),'position');
growShiftAx = [-0.65,-0.3,-0.2,0];%grow width, grow height, shift X, shift Y
haxPosA = [haxPosA(1)-haxPosA(3)*growShiftAx(1)/2,haxPosA(2)-haxPosA(4)*growShiftAx(2)/2,...
    haxPosA(3)*(1+growShiftAx(1)),haxPosA(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(1),'position',haxPosA,'ylim',[ydrop 1],'xlim',[0 8],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')

hax(2) = axes;
haxPosB = get(hax(2),'position');
growShiftAx = [-0.65,-0.3,0.25,0];%grow width, grow height, shift X, shift Y
haxPosB = [haxPosB(1)-haxPosB(3)*growShiftAx(1)/2,haxPosB(2)-haxPosB(4)*growShiftAx(2)/2,...
    haxPosB(3)*(1+growShiftAx(1)),haxPosB(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(2),'position',haxPosB,'ylim',[ydrop 1],'xlim',[0 8],'xtick',[],...
    'xticklabel',[],'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')

rangeB = 360;
hplot = zeros(2,1);
for iterD = 1:2
    
    for iterP = 1:2
        if iterP == 1
            plotY = backingDataCell{iterD};
            if dataOp == 1
                rangeA = 30;
            else
                rangeA = 60;
            end
            plotY = (plotY)/rangeA+0.5;
            
        else
            plotY = turningDataCell{iterD};
            plotY = (plotY)/rangeB;
        end
        if iterD == 1
            colr = 'r';
        else
            colr = 'k';
        end
        if dataOp == 1
            notchStr = 'on';
        else
            notchStr = 'off';
        end
        hbox = boxplot(hax(iterP),plotY,'positions',[1 3.5 6]+(iterD-1),'boxstyle','outline',...
            'colors',colr,'notch',notchStr,'widths',0.75);
        set(hbox(7,:),'marker','.','markeredgecolor',colr,'markersize',12)
        if iterD == 1
            text(1.5,-0.12,{'Lights','Off'},'rotation',0,'horizontalalignment','center','fontsize',12,'parent',hax(iterP))
            text(4,-0.12,{'Lights','On'},'rotation',0,'horizontalalignment','center','fontsize',12,'parent',hax(iterP))
            text(6.5,-0.12,{'Lights','Off'},'rotation',0,'horizontalalignment','center','fontsize',12,'parent',hax(iterP))
        end
    end
end
set(hax(1),'position',haxPosA,'ylim',[ydrop 1],'xlim',[0 8],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')

set(hax(2),'position',haxPosB,'ylim',[ydrop 1],'xlim',[0 8],'xtick',[],...
    'xticklabel',[],'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')

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

yticks = linspace(0,1,7);
botlabel = round(linspace(-rangeA/2,rangeA/2,7));
yticklabel = cellstr(num2str(botlabel(1:end)'));
skipvals = [1 3 5 7];
yticklabel(skipvals) = repmat({''},numel(skipvals),1);
for iterT = 1:numel(yticks)
    text(-0.5,yticks(iterT),strtrim(yticklabel{iterT}),'rotation',0,...
        'horizontalalignment','right','fontsize',14,'parent',hax(1))
end
yplot = cat(1,yticks,yticks,NaN(1,numel(yticks)));
xplot = cat(1,zeros(1,numel(yticks))-0.10,zeros(1,numel(yticks))+0.1,zeros(1,numel(yticks)));
plot(xplot,yplot,'k','linewidth',2,'parent',hax(1))
plot([0.1 0.1],[yticks(1) yticks(end)],'k','linewidth',2,'parent',hax(1))

yticks = linspace(0,1,5);
toplabel = (linspace(0,rangeB,5));
yticklabel = cellstr(num2str(toplabel(:)));
skipvals = [2 4];
yticklabel(skipvals) = repmat({''},numel(skipvals),1);
for iterT = 1:numel(yticks)
    text(-0.5,yticks(iterT),strtrim(yticklabel{iterT}),'rotation',0,...
        'horizontalalignment','right','fontsize',14,'parent',hax(2))
end
yplot = cat(1,yticks,yticks,NaN(1,numel(yticks)));
xplot = cat(1,zeros(1,numel(yticks))-0.10,zeros(1,numel(yticks))+0.1,zeros(1,numel(yticks)));
plot(xplot,yplot,'k','linewidth',2,'parent',hax(2))
plot([0.1 0.1],[yticks(1) yticks(end)],'k','linewidth',2,'parent',hax(2))

plot([0 8],[0.5 0.5],'k','linewidth',1,'parent',hax(1))

text(-2.5,0.5,{'Angular Velocity';'(deg/sec)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(2))
text(-2.5,0.5,{'Speed';'(mm/sec)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(1))
text(4,1.1,'Turning','rotation',0,'horizontalalignment','center','fontsize',14,'parent',hax(2))
text(4,1.1,'Walking','rotation',0,'horizontalalignment','center','fontsize',14,'parent',hax(1))

% hleg = legend(hplot,labelB);
% hlegpos = get(hleg,'pos');
% hlegpos(1:2) = [0.85 0.9];
% set(hleg,'pos',hlegpos)
if dataOp == 1
    pdfName = 'bilateral_LC16_backup_and_turn_quantification_summary';
elseif dataOp == 2
    pdfName = 'LC16vCtrl_revisions_firstTrial_summary';
else
    pdfName = 'LC16vCtrl_revisions_lastTrial_summary';
end
%%
plotMode = 3;
close all
ydrop = -0.2;
figure
set(gcf,'pos',[2137 264 1191 633],'color','w')
hax = axes;
haxPosA = get(hax(1),'position');
growShiftAx = [-0.6,-0.6,0.1,0.2];%grow width, grow height, shift X, shift Y
haxPosA = [haxPosA(1)-haxPosA(3)*growShiftAx(1)/2,haxPosA(2)-haxPosA(4)*growShiftAx(2)/2,...
    haxPosA(3)*(1+growShiftAx(1)),haxPosA(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(1),'position',haxPosA,'ylim',[ydrop 1],'xlim',[-3.5 150],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,...
    'color','none','xcolor','none','ycolor','none')

hax(2) = axes;
haxPosB = get(hax(2),'position');
growShiftAx = [-0.6,-0.6,0.1,-0.2];%grow width, grow height, shift X, shift Y
haxPosB = [haxPosB(1)-haxPosB(3)*growShiftAx(1)/2,haxPosB(2)-haxPosB(4)*growShiftAx(2)/2,...
    haxPosB(3)*(1+growShiftAx(1)),haxPosB(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
ydrop = -0.1;
set(hax(2),'position',haxPosB,'ylim',[ydrop 1],'xlim',[-3 150],'xtick',[],...
    'xticklabel',[],'nextplot','add','box','off','fontsize',14,...
    'color','none','xcolor','none','ycolor','none')
xticks = (0:25:150);
xticklabel = cellstr(num2str((xticks'-50)/50));
if plotMode == 3
    for iterT = 1:numel(xticks)
        text(xticks(iterT),-0.2,strtrim(xticklabel{iterT}),'rotation',0,'horizontalalignment','center','fontsize',14)
    end
    
    xplot = cat(1,xticks,xticks,NaN(1,numel(xticks)));
    yplot = cat(1,zeros(1,numel(xticks))+ydrop/2,zeros(1,numel(xticks))+ydrop,zeros(1,numel(xticks)));
    
    plot(xplot,yplot,'k','linewidth',2,'parent',hax(2))
    plot([0 150],[ydrop ydrop]/2,'k','linewidth',2,'parent',hax(2))
end
rangeB = 360;
xshift = 0;
smoothLeg = 2;
for iterD = 1:2
    for iterP = 1:2
        if plotMode == 1 && iterD == 2
            continue
        elseif plotMode == 2 && iterD == 1
            continue
        end
        if iterP == 1
            plotY = backingDataCell{iterD};
            if dataOp == 1
                rangeA = 30;
            else
                rangeA = 60;
            end
            plotY = (plotY)/rangeA+0.5;
        else
            plotY = turningDataCell{iterD};
            plotY = (plotY)/rangeB;
        end
        plotYm = median(plotY,1);
        if smoothLeg > 0
            plotYm = smooth(plotYm,smoothLeg*2+1);
            plotYm(1:smoothLeg) = NaN;
            plotYm(end-smoothLeg:end) = NaN;
        end
        plotXm = (1:numel(plotYm))+xshift;
        if iterD == 1
            colr = 'r';
        else
            if plotMode == 3
                colr = 'k';
            else
                colr = 'g';
            end
        end
        hplot(iterD) = plot(plotXm,plotYm,'color',colr,'linewidth',5,'parent',hax(iterP));
        
        q1 = prctile(plotY,25);  q3 = prctile(plotY,75);
        err = 1.57*(q3-q1)/sqrt(size(plotY,1));
        if smoothLeg > 0
            q1 = smooth(q1,smoothLeg*2+1)';
            q1(1:smoothLeg) = [];
            q1(end-smoothLeg:end) = [];
            q3 = smooth(q3,smoothLeg*2+1)';
            q3(1:smoothLeg) = [];
            q3(end-smoothLeg:end) = [];
        end
        
%         q1 = plotYm-std(plotY).^2;  q3 = plotYm+std(plotY).^2;
        plotY = [q1 fliplr(q3) q1(1)];
%         plotY = [plotYm+err;plotYm-err;NaN(size(q1))];
        plotX = [(1:numel(q1)) fliplr(1:numel(q1)) 1]+smoothLeg;
        plot(plotX',plotY','color',colr,'linewidth',1,'parent',hax(iterP),'linestyle','-')
        if plotMode ~= 3
            plot([1 1 NaN 150 150],[q1(1) q3(1) NaN q1(150) q3(150)],'color',colr,'linewidth',1,'parent',hax(iterP))
        end
    end
end
%
yticks = linspace(0,1,7);
botlabel = round(linspace(-rangeA/2,rangeA/2,7));
yticklabel = cellstr(num2str(botlabel(1:end)'));
skipvals = [1 3 5 7];
yticklabel(skipvals) = repmat({''},numel(skipvals),1);
if plotMode == 3
    for iterT = 1:numel(yticks)
        text(-5,yticks(iterT),strtrim(yticklabel{iterT}),'rotation',0,...
            'horizontalalignment','right','fontsize',14,'parent',hax(1))
    end
    yplot = cat(1,yticks,yticks,NaN(1,numel(yticks)));
    xplot = cat(1,zeros(1,numel(yticks))-10,zeros(1,numel(yticks))-1,zeros(1,numel(yticks)));
    
    plot(xplot,yplot,'k','linewidth',2,'parent',hax(1))
    plot([-1 -1],[yticks(1) yticks(end)],'k','linewidth',2,'parent',hax(1))
    
    yticks = linspace(0,1,5);
    toplabel = (linspace(0,rangeB,5));
    yticklabel = cellstr(num2str(toplabel(:)));
    skipvals = [2 4];
    yticklabel(skipvals) = repmat({''},numel(skipvals),1);
    for iterT = 1:numel(yticks)
        text(-5,yticks(iterT),strtrim(yticklabel{iterT}),'rotation',0,...
            'horizontalalignment','right','fontsize',14,'parent',hax(2))
    end
    yplot = cat(1,yticks,yticks,NaN(1,numel(yticks)));
    xplot = cat(1,zeros(1,numel(yticks))-10,zeros(1,numel(yticks))-1,zeros(1,numel(yticks)));
    
    plot(xplot,yplot,'k','linewidth',2,'parent',hax(2))
    plot([-1 -1],[yticks(1) yticks(end)],'k','linewidth',2,'parent',hax(2))
    plot([0 150],[0.5 0.5],'k','linewidth',2,'parent',hax(1))
    
    text(-40,0.5,{'Angular Velocity';'(deg/sec)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(2))
    text(-40,0.5,{'Speed';'(mm/sec)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(1))
    text(75,-0.4,'Time (sec)','rotation',0,'horizontalalignment','center','fontsize',14,'parent',hax(2))
    legend(hplot,labelB)
end
trialInfo = {'Trial count'
    ['LC16: ' num2str(size(backingDataCell{1},1))]
    ['Ctrl: ' num2str(size(backingDataCell{2},1))]};
text(-155,0.5,trialInfo,'rotation',0,'horizontalalignment','left','fontsize',14,'parent',hax(2))
flyInfo = {'Fly count'
    ['LC16: ' num2str(flyCounts(1))]
    ['Ctrl: ' num2str(flyCounts(2))]};
text(-155,1.5,flyInfo,'rotation',0,'horizontalalignment','left','fontsize',14,'parent',hax(2))
if plotMode == 1
    pdfName = 'bilateral_LC16_backup_and_turn_quantification_timeseriesA';
elseif plotMode == 2
    pdfName = 'bilateral_LC16_backup_and_turn_quantification_timeseriesB';
else
    pdfName = 'bilateral_LC16_backup_and_turn_quantification_timeseriesC';
end
if dataOp == 2
    pdfName = 'LC16vCtrl_revisions_firstTrial';
elseif dataOp == 3
    pdfName = 'LC16vCtrl_revisions_lastTrial';
end
%%
% writePath = fullfile('Y:\Ming_RubinLab\LC16_bilateral_activation_figures',[pdfName '.pdf']);
% export_fig(writePath,'-nocrop')
writePath = fullfile('Y:\Ming_RubinLab\LC16_bilateral_activation_figures',[pdfName '.eps']);
export_fig(writePath,'-nocrop')