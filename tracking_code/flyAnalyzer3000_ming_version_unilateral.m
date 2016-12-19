% clear all
% close all
% clc
op_sys = system_dependent('getos');
if strfind(op_sys,'Microsoft Windows 7')
    dm11Dir = [filesep filesep 'dm11' filesep 'cardlab'];
else
    dm11Dir = [filesep 'Volumes' filesep 'cardlab'];
end
dataDir = fullfile(dm11Dir,'CARD LAB SHARE','Ming_temp','unilateral_LC16','all_examples_crop');
dataFiles = dir(fullfile(dataDir,'*filteredData.mat'));
dataFiles = {dataFiles(:).name};
dataLoading = cellfun(@(x) load(fullfile(dataDir,x)),dataFiles,'uniformoutput',false);
dataLoading = cellfun(@(x) x.analyzer_record,dataLoading,'uniformoutput',false);
dataSet = cat(1,dataLoading{:});
dataCt = size(dataSet,1);
initPosPath = fullfile(dataDir,'initial_positions.mat');
load(initPosPath)
exceldir = 'Z:\CARD LAB SHARE\Ming_temp\unilateral_LC16';
excelname = 'LC16_unilateral_videos_notes.xlsx';
exceltable = readtable(fullfile(exceldir,excelname),'Sheet','display_video','ReadRowNames',true);
rowNames = exceltable.Properties.RowNames;
% rowNames = cellfun(@(x) x(1:end-4),rowNames,'uniformoutput',false);
exceltable.Properties.RowNames = rowNames;
dataSet = [initPos exceltable dataSet];

patternTable = table2cell(readtable(fullfile(exceldir,excelname),'Sheet','pattern'));
keyTable = table2cell(readtable(fullfile(exceldir,excelname),'Sheet','key'));
patternCt = size(patternTable,1);
%%
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
    ydata = turning;
    xdata = (rand(numel(ydata),1)-0.5)/2+iterP;
    plot(xdata,ydata,'.','markersize',20,'color',colrs{grpRef})
    hold on
    
    dataCell{grpRef,1} = [dataCell{grpRef,1};forwardBack];
    dataCell{grpRef,2} = [dataCell{grpRef,2};sideways];
    dataCell{grpRef,3} = [dataCell{grpRef,3};turning];
    subNames = subSet.Properties.RowNames;
    exceltable.distance(subNames) = round(forwardBack*10)/10;
    exceltable.turning(subNames) = round(turning);
end
% writetable(exceltable,fullfile(exceldir,excelname),'Sheet','display_video','WriteRowNames',true);
%%

funCell = {@(x) x(x > 0),@(x) x(x < 0)};
funrefs = [1 2 1 2 1 2];
datarefs = [2 2 1 1 3 3];
dataCt = 6;
px = cell(dataCt,1);
py = cell(dataCt,1);
eYall = zeros(3,dataCt); mYall = zeros(3,dataCt);

for iterL = 1:dataCt
    
    ydata = dataCell{datarefs(iterL),3};
    ydata = abs(funCell{funrefs(iterL)}(ydata));
    x = (rand(numel(ydata),1)-0.5)/2+iterL;
    px{iterL} = x;
    py{iterL} = ydata;
    
    %standard error of the median
    q1 = prctile(ydata,25); q2 = prctile(ydata,50);
    q3 = prctile(ydata,75); n = numel(ydata);
    eBar = ([q2 - 1.57*(q3 - q1)/sqrt(n);q2 + 1.57*(q3 - q1)/sqrt(n)]);
    mBar = ([q2;q2]);
    
    eYall(:,iterL) = [eBar' NaN]';
    mYall(:,iterL) = [mBar' NaN]';
end
figure
backC = [1 1 1];
set(gcf,'units','normalized')
set(gcf,'position',[  0.1620    0.4150    0.4974    0.4750])
figPos = get(gcf,'position');
set(gcf,'color',backC,'position',figPos)
hax = zeros(1,dataCt);
hax(1) = axes;
haxPos = get(hax(1),'position');
growShiftAx = [-0.2,-0.1,0.1,0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on
x = cat(1,px{:});
ydata = cat(1,py{:});
plot(x,ydata,'.','color','k','markersize',18)
% hold on

% eYall(eYall < 0) = 0; eYall(eYall > 1) = 1;
medianX = repmat((1:dataCt),3,1);
medianX = medianX+repmat([-.1;.1;0],1,dataCt);
errorX = repmat((1:dataCt),3,1);
% plot([medianX;errorX],[mYall;eYall],'color','k','linewidth',3);

ytickpos = [0 180 360];
xtickpos = [1 2 3 4 5 6];
fontC = [0 0 0];
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[ytickpos(1) ytickpos(end)],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
xtickval = {'Left','Right','Contra','Ipsi','Left','Right'};

set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.07];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),(xtickval{iterL}),'rotation',20,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end
ytickval = cellstr(num2str(ytickpos'));
for iterL = 1:numel(ytickval)
    text(0,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end

pdfName = 'LC16_unilateral_activation_turning.pdf';
text(0,400-5,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
%
addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
writePath = fullfile(exceldir,pdfName);
% export_fig(writePath)

%%
plot([0.6 3+0.5],[0 0],'color','k','linewidth',1)
hold on
datarefs = [2 1 3];
for iterL = 1:3
    ydata = dataCell{datarefs(iterL),1};
    xdata = (rand(numel(ydata),1)-0.5)/2+iterL;
    plot(xdata,ydata,'.','color','k','markersize',18)
end
set(gcf,'units','normalized')
set(gcf,'color',backC)
set(gcf,'position',[0.3932    0.4433    0.2521    0.4708])
ytickpos = [-20 0 20];
xtickpos = [1 2 3];
fontC = [0 0 0];
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[ytickpos(1) ytickpos(end)],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
xtickval = {'Bilateral','Unilateral','Empty'};

set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.07];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),(xtickval{iterL}),'rotation',20,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end
ytickval = cellstr(num2str(ytickpos'));
for iterL = 1:numel(ytickval)
    text(0,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end


pdfName = 'LC16_unilateral_activation_backing.pdf';
text(-0.75,23,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
%
addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
writePath = fullfile(exceldir,pdfName);
% export_fig(writePath)


%%
return
xshift = 0;
xspace = 30;
xtickcell = cell(1,2);
xlinebase = (0:25:150);
ylines = zeros(size(xlinebase));
ylines = [ylines;ylines+2;ylines*NaN];
colrs = {'b','g','r'};
for iterD = 1:2
    dists = sqrt(sum([(hPosInit(:,1)-tPosInit(:,1)).^2,(hPosInit(:,2)-tPosInit(:,2)).^2],2));
    avgL = mean(dists);
    %assuming 2.5mm avg fly length
    pix2mm = @(x) x*(2.5/avgL);
    %50 fps recording rate
    frm2sec = @(x) x*50;
    
    dataFiles = {dataFiles(:).name};
    dataLoading = cellfun(@(x) load(fullfile(dataDir,x)),dataFiles,'uniformoutput',false);
    dataSet = cat(1,dataLoading{:});
    dataCt = numel(dataSet);
%     return
    sideways = cell(dataCt,1);
    forwardBack = cell(dataCt,1);
    turning = cell(dataCt,1);
    origFilterWin = 19;
    part2cut = (origFilterWin-1)/2;
    for iterV = 1:dataCt
        dataAccel = dataSet(iterV).filtered_accel;
        dataVel = dataSet(iterV).filtered_vel;
        dataPos = dataSet(iterV).filtered_pos;
        dataVel(1:part2cut,:) = [];
        dataVel(end-part2cut+1:end,:) = [];
        dataPos(1:part2cut,:) = [];
        dataPos(end-part2cut+1:end,:) = [];
        
        [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
        [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
        relXvel = frm2sec(pix2mm(relXvel));
        relYvel = frm2sec(pix2mm(relYvel));
        forwardBack{iterV} = smooth(relXvel,25)';
        %     sideways{iterV} = smooth(abs(relYvel),45)';
        turning{iterV} = frm2sec(smooth(abs(dataVel(:,3)),25)');
        %     turning{iterV} = abs(dataVel(:,3))';
    end
    %
    
%     xlines = repmat(xlinebase,3,1)+xshift;
    xlines = repmat(xlinebase,3,1);
    xtickcell{iterD} = unique(xlines);
    plot(xlines,ylines,'color','k','linewidth',2)
    hold on
%     plot([0 150]+xshift,[0.5 0.5],'color',[.5 .5 .5],'linewidth',2)
%     plot([0 150]+xshift,[1 1],'color',[.5 .5 .5],'linewidth',2)
    plot([0 150],[0.5 0.5],'color',[.5 .5 .5],'linewidth',2)
    plot([0 150],[1 1],'color',[.5 .5 .5],'linewidth',2)
    
    plotY = cat(1,forwardBack{:});
    % plotY = cat(1,sideways{:});
    % rangeA = range(plotY(:))
    rangeA = 120;
    plotY = (plotY)/rangeA+0.5;
    plotYm = mean(plotY,1);
    plotXm = (1:numel(plotYm))+part2cut+xshift;
    colr = colrs{iterD};
    plot(plotXm,plotYm,'color',colr,'linewidth',5)
    err = std(plotY);
    plotY = [err;-err]+[plotYm;plotYm];
    plotY = [err;-err;NaN(size(err))]+[plotYm;plotYm;plotYm];
    plotX = repmat((1:size(plotY,2)),size(plotY,1),1)+part2cut+xshift;
%     plot(plotX',plotY','color',colr,'linewidth',2)
    plot(plotX(:),plotY(:),'color',colr,'linewidth',1)
    
    plotYb = cat(1,turning{:});
    % rangeB = range(plotYb(:))
    rangeB = pi*3;
    plotYb = (plotYb)/rangeB+1;
    plotYmB = mean(plotYb,1);
    plot(plotXm,plotYmB,'color',colr,'linewidth',5)
    err = std(plotYb);
%     plotYb = [err;-err]+[plotYmB;plotYmB];
    plotYb = [err;-err;NaN(size(err))]+[plotYmB;plotYmB;plotYmB];
    plotX = repmat((1:size(plotYb,2)),size(plotYb,1),1)+part2cut+xshift;
%     plot(plotX',plotYb','color',colr,'linewidth',2)
    plot(plotX(:),plotYb(:),'color',colr,'linewidth',1)
    xshift = 0.5;
%     xshift = max(xlines(:))+xshift+xspace;
end
%
oldyticks = get(gca,'ytick');
newyticks = linspace(min(oldyticks),max(oldyticks),13);
set(gca,'ytick',newyticks)
botlabel = round(linspace(-rangeA/2,rangeA/2,7));
toplabel = (linspace(0,rangeB,7));
yticklabel = cellstr(num2str([botlabel(1:end-1)';toplabel(:)]));
yticklabel(1:end) = {'','-40','','0','','40','0','','1','','2','','3'}';
set(gca,'yticklabel',yticklabel)

% xticks = cat(1,xtickcell{:});
xticks = xtickcell{1};
set(gca,'xtick',xticks)
xticklabel = repmat(cellstr(num2str(xlinebase(:))),1,1);
set(gca,'xticklabel',xticklabel)

text(-20,1.5,'Angular Velocity (times  \pi)','rotation',90,'horizontalalignment','center','fontsize',14)
text(-20,0.5,'Speed (mm/sec)','rotation',90,'horizontalalignment','center','fontsize',14)
%
% set(gca,'ticklength',[0 0])
set(gca,'box','off','fontsize',14)
%
text(-10,0.2,'backward','rotation',90,'horizontalalignment','center','fontsize',14)
text(-10,0.8,'forward','rotation',90,'horizontalalignment','center','fontsize',14)
set(gcf,'color','w','position',[2137 264 1191 633])
%%
addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-5735e6d')
pdfName = 'matlab_figure_backup_and_turn_quantification.pdf';
writePath = fullfile('Z:\CARD LAB SHARE\Ming_temp',pdfName);
export_fig(writePath)