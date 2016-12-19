clearvars -except unilateralTable

op_sys = system_dependent('getos');
if strfind(op_sys,'Microsoft Windows 7')
    dm11Dir = [filesep filesep 'dm11' filesep 'cardlab'];
else
    dm11Dir = [filesep 'Volumes' filesep 'cardlab'];
end
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
flyCtTaly = zeros(size(keyTable,1),1);
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
    rad2deg = @(x) x.*(180/pi);
    pixL = 0.11;
    %assuming 2.5mm avg fly length
    pix2mm = @(x) x*pixL;
    %50 fps recording rate
    frm2sec = @(x) x*50;
    
    dataCt = size(subSet,1);
%     return
    sidewaysB = zeros(dataCt,1);
    forwardBackB = zeros(dataCt,1);
    turningB = zeros(dataCt,1);
    origFilterWin = 9;
    part2cut = (origFilterWin-1)/2+1;
    for iterV = 1:dataCt
        
        posArray = subSet.filtered_pos{iterV};
        xPos = posArray(:,1);
        yPos = posArray(:,2);
        thetaFilt = posArray(:,3);
%         thetaRaw = subSet.bot_theta{iterV};
%         posRaw = subSet.bot_centroid{iterV};
        
%         dataPos = [posRaw thetaRaw];
%         dataVel = cat(1,[0 0 0],diff(dataPos));
        dataVel = subSet.filtered_vel{iterV};
        dataPos = subSet.filtered_pos{iterV};
        dataVel(1:part2cut-1,:) = NaN;
        dataVel(end-part2cut:end,:) = NaN;
        dataPos(1:part2cut-1,:) = NaN;
        dataPos(end-part2cut:end,:) = NaN;
        [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
        [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
        relXvel = (pix2mm(relXvel));
        relYvel = (pix2mm(relYvel));
        
        forwardBack = (relXvel);
        goodndcs = find(abs(forwardBack) < 1);
        forwardBack = interp1(goodndcs,forwardBack(goodndcs),(1:numel(forwardBack)),'linear','extrap');
        forwardBack = frm2sec(forwardBack);
        
        sideways = (relYvel)*patternTable{iterP,3};
        sideways(abs(sideways) > 2) = 0;
        
        turning = rad2deg(dataVel(:,3))*patternTable{iterP,3};
        goodndcs = find(abs(turning) < 30);
        turning = interp1(goodndcs,turning(goodndcs),(1:numel(turning)),'linear','extrap');
        turning = frm2sec(turning);
        
%         instanceList = table2cell(exceltable(exptRef,:));
        
        turning(forwardBack > 60) = 0;
        turning(forwardBack < -30) = 0;
        forwardBack(forwardBack > 60) = 0;
        forwardBack(forwardBack < -30) = 0;
        
        
        
        
%         dataVel = subSet.filtered_vel{iterV};
%         dataPos = subSet.filtered_pos{iterV};
%         dataVel(1:part2cut-1,:) = [];
%         dataVel(end-part2cut:end,:) = [];
%         dataPos(1:part2cut-1,:) = [];
%         dataPos(end-part2cut:end,:) = [];
%         [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
%         [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
%         relXvel = (pix2mm(relXvel));
%         relYvel = (pix2mm(relYvel));
        forwardBackB(iterV) = mean(forwardBack);
        sidewaysB(iterV) = mean(sideways)*(-1);
        turningB(iterV) = mean(turning)*(-1);
        
%         forwardBack(iterV) = sum(smooth(relXvel,25));
%         sideways(iterV) = sum(smooth(abs(relYvel),45));
%         turning(iterV) = (sum((smooth((dataVel(:,3)/(pi/180)),25))));
    end
%     ydata = turningB;
%     xdata = (rand(numel(ydata),1)-0.5)/2+iterP;
%     plot(xdata,ydata,'.','markersize',20,'color',colrs{grpRef})
%     hold on
    
    dataCell{grpRef,1} = [dataCell{grpRef,1};forwardBackB];
    dataCell{grpRef,2} = [dataCell{grpRef,2};sidewaysB];
    dataCell{grpRef,3} = [dataCell{grpRef,3};turningB];
    subNames = subSet.Properties.RowNames;
    exceltable.distance(subNames) = round(forwardBackB*10)/10;
    exceltable.turning(subNames) = round(turningB);
    flyCtTaly(grpRef) = flyCtTaly(grpRef)+1;
end
for iterL = 1:size(keyTable,1)
    if ~exist('unilateralTable','var')
        unilateralTable = table(cell(6,1),cell(6,1),'VariableNames',{'CellType','data'},...
            'RowNames',{'LC16_Unilateral','LC16_Bilateral','LC16_Empty',...
            'LC10_Unilateral','LC10_Bilateral','LC10_Empty'});
    end
    tableID = cat(2,'LC16_',keyTable{iterL,2});
    unilateralTable.flyCt(tableID) = flyCtTaly(iterL);
    unilateralTable.data(tableID) = dataCell(iterL,3);
    unilateralTable.CellType(tableID) = {'LC16'};
end
% writetable(exceltable,fullfile(exceldir,excelname),'Sheet','display_video','WriteRowNames',true);
%%
plotOps = {'walking','turning'};
plotChoice = 1;
saveFig = 0;
if plotChoice == 1
    datarefs = [3 1 2];
else
    datarefs = [3 1 2];
end
figure
backC = [1 1 1];
set(gcf,'units','normalized')
set(gcf,'position',[  0.1620    0.4150    0.4974    0.4750])
figPos = get(gcf,'position');
set(gcf,'color',backC,'position',figPos)
% hax = zeros(1,dataCt);
hax(1) = axes;
haxPos = get(hax(1),'position');
growShiftAx = [-0.75,-0.4,0.1,-0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on
if strcmp(plotOps{plotChoice},'walking')
    xtickpos = [-20 0 20];
else
    xtickpos = [-360 0 360];
    
end
ydrop = -0.1;
ytickpos = [1 2 3]-0.5;
fontC = [0 0 0];
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[-0.2 numel(datarefs)],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','w','xcolor','w')

plot([0 0],[0 numel(datarefs)],'color',[0 0 0],'linewidth',2)
xplot = cat(1,xtickpos,xtickpos,NaN(1,numel(xtickpos)));
yplot = cat(1,zeros(1,numel(xtickpos))+ydrop/2,zeros(1,numel(xtickpos))+ydrop,zeros(1,numel(xtickpos)));
plot(xplot,yplot,'k','linewidth',2,'parent',hax(1))
plot([xtickpos(1) xtickpos(end)],[ydrop ydrop]/2,'k','linewidth',2,'parent',hax(1))
fontSz = 12;
kruskData = cell(numel(datarefs),1);
kruskRefs = cell(numel(datarefs),1);
for iterL = 1:numel(datarefs)
    if strcmp(plotOps{plotChoice},'walking')
        fitdata = dataCell{datarefs(iterL),1};
    else
        fitdata = dataCell{datarefs(iterL),3};
    end
    ydata = fitdata;
    xdata = (rand(numel(ydata),1)-0.5)/2+iterL-0.5;
    
    
    hbox = boxplot(hax(1),ydata,'positions',iterL-0.5,'boxstyle','outline',...
        'colors','k','notch','on','widths',0.35,'orientation','horizontal');
    set(hbox(7,:),'marker','.','markeredgecolor','k','markersize',12)
    xdata = get(hbox(5),'xdata');
    ydata = get(hbox(5),'ydata');
    hfill = fill(xdata,ydata,'w');
    delete(hbox)
    ydata = fitdata;
    hbox = boxplot(hax(1),ydata,'positions',iterL-0.5,'boxstyle','outline',...
        'colors','k','notch','on','widths',0.35,'orientation','horizontal');
    set(hbox(7,:),'marker','.','markeredgecolor','k','markersize',12)
    
    text(xtickpos(1)-(1.5)*range(xtickpos),ytickpos(iterL),['n = ' num2str(numel(ydata))],'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',round(fontSz*0.8));
    kruskData{iterL} = fitdata(:);
    kruskRefs{iterL} = zeros(numel(fitdata),1)+iterL;
end

kruskRefs = cat(1,kruskRefs{:});
kruskDataB = cat(1,kruskData{:});
% kruskGroup = (kruskRefs == 1 | kruskRefs == 2);
% [pVal,tabl,stats] = kruskalwallis(kruskDataB(kruskGroup),kruskRefs(kruskGroup),'off')
kruskGroup = (kruskRefs == 2 | kruskRefs == 3);
[pVal,tabl,stats] = kruskalwallis(kruskDataB(kruskGroup),kruskRefs(kruskGroup),'off')

set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[-0.2 numel(datarefs)],'tickdir','out','xtick',xtickpos,...
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
xticklabels = {'Backward',{'Velocity';'(mm/sec)'},'Forward'};
set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) -0.1];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),cat(1,{strtrim(xtickval{iterL})},xticklabels{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',fontSz,'verticalalignment','top');
end
ytickval = {'Unilateral','Bilateral','Empty'};
for iterL = 1:numel(datarefs)
    text(xtickpos(1)-(0.5)*range(xtickpos),ytickpos(iterL),(ytickval{datarefs(iterL)}),'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',fontSz);
end
if strcmp(plotOps{plotChoice},'walking')
    pdfName = 'LC16_unilateral_activation_walking.pdf';
else
    pdfName = 'LC16_unilateral_activation_turning.pdf';
end
text(xtickpos(1),numel(datarefs)+.5,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',fontSz);
if saveFig
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab',pdfName);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName(1:end-4) '.eps']);
    export_fig(writePath,'-nocrop')
end
