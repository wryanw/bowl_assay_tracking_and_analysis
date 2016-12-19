%% unilateral comparrison graph LC16 and LC10
rng('default');
rng(19);
ipsiVlr = 2;
saveFig = 1;
datarefs = [6 4 5 3 1 2];
figure
backC = [1 1 1];
set(gcf,'units','normalized')
set(gcf,'position',[1.1995    0.0783    0.5292    0.7958])
figPos = get(gcf,'position');
set(gcf,'color',backC,'position',figPos)
% hax = zeros(1,dataCt);
hax(1) = axes;
haxPos = get(hax(1),'position');
growShiftAx = [-0.5,-0.6,0,-0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on

% xtickpos = [-360 0 360];
xtickpos = [-270 0 270];

ydrop = -0.2;
ytickpos = [1 2 3 4 5 6]-0.5;
ytickpos(4:end) = ytickpos(4:end)+0.5;
fontC = [0 0 0];
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+range(xtickpos)*0.1],...
    'ticklength',[0.02 0.01],'YLim',[-0.2 max(ytickpos)+1],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','w','xcolor','w')

plot([0 0],[0 max(ytickpos)+0.5],'color',[0 0 0],'linewidth',2)
xplot = cat(1,xtickpos,xtickpos,NaN(1,numel(xtickpos)));
yplot = cat(1,zeros(1,numel(xtickpos))+ydrop/2,zeros(1,numel(xtickpos))+ydrop,zeros(1,numel(xtickpos)));
plot(xplot,yplot,'k','linewidth',2,'parent',hax(1))
plot([xtickpos(1) xtickpos(end)],[ydrop ydrop]/2,'k','linewidth',2,'parent',hax(1))
dataLabels = unilateralTable.Properties.RowNames;

specColrOps = {[255 255 255]
%     [152,78,163] %purple
    [255,127,0] %orange
    [200 200 200]};
%     [228,26,28] %red
% [55,126,184]}; %blue
%     [77,175,74] %green

hplot = zeros(numel(specColrOps),1);
% specColrOps = cellfun(@(x) x/255,specColrOps,'uniformoutput',false);
specColrOps = cellfun(@(x) double(uint8(x.^1))/255,specColrOps,'uniformoutput',false);
legendTable = table(specColrOps(1:end),hplot(:),'RowNames',fliplr({'Bilateral Labeling','Unilateral Labeling','No Labeling'}),...
    'VariableNames',{'specColrOps','hplot'});
colrRefs = [1 2 3 1 2 3];
for iterL = 1:numel(datarefs)
    fitdata = unilateralTable.data{datarefs(iterL)};
    if isempty(strfind(dataLabels{datarefs(iterL)},'Unilateral'))
        if ipsiVlr == 1
            yflip = repmat([1 -1],1,round(numel(fitdata)/2));
            fitdata = fitdata.*(yflip(1:numel(fitdata)))';
        end
    end
    ydata = fitdata;
    colr = specColrOps{colrRefs(iterL)};
    xpos = ytickpos(iterL);
    hbox = boxplot(hax(1),ydata,'positions',xpos,'boxstyle','outline',...
        'colors','k','notch','on','widths',0.6,'orientation','horizontal');
    set(hbox(7,:),'marker','.','markeredgecolor','k','markersize',12)
    xdata = get(hbox(5),'xdata');
    ydata = get(hbox(5),'ydata');
    hfill = fill(xdata,ydata,colr);
    delete(hbox)
    ydata = fitdata;
    hbox = boxplot(hax(1),ydata,'positions',xpos,'boxstyle','outline',...
        'colors','k','notch','on','widths',0.6,'orientation','horizontal');
    set(hbox(7,:),'marker','.','markeredgecolor','k','markersize',12)
    for iterH = 1:numel(hbox)
        xvals = get(hbox(iterH),'xdata');
        if max(xvals(:)) > max(xtickpos)
            xvals(xvals > max(xtickpos)) = max(xtickpos);
            set(hbox(iterH),'xdata',xvals)
        end
    end
    text(xtickpos(1)-(0.3)*range(xtickpos),ytickpos(iterL)-0.3,['data ct = ' num2str(numel(ydata))],'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',10);
    text(xtickpos(1)-(0.3)*range(xtickpos),ytickpos(iterL)-0,['fly ct = ',...
        num2str(unilateralTable.flyCt(datarefs(iterL)))],'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',10);
    text(xtickpos(1)-(0.3)*range(xtickpos),ytickpos(iterL)+0.3,(dataLabels{datarefs(iterL)}),'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',10,'interpreter','none');
    if iterL < 4
        legendTable.hplot(iterL) = hfill(1);
    end
end
legendTable = flipud(legendTable);
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+range(xtickpos)*0.1],...
    'ticklength',[0.02 0.01],'YLim',[-0.2 max(ytickpos)+1],'tickdir','out','xtick',xtickpos,...
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


hleg = legend(legendTable.hplot,legendTable.Properties.RowNames,'interpreter','none');
hlegpos = get(hleg,'pos');
hlegpos(1:2) = [0.28 0.65];
hlegpos(4) = hlegpos(4)*1.5;
set(hleg,'pos',hlegpos,'box','off','fontsize',14)

xPts = max(xtickpos)+range(xtickpos)*[0.04 0.08 0.14];
xdata = [xPts(1) xPts(2) xPts(2) xPts(1)];
ydata = repmat([ytickpos(1)-0.4 ytickpos(3)+0.4],2,1);
plot(xdata(:),ydata(:),'color','k','linewidth',2)
text(xPts(3),ytickpos(2),'LC10','rotation',-90,...
        'color',fontC,'horizontalalignment','center','fontsize',16,'interpreter','none');
xdata = [xPts(1) xPts(2) xPts(2) xPts(1)];
ydata = repmat([ytickpos(4)-0.4 ytickpos(6)+0.4],2,1);
plot(xdata(:),ydata(:),'color','k','linewidth',2)
text(xPts(3),ytickpos(5),'LC16','rotation',-90,...
        'color',fontC,'horizontalalignment','center','fontsize',16,'interpreter','none');

xtickval = strtrim(cellstr(num2str((xtickpos)')));
set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) -0.2];
for iterL = 1:numel(xtickval)
    if iterL == 1
        if ipsiVlr == 1
            xlabel = {[xtickval{iterL} char(176)];'Contralateral'};
        else
%             xlabel = {[xtickval{iterL} char(176)];'Left'};
            xlabel = {[xtickval{iterL} char(176)]};
%             textborder(xtickpos(iterL),xlabelPos(2),{'';'';'';'Away from unilateral';'activation'},...
%                 mean(cat(1,specColrOps{2},[1 1 1]*0.5)),[0 0 0]+0.8,'rotation',0,...
%                 'horizontalalignment','center','fontsize',16,...
%                 'verticalalignment','top');
            text(xtickpos(iterL)-mean(xtickpos(1:2))/2,xlabelPos(2),{'';'';'Away from unilateral';'activation'},...
                'color',mean(cat(1,specColrOps{2},[1 1 1]*0.5)),'rotation',0,...
                'horizontalalignment','center','fontsize',16,...
                'verticalalignment','top','fontweight','bold');
        end
    elseif iterL == numel(xtickval)
        if ipsiVlr == 1
            xlabel = {[xtickval{iterL} char(176)];'Ipsilateral'};
        else
%             xlabel = {[xtickval{iterL} char(176)];'Right'};
            xlabel = {[xtickval{iterL} char(176)]};
%             textborder(xtickpos(iterL),xlabelPos(2),{'';'';'';'Toward unilateral';'activation'},...
%                 specColrOps{2}.^1,[0 0 0]+0.8,'rotation',0,...
%                 'horizontalalignment','center','fontsize',16,...
%                 'verticalalignment','top');
            text(xtickpos(iterL)+mean(xtickpos(1:2))/2,xlabelPos(2),{'';'';'Toward unilateral';'activation'},...
                'color',mean(cat(1,specColrOps{2},[1 1 1]*0.5)),'rotation',0,...
                'horizontalalignment','center','fontsize',16,...
                'verticalalignment','top','fontweight','bold');
%             textborder(x, y, string, text_color, border_color, varargin)
        end
    else
        xlabel = [xtickval{iterL} char(176)];
    end
    text(xtickpos(iterL),xlabelPos(2),xlabel,'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',16,...
        'verticalalignment','top');
end
if ipsiVlr == 2
%     plot([-200 -110],[4.6 5],'color','k','linewidth',1)
%     text(-200,4.6,'Away from activation','rotation',0,...
%         'color',fontC,'horizontalalignment','center','fontsize',14,...
%         'verticalalignment','top');
%     plot([150 230],[1.5 1],'color','k','linewidth',1)
%     text(230,1,'Toward activation','rotation',0,...
%         'color',fontC,'horizontalalignment','center','fontsize',14,...
%         'verticalalignment','top');
end
text(median(xtickpos),xlabelPos(2)-0.75,{'Angular velocity (deg/sec)'},'rotation',0,...
    'color',fontC,'horizontalalignment','center','fontsize',14);
if ipsiVlr == 1
    pdfName = 'LC16vLC10_unilateral_turning_ipsilateralVScontralateral.pdf';
else
    pdfName = 'LC16vLC10_unilateral_turning_leftVSright.pdf';
end
text(xtickpos(1),numel(datarefs)+4,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
if saveFig
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab',pdfName);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName(1:end-4) '.eps']);
    export_fig(writePath,'-nocrop')
end

%%

% clearvars
axesObj = findobj(gcf,'Type','Axes');
for iterA = 1:numel(axesObj)
    plotObj = findobj(axesObj(iterA),'Type','Line');
    set(plotObj(1),'LineWidth',8)
    set(plotObj(2),'LineWidth',8)
    set(plotObj(3),'MarkerSize',20,'marker','.','markeredgecolor','auto','markerfacecolor','none','color',[113 171 77]/255)
    set(plotObj(4),'MarkerSize',35,'marker','.','markeredgecolor','auto','markerfacecolor','none','color','k')
    set(plotObj(5),'LineWidth',3,'color',[.3 .3 .3])
    set(plotObj(6),'LineWidth',3,'color',[.3 .3 .3])
%     get(plotObj(5))
end
textObj = findobj(gcf,'Type','Text');
for iterT = 1:numel(textObj)
    set(textObj(iterT),'fontsize',20,'color',[.3 .3 .3],'fontweight','bold')
end
%%
axesObj = findobj(gcf,'Type','Axes');
for iterA = 1:numel(axesObj)
    if iterA == 1 || iterA == 3
        plotC = cA{1};
    else
        plotC = cA{2};
        set(axesObj(iterA),'position',get(axesObj(iterA-1),'position'))
    end
    plotObj = findobj(axesObj(iterA),'Type','Line');
    set(plotObj(1),'LineWidth',8,'color',plotC)
    set(plotObj(2),'LineWidth',3,'linestyle','--','color',plotC)
    set(plotObj(3),'LineWidth',3,'color',[.3 .3 .3])
    set(plotObj(4),'LineWidth',3,'color',[.3 .3 .3])
%     get(plotObj(3))
end
textObj = findobj(gcf,'Type','Text');
for iterT = 1:numel(textObj)
    set(textObj(iterT),'fontsize',20,'color',[.3 .3 .3],'fontweight','bold')
end

%% jitter plots for pez forward, back, turn
graphTable = makeGraphingTable;
%%

clearvars -except graphTable

% graphName = 'relMotion_walking_Vel_FullDur';
graphName = 'relMotion_turning_Vel_FullDur';
makeGraphOptionsStruct([],graphName)
optionsPath = 'Y:\Data_pez3000_analyzed\WRW_graphing_variables\graphOptions.mat';
%%%%%
makeGraphOptionsStruct([],graphName)
[plotTable,excelTable] = addPlotData(graphTable);
%%%%%
graphOptionsLoading = load(optionsPath);
graphOptions = graphOptionsLoading.graphOptions;
excelPath = graphOptions.excelPath;
sheetName2plot = graphOptions.sheetName;

saveFig = 1;
write2excel = 1;
ctrls = {'empty-ctrl'
    'SS02638-split-ctrl-A'
    'SS02638-split-ctrl-B'
    '20B-split-ctrl-A'
    '20B-split-ctrl-B'
    '46B-split-ctrl-A'
    '46B-split-ctrl-B'
    '77B-split-ctrl-A'
    '77B-split-ctrl-B'};
ctrlData = cell2mat(plotTable.returnData(ctrls));
ctrlData(isnan(ctrlData)) = [];
figure
hmKeep = logical(excelTable.keep_heatmap);
excelNames = excelTable.Properties.RowNames;
dataCellB = plotTable.returnData(excelNames(hmKeep));
dataorder = excelTable.order_heatmap(hmKeep);
labelB = excelTable.cellType(hmKeep);
% [~,dataorder] = sort(sortdata);
if ~isempty(strfind(graphName,'turn'))
    ctrlData = ctrlData/1000;
    ctrlData(ctrlData > 1.5) = [];
    ystep = 0.5;
    ylines = ystep*(0:3);
    yoff = -ystep*0.25;
else
    ctrlData = ctrlData*1000;
    ystep = 15;
    ylines = [-ystep 0 ystep];
    yoff = -ystep*1.25;
end
posThresh = prctile(ctrlData,97.5);
negThresh = prctile(ctrlData,2.5);
for iterY = 1:numel(ylines)
    plot([0 numel(dataorder)]+0.5,[ylines(iterY) ylines(iterY)],'linewidth',1,'color',[1 1 1]*0.7)
    text(0,ylines(iterY),num2str(ylines(iterY)),'fontsize',14,'horizontalalignment','right')
    hold on
end
plot([0 numel(dataorder)]+0.5,[0 0],'linewidth',1.5,'color',[1 1 1]*0)
for iterD = 1:numel(dataorder)
    dataref = dataorder(iterD);
    ydata = dataCellB{iterD};
    if ~isempty(strfind(graphName,'turn'))
        ydata = ydata/1000;
    else
        ydata = ydata*1000;
    end
    xdata = (rand(numel(ydata),1)-0.5)*0.5+dataref;
    ydata(ydata > max(ylines)) = NaN;
    ydata(ydata < min(ylines)) = NaN;
    ydata = ydata(:);
    xdata = xdata(:);
    plot([dataref dataref],[yoff*0.95 max(ylines)],'linewidth',1,'color',[1 1 1]*0.7)
    if numel(ydata) > 15 && ~strcmp(labelB{iterD},'LPLC2')
        excelPos = strcmp(excelTable.cellType(hmKeep),labelB{iterD});
        if ~isempty(strfind(graphName,'turning'))
            ydata = abs(ydata);
            midDataX = xdata(ydata < posThresh);
            midDataY = ydata(ydata < posThresh);
            freqTurn = sum(ydata > posThresh)/numel(ydata);
            excelTable.FlyPez_turning_prct(excelPos) = freqTurn;
        else
            midDataX = xdata(ydata < posThresh & ydata > negThresh);
            midDataY = ydata(ydata < posThresh & ydata > negThresh);
            freqBack = sum(ydata > posThresh)/numel(ydata);
            freqFwrd = sum(ydata < negThresh)/numel(ydata);
            excelTable.FlyPez_progressing_prct(excelPos) = freqFwrd;
            excelTable.FlyPez_regressing_prct(excelPos) = freqBack;
        end
        plot(xdata(:),ydata(:),'.','color','k')
        plot(midDataX(:),midDataY(:),'.','color',[0 0 0]+0.5)
    else
%         disp(labelB{iterD})
    end
    
%     labelA{dataref} = regexprep(labelA{dataref},'_male','');
%     labelA{dataref} = regexprep(labelA{dataref},'_attP18','');
%     labelA{dataref} = regexprep(labelA{dataref},'_attP40','');
%     labelA{dataref} = regexprep(labelA{dataref},'_attP2','');
%     text(dataref,xoff,labelA{dataref},'rotation',45,'horizontalalignment','left','interpreter','none')
    if ~isempty(strfind(graphName,'turn'))
        labely = yoff*1.5;
    else
        labely = yoff;
    end
    text(dataref,labely,labelB{iterD},'rotation',45,'horizontalalignment','right','interpreter','none','fontsize',14)
end
ylim = get(gca,'ylim');
ylim(1) = yoff*1.5;
if ~isempty(strfind(graphName,'turn'))
    text(0,ylim(2)+diff(ylim)*0.2,'FlyPez','fontsize',18)
    text(-1,median(ylines),'Turning  (deg/ms)','rotation',90,'fontsize',14,'horizontalalignment','center')
    plot([0 numel(dataorder)]+0.5,repmat(posThresh,1,2),'color',[0.8 0 0],'linewidth',1,'linestyle','--')
else
    text(0,ylim(2)+diff(ylim)*0.2,'FlyPez','fontsize',18)
    text(-1,median(ylines),'Walking  (mm/ms)','rotation',90,'fontsize',14,'horizontalalignment','center')
    plot([0 numel(dataorder)]+0.5,repmat(posThresh,1,2),'color',[0.8 0 0],'linewidth',1,'linestyle','--')
    plot([0 numel(dataorder)]+0.5,repmat(negThresh,1,2),'color',[0.8 0 0],'linewidth',1,'linestyle','--')
end

haxPos = get(gca,'position');
if ~isempty(strfind(graphName,'turn'))
    growShiftAx = [-0,-0.4,0,-0.1];%grow width, grow height, shift X, shift Y
else
    growShiftAx = [-0,-0.1,0,-0.1];%grow width, grow height, shift X, shift Y
end
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(gca,'ylim',ylim,'box','off','xcolor','w','ycolor','w','pos',haxPos)
set(gcf,'color','w')
set(gcf,'position',[101 278 1750 440])
if saveFig
    if ~isempty(strfind(graphName,'turn'))
        pdfName = 'Turning_pez_assay_jitterPlot.pdf';
    else
        pdfName = 'Longitudinal_motion_pez_assay_jitterPlot.pdf';
    end
    writePath = fullfile('Y:\Ming_RubinLab',pdfName);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName(1:end-4) '.eps']);
    export_fig(writePath,'-nocrop')
end

if write2excel == 1
    try
        writetable(excelTable,excelPath,'Sheet',sheetName2plot,'WriteRowNames',true)
    catch
        disp('could not write to excel file')
    end
end

%% box plots for jumping AND backing figures
graphTable = makeGraphingTable;
%%
rng('default');
rng(19);

clearvars -except graphTable

optionsPath = 'Z:\Data_pez3000_analyzed\WRW_graphing_variables\graphOptions.mat';
graphOptionsLoading = load(optionsPath);
graphOptions = graphOptionsLoading.graphOptions;
excelPath = graphOptions.excelPath;
[~,excelName] = fileparts(excelPath);
if strcmp(excelName,'backing_on_pez')
    graphName = 'relMotion_walking_Dist_FullDur';
else
    graphName = 'wing2fot_histogram';
end
makeGraphOptionsStruct([],graphName)
%%%%%
[plotTable,excelTable] = addPlotData(graphTable);
%%%%%
close all
saveFig = 2;
if ~isempty(strfind(graphName,'relMotion_walking'))
    growShiftAx = [-0.5,-0.45,0,0.1];%grow width, grow height, shift X, shift Y
    colors = cat(2,{[0 0 0],[1 0 0],[0 0 0]},repmat({[0 0 1]},1,5));
    ymax = 2.5;
    ymin = -2.5;
    pdfName = 'backing_on_pez';
else
    growShiftAx = [-0.5,-0.45,0,0.1];%grow width, grow height, shift X, shift Y
    colors = {[0 0 0],[0 0 0],[0 0 0]};
    colors = cat(2,colors,colors,colors);
    ymax = 7;
    ymin = -0.1;
    pdfName = 'takeoff_sequence_duration';
end
[~,dataOrder] = sort(plotTable.order);
yrange = ymax-ymin;
plotTable = plotTable(dataOrder,:);
% colors = colors(dataOrder);
dataCell = plotTable.returnData;
dataCt = numel(dataCell);
botlabel = plotTable.plotID;

figure
set(gcf,'pos',[2062 264 866 633],'color','w')
hax = axes;
haxPosA = get(hax(1),'position');
haxPosA = [haxPosA(1)-haxPosA(3)*growShiftAx(1)/2,haxPosA(2)-haxPosA(4)*growShiftAx(2)/2,...
    haxPosA(3)*(1+growShiftAx(1)),haxPosA(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(1),'position',haxPosA,'ylim',[ymin ymax+yrange*0.1],'xlim',[0 dataCt],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,'xcolor','none','ycolor','none')
if ~isempty(strfind(graphName,'relMotion_walking'))
    plot([0 dataCt],[0 0],'k','linewidth',2,'parent',hax(1),'color',[1 1 1]*0.5)
end
xticks = (1:dataCt)-0.3;
kruskData = cell(dataCt,1);
kruskRefs = cell(dataCt,1);
boxVals = [];
colr = [0 0 0];
for iterD = 1:dataCt
    plotY = dataCell{iterD};
    maxCt = 130;
    if numel(plotY) > maxCt
        permVec = randperm(numel(plotY));
        plotY = plotY(permVec(1:maxCt));
    end
    boxVals = cat(1,boxVals,[plotY(:) zeros(numel(plotY),1)+iterD]);
    text(xticks(iterD),ymin-yrange*0.15,['n = ' num2str(numel(plotY))],'rotation',30,...
        'horizontalalignment','right','fontsize',14,'parent',hax(1),'interpreter','none')
    text(xticks(iterD),ymin-yrange*0.3,strtrim(botlabel{iterD}),'rotation',30,...
        'horizontalalignment','right','fontsize',14,'parent',hax(1),'interpreter','none')
    kruskData{iterD} = dataCell{iterD}(:);
    kruskRefs{iterD} = zeros(numel(dataCell{iterD}),1)+iterD;
end
hbox = boxplot(hax(1),boxVals(:,1),boxVals(:,2),'positions',xticks,'boxstyle','outline',...
    'colors',colr,'notch','on','widths',0.5);
set(hbox(end,:),'marker','.','markeredgecolor',colr,'markersize',12)
set(hax(1),'position',haxPosA,'ylim',[ymin ymax+yrange*0.1],'xlim',[0 dataCt],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,'xcolor','none','ycolor','none','color','none')
kruskRefs = cat(1,kruskRefs{:});
kruskDataB = cat(1,kruskData{:});
kruskGroup = (kruskRefs == 1 | kruskRefs == 2);
[pVal,tabl,stats] = kruskalwallis(kruskDataB(kruskGroup),kruskRefs(kruskGroup),'off');
% line([xticks(1) xticks(2)]+[0.1 -0.1],[ymax ymax]+yrange*0.1,'color','k','linewidth',2)
% text(mean([xticks(1) xticks(2)]),ymax+yrange*0.175,['p = ' num2str(pVal,2)],'rotation',0,...
%     'horizontalalignment','center','fontsize',14,'parent',hax(1),'interpreter','none')

kruskGroup = (kruskRefs == 2 | kruskRefs == 3);
[pVal,tabl,stats] = kruskalwallis(kruskDataB(kruskGroup),kruskRefs(kruskGroup),'off');
% line([xticks(2) xticks(3)]+[0.1 -0.1],[ymax ymax]+yrange*0.1,'color','k','linewidth',2)
% text(mean([xticks(2) xticks(3)]),ymax+yrange*0.175,['p = ' num2str(pVal,2)],'rotation',0,...
%     'horizontalalignment','center','fontsize',14,'parent',hax(1),'interpreter','none')
% text(mean([xticks(1) xticks(3)]),ymax+yrange*0.3,'Kruskal-Wallis nonparametric ANOVA test','rotation',0,...
%     'horizontalalignment','center','fontsize',14,'parent',hax(1),'interpreter','none')

linewidth = 2.5;
obj = findobj(gcf,'Tag','Upper Adjacent Value');
set(obj,'visible','off')
obj = findobj(gcf,'Tag','Lower Adjacent Value');
set(obj,'visible','off')
obj = findobj(gcf,'Tag','Upper Whisker');
set(obj,'LineStyle','-','LineWidth',linewidth)%,'color','k');
obj = findobj(gcf,'Tag','Lower Whisker');
set(obj,'LineStyle','-','LineWidth',linewidth)%,'color','k');
obj = findobj(gcf,'Tag','Box');
set(obj,'LineWidth',linewidth)
obj = findobj(gcf,'Tag','Median');
set(obj,'LineWidth',linewidth)

if ~isempty(strfind(graphName,'relMotion_walking'))
    text(-1.4,yrange/2+ymin,'Distance (mm)','rotation',90,...
        'horizontalalignment','center','fontsize',14,'parent',hax(1))
    text(-0.9,yrange*0.75+ymin,'Forward','rotation',90,...
        'horizontalalignment','center','fontsize',14,'parent',hax(1))
    text(-0.9,yrange*0.25+ymin,'Backward','rotation',90,...
        'horizontalalignment','center','fontsize',14,'parent',hax(1))
%     annotation('arrow',[0.24 0.24],[0.73 0.83]-0.075)
%     annotation('arrow',[0.24 0.24],[0.61 0.51]-0.075)
    ytickval = [ymin 0 ymax];
    yticks = ytickval;
else
    text(-0.9,yrange/2.5,{'Takeoff sequence','duration (ms)'},'rotation',90,...
        'horizontalalignment','center','fontsize',14,'parent',hax(1))
    ytickval = unique(round(logspace(log10(1),log10(300),5)/5)*5);
    ytickval = [1 2 ytickval(2:end)];
    yticks = log(ytickval);
end
for iterT = 1:numel(yticks)
    text(-0.25,yticks(iterT),num2str(ytickval(iterT)),'rotation',0,...
        'horizontalalignment','right','fontsize',14,'parent',hax(1),'interpreter','none')
end

yplot = cat(1,yticks,yticks,NaN(1,numel(yticks)));
xplot = cat(1,zeros(1,numel(yticks))-0.10,zeros(1,numel(yticks))+0.1,zeros(1,numel(yticks)));
plot(xplot,yplot,'k','linewidth',2,'parent',hax(1))
plot([0.1 0.1],[yticks(1) yticks(end)],'k','linewidth',2,'parent',hax(1))

% plot([0 8],[0.5 0.5],'k','linewidth',1,'parent',hax(1))

% hleg = legend(hplot,labelB);
% hlegpos = get(hleg,'pos');
% hlegpos(1:2) = [0.85 0.9];
% set(hleg,'pos',hlegpos)

if saveFig
    %%
    pdfName = [pdfName '_revisions']
%     writePath = fullfile('Z:\Ming_RubinLab',[pdfName '.pdf']);
%     export_fig(writePath,'-nocrop')
    writePath = fullfile('Z:\Ming_RubinLab',[pdfName '.eps']);
    export_fig(writePath,'-nocrop')
end

%% wing raising

excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2.xlsx';
excelPath = fullfile(excelDir,excelName);
dataTable = readtable(excelPath,'Sheet','wing_extension_scoring');
%%
saveFig = 1;
% close all
hfig = figure;
set(gcf,'pos',[2062 264 766 633],'color','w')
hax = zeros(4,1);
growShiftAxCell = {[-0.4,-0.35,0.0,0.2] %grow width, grow height, shift X, shift Y
    [-0.4,-0.35,0.0,-0.1]
    [-0.45,-0.7,0.1,0.2] %grow width, grow height, shift X, shift Y
    [-0.45,-0.7,0.1,-0.1]};
legAx = zeros(3,1);
for i = 1:2
% for i = 3:4
    hax(i) = axes;
    growShiftAx = growShiftAxCell{i};
    haxPosA = get(hax(i),'position');
    haxPosA = [haxPosA(1)-haxPosA(3)*growShiftAx(1)/2,haxPosA(2)-haxPosA(4)*growShiftAx(2)/2,...
        haxPosA(3)*(1+growShiftAx(1)),haxPosA(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
    axisRef = dataTable.axisRef == i;
    dataCt = sum(axisRef);
    ymin = -10;
    ymax = 100;
    set(hax(i),'position',haxPosA,'ylim',[-100 150],'xlim',[0 dataCt],'xtick',[],...
        'nextplot','add','box','off','fontsize',14,'xcolor','none','ycolor','none','color','none')
    plot([0 dataCt],[ymin ymin]+1,'linewidth',2,'color','k')
    plot([0 0],[ymin+1 ymax],'linewidth',2,'color','k')
    xplot = repmat([0 dataCt NaN]',3,1);
    yplot = repmat([0 50 100],3,1);
    plot(xplot(:),yplot(:),'linewidth',1,'color',[1 1 1]*0.5)
    plotData = dataTable.x_wingext(axisRef);
    plotCt = dataTable.trial(axisRef);
    for j = 1:3
        text(-0.25,(j-1)*50,num2str((j-1)*50),'fontsize',16,'color','k','parent',hax(i),'horizontalalignment','right')
    end
    ydata = dataTable.x_wingext(axisRef)*100;
    xdata = (1:dataCt)-0.5;
%     ydata = [ydata';ydata';NaN(1,size(ydata,1))];
%     plot([xdata-0.2;xdata+0.2;xdata],ydata,'color','k','linewidth',2);
    plot(xdata,ydata,'.','markersize',38,'color','k')
    trialCt = dataTable.trial(axisRef);
    behaviorCt = dataTable.x_wingext(axisRef).*dataTable.trial(axisRef);
    alpha = 0.05;
    [~,pci] = binofit(behaviorCt,trialCt,alpha);
    errorY = pci;
    errorY = [errorY';NaN(1,size(errorY,1))]*100;
    plot(repmat(xdata,3,1),errorY,'color','k','linewidth',2);
    gender = dataTable.gender(axisRef);
    for iterL = 1:numel(ydata)
        if strcmp(gender{iterL},'mixed')
            hplot = plot(xdata(iterL),ydata(iterL),'.','markersize',28,'parent',hax(i),'color',[1 1 1]*0.6);
            legAx(3) = hplot;
        elseif strcmp(gender{iterL},'female')
            hplot = plot(xdata(iterL),ydata(iterL),'.','markersize',28,'parent',hax(i),'color',[1 0.4 0.4]*0.8);
            legAx(2) = hplot;
        else
            hplot = plot(xdata(iterL),ydata(iterL),'.','markersize',28,'parent',hax(i),'color',[0.4 0.4 1]*0.8);
            legAx(1) = hplot;
        end
    end
    plotLabel = [dataTable.unilateral_bilateral{i} ' %'];
    text(-1.5,50,{plotLabel;'(Arena)'},'parent',hax(i),'color','k',...
        'fontsize',16,'rotation',90,'interpreter','none','horizontalalignment','center');
    plot([1 1],[-75 100],'--','color',[0 0 0]+0.5,'linewidth',1,'parent',hax(i))
    plot([3 3],[-75 100],'--','color',[0 0 0]+0.5,'linewidth',1,'parent',hax(i))
    plot([5 5],[-75 100],'--','color',[0 0 0]+0.5,'linewidth',1,'parent',hax(i))
end
plot([xdata(end)-0.6 xdata(end)+1],[50 50],'linewidth',4,'color','w','parent',hax(i-1))
hleg = legend(legAx,{'Male','Female','Mixed'});
hlegpos = get(hleg,'pos');
hlegpos(1:2) = [0.66 0.72];
set(hleg,'pos',hlegpos,'box','off','color','w')
cellTypes = dataTable.driver_cellType(axisRef);
gender = dataTable.gender(axisRef);
for k = 1:dataCt
    if k == 1
        labx = k-0.5;
    elseif k == 2 || k == 4 || k == 6
        labx = k;
    else
        continue
    end
%     text(labx,-30,cat(2,cellTypes{k},' ',gender{k}),'parent',hax(i),'color','k',...
%         'fontsize',16,'rotation',45,'interpreter','none','horizontalalignment','right');
    text(labx,-100,cellTypes{k},'parent',hax(i),'color','k',...
        'fontsize',16,'rotation',45,'interpreter','none','horizontalalignment','right');
end

if hax(1) == 0
    pdfName = 'wing_extension_AvsD';
    htitle = text(dataCt/2,103,'Wing extension','interpreter','none','parent',hax(3),...
        'horizontalalignment','center','fontsize',18);
else
    pdfName = 'wing_extension';
    htitle = text(dataCt/2,110,'Wing extension','interpreter','none','parent',hax(1),...
        'horizontalalignment','center','fontsize',18);
end
tpos = get(htitle,'pos');
tpos(2) = tpos(2)*1.2;
set(htitle,'pos',tpos)
if saveFig
    %%
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.pdf']);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.eps']);
    export_fig(writePath,'-nocrop')
end
