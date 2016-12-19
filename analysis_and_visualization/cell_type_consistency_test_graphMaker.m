clear all
excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2.xlsx';
sheetName = 'cellType_specific_tests';
% sheetName = 'LC10_cellType_specific_tests';
excelTableA = readtable(fullfile(excelDir,excelName),'Sheet',sheetName,'readrownames',true);
excelTableA = excelTableA(logical(excelTableA.Keepers),:);
%%
saveFig = 1;
multChoice = 1;
multOps = {'LC6','LC16'};

% close all
ydrop = -0.1;
figure
set(gcf,'pos',[2137 264 1191 633],'color','w')
hax = axes;
haxPosA = get(hax(1),'position');
growShiftAx = [-0.8,-0.3,-0.25,-0.15];%grow width, grow height, shift X, shift Y
haxPosA = [haxPosA(1)-haxPosA(3)*growShiftAx(1)/2,haxPosA(2)-haxPosA(4)*growShiftAx(2)/2,...
    haxPosA(3)*(1+growShiftAx(1)),haxPosA(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(1),'position',haxPosA,'ylim',[ydrop 1],'xlim',[0 2],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')

hax(2) = axes;
haxPosB = get(hax(2),'position');
growShiftAx = [-0.8,-0.3,0.05,-0.15];%grow width, grow height, shift X, shift Y
haxPosB = [haxPosB(1)-haxPosB(3)*growShiftAx(1)/2,haxPosB(2)-haxPosB(4)*growShiftAx(2)/2,...
    haxPosB(3)*(1+growShiftAx(1)),haxPosB(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(2),'position',haxPosB,'ylim',[ydrop 1],'xlim',[0 2],'xtick',[],...
    'xticklabel',[],'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')

hax(3) = axes;
haxPosB = get(hax(3),'position');
growShiftAx = [-0.8,-0.3,0.35,-0.15];%grow width, grow height, shift X, shift Y
haxPosB = [haxPosB(1)-haxPosB(3)*growShiftAx(1)/2,haxPosB(2)-haxPosB(4)*growShiftAx(2)/2,...
    haxPosB(3)*(1+growShiftAx(1)),haxPosB(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(3),'position',haxPosB,'ylim',[ydrop 1],'xlim',[0 2],'xtick',[],...
    'xticklabel',[],'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')

groupOps = unique(excelTableA.plotGroup);
posRefs = strcmp(excelTableA.clusterRef,'pos');
specStrOps = unique(excelTableA.subGroups);
colorOrder = [1 2 4 5 3];
labelOrder = [2 3 1 5 4];
specColrOps = {[228,26,28]
    [120 120 120]
    [55,126,184]
    [77,175,74]
%     [152,78,163]
    [255,127,0]};%orrange
specColrOps = specColrOps(colorOrder);
hplot = zeros(numel(specStrOps),1);

% specColrOps = cellfun(@(x) double(uint8(x))/255,specColrOps,'uniformoutput',false);
specColrOps = cellfun(@(x) double(uint8(x.^1.05))/255,specColrOps,'uniformoutput',false);
% specColrOps = cellfun(@(x) double(uint8(x*1.35))/255,specColrOps,'uniformoutput',false);
specColrOps = flipud(specColrOps);
legendTable = table(specColrOps(:),hplot(:),'RowNames',specStrOps(:),'VariableNames',{'specColrOps','hplot'});
plotStrOps = {'Jumping','regressing','turning'};
groupOps = repmat(multOps(multChoice),1,3);
for iterD = 1:numel(plotStrOps)
    plot([0 2],[ydrop ydrop],'linewidth',2,'color','k','parent',hax(iterD))
    plot([0 0],[ydrop 1],'linewidth',2,'color','k','parent',hax(iterD))
    xplot = repmat([0 2 NaN]',3,1);
    yplot = repmat([0 0.5 1],3,1);
    plot(xplot(:),yplot(:),'linewidth',1,'color',[1 1 1]*0.5,'parent',hax(iterD))

    yticks = linspace(0,1,5);
    botlabel = round(linspace(0,100,5));
    yticklabel = cellstr(num2str(botlabel(1:end)'));
    skipvals = [2 4];
    yticklabel(skipvals) = repmat({''},numel(skipvals),1);
    for iterT = 1:numel(yticks)
        text(-0.1,yticks(iterT),strtrim(yticklabel{iterT}),'rotation',0,...
            'horizontalalignment','right','fontsize',14,'parent',hax(iterD))
    end
%     yplot = cat(1,yticks,yticks,NaN(1,numel(yticks)));
%     xplot = cat(1,zeros(1,numel(yticks))-0,zeros(1,numel(yticks))+8,zeros(1,numel(yticks)));
%     plot(xplot,yplot,'color',[.5 .5 .5],'linewidth',1,'parent',hax(iterD))
    if ~max(strcmp(plotStrOps{iterD},'Jumping'))
        excelTable = excelTableA(excelTableA.PezMotionDataCt > 5,:);
    else
        excelTable = excelTableA;
    end
    posRefs = strcmp(excelTable.clusterRef,'pos');
    specStrOps = unique(excelTable.subGroups);
    plotRefsAll = strcmp(excelTable.plotGroup,groupOps{iterD});
    legStr = {};
    for iterP = 1:2
        if iterP == 1
            plotRefs = plotRefsAll & ~posRefs;
        else
            plotRefs = plotRefsAll & posRefs;
        end
        data = excelTable.(plotStrOps{iterD});
        dataCt = excelTable.PezMotionDataCt(plotRefs);
        
        ydata = data(plotRefs);
        xdata = (1:numel(ydata));
        xdata = ((xdata-min(xdata))/range(xdata)-0.5)*0.4+0.5+(iterP-1);
        
%         plot(xdata,ydata,'.','markersize',22,'parent',hax(iterD),'color','k')
        for iterN = fliplr(1:numel(specStrOps))
            specRefs = plotRefs & strcmp(excelTable.subGroups,specStrOps{iterN});
            specRefs(~plotRefs) = [];
            yplot = ydata(specRefs);
            xplot = xdata(specRefs);
            hplotA = [];
            for iterS = 1:numel(yplot)
                hplotPreA = plot(xplot(iterS),yplot(iterS),'.','markersize',48,'parent',hax(iterD),'color','k');
                hplotPre = plot(xplot(iterS),yplot(iterS),'.','markersize',38,'parent',hax(iterD),'color',specColrOps{iterN});
                if iterS == 1
                    hplotA = hplotPre;
                end
            end
            if ~isempty(hplotA)
                legendTable.hplot(specStrOps{iterN}) = hplotA(1);
            end
        end
        
        ydata = data(plotRefs);
        meanval = mean(ydata);
        stdval = std(ydata);
        ydata = [meanval meanval NaN meanval+stdval meanval-stdval];
        ydata(ydata > 1) = 1;
        ydata(ydata < 0) = 0;
        xdata = [0.1 0.9 0 0.5 0.5]+iterP-1;
        plot(xdata,ydata,'k','linewidth',5,'parent',hax(iterD))
        
    end
    text(1,1.2,groupOps{iterD},'rotation',0,'horizontalalignment','center','fontsize',14,'parent',hax(iterD))
end
if strcmp(multOps{multChoice},'LC16')
    text(-0.75,0.5,{'Regression %';'(Arena)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(2))
    text(-0.75,0.5,{'Turning %';'(Arena)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(3))
else
    text(-0.75,0.5,{'Regression %';'(Single-fly)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(2))
    text(-0.75,0.5,{'Turning %';'(Single-fly)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(3))
end
text(-0.75,0.5,{'Jumping %';'(Single-fly)'},'rotation',90,'horizontalalignment','center','fontsize',14,'parent',hax(1))
legNames = legendTable.Properties.RowNames;
legNames = legNames(labelOrder);
legHplot = legendTable.hplot;
legHplot = legHplot(labelOrder);
hleg = legend(legHplot,legNames,'interpreter','none');
hlegpos = get(hleg,'pos');
hlegpos(1:2) = [0.8 0.78];
set(hleg,'pos',hlegpos,'box','off')
if strcmp(multOps{multChoice},'LC16')
    pdfName = 'LC16_multiCompare';
else
    pdfName = 'LC6_multiCompare';
end

if saveFig
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab\cell_type_consistency_tests',[pdfName '.pdf']);
    export_fig(writePath)
    writePath = fullfile('Y:\Ming_RubinLab\cell_type_consistency_tests',[pdfName '.eps']);
    export_fig(writePath)
end

%%
clc
excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2.xlsx';
sheetName = 'cellType_specific_tests';
excelTableA = readtable(fullfile(excelDir,excelName),'Sheet',sheetName,'readrownames',true);
excelTableA = excelTableA(logical(excelTableA.Keepers),:);

lineRefs = {'OL0077B_blind','OL0077B'
    'OL0046B_blind','OL0046B'
    'OL0046B_blind','OL0046B'};
assayRefs = {'Jumping','PezJumperDataCt'
    'regressing','BowlDataCt'
    'turning','BowlDataCt'};
lineList = excelTableA.pezRefs;
for iterP = 1:size(lineRefs,1)
    disp([lineRefs{iterP,1} ' -vs- ' lineRefs{iterP,2} ' --- ' assayRefs{iterP,1}]);
    tableRef1 = find(strcmp(lineList,lineRefs{iterP,1}),1,'first');
    tableRef2 = find(strcmp(lineList,lineRefs{iterP,2}),1,'first');
    n1 = excelTableA.(assayRefs{iterP,2})(tableRef1);
    n2 = excelTableA.(assayRefs{iterP,2})(tableRef2);
    x1 = round(excelTableA.(assayRefs{iterP,1})(tableRef1)*n1);
    x2 = round(excelTableA.(assayRefs{iterP,1})(tableRef2)*n2);
    [zvalue,pvalue] = binoPval(x1,x2,n1,n2)
end
%%
clear all
excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2.xlsx';
% sheetName = 'cellType_specific_tests';
sheetName = 'LC10_cellType_specific_tests';
excelTable = readtable(fullfile(excelDir,excelName),'Sheet',sheetName,'readrownames',true);
excelTable = excelTable(excelTable.Order,:);
groups = unique(excelTable.subGroups);
multChoice = 1;
multOps = {'experimentals','controls'};
if ~strcmp(multOps{multChoice},'controls')
    excelTable = excelTable(strcmp(excelTable.subGroups,'unique splits'),:);
else
    excelTable = excelTable(~strcmp(excelTable.subGroups,'unique splits'),:);
end

%%


totalVids = excelTable.TotalFlies;
totalReaching = excelTable.reaching;

xCt = numel(totalVids);

close all

ydrop = -0.1;
figure
set(gcf,'pos',[2137 264 1191 633],'color','w')
hax = axes;
haxPosA = get(hax(1),'position');
if ~strcmp(multOps{multChoice},'controls')
    growShiftAx = [-0.1,-0.4,0,0.1];%grow width, grow height, shift X, shift Y
else
    growShiftAx = [-0.75,-0.4,0,0.1];%grow width, grow height, shift X, shift Y
end
haxPosA = [haxPosA(1)-haxPosA(3)*growShiftAx(1)/2,haxPosA(2)-haxPosA(4)*growShiftAx(2)/2,...
    haxPosA(3)*(1+growShiftAx(1)),haxPosA(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(1),'position',haxPosA,'ylim',[ydrop 1],'xlim',[0 xCt+1],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,'xcolor','w','ycolor','w')
yticks = linspace(0,1,5);
botlabel = round(linspace(0,100,5));
yticklabel = cellstr(num2str(botlabel(1:end)'));
skipvals = [2 4];
yticklabel(skipvals) = repmat({''},numel(skipvals),1);
for iterT = 1:numel(yticks)
    text(-0.1,yticks(iterT),strtrim(yticklabel{iterT}),'rotation',0,...
        'horizontalalignment','right','fontsize',14,'parent',hax(1))
end
yplot = cat(1,yticks,yticks,NaN(1,numel(yticks)));
xplot = cat(1,zeros(1,numel(yticks))-0,zeros(1,numel(yticks))+xCt+1,zeros(1,numel(yticks)));
plot(xplot,yplot,'color',[.5 .5 .5],'linewidth',1,'parent',hax(1))

if ~strcmp(multOps{multChoice},'controls')
    specColrOps = {[0 0 0]};
    %     [152,78,163]};
else
    specColrOps = {[228,26,28]
        [55,126,184]
        [77,175,74]
        [255,127,0]};
end
specColrOps = cellfun(@(x) x/255,specColrOps,'uniformoutput',false);
specColrOps = flipud(specColrOps);
groups = unique(excelTable.subGroups);
hplot = zeros(numel(groups),1);
legendTable = table(specColrOps(:),hplot(:),'RowNames',groups(:),'VariableNames',{'specColrOps','hplot'});

xVals = (1:xCt);
xGridPlot = repmat(xVals,3,1);
yMin = 0;
yMax = 100;
plotC = [0 0 0];
barData = totalReaching*100;
yData = [barData(:)';NaN(1,xCt)];
xPlot = xGridPlot(1:2,:);
for iterP = 1:xCt
    xPlot = iterP;
    yData = totalReaching(iterP)*100;
    plotC = specColrOps{strcmp(groups,excelTable.subGroups{iterP})};
    plot(xPlot,(yData-yMin)/(yMax-yMin),'color',plotC,...
        'marker','.','linestyle','none','markersize',30);
    alpha = 0.05;
    [~,pci] = binofit(round(totalReaching.*totalVids),totalVids,alpha);
    errorY = pci;
    errorY = [errorY';NaN(1,size(errorY,1))]*100;
    hplot = plot(xGridPlot(:,iterP),(errorY(:,iterP)-yMin)/(yMax-yMin),'color',plotC,'linewidth',2);
    
    if ~isempty(hplot)
        legendTable.hplot(excelTable.subGroups{iterP}) = hplot(1);
    end
end

xticks = (1:xCt);
xticklabel = excelTable.Properties.RowNames;
for iterT = 1:numel(xticks)
    text(xticks(iterT),-0.1,strtrim(xticklabel{iterT}),'rotation',45,...
        'horizontalalignment','right','fontsize',14,'parent',hax(1),'interpreter','none')
end
if strcmp(multOps{multChoice},'controls')
    hleg = legend(legendTable.hplot,legendTable.Properties.RowNames,'interpreter','none');
    hlegpos = get(hleg,'pos');
    hlegpos(1:2) = [0.8 0.78];
    set(hleg,'pos',hlegpos)
end

if ~strcmp(multOps{multChoice},'controls')
    pdfName = 'LC10_multiple_split_comparison';
else
    pdfName = 'LC10_multiple_split_comparison_controls';
end
text(4,1.2,pdfName,'rotation',0,'horizontalalignment','center','fontsize',14,'parent',hax(1),'interpreter','none')

if saveFig
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab\cell_type_consistency_tests',[pdfName '.pdf']);
    export_fig(writePath)
    writePath = fullfile('Y:\Ming_RubinLab\cell_type_consistency_tests',[pdfName '.eps']);
    export_fig(writePath)
end