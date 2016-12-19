%% New LC10 multi-compare

clear all
excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2.xlsx';
sheetName = 'LC10_cellType_specific_tests';
excelTable = readtable(fullfile(excelDir,excelName),'Sheet',sheetName,'readrownames',true);
%%
rng('default');
rng(22);
saveFig = 1;
close all
figure
set(gcf,'pos',[2137 264 1191 633],'color','w')
hax = axes;
ymin = -10;
ymax = 100;
dataCt = 3;
haxPosA = get(hax(1),'position');
growShiftAx = [-0.65,-0.3,-0.15,0.05];%grow width, grow height, shift X, shift Y
haxPosA = [haxPosA(1)-haxPosA(3)*growShiftAx(1)/2,haxPosA(2)-haxPosA(4)*growShiftAx(2)/2,...
    haxPosA(3)*(1+growShiftAx(1)),haxPosA(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(1),'position',haxPosA,'ylim',[ymin ymax],'xlim',[0 dataCt],'xtick',[],...
    'nextplot','add','box','off','fontsize',14,'xcolor','none','ycolor','none','color','none')
plot([0 dataCt],[ymin ymin]+1,'linewidth',2,'color','k')
plot([0 0],[ymin+1 ymax],'linewidth',2,'color','k')
xplot = repmat([0 dataCt NaN]',3,1);
yplot = repmat([0 50 100],3,1);
plot(xplot(:),yplot(:),'linewidth',1,'color',[1 1 1]*0.5)

specStrOps = unique(excelTable.subGroups(excelTable.color > 1));
specStrOps = cat(1,specStrOps(:),{'unique splits'});
specColrOps = {[120 120 120]
    [152,78,163] %purple
    [77,175,74] %green
    [228,26,28] %red
    [255,127,0]}; %orange
%     [55,126,184]}; %blue
hplot = zeros(numel(specColrOps),1);
% specColrOps = cellfun(@(x) x/255,specColrOps,'uniformoutput',false);
specColrOps = cellfun(@(x) double(uint8(x.^1.05))/255,specColrOps,'uniformoutput',false);
legendTable = table(specColrOps(1:end),hplot(:),'RowNames',specStrOps(:),'VariableNames',{'specColrOps','hplot'});
groupAD = logical(excelTable.groupA) & logical(excelTable.groupD) | strcmp(excelTable.subGroups,'blind');
groupBC = logical(excelTable.groupB) & logical(excelTable.groupC);
groupBC = groupBC & ~(logical(excelTable.groupA) | logical(excelTable.groupD));
groupCtrl = ~(strcmp(excelTable.subGroups,'split GAL4 lines') | strcmp(excelTable.subGroups,'blind'));

iterD = 1;
xlabels = {'a and/or d','Neither a nor d','Negative Controls'};
for iterP = 1:3
    groupVec = excelTable.plotGroup == iterP;
    colorVec = excelTable.color(groupVec);
    ydata = excelTable.reaching(groupVec)*100;
    dataLabels = excelTable.subGroups(groupVec);
    xdata = randperm(numel(ydata));
    xdata = ((xdata-min(xdata))/range(xdata)-0.5)*0.4+0.5+(iterP-1);
    for iterS = 1:numel(ydata)
        if strcmp(dataLabels{iterS},'a only') || strcmp(dataLabels{iterS},'d only')
            plot([xdata(iterS) xdata(iterS)+0.05],[ydata(iterS) ydata(iterS)-10],'color','k','linewidth',2)
            text(xdata(iterS)+0.1,ydata(iterS)-10,dataLabels{iterS},'rotation',0,'horizontalalignment','left',...
                'fontsize',16,'parent',hax(iterD),'interpreter','none')
        end
        hplotPreA = plot(xdata(iterS),ydata(iterS),'.','markersize',48,'parent',hax(iterD),'color','k');
        hplotA = plot(xdata(iterS),ydata(iterS),'.','markersize',38,'parent',hax(iterD),'color',specColrOps{colorVec(iterS)});
        if ~isempty(hplotA) && colorVec(iterS) > 0
            strVecRef = find(colorVec(iterS) == excelTable.color);
            legendTable.hplot(excelTable.subGroups{strVecRef(1)}) = hplotA(1);
        end
    end
    meanval = mean(ydata);
    stdval = std(ydata);
    ydata = [meanval meanval NaN meanval+stdval meanval-stdval];
    xdata = [0.2 0.8 0 0.5 0.5]+iterP-1;
    plot(xdata,ydata,'k','linewidth',4,'parent',hax(iterD))
    
    text(iterP-0.5,-15,xlabels{iterP},'rotation',45,...
        'horizontalalignment','right','fontsize',16,'parent',hax(iterD),'interpreter','none')
end


legendTable(legendTable.hplot == 0,:) = [];
hleg = legend(legendTable.hplot,legendTable.Properties.RowNames,'interpreter','none');
hlegpos = get(hleg,'pos');
hlegpos(1:2) = [0.37 0.67];
set(hleg,'pos',hlegpos,'box','off')

iterD = 2;
hax(2) = axes;
haxPosB = get(hax(2),'position');
iterP = 4;
groupVec = find(excelTable.plotGroup == iterP);
sortVec = excelTable.Order(groupVec);
[~,sortVec] = sort(sortVec);
groupVec = groupVec(sortVec);
ydata = excelTable.reaching(groupVec)*100;
dataCt = numel(ydata);
growShiftAx = [-0.7,-0.3,0.2,0.05];%grow width, grow height, shift X, shift Y
haxPosB = [haxPosB(1)-haxPosB(3)*growShiftAx(1)/2,haxPosB(2)-haxPosB(4)*growShiftAx(2)/2,...
    haxPosB(3)*(1+growShiftAx(1)),haxPosB(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
set(hax(2),'position',haxPosB,'ylim',[ymin ymax],'xlim',[0 dataCt],'xtick',[],...
    'xticklabel',[],'nextplot','add','box','off','fontsize',14,'xcolor','none','ycolor','none','color','none')
plot([0 dataCt],[ymin ymin]+1,'linewidth',2,'color','k')
plot([0 0],[ymin+1 ymax],'linewidth',2,'color','k')
xplot = repmat([0 dataCt NaN]',3,1);
yplot = repmat([0 50 100],3,1);
plot(xplot(:),yplot(:),'linewidth',1,'color',[1 1 1]*0.5)
yplot = repmat([ymin ymax NaN]',3,1);
xplot = repmat([2 4 4],3,1);
plot(xplot(:),yplot(:),'--','linewidth',1,'color',[1 1 1]*0.5)

xdata = (1:numel(ydata))-0.5+repmat([0.1 -0.1],1,3);
plot(xdata,ydata,'.','markersize',48,'parent',hax(iterD),'color','k');
xlabels = excelTable(groupVec,:).Properties.RowNames;
for iterL = 1:numel(xlabels)
    text(xdata(iterL),-15,xlabels{iterL},'rotation',45,...
        'horizontalalignment','right','fontsize',16,'parent',hax(iterD),'interpreter','none')
end
alpha = 0.05;
totalVids = excelTable.TotalFlies(groupVec);
[~,pci] = binofit(round((ydata/100).*totalVids),totalVids,alpha);
errorY = pci;
errorY = [errorY';NaN(1,size(errorY,1))]*100;
hplot = plot(repmat(xdata,3,1),errorY,'color','k','linewidth',2);
plot(xdata,ydata,'.','markersize',38,'parent',hax(iterD),'color',[1 1 1]*0.8);
yticks = [0 50 100];
for iterT = 1:numel(yticks)
    text(-0.25,yticks(iterT),num2str(yticks(iterT)),'rotation',0,...
        'horizontalalignment','right','fontsize',16,'parent',hax(1),'interpreter','none')
end
text(-0.75,50,{'Reaching %';'(Single-fly)'},'rotation',90,'horizontalalignment','center','fontsize',16,'parent',hax(1))

pdfName = 'LC10_multiCompare_newGroupings';
text(-0.75,120,pdfName,'fontsize',16,'parent',hax(1),'interpreter','none')

if saveFig
    %%
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.pdf']);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.eps']);
    export_fig(writePath,'-nocrop')
end
