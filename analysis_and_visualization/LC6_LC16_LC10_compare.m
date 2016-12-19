clearvars

excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW.xlsx';
dataExcelTable = readtable(fullfile(excelDir,excelName),'Sheet','mult_compare');

%%
close all
graphOps = 1;
ygrid = [0 25 50 75 100];
if graphOps == 1
    pdfName = 'Average_backward_walking_by_group';
    axislabel = 'percent';
    rearranger = [5 6 2 3 4];
elseif graphOps == 2
    pdfName = 'Average_jump_rates_by_group';
    axislabel = 'percent';
%     rearranger = [4 5 2 3];
    rearranger = [5 6 2 3 4];
else
    pdfName = 'Average_reaching_by_group';
    axislabel = 'percent';
    rearranger = [5 6 2 3 4];
end
grpLabel = unique(dataExcelTable.CellTypes);
grpLabel = grpLabel(rearranger);
grpCt = numel(grpLabel);

xplot = [0 grpCt NaN];
plot(repmat(xplot,numel(ygrid),1)',repmat(ygrid',1,numel(xplot))','color',[1 1 1]*0.5)
hold on

for iterG = 1:grpCt
    if graphOps == 1
        penetrance = (dataExcelTable.Backing(strcmp(dataExcelTable.CellTypes,grpLabel{iterG})));
        ydata = penetrance.*100;
    elseif graphOps == 2
        penetrance = (dataExcelTable.Jumping(strcmp(dataExcelTable.CellTypes,grpLabel{iterG})));
        ydata = penetrance.*100;
    else
        penetrance = (dataExcelTable.reaching(strcmp(dataExcelTable.CellTypes,grpLabel{iterG})));
        ydata = penetrance.*100;
    end
    ydata(isnan(ydata)) = [];
    meanY = mean(ydata);
    stdY = std(ydata);
    if round(numel(ydata)/2)*2 ~= numel(ydata)
        ydata = cat(1,ydata(:),NaN);
    end
    
    spreadVal = 0.4;
    xdata = linspace(-spreadVal,spreadVal,numel(ydata))+iterG-0.5;
    plot([xdata(1) xdata(end)],[meanY meanY],'linewidth',2,'color','k')
    plot([mean(xdata) mean(xdata)],[meanY-stdY meanY+stdY],'linewidth',2,'color','k')
    spreadVal = 0.2;
    xdata = linspace(-spreadVal,spreadVal,numel(ydata))+iterG-0.5;
    plot(xdata,ydata,'.','markersize',22,'color','k')
    
end
axis normal, box off, axis off, hold on
fpos = get(gcf,'position');
set(gcf,'color',[1 1 1],'position',[2234 501 660 520])
apos = get(gca,'position');
set(gca,'position',([0.3 0.3 0.5 0.5]),'ylim',[min(ygrid) max(ygrid)])
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');

rowLabel = cellstr(num2str(ygrid'));
rowCt = numel(rowLabel);
for iterT = 1:rowCt
    text(-0.5,ygrid(iterT),rowLabel{iterT},'HorizontalAlignment','right',...
        'fontsize',12)
end

colCt = grpCt;
colPos = (1:colCt)-0.5;
for iterT = 1:colCt
    label = (grpLabel{iterT});
    text(colPos(iterT),-range(ygrid)*.15+min(ygrid),label,'HorizontalAlignment','right',...
        'fontsize',15,'rotation',45,'interpreter','none')
end
text(-1.5,mean(ygrid),axislabel,'rotation',90,'fontsize',15,'HorizontalAlignment','center')
text(-1.5,range(ygrid)*1.25+min(ygrid),pdfName,'rotation',0,'fontsize',15,'HorizontalAlignment','left',...
    'interpreter','none')
%%
addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
writeDir = 'Y:\Ming_RubinLab';
writePath = fullfile(writeDir,[pdfName '.pdf']);
export_fig(writePath)
writePath = fullfile(writeDir,[pdfName '.eps']);
export_fig(writePath)