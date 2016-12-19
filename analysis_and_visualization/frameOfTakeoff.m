close all
excelDir = 'Y:\Ming_RubinLab';
excelName = 'mingPaperFinalList.xlsx';
sheetName = 'dataLookup';
lookExcelTable = readtable(fullfile(excelDir,excelName),'Sheet',sheetName);
sheetName = 'dataTable';
dataExcelTable = readtable(fullfile(excelDir,excelName),'Sheet',sheetName);
%%
maxCt = max(graphDataStruct.dataCount);
varNames = dataExcelTable.Properties.VariableNames;
blankData = repmat({NaN(maxCt,1)},1,numel(varNames));
dataExcelTable = table(blankData{:},'VariableNames',varNames);

labelOps = lookExcelTable.groups;
for iterD = 1:dataCt
    dataID = graphDataStruct.dataID{iterD};
    dataVals = graphDataStruct.returnData{iterD};
    dataVals = cat(1,dataVals(:),NaN(maxCt-numel(dataVals),1));
    datandx = strcmp(dataID,labelOps);
    if max(datandx)
        dataRowName = lookExcelTable.label(datandx);
        dataExcelTable.(dataRowName{1}) = dataVals;
        %             num2str(jumpCount)
    end
end
%%
close all
clearvars
excelDir = 'Y:\Ming_RubinLab';
excelName = 'mingPaperFinalList.xlsx';
sheetName = 'dataTable';
dataExcelTable = readtable(fullfile(excelDir,excelName),'Sheet',sheetName);
sheetName = 'dataLookup';
lookExcelTable = readtable(fullfile(excelDir,excelName),'Sheet',sheetName);
%%
grpLabel = lookExcelTable.groups;
rearranger = [2 3 4 5];
grpLabel = grpLabel(rearranger);
grpCt = numel(grpLabel);
ygrid = [0 50 100 150];
pdfName = 'Time_from_photoactivation_to_end_of_jump';
axislabel = 'ms';
xplot = [0 grpCt NaN];
plot(repmat(xplot,numel(ygrid),1)',repmat(ygrid',1,numel(xplot))','color',[1 1 1]*0.5)
hold on

for iterG = 1:grpCt
    
    ydata = dataExcelTable(:,strcmp(lookExcelTable.groups,grpLabel{iterG}));
    ydata = cell2mat(table2cell(ydata));
    ydata(isnan(ydata)) = [];
    q1 = prctile(ydata,25); q2 = prctile(ydata,50);
    q3 = prctile(ydata,75); n = numel(ydata);
    err = 1.57*(q3-q1)/sqrt(n);
    meanY = q2;
    stdY = err;
    
    spreadVal = 0.4;
    xdata = linspace(-spreadVal,spreadVal,numel(ydata))+iterG-0.5;
    plot([xdata(1) xdata(end)],[meanY meanY],'linewidth',2,'color','k')
    plot([min(xdata) min(xdata)],[meanY-stdY meanY+stdY],'linewidth',2,'color','k')
    spreadVal = 0.2;
    xdata = linspace(-spreadVal,spreadVal,numel(ydata))+iterG-0.5;
    plot(xdata,ydata,'.','markersize',12,'color','k')
    
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
    label = lookExcelTable.label(strcmp(lookExcelTable.groups,grpLabel{iterT}));
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