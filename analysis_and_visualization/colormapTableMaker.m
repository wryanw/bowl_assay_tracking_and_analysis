refsMode = 3;
% 1 - save only
% 2 - load

excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2_swapped.xlsx';
excelTable = readtable(fullfile(excelDir,excelName),'Sheet','heatmap_compiled');
prctVal = 97.7;%%%%%%%%%%%%%%%%%%%%%%%
%%
figure
rowNames = excelTable.CellTypes;
rowCt = numel(rowNames);
% dataRefs = {'Jumping_heatmap_value','Turning_heatmap_value','Backing_heatmap_value','Reaching_heatmap'};
dataRefs = excelTable.Properties.VariableNames;
dataRefs = dataRefs(2:end-1);
% dataRefs = {'Jumping_heatmap_value','Turning_median','Backing_median','Reaching_heatmap'};
colCt = numel(dataRefs);

rescale = 0;
dims = 50;
stretch = 0.7;%xy - less than one make wider
% colorstr = 'bone';
% colorstr = 'jet';
% colorstr = 'hot';
colorstr = 'gray';
colorfun = str2func(colorstr);
binCt = 3;
resultIm = [];
padmag = 2;
% remapper = (colorfun(256).*255);
% remapper = (remapper/(255/binCt))*(255/binCt);
remapper = uint8(repmat(linspace(-0.5,binCt+0.5,256)',1,3))*(255/binCt);
remapper = double(remapper).^0.8;
remapper = (remapper-min(remapper(:)))/range(remapper(:));
remapper = uint8(remapper*255);
for iterC = 1:5
    iterData = excelTable.(dataRefs{iterC});
    iterData = abs(iterData-1);
    if rescale == 1
        iterData = (iterData-min(iterData(:)))/range(iterData(:));
    end
    rowIm = [];
    for iterR = 1:rowCt
        intIm = repmat(uint8(iterData(iterR)*255),dims,round(dims/stretch));
        mapIm = uint8(repmat(zeros(dims,round(dims/stretch)),[1 1 3]));
        for iterM = 1:3
            mapIm(:,:,iterM) = intlut(intIm,remapper(:,iterM));
        end
        mapIm(end-padmag+1:end,:,:) = 0;
        mapIm(:,end-padmag+1:end,:) = 0;
        mapIm(1:padmag,:,:) = 0;
        mapIm(:,1:padmag,:) = 0;
        rowIm = [rowIm;mapIm];
        
    end
    resultIm = [resultIm rowIm];
end
resultIm = padarray(resultIm,[padmag padmag],0);
spaceIm = rowIm*0+255;
spaceIm(spaceIm ~= 255) = 255;
spaceIm = padarray(spaceIm,[padmag padmag],255);
resultIm = [resultIm spaceIm];
resultImB = resultIm;
resultIm = [];
for iterC = 6:10
    iterData = excelTable.(dataRefs{iterC});
    iterData = abs(iterData-1);
    if rescale == 1
        iterData = (iterData-min(iterData(:)))/range(iterData(:));
    end
    rowIm = [];
    for iterR = 1:rowCt
        intIm = repmat(uint8(iterData(iterR)*255),dims,round(dims/stretch));
        mapIm = uint8(repmat(zeros(dims,round(dims/stretch)),[1 1 3]));
        for iterM = 1:3
            mapIm(:,:,iterM) = intlut(intIm,remapper(:,iterM));
        end
        if excelTable.TotalFlies_FlyPez(iterR) < 15 && iterC > 7
            mapIm = mapIm*255;
            diagNdcs = (linspace(1,size(mapIm,1),size(mapIm,2)));
            diagVal = 0.5;
            for iterDiag = 1:size(mapIm,2)
                rowDiag = round(diagNdcs(iterDiag)-diagVal:diagNdcs(iterDiag)+diagVal);
                colDiag = round(iterDiag-diagVal:iterDiag+diagVal);
                rowDiag(rowDiag < 1) = 1;
                colDiag(colDiag < 1) = 1;
                rowDiag(rowDiag > size(mapIm,1)) = size(mapIm,1);
                colDiag(colDiag > size(mapIm,2)) = size(mapIm,2);
                mapIm(round(rowDiag),round(colDiag),:) = 150;
                mapImB = imresize(mapIm,2);
                mapIm = imresize(mapImB,size(mapIm(:,:,1)));
            end
        end
        mapIm(end-padmag+1:end,:,:) = 0;
        mapIm(:,end-padmag+1:end,:) = 0;
        mapIm(1:padmag,:,:) = 0;
        mapIm(:,1:padmag,:) = 0;
        rowIm = [rowIm;mapIm];
        
    end
    resultIm = [resultIm rowIm];
end
resultIm = padarray(resultIm,[padmag padmag],0);
resultIm = [resultImB resultIm];

binCt = binCt+1;
iterData = fliplr(linspace(0,1,binCt));
rowIm = [];
for iterR = 1:binCt
    intIm = repmat(uint8(iterData(iterR)*255),dims,round(dims/stretch));
    mapIm = uint8(repmat(zeros(dims,round(dims/stretch)),[1 1 3]));
    for iterM = 1:3
        mapIm(:,:,iterM) = intlut(intIm,remapper(:,iterM));
    end
    rowIm = [rowIm mapIm];
end
rowIm(end-padmag+1:end,:,:) = 0;
rowIm(:,end-padmag+1:end,:) = 0;
rowIm(1:padmag,:,:) = 0;
rowIm(:,1:padmag,:) = 0;
rowIm = uint8([rowIm+255;rowIm]);
rowIm = [zeros(size(rowIm,1),size(resultIm,2)-size(rowIm,2),3)+255 rowIm];
resultIm = [resultIm;rowIm];
resultIm = padarray(resultIm,[15 15],255);

imshow(resultIm)

% axis normal, box off, axis off, hold on
hold on
fpos = get(gcf,'position');
set(gcf,'color',[1 1 1],'position',round(fpos.*[0.8 0.8 1.5 1]))
apos = get(gca,'position');
set(gca,'position',(apos.*[1.2 0.5 0.9 0.8]))
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
rowPos = linspace(ylim(1),ylim(2)*0.914,rowCt+1);
rowPos = (rowPos(1:end-1)+(diff(rowPos)/2))*0.97+0.015*ylim(2);
for iterT = 1:rowCt
    text(-5,rowPos(iterT),rowNames{iterT},'HorizontalAlignment','right',...
        'fontsize',12)
    %     text(xlim(2)+40,rowPos(iterT),num2str(excelTable.Total_flies(iterT)),'HorizontalAlignment','right',...
    %         'fontsize',12)
end

dataRefs = [dataRefs(1:5) {''} dataRefs(6:10)];
colCt = numel(dataRefs);
colPos = linspace(xlim(1),xlim(2),colCt+1);
colPos = (colPos+mean(diff(colPos)/2))*0.97;
for iterT = 1:numel(dataRefs)
    label = strsplit(dataRefs{iterT},'_');
    text(colPos(iterT),-15,label{1},'HorizontalAlignment','left',...
        'fontsize',15,'rotation',45)
end

colormapNames = round(linspace(0,100,binCt+1));
colormapNames = (colormapNames);
colorPos = linspace(xlim(1),xlim(2),colCt+1);
colorPos = colorPos(end-binCt:end)*0.97;
for iterT = 1:binCt+1
    text(colorPos(iterT),ylim(2)+5,num2str(colormapNames(iterT)),'HorizontalAlignment','left',...
        'fontsize',12)
    %     text(xlim(2)+40,rowPos(iterT),num2str(excelTable.Total_flies(iterT)),'HorizontalAlignment','right',...
    %         'fontsize',12)
end


if max(cellfun(@(x) ~isempty(strfind(x,'median')),dataRefs))
    pdfName = 'colormapTable_medians';
else
    %     pdfName = 'colormapTable_logPvalues';
    pdfName = 'colormapTable_percent';
end
pdfName = cat(2,pdfName,'_FlyPez_',regexprep(num2str(prctVal),'\.','pt'),'prctile_500ms_',num2str(binCt),'bins');
if rescale == 1
%     pdfName = [pdfName '_rescaled'];
else
%     pdfName = [pdfName '_rawValues'];
end
text(-120,-300,pdfName,'interpreter','none','fontsize',20);
writePath = fullfile(excelDir,[pdfName '.pdf']);
export_fig(writePath,'-nocrop')
writePath = fullfile(excelDir,[pdfName '.eps']);
export_fig(writePath,'-nocrop')