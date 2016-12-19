clear all
close all
clc
dataDirs = {'Y:\Ming_RubinLab\LC16_CsChrimson';'Y:\Ming_RubinLab\pBDPGAL4U_CsChrimson'};
xshift = 0;
xspace = 30;
xtickcell = cell(1,2);
xlinebase = (0:25:150);
ylines = zeros(size(xlinebase));
ylines = [ylines;ylines+2;ylines*NaN];
sheetNames = {'LC16','pBDPGAL4'};
for iterD = 1:2
    dataFiles = dir(fullfile(dataDirs{iterD},'*filteredData.mat'));
    
    initPosPath = fullfile(dataDirs{iterD},'initial_positions.mat');
    load(initPosPath)
    hPosInit = initPos.headPos;
    tPosInit = initPos.tailPos;
    dists = sqrt(sum([(hPosInit(:,1)-tPosInit(:,1)).^2,(hPosInit(:,2)-tPosInit(:,2)).^2],2));
    avgL = mean(dists);
    %assuming 2.5mm avg fly length
    pix2mm = @(x) x*(2.5/avgL);
    %50 fps recording rate
    frm2sec = @(x) x*50;
    pixL = 0.11;
    dataFiles = {dataFiles(:).name};
    dataLoading = cellfun(@(x) load(fullfile(dataDirs{iterD},x)),dataFiles,'uniformoutput',false);
    dataCt = numel(dataLoading);
    rowNames = cellfun(@(x) regexprep(x,'_filteredData.mat',''),dataFiles,'uniformoutput',false)';
    dataTable = table(cell(dataCt,1),cell(dataCt,1),cell(dataCt,1),'RowNames',...
        rowNames,'VariableNames',{'filtered_pos','filtered_vel','filtered_accel'});
    for i = 1:numel(rowNames)
        if isfield(dataLoading{i},'diffWin')
            dataTable.filtered_pos(rowNames{i}) = {dataLoading{i}.filtered_pos};
            dataTable.filtered_vel(rowNames{i}) = {dataLoading{i}.filtered_vel};
            dataTable.filtered_accel(rowNames{i}) = {dataLoading{i}.filtered_accel};
        else
            dataTable(rowNames{i},:) = dataLoading{i}.analyzer_record;
        end
        
    end

%     dataSet = cat(1,dataLoading{:});
%     dataCt = numel(dataSet);
    sideways = cell(dataCt,1);
    forwardBack = cell(dataCt,1);
    
    turning = cell(dataCt,1);
    origFilterWin = 19;
    part2cut = (origFilterWin-1)/2+1;
    excelData = zeros(dataCt,6);
    for iterV = 1:dataCt
        videoID = rowNames{iterV};
        dataAccel = dataTable.filtered_accel{videoID};
        dataVel = dataTable.filtered_vel{videoID};
        dataPos = dataTable.filtered_pos{videoID};
        dataVel(1:part2cut-1,:) = [];
        dataVel(end-part2cut:end,:) = [];
        dataPos(1:part2cut-1,:) = [];
        dataPos(end-part2cut:end,:) = [];
        
        [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
        [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
        relXvel = frm2sec(pix2mm(relXvel));
        relYvel = frm2sec(pix2mm(relYvel));
        forwardBack{iterV} = smooth(relXvel,25)';
        %     sideways{iterV} = smooth(abs(relYvel),45)';
        turning{iterV} = frm2sec(smooth(abs(dataVel(:,3)),25)');
        %     turning{iterV} = abs(dataVel(:,3))';
        vecFrag = (numel(relXvel)-50)/2;%because one second 'on' in middle (50 frames)
        excelData(iterV,1) = sum(forwardBack{iterV}(1:vecFrag));
        excelData(iterV,4) = sum(turning{iterV}(1:vecFrag));
        excelData(iterV,2) = sum(forwardBack{iterV}(vecFrag+1:end-vecFrag));
        excelData(iterV,5) = sum(turning{iterV}(vecFrag+1:end-vecFrag));
        excelData(iterV,3) = sum(forwardBack{iterV}(end-vecFrag+1:end));
        excelData(iterV,6) = sum(turning{iterV}(end-vecFrag+1:end));
    end
    rowNames = cellfun(@(x) x(1:end-17),dataFiles,'uniformoutput',false);
    varNames = {'distance_before','distance_during','distance_after',...
        'turning_before','turning_during','turning_after'};
    excelData = round(excelData);
    excelData = table(excelData(:,1),excelData(:,2),excelData(:,3),...
        excelData(:,4),excelData(:,5),excelData(:,6),...
        'RowNames',rowNames,'VariableNames',varNames);
%     writetable(excelData,'Z:\CARD LAB SHARE\Ming_temp\LC16_Bilateral_activation.xlsx',...
%         'Sheet',sheetNames{iterD},'WriteRowNames',true);
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
    if iterD == 1
        colr = 'b';
    else
        colr = 'g';
    end
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