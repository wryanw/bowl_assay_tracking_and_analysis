clear all
close all
clc
archDir = '\\tier2\card';
destdir = fullfile(archDir,'Ming_RubinLab');
dataDir{1} = fullfile(archDir,'Ming_RubinLab','LC16_CsChrimson');
dataDir{2} = fullfile(archDir,'Ming_RubinLab','pBDPGAL4U_CsChrimson');
dataTable = cell(2,1);
for iterLoad = 1:2
    vidList = dir(fullfile(dataDir{iterLoad},'*.avi'));
    vidList = {vidList(:).name};
    dataFiles = dir(fullfile(dataDir{iterLoad},'*filteredData.mat'));
    dataFiles = {dataFiles(:).name};
    dataLoading = cellfun(@(x) load(fullfile(dataDir{iterLoad},x)),dataFiles,'uniformoutput',false);
    rowNames = cellfun(@(x) regexprep(x,'_filteredData.mat',''),dataFiles,'uniformoutput',false)';
    dataCt = numel(rowNames);
    dataTable{iterLoad} = table(cell(dataCt,1),cell(dataCt,1),cell(dataCt,1),'RowNames',...
        rowNames,'VariableNames',{'filtered_pos','filtered_vel','filtered_accel'});
    for i = 1:numel(rowNames)
        if isfield(dataLoading{i},'diffWin')
            dataTable{iterLoad}.filtered_pos(rowNames{i}) = {dataLoading{i}.filtered_pos};
            dataTable{iterLoad}.filtered_vel(rowNames{i}) = {dataLoading{i}.filtered_vel};
            dataTable{iterLoad}.filtered_accel(rowNames{i}) = {dataLoading{i}.filtered_accel};
        else
            dataTable{iterLoad}(rowNames{i},:) = dataLoading{i}.analyzer_record;
        end
    end
end
%%
dataCell = cell(2,6);
colrs = {'r','g','b'};
for grpRef = 1:2
    pixL = 0.11;
    %assuming 2.5mm avg fly length
%     pix2mm = @(x) x*(2.5/avgL);
    pix2mm = @(x) x*pixL;
    %50 fps recording rate
    frm2sec = @(x) x*50;
    dataCt = size(dataTable{grpRef},1);
%     return
    origFilterWin = 19;
    part2cut = (origFilterWin-1)/2+1;
    excelData = zeros(dataCt,6);
    for iterV = 1:dataCt
        dataVel = dataTable{grpRef}.filtered_vel{iterV};
        dataPos = dataTable{grpRef}.filtered_pos{iterV};
        dataVel(1:part2cut-1,:) = [];
        dataVel(end-part2cut:end,:) = [];
        dataPos(1:part2cut-1,:) = [];
        dataPos(end-part2cut:end,:) = [];
        [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
        [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
        relXvel = (pix2mm(relXvel));
        relYvel = (pix2mm(relYvel));
        forwardBack = (relXvel);
        sideways = (relYvel);
        turning = (abs(dataVel(:,3)/(pi/180)));
        
        vecFrag = (numel(relXvel)-50)/2;%because one second 'on' in middle (50 frames)
        compFiltWin = 4;
        excelData(iterV,1) = nansum(forwardBack(1:vecFrag));
        excelData(iterV,4) = nansum(turning(1:vecFrag));
        excelData(iterV,2) = nansum(forwardBack(vecFrag+1+compFiltWin:end-vecFrag-compFiltWin));
        excelData(iterV,5) = nansum(turning(vecFrag+1+compFiltWin:end-vecFrag-compFiltWin));
        excelData(iterV,3) = nansum(forwardBack(end-vecFrag+1:end));
        excelData(iterV,6) = nansum(turning(end-vecFrag+1:end));
        
%         forwardBack(iterV) = sum(smooth(relXvel,25));
%         sideways(iterV) = sum(smooth(abs(relYvel),45));
%         turning(iterV) = (sum((smooth((dataVel(:,3)/(pi/180)),25))));
    end
%     ydata = turning;
%     xdata = (rand(numel(ydata),1)-0.5)/2+grpRef;
%     plot(xdata,ydata,'.','markersize',20,'color',colrs{grpRef})
%     hold on
    
    dataCell{grpRef,1} = excelData(:,1);
    dataCell{grpRef,2} = excelData(:,2);
    dataCell{grpRef,3} = excelData(:,3);
    dataCell{grpRef,4} = excelData(:,4);
    dataCell{grpRef,5} = excelData(:,5);
    dataCell{grpRef,6} = excelData(:,6);
end
% writetable(exceltable,fullfile(exceldir,excelname),'Sheet','display_video','WriteRowNames',true);
%% Figure for publishing !!!!!!!!!!!!!!!!!!!!!!!!!!!
saveFig = 1;
plotMode = 1; % 1 - show dots, 2 - hide dots
dataMode = 2;
% 1 - forward v back
% 2 - turning

datarefs = [1 3 5 2 4 6];
if dataMode == 2
    datarefs = datarefs+6;
end

dataCt = 6;
px = cell(dataCt,1);
py = cell(dataCt,1);
eYall = zeros(3,dataCt); mYall = zeros(3,dataCt);

for iterL = 1:dataCt
    
    ydata = dataCell{datarefs(iterL)};
    x = (rand(numel(ydata),1)-0.5)*0.25+iterL;%+datarefs(iterL);
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
growShiftAx = [-0.2,-0.2,0.1,-0.05];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on
fontC = [0 0 0];
if dataMode == 1
    pdfName = 'LC16_bilateral_activation_forwardVSback_scatter.pdf';
else
    pdfName = 'LC16_bilateral_activation_turning_scatter.pdf';
end

x = cat(1,px{:});
ydata = cat(1,py{:});
if plotMode == 1
    plot(x,ydata,'.','color','k','markersize',18)
end
% hold on

medianX = repmat((1:dataCt),3,1);
medianX = medianX+repmat([-.3;.3;0],1,dataCt);
errorX = repmat((1:dataCt),3,1)-.3;
if plotMode == 1
    plot([medianX;errorX],[mYall;eYall],'color','k','linewidth',4);
else
    plot(medianX(3,:),mYall(1,:),'.','color','k','markersize',36);
    plot(errorX,eYall,'color','k','linewidth',2);
end

if dataMode == 1
    ytickpos = [-30 0 30];
    plot(get(gca,'xlim'),[0 0],'color','k','linewidth',1)
else
    ytickpos = [0 180 360];
    plot(get(gca,'xlim'),[0 0],'color','k','linewidth',1)
    plot(get(gca,'xlim'),[180 180],'color','k','linewidth',1)
    plot(get(gca,'xlim'),[360 360],'color','k','linewidth',1)
end
xtickpos = [1 2 3 4 5 6];
text(0,max(ytickpos)+0.2*range(ytickpos),pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[ytickpos(1) ytickpos(end)],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
xtickval = {'Before','During','After','Before','During','After'};

set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.07];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),(xtickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',18);
end
ytickval = cellstr(num2str(ytickpos'));
for iterL = 1:numel(ytickval)
    text(0,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end

if saveFig
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile(destdir,pdfName);
    export_fig(writePath)
    writePath = fullfile(destdir,[pdfName(1:end-4) '.eps']);
    export_fig(writePath)
end
%%
figure
datarefs = [1 3 5 2 4 6]+6;
backC = [1 1 1];
set(gcf,'units','normalized')
set(gcf,'position',[  0.1620    0.4150    0.4974    0.4750])
figPos = get(gcf,'position');
set(gcf,'color',backC,'position',figPos)
hax = zeros(1,dataCt);
hax(1) = axes;
haxPos = get(hax(1),'position');
growShiftAx = [-0.2,-0.2,0.1,-0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on
xtickpos = [-360 0 360];
ytickpos = [1 2 3]-0.5;
fontC = [0 0 0];

for iterL = 1:dataCt
    fitdata = dataCell{datarefs(iterL)};
    conv = false;
    iterations = 3;
    options = statset('MaxIter',500);
    fitX = linspace(-360,360,1000);
    while ~conv && iterations > 1
        iterations = iterations - 1;     int_cond = [];
        for iterI = 1:(iterations - 1);
            int_cond = [int_cond;ones(round(length(fitdata)/iterations),1)*iterI];
        end
        if isempty(iterI), iterI = 0; end
        int_cond = [int_cond;ones((length(fitdata) - length(int_cond)),1)*(iterI+1)];
        try
            obj = gmdistribution.fit(fitdata,iterations,'Options',options,'Start',int_cond); %test to see if 3 seperate peeks
            conv = obj.Converged;
        catch ME
            getReport(ME)
        end
    end
    %         warning(orig_state)
    if exist('obj','var')
        fitY = pdf(obj,fitX');
        fitY = (fitY-min(fitY(:)))/range(fitY)*0.9;
        plot(fitX,fitY+iterL-1,'color',[0 0 0],'LineWidth',3);
        clearvars obj
    end
end

set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[0 dataCt],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
xtickval = cellstr(num2str(xtickpos'));
set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) -0.5];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),(xtickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end
plot([xtickpos(2) xtickpos(2)]+10,[0 dataCt],'color',[0 0 0],'linewidth',2)
ytickval = {'Bilateral','Unilateral','Empty'};
for iterL = 1:numel(ytickval)
    text(-400,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',20);
end
pdfName = 'LC16_bilateral_activation_turning.pdf';
text(-300,dataCt+.5,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
%
addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
writePath = fullfile(destdir,pdfName);
% export_fig(writePath)

%%
plot([0.6 3+0.5],[0 0],'color','k','linewidth',1)
hold on
datarefs = [2 1 3];
for iterL = 1:3
%     ydata = dataCell{datarefs(iterL),1};
%     xdata = (rand(numel(ydata),1)-0.5)/2+iterL;
%     plot(xdata,ydata,'.','color','k','markersize',18)
    fitdata = dataCell{datarefs(iterL),1};
    conv = false;
    iterations = 3;
    options = statset('MaxIter',500);
    fitX = linspace(-360,360,1000);
    while ~conv && iterations > 1
        iterations = iterations - 1;     int_cond = [];
        for iterI = 1:(iterations - 1);
            int_cond = [int_cond;ones(round(length(fitdata)/iterations),1)*iterI];
        end
        if isempty(iterI), iterI = 0; end
        int_cond = [int_cond;ones((length(fitdata) - length(int_cond)),1)*(iterI+1)];
        try
            obj = gmdistribution.fit(fitdata,iterations,'Options',options,'Start',int_cond); %test to see if 3 seperate peeks
            conv = obj.Converged;
        catch ME
            getReport(ME)
        end
    end
    %         warning(orig_state)
    if exist('obj','var')
        fitY = pdf(obj,fitX');
        fitY = fitY/range(fitY);
        plot(fitY+iterL-0.5,fitX,'color',[0 0 0],'LineWidth',3);
%         plot(xdata,ydata,'color','k','markersize',18)
        clearvars obj
    end
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
    text(xtickpos(iterL),xlabelPos(2),(xtickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',16);
end
ytickval = cellstr(num2str(ytickpos'));
for iterL = 1:numel(ytickval)
    text(0,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end


pdfName = 'LC16_bilateral_activation_backing.pdf';
text(-0.75,26,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
%
addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
writePath = fullfile(destdir,pdfName);
% export_fig(writePath)


