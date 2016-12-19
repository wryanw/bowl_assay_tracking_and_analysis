clear all
close all
clc
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
datarefs = [3 1 2];
figure
backC = [1 1 1];
set(gcf,'units','normalized')
set(gcf,'position',[  0.1620    0.4150    0.4974    0.4750])
figPos = get(gcf,'position');
set(gcf,'color',backC,'position',figPos)
% hax = zeros(1,dataCt);
hax(1) = axes;
haxPos = get(hax(1),'position');
growShiftAx = [-0.2,-0.2,0.1,-0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on
xtickpos = [-360 0 360];
ytickpos = [1 2 3]-0.5;
fontC = [0 0 0];
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[0 3],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
plot([0 0],[0 3],'color',[0 0 0],'linewidth',2)
for iterL = 1:3
    fitdata = dataCell{datarefs(iterL),3};
    x = fitdata;
    pdf_evmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
        p*evpdf(x,mu1,sigma1) + (1-p)*evpdf(x,mu2,sigma2);
    
    pStart = .5;
    muStart = quantile(x,[.25 .75]);
    sigmaStart = sqrt(var(x)-0.25*diff(muStart).^2);
    start = [pStart muStart sigmaStart sigmaStart];
    
    lb = [0 -Inf -Inf 0 0];
    ub = [1 Inf Inf Inf Inf];
    
    options = statset('MaxIter',500,'MaxFunEvals',600);
    try
        paramEsts = mle(x,'pdf',pdf_evmixture,'start',start,'options',options);
        %     paramEsts = mle(x,'pdf',pdf_evmixture,'start',start,'lower',lb,'upper',ub,'options',options);
        
        evlike([paramEsts(2) paramEsts(4)],x)
        evlike([paramEsts(3) paramEsts(5)],x)
        
        xgrid = linspace(1.1*min(x),1.1*max(x),200);
        pdfgrid = pdf_evmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
    catch
        params = evfit(x);
        pdfgrid = evpdf(sort(x),params(1),params(2));
        xgrid = sort(x);
    end
    
    ydata = fitdata;
    xdata = (rand(numel(ydata),1)-0.5)/2+iterL-0.5;
    plot(ydata,xdata,'.','color',[0 0 0],'markersize',16)
    fitY = pdfgrid;
    fitY = fitY/range(fitY)*0.8;
    plot(xgrid,fitY+iterL-1,'color',[0 0 0],'LineWidth',3);
    continue
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
        ydata = fitdata;
        xdata = (rand(numel(ydata),1)-0.5)/2+iterL-0.5;
        plot(ydata,xdata,'.','color',[0 0 0],'markersize',16)
        fitY = pdf(obj,fitX');
        fitY = fitY/range(fitY);
        plot(fitX,fitY+iterL-1,'color',[0 0 0],'LineWidth',3);
        clearvars obj
    end
end
xtickval = cellstr(num2str(xtickpos'));
set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) -0.5];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),(xtickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end
ytickval = {'Unilateral','Bilateral','Empty'};
for iterL = 1:numel(ytickval)
    text(-400,ytickpos(iterL),(ytickval{datarefs(iterL)}),'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',20);
end
pdfName = 'LC16_unilateral_activation_turning.pdf';
text(-300,3.5,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
%
% addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
% writePath = fullfile(exceldir,pdfName);
% export_fig(writePath)
% writePath = fullfile(exceldir,[pdfName(1:end-4) '.eps']);
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
        ydata = fitdata;
        xdata = (rand(numel(ydata),1)-0.5)/2+iterL;
        plot(xdata,ydata,'.','color',[0 0 0],'markersize',16)
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
xtickval = {'Unilateral','Bilateral','Empty'};
set(get(gca,'xlabel'),'units','normalized')
xlabelPos = get(get(gca,'xlabel'),'position');
xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.07];
for iterL = 1:numel(xtickval)
    text(xtickpos(iterL),xlabelPos(2),(xtickval{datarefs(iterL)}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',16);
end
ytickval = cellstr(num2str(ytickpos'));
for iterL = 1:numel(ytickval)
    text(0,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','center','fontsize',20);
end


pdfName = 'LC16_unilateral_activation_backing.pdf';
text(-0.75,26,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
%
% addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
% writePath = fullfile(exceldir,pdfName);
% export_fig(writePath)
% writePath = fullfile(exceldir,[pdfName(1:end-4) '.eps']);
% export_fig(writePath)

