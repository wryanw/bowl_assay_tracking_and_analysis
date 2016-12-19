
clear all
close all
clc

op_sys = system_dependent('getos');
if strfind(op_sys,'Microsoft Windows 7')
    dm11Dir = [filesep filesep 'tier2' filesep 'card'];
else
    dm11Dir = [filesep 'Volumes' filesep 'cardlab'];
end
exceldir = fullfile(dm11Dir,'Ming_RubinLab','LC16_bilateral_quantification');
excelname = 'LC16_Bilateral_activation.xlsx';
exceltableA = readtable(fullfile(exceldir,excelname),'Sheet','LC16','ReadRowNames',true);
rowNamesA = exceltableA.Properties.RowNames;
exceltableB = readtable(fullfile(exceldir,excelname),'Sheet','pBDPGAL4','ReadRowNames',true);
rowNamesB = exceltableB.Properties.RowNames;

%%
figure
backC = [1 1 1];
set(gcf,'units','normalized')
set(gcf,'position',[  0.1620    0.4150    0.4974    0.4750])
figPos = get(gcf,'position');
set(gcf,'color',backC,'position',figPos)
hax = axes;
haxPos = get(hax(1),'position');
growShiftAx = [-0.2,-0.2,0.1,-0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
hold on
xtickpos = [0 360];
ytickpos = [1 2 3]-0.5;
fontC = [0 0 0];
set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[0 3],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
plot([0 0],[0 3],'color',[0 0 0],'linewidth',2)
dataCell = table2array(exceltableA)
for iterL = 1:3
    fitdata = dataCell(:,iterL+3);
    conv = false;
    iterations = 2;
    options = statset('MaxIter',500);
    fitX = linspace(0,360,1000);
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
%         fitY = fitY/0.05;
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
ytickval = {'Bilateral','Unilateral','Empty'};
for iterL = 1:numel(ytickval)
    text(-400,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',20);
end
pdfName = 'LC16_bilateral_activation_turning.pdf';
text(-300,3.5,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
%
addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
writePath = fullfile(exceldir,pdfName);
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
% addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
% writePath = fullfile(exceldir,pdfName);
% export_fig(writePath)


