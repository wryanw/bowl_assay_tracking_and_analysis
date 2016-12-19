clearvars -except unilateralTable
% dataDir = 'D:\rai_tracking\Rai2Ryan_all';
dataDir = 'Y:\Ming_RubinLab\LC10_quantification_newAvi';
dataList = dir(fullfile(dataDir,'*.mat'));
dataList = {dataList(:).name}';

vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileNames = fileNames(1:end-3);
fileCt = numel(fileNames);
% exceldir = 'D:\rai_tracking\Rai2Ryan_all';
% excelname = 'Stochastic_MDN.xlsx';
exceldir = 'Y:\Ming_RubinLab\LC10_quantification_newAvi';
excelname = 'LC10_unilateral_labeling_050715_Ming.xlsx';
exceltable = readtable(fullfile(exceldir,excelname),'Sheet','Sheet2','ReadRowNames',true);
patternOps = {'/'
    'Both'
    'Left'
    'Left > Right'
    'Left >> Right'
    'Right'
    'Right > Left'
    'Right >>Left'
    'no LC10'};
patternCode = {'NA'
    'B'
    'L'
    'L'
    'L'
    'R'
    'R'
    'R'
    'E'};
exceltable.pattern_Nern = strtrim(exceltable.pattern_Nern);
for iterC = 1:numel(patternOps)
    exceltable.pattern_Nern(strcmp(exceltable.pattern_Nern,patternOps{iterC})) = patternCode(iterC);
end
pattern2flip = ones(numel(exceltable.pattern_Nern),1);
pattern2flip(strcmp(exceltable.pattern_Nern,'L')) = -1;

%
showIndicatorPlots = 0;
showThetaScroll = 0;
sumPlot = 1;


plotMode = 2;
if plotMode == 1
    grpList = strtrim(exceltable.pattern_8categories);
    nametag = '_8categories';
elseif plotMode == 2
    grpList = strtrim(exceltable.pattern_4categories);
    nametag = '_4categories';
elseif plotMode == 3
    grpList = strtrim(exceltable.pattern_4categories_others);
    nametag = '_4categories_others';
else
    grpList = strtrim(exceltable.pattern_3categories);
    nametag = '_3categories';
end
grpOps = unique(grpList);
grpOps(strcmp(grpOps,'')) = [];
grpOps(strcmp(grpOps,'Others')) = [];
% grpOps = {'E','B','U'};
exptList = exceltable.Properties.RowNames;
grpCt = numel(grpOps);
dataCell = cell(grpCt,3);
flyCtTaly = zeros(grpCt,1);
colrs = {'r','g','b'};
for iterD = 1:fileCt
    videoID = fileNames{iterD};
    vidRefPart = videoID(strfind(videoID,'_v')+2:end);
    vidRef = vidRefPart(1:strfind(vidRefPart,'_')-1);
    metaTablePath = fullfile(dataDir,[videoID 'metaInformation.mat']);
    metaTable = load(metaTablePath);
    metaTable = metaTable.metaTable;
    indicatorPath = fullfile(dataDir,[videoID 'indicator.mat']);
    indicatorVars = load(indicatorPath);
    indicatorVector = indicatorVars.indicatorVector;
    deltaLights = indicatorVars.deltaLights;
    indicatorPos = indicatorVars.indicatorPos;
    lightsOn = find(deltaLights == 1);
    lightsOff = find(deltaLights == 0);
    if lightsOn(1) > lightsOff(1)
        lightsOff(1) = [];
    end
    if numel(lightsOn) > numel(lightsOff)
        lightsOn = lightsOn(1:numel(lightsOff));
    elseif numel(lightsOn) < numel(lightsOff)
        lightsOff = lightsOff(1:numel(lightsOn));
    end
    %     lightsOn = lightsOn(1:5);%restricts experiment count to 5
    %     lightsOff = lightsOff(1:5);
    
    onDur = abs(min(lightsOn-lightsOff));
    onDur = round(onDur/2);
    offDur = abs(min(lightsOff-lightsOn));
    if showIndicatorPlots
        disp(videoID)
        plot(deltaLights)
        hold all
        indicatorVector = (indicatorVector-min(indicatorVector(:)))/range(indicatorVector);
        plot(indicatorVector)
        uiwait(gcf)
        continue
    end
    flyCt = size(metaTable,1);
    rowNames = metaTable.Properties.RowNames;
    for iterF = 1:flyCt
        flyID = rowNames{iterF};
        flyRect = str2double(strsplit(flyID,'_'));
        flyRectY = [flyRect(4) flyRect(2)];
        otherID = rowNames{abs(iterF-3)};
        flyRect = str2double(strsplit(otherID,'_'));
        otherRectY = [flyRect(4) flyRect(2)];
        if max(flyRectY) < min(otherRectY)
            exptName = [vidRef 'a'];
        else
            exptName = [vidRef 'b'];
        end
        exptRef = find(strcmp(exptList,exptName));
        if isempty(exptRef)
            continue
        end
        grpRef = find(strcmp(grpOps,grpList{exptRef}));
        if isempty(grpRef)
            continue
        end
        
        posArray = metaTable.filtered_pos{flyID};
        xPos = posArray(:,1);
        yPos = posArray(:,2);
        thetaFilt = posArray(:,3);
        thetaRaw = metaTable.bot_theta{flyID};
        posRaw = metaTable.bot_centroid{flyID};
        
        hPosInit = metaTable.headPos(iterF,:);
        tPosInit = metaTable.tailPos(iterF,:);
        dists = sqrt(sum([(hPosInit(:,1)-tPosInit(:,1)).^2,(hPosInit(:,2)-tPosInit(:,2)).^2],2));
        avgL = mean(dists);
        pixL = 0.11;
        rad2deg = @(x) x.*(180/pi);
        %assuming 2.5mm avg fly length
        %     pix2mm = @(x) x*(2.5/avgL);
        pix2mm = @(x) x*pixL;
        %50 fps recording rate
        frm2sec = @(x) x*50;
        origFilterWin = 9;
        part2cut = (origFilterWin-1)/2+1;
%         dataPos = [posRaw thetaRaw];
%         dataVel = cat(1,[0 0 0],diff(dataPos));
        dataVel = metaTable.filtered_vel{flyID};
        dataPos = metaTable.filtered_pos{flyID};
        dataVel(1:part2cut-1,:) = NaN;
        dataVel(end-part2cut:end,:) = NaN;
        dataPos(1:part2cut-1,:) = NaN;
        dataPos(end-part2cut:end,:) = NaN;
        [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
        [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
        relXvel = (pix2mm(relXvel));
        relYvel = (pix2mm(relYvel));
        
        forwardBack = (relXvel);
        goodndcs = find(abs(forwardBack) < 1);
        forwardBack = interp1(goodndcs,forwardBack(goodndcs),(1:numel(forwardBack)),'linear','extrap');
        forwardBack = frm2sec(forwardBack);
        
        sideways = (relYvel)*pattern2flip(exptRef)*(-1);
        sideways(abs(sideways) > 2) = 0;
        
        turning = rad2deg(dataVel(:,3))*pattern2flip(exptRef);
        goodndcs = find(abs(turning) < 30);
        turning = interp1(goodndcs,turning(goodndcs),(1:numel(turning)),'linear','extrap');
        turning = frm2sec(turning)*(-1);
        
        instanceList = table2cell(exceltable(exptRef,:));
        
        turning(forwardBack > 60) = 0;
        turning(forwardBack < -30) = 0;
        forwardBack(forwardBack > 60) = 0;
        forwardBack(forwardBack < -30) = 0;
        
        for iterL = 1:numel(lightsOn)
            if iterL > numel(lightsOff)
                continue
            end
            if isnan(instanceList{iterL})
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'NA')
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'j')
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'c')
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'/')
                continue
            end
            dataBegin = lightsOn(iterL);
            dataEnd = dataBegin+onDur;
            if dataEnd > numel(forwardBack)
                continue
            end
            if sumPlot
                dataCell{grpRef,1} = [dataCell{grpRef,1};mean(forwardBack(dataBegin:dataEnd))];
                dataCell{grpRef,2} = [dataCell{grpRef,2};mean(sideways(dataBegin:dataEnd))];
                dataCell{grpRef,3} = [dataCell{grpRef,3};mean(turning(dataBegin:dataEnd))];
            else
                dataEnd = round((dataEnd-dataBegin)/2)+dataEnd;
                if dataEnd > numel(forwardBack)
                    continue
                end
                dataCell{grpRef,1} = [dataCell{grpRef,1};(forwardBack(dataBegin:dataEnd))'];
                dataCell{grpRef,2} = [dataCell{grpRef,2};(sideways(dataBegin:dataEnd))'];
                dataCell{grpRef,3} = [dataCell{grpRef,3};(turning(dataBegin:dataEnd))'];
            end
            
            %             subNames = metaTable.Properties.RowNames;
            %             exceltable.distance(subNames) = round(forwardBack*10)/10;
            %             exceltable.turning(subNames) = round(turning);
        end
        flyCtTaly(grpRef) = flyCtTaly(grpRef)+1;
        if showThetaScroll
            thetaOmit = abs(diff(thetaRaw)) > pi/8;
            thetaRaw(thetaOmit) = NaN;
            thetaInterp = interp1(find(~thetaOmit),thetaRaw(~thetaOmit),(1:numel(thetaRaw)),'pchip');
            thetaInput = thetaRaw;
            startCt = 1000;
            scrollVal = 1;
            scrollCt = (numel(thetaInput)-startCt)/scrollVal;
            for iterSc = 1:scrollCt-1
                dataRange = (iterSc*scrollVal:1000+iterSc*scrollVal);
                [xT,yT] = pol2cart(thetaInput(dataRange),(1:numel(dataRange))');
                %         xT(thetaOmit) = NaN;
                %             xT = xT(dataRange);
                %             yT = yT(dataRange);
                if iterSc == 1
                    %                 hP = plot(xT,yT,'.','markersize',4);
                    hP = plot(xT,yT);
                else
                    set(hP,'xdata',xT,'ydata',yT)
                    set(gca,'YLim',[-startCt,startCt],'XLim',[-startCt,startCt])
                end
                drawnow
            end
        end
    end
    
end

dataLabels = grpOps;

for iterL = 1:size(dataLabels,1)
%     ydata = dataCell{iterL,3};
%     xdata = (rand(numel(ydata),1)-0.5)/2+iterL;
%     plot(xdata,ydata,'.','markersize',20,'color',colrs{iterL})
%     hold on
    if ~exist('unilateralTable','var')
        unilateralTable = table(cell(6,1),cell(6,1),zeros(6,1),'VariableNames',{'CellType','data','flyCt'},...
            'RowNames',{'LC16_Unilateral','LC16_Bilateral','LC16_Empty',...
            'LC10_Unilateral','LC10_Bilateral','LC10_Empty'});
    end
    tableID = cat(2,'LC10_',dataLabels{iterL});
    unilateralTable.data(tableID) = dataCell(iterL,3);
    unilateralTable.flyCt(tableID) = flyCtTaly(iterL);
    unilateralTable.CellType(tableID) = {'LC10'};
end
%%
if sumPlot
    datarefs = [2 3 1];
    saveFigs = 0;
    ytickval = dataLabels;
    figure
    backC = [1 1 1];
    set(gcf,'units','normalized')
    set(gcf,'position',[  0.1620    0.4150    0.4974    0.4750])
    figPos = get(gcf,'position');
    set(gcf,'color',backC,'position',figPos)
    if plotMode > 1
        haxB = axes;
    end
    haxA = axes;
    haxPos = get(haxA,'position');
    growShiftAx = [-0.2,-0.2,0.0,-0.02];%grow width, grow height, shift X, shift Y
    haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
        haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
    xtickpos = [-360 0 360];
    ytickpos = (1:numel(datarefs))-0.5;
    fontC = [0 0 0];
    set(haxA,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
        'ticklength',[0.02 0.01],'YLim',[0 numel(datarefs)],'tickdir','out','xtick',xtickpos,...
        'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w','nextplot','add')
    plot([0 0],[0 numel(datarefs)],'color',[0 0 0],'linewidth',2,'parent',haxA)
    ydrop = -0.1;
    ytickpos = [1 2 3]-0.5;
    fontC = [0 0 0];
    set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
        'ticklength',[0.02 0.01],'YLim',[-0.2 numel(datarefs)],'tickdir','out','xtick',xtickpos,...
        'yticklabel',[],'ytick',ytickpos,'ycolor','w','xcolor','w')
    
    plot([0 0],[0 numel(datarefs)],'color',[0 0 0],'linewidth',2)
    xplot = cat(1,xtickpos,xtickpos,NaN(1,numel(xtickpos)));
    yplot = cat(1,zeros(1,numel(xtickpos))+ydrop/2,zeros(1,numel(xtickpos))+ydrop,zeros(1,numel(xtickpos)));
    plot(xplot,yplot,'k','linewidth',2,'parent',haxA(1))
    plot([xtickpos(1) xtickpos(end)],[ydrop ydrop]/2,'k','linewidth',2,'parent',haxA(1))
    
    if plotMode > 1
        growShiftAx = [-0.8,0,0.4,0];%grow width, grow height, shift X, shift Y
        haxPos2 = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
            haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
        
        set(haxB,'position',haxPos2,'color','none','box','off','Xlim',[0.5 3.5],...
            'YLim',[0 numel(datarefs)],'tickdir','out','xtick',[1 2 3],'nextplot','add',...
            'ytick',ytickpos,'ycolor','w','xcolor','w')
        
    end
    
    for iterL = 1:numel(datarefs)
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
        %         try
        %             paramEsts = mle(x,'pdf',pdf_evmixture,'start',start,'options',options);
        %             %     paramEsts = mle(x,'pdf',pdf_evmixture,'start',start,'lower',lb,'upper',ub,'options',options);
        %             evlike([paramEsts(2) paramEsts(4)],x);
        %             evlike([paramEsts(3) paramEsts(5)],x);
        %             fitX = linspace(1.1*min(x),1.1*max(x),200);
        %             pdfgrid = pdf_evmixture(fitX,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
        %         catch
        fitX = linspace(1.1*min(x),1.1*max(x),200);
        paramEsts = evfit(x);
        evlike([paramEsts(1) paramEsts(2)],x);
        pdfgrid = evpdf(fitX,paramEsts(1),paramEsts(2));
        %         end
        
        
        ydata = fitdata;
        xdata = (rand(numel(ydata),1)-0.25)/2+iterL-0.5;
        hbox = boxplot(gca,ydata,'positions',iterL-0.5,'boxstyle','outline',...
            'colors','k','notch','on','widths',0.75,'orientation','horizontal');
        set(hbox(7,:),'marker','.','markeredgecolor','k','markersize',12)
        continue
        plot(ydata,xdata,'.','color',[0 0 0],'markersize',16,'parent',haxA)
        fitY = pdfgrid;
        fitY = (fitY-min(fitY(:)))/(range(fitY)*1.1);
        plot(fitX,fitY+iterL-1,'color',[0 0 0],'LineWidth',3,'parent',haxA);
        clearvars obj
        medData = median(fitdata);
        plot([medData medData],[iterL-1 iterL-.2],'color',[.8 0 0],'linewidth',2,'parent',haxA)
        err = 1.57*(prctile(fitdata,75)-prctile(fitdata,25))/sqrt(numel(fitdata));
        plot([medData-err medData+err],[iterL-1 iterL-1],'color',[.8 0 0],'linewidth',2,'parent',haxA)
        
    end
    set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
        'ticklength',[0.02 0.01],'YLim',[-0.2 numel(datarefs)],'tickdir','out','xtick',xtickpos,...
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
    
    xtickval = cellstr(num2str(xtickpos'));
    set(get(haxA,'xlabel'),'units','normalized')
    xlabelPos = get(get(haxA,'xlabel'),'position');
    xlabelPos(1:2) = [mean(xtickpos) -0.5];
    for iterL = 1:numel(xtickval)
        text(xtickpos(iterL),xlabelPos(2),(xtickval{iterL}),'rotation',0,...
            'color',fontC,'horizontalalignment','center','fontsize',20,'parent',haxA);
    end
    
    for iterL = 1:numel(datarefs)
        text(-400,ytickpos(iterL),(ytickval{datarefs(iterL)}),'rotation',0,...
            'color',fontC,'horizontalalignment','right','fontsize',20,'parent',haxA);
    end
    pdfName = ['LC10_unilateral_activation_turning' nametag];
    text(-500,numel(datarefs)+0.25,pdfName,...
        'horizontalalignment','left','interpreter','none',...
        'rotation',0,'color',fontC,'fontsize',20,'parent',haxA);
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    if saveFigs
        %%
        writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.pdf']);
        export_fig(writePath)
        writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.eps']);
        export_fig(writePath)
    end
    %%
    plot([0.6 3+0.5],[0 0],'color','k','linewidth',1)
    hold on
    datarefs = (1:grpCt);
    for iterL = 1:grpCt
        %     ydata = dataCell{datarefs(iterL),1};
        %     xdata = (rand(numel(ydata),1)-0.5)/2+iterL;
        %     plot(xdata,ydata,'.','color','k','markersize',18)
        fitdata = dataCell{datarefs(iterL),1};
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
        %         try
        %             paramEsts = mle(x,'pdf',pdf_evmixture,'start',start,'options',options);
        %             %     paramEsts = mle(x,'pdf',pdf_evmixture,'start',start,'lower',lb,'upper',ub,'options',options);
        %             evlike([paramEsts(2) paramEsts(4)],x);
        %             evlike([paramEsts(3) paramEsts(5)],x);
        %             fitX = linspace(1.1*min(x),1.1*max(x),200);
        %             pdfgrid = pdf_evmixture(fitX,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
        %         catch
        fitX = linspace(1.1*min(x),1.1*max(x),200);
        paramEsts = evfit(x);
        evlike([paramEsts(1) paramEsts(2)],x);
        pdfgrid = evpdf(fitX,paramEsts(1),paramEsts(2));
        %         end
        ydata = fitdata;
        xdata = (rand(numel(ydata),1)-0.5)/2+iterL;
        plot(xdata,ydata,'.','color',[0 0 0],'markersize',16)
        fitY = pdfgrid;
        fitY = (fitY-min(fitY(:)))/(range(fitY)*1.1);
        plot(fitY+iterL-0.5,fitX,'color',[0 0 0],'LineWidth',3);
        
        medData = median(fitdata);
        plot([iterL-0.4 iterL],[medData medData],'color',[.8 0 0],'linewidth',2)
        err = 1.57*(prctile(fitdata,75)-prctile(fitdata,25))/sqrt(numel(fitdata));
        plot([iterL-0.4 iterL-0.4],[medData-err medData+err],'color',[.8 0 0],'linewidth',2)
        
    end
    set(gcf,'units','normalized')
    set(gcf,'color',backC)
    set(gcf,'position',[0.332    0.4433    0.321    0.4708])
    ytickpos = [-5 0 10];
    xtickpos = (1:grpCt);
    fontC = [0 0 0];
    set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
        'ticklength',[0.02 0.01],'tickdir','out','xtick',xtickpos,'YLim',[ytickpos(1) ytickpos(end)],...
        'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
    xtickval = dataLabels;
    set(get(gca,'xlabel'),'units','normalized')
    xlabelPos = get(get(gca,'xlabel'),'position');
    xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.07];
    for iterL = 1:numel(xtickval)
        text(xtickpos(iterL),xlabelPos(2),(xtickval{datarefs(iterL)}),'rotation',45,...
            'color',fontC,'horizontalalignment','center','fontsize',16);
    end
    ytickval = cellstr(num2str(ytickpos'));
    for iterL = 1:numel(ytickval)
        text(0,ytickpos(iterL),(ytickval{iterL}),'rotation',0,...
            'color',fontC,'horizontalalignment','center','fontsize',20);
    end
    
    
    pdfName = ['LC10_unilateral_activation_backing' nametag];
    text(-.5,26,pdfName,...
        'horizontalalignment','left','interpreter','none',...
        'rotation',0,'color',fontC,'fontsize',20);
    %
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    if saveFigs
        writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.pdf']);
        export_fig(writePath)
        writePath = fullfile('Y:\Ming_RubinLab',[pdfName '.eps']);
        export_fig(writePath)
    end
else
    %%
    xshift = 0;
    xspace = 30;
    xtickcell = cell(1,numel(grpOps));
    xlinebase = (0:25:150);
    ylines = zeros(size(xlinebase));
    ylines = [ylines;ylines+2;ylines*NaN];
    colrOps = {'r','g','b'};
    for iterD = 1:numel(grpOps)
        %     xlines = repmat(xlinebase,3,1)+xshift;
        xlines = repmat(xlinebase,3,1);
        xtickcell{iterD} = unique(xlines);
        plot(xlines,ylines,'color','k','linewidth',2)
        hold on
        %     plot([0 150]+xshift,[0.5 0.5],'color',[.5 .5 .5],'linewidth',2)
        %     plot([0 150]+xshift,[1 1],'color',[.5 .5 .5],'linewidth',2)
        plot([0 150],[0.5 0.5],'color',[.5 .5 .5],'linewidth',2)
        plot([0 150],[1 1],'color',[.5 .5 .5],'linewidth',2)
        plotY = dataCell{iterD,1}';
        % plotY = cat(1,sideways{:});
        % rangeA = range(plotY(:))
        rangeA = 120;
        plotY = (plotY)/rangeA+0.5;
        plotYm = mean(plotY,1);
        plotXm = (1:numel(plotYm))+part2cut+xshift;
        colr = colrOps{iterD};
        plot(plotXm,plotYm,'color',colr,'linewidth',5)
        err = std(plotY);
        plotY = [err;-err]+[plotYm;plotYm];
        plotY = [err;-err;NaN(size(err))]+[plotYm;plotYm;plotYm];
        plotX = repmat((1:size(plotY,2)),size(plotY,1),1)+part2cut+xshift;
        %     plot(plotX',plotY','color',colr,'linewidth',2)
        plot(plotX(:),plotY(:),'color',colr,'linewidth',1)
        
        plotYb = dataCell{iterD,3}';
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
end

%%
for iterD = 1:3
    yOff = iterD*2;
    ydata = dataCell{iterD,1}'+yOff;
    xdata = repmat(1:size(ydata,1),size(ydata,2),1)';
    plot(xdata,ydata,'color',colrOps{iterD})
    hold on
    plot([1 size(ydata,1)],[yOff yOff],'k')
end