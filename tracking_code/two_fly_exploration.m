clearvars
dataDir = 'D:\rai_tracking\Rai2Ryan_all';
% dataDir = 'Y:\Ming_RubinLab\LC10_quantification_newAvi';
dataList = dir(fullfile(dataDir,'*.mat'));
dataList = {dataList(:).name}';

vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileNames = fileNames(1:end-3);
fileCt = numel(fileNames);
exceldir = 'D:\rai_tracking\Rai2Ryan_all';
excelname = 'Stochastic_MDN-1.xlsx';
exceltable = readtable(fullfile(exceldir,excelname),'Sheet','MDN_manual','ReadRowNames',true);
exceltableB = readtable(fullfile(exceldir,excelname),'Sheet','Discard trials','ReadRowNames',true);
%%
showIndicatorPlots = 0;
showThetaScroll = 0;
sumPlot = 1;

exceltable.pattern = strtrim(exceltable.pattern);
pattern2flip = ones(numel(exceltable.pattern),1);
pattern2flip(strcmp(exceltable.pattern,'L')) = -1;
exceltable.pattern(strcmp(exceltable.pattern,'L')) = {'U'};
exceltable.pattern(strcmp(exceltable.pattern,'R')) = {'U'};
grpOps = unique(exceltable.pattern);
exptList = exceltable.Properties.RowNames;
exptListB = exceltable.Properties.RowNames;

dataCell = cell(3,3);
colrs = {'r','g','b'};
tally = [];
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
    
    onDur = abs(min(lightsOn-lightsOff));
    onOffset = 10;
    onDur = 20;
    
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
        tally = [tally;mean(flyRectY)];
        if mean(flyRectY) < 500
            exptName = [vidRef 'A'];
        else
            exptName = [vidRef 'B'];
        end
        exptRef = find(strcmp(exptList,exptName));
        exptRefB = find(strcmp(exptListB,exptName));
        if isempty(exptRef)
            continue
        end
        grpRef = find(strcmp(grpOps,strtrim(exceltable.pattern{exptName})));
        
        posArray = metaTable.filtered_pos{flyID};
        xPos = posArray(:,1);
        yPos = posArray(:,2);
        thetaFilt = posArray(:,3);
        thetaRaw = metaTable.bot_theta{flyID};
        
        hPosInit = metaTable.headPos(iterF,:);
        tPosInit = metaTable.tailPos(iterF,:);
        dists = sqrt(sum([(hPosInit(:,1)-tPosInit(:,1)).^2,(hPosInit(:,2)-tPosInit(:,2)).^2],2));
        avgL = mean(dists);
        pixL = 0.11;
        %assuming 2.5mm avg fly length
        %     pix2mm = @(x) x*(2.5/avgL);
        pix2mm = @(x) x*pixL;
        %50 fps recording rate
        frm2sec = @(x) x*50;
        origFilterWin = 9;
        part2cut = (origFilterWin-1)/2+1;
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
        sideways = (relYvel)*pattern2flip(exptRef);
        turning = (dataVel(:,3)/(pi/180))*pattern2flip(exptRef);
        forwardBack(abs(forwardBack) > 2) = 0;
        sideways(abs(sideways) > 2) = 0;
        turning(abs(turning) > 15) = 0;
        instanceList = table2cell(exceltable(exptName,:));
        instanceListB = table2cell(exceltableB(exptName,:));
        sum(~strcmp(strtrim(instanceListB(1:numel(lightsOn))),''))
        for iterL = 1:numel(lightsOn)
            if iterL > numel(lightsOff)
                continue
            end
            if ~strcmp(strtrim(instanceListB{iterL}),'')
%                 instanceListB{iterL}
%                 disp(exptName)
                continue
            end
            if strcmp(strtrim(instanceList{iterL}),'NA')
                continue
            end
            
            dataBegin = lightsOn(iterL)+onOffset;
            dataEnd = dataBegin+onDur;
            if dataEnd > numel(forwardBack)
                continue
            end
            if sumPlot
                dataCell{grpRef,1} = [dataCell{grpRef,1};sum(forwardBack(dataBegin:dataEnd))];
                dataCell{grpRef,2} = [dataCell{grpRef,2};sum(sideways(dataBegin:dataEnd))];
                dataCell{grpRef,3} = [dataCell{grpRef,3};sum(turning(dataBegin:dataEnd))];
            else
%                 dataEnd = round((dataEnd-dataBegin)/2)+dataEnd;
%                 if dataEnd > numel(forwardBack)
%                     continue
%                 end
                dataCell{grpRef,1} = [dataCell{grpRef,1};(forwardBack(dataBegin:dataEnd))'];
                dataCell{grpRef,2} = [dataCell{grpRef,2};(sideways(dataBegin:dataEnd))'];
                dataCell{grpRef,3} = [dataCell{grpRef,3};(turning(dataBegin:dataEnd))'];
            end
            %             subNames = metaTable.Properties.RowNames;
            %             exceltable.distance(subNames) = round(forwardBack*10)/10;
            %             exceltable.turning(subNames) = round(turning);
        end
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

dataLabels = {'Bilateral','Empty','Unilateral'};
% dataLabels = grpOps;
%%
if sumPlot
    datarefs = [3 1 2];
    ytickval = dataLabels;
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
        conv = false;
        iterations = 2;
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
            fitY = (fitY-min(fitY(:)))/(range(fitY)*1.1);
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
    
    for iterL = 1:numel(ytickval)
        text(-400,ytickpos(iterL),(ytickval{datarefs(iterL)}),'rotation',0,...
            'color',fontC,'horizontalalignment','right','fontsize',20);
    end
    pdfName = 'MDN_unilateral_activation_turning.pdf';
    text(-300,3.5,pdfName(1:end-4),...
        'horizontalalignment','left','interpreter','none',...
        'rotation',0,'color',fontC,'fontsize',20);
    %
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile(exceldir,pdfName);
    export_fig(writePath)
    writePath = fullfile(exceldir,[pdfName(1:end-4) '.eps']);
    export_fig(writePath)
    
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
        iterations = 2;
        options = statset('MaxIter',500);
        ytickpos = [-10 0 10];
        fitX = linspace(ytickpos(1),ytickpos(end),1000);
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
            fitY = (fitY-min(fitY(:)))/(range(fitY)*1.1);
            plot(fitY+iterL-0.5,fitX,'color',[0 0 0],'LineWidth',3);
            %         plot(xdata,ydata,'color','k','markersize',18)
            clearvars obj
        end
    end
    set(gcf,'units','normalized')
    set(gcf,'color',backC)
    set(gcf,'position',[0.3932    0.4433    0.2521    0.4708])
    
    xtickpos = [1 2 3];
    fontC = [0 0 0];
    set(gca,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
        'ticklength',[0.02 0.01],'tickdir','out','xtick',xtickpos,'YLim',[ytickpos(1) ytickpos(end)],...
        'yticklabel',[],'ytick',ytickpos,'ycolor','k','xcolor','w')
    xtickval = dataLabels;
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
    
    
    pdfName = 'MDN_unilateral_activation_backing.pdf';
    text(-0.75,12,pdfName(1:end-4),...
        'horizontalalignment','left','interpreter','none',...
        'rotation',0,'color',fontC,'fontsize',20);
    %
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile(exceldir,pdfName);
    export_fig(writePath)
    writePath = fullfile(exceldir,[pdfName(1:end-4) '.eps']);
    export_fig(writePath)
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
        plotYb = (plotYb)/rangeB;
        plotYmB = mean(plotYb,1);
        plot(plotXm,plotYmB+iterD,'color',colr,'linewidth',5)
        err = std(plotYb);
        %     plotYb = [err;-err]+[plotYmB;plotYmB];
        plotYb = [err;-err;NaN(size(err))]+[plotYmB;plotYmB;plotYmB];
        plotX = repmat((1:size(plotYb,2)),size(plotYb,1),1)+part2cut+xshift;
        %     plot(plotX',plotYb','color',colr,'linewidth',2)
        plot(plotX(:),plotYb(:)+iterD,'color',colr,'linewidth',1)
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
grphRef = 1;
for iterD = 1:3
    if grphRef == 3
        yOff = iterD*30;
    else
        yOff = iterD*2;
    end
    ydata = dataCell{iterD,grphRef}'+yOff;
    xdata = repmat(1:size(ydata,1),size(ydata,2),1)';
    plot(xdata,ydata,'color',colrOps{iterD})
    hold on
    plot([1 size(ydata,1)],[yOff yOff],'k')
end

%%

grphRef = 1;
for iterD = 1:3
    if grphRef == 3
        yOff = iterD*30;
    else
        yOff = iterD*2;
    end
    ydata = median(dataCell{iterD,grphRef});
    xdata = iterD;
    plot(xdata,ydata,'.','color',colrOps{iterD},'markersize',22)
    hold all
end
