clear all

dataDir = fullfile('Y:\Ming_RubinLab','LC16_CsChrimson');
vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
dataFiles = dir(fullfile(dataDir,'*filteredData.mat'));
dataFiles = {dataFiles(:).name};
dataLoading = cellfun(@(x) load(fullfile(dataDir,x)),dataFiles,'uniformoutput',false);
rowNames = cellfun(@(x) regexprep(x,'_filteredData.mat',''),dataFiles,'uniformoutput',false)';
dataCt = numel(rowNames);
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
%%

pixL = 0.11;
pix2mm = @(x) x*pixL;
%50 fps recording rate
frm2sec = @(x) x*50;
for iterV = [2]
% for iterV = [2 4 5 6 9 13 15 18]
%     for iterV = 1:numel(rowNames)
    % vidName = 'OL0046B_UAS-Chrimson_invert_substack1724-1873_fly12.avi';
    % vidName = vidList{iterV}
    videoID = rowNames{iterV};
    vidName = [videoID '.avi']
    iterV
    tracker_data_path = fullfile(dataDir,[videoID '.mat']);
    load(tracker_data_path)
    dataAccel = dataTable.filtered_accel{videoID};
    dataVel = dataTable.filtered_vel{videoID};
    dataPos = dataTable.filtered_pos{videoID};
%     dataPos = [tracker_record.bot_centroid{videoID} tracker_record.bot_theta{videoID}];
    [velT,velR] = cart2pol(dataVel(:,1),-dataVel(:,2));
    [relXvel,relYvel] = pol2cart(velT-dataPos(:,3),velR);
    relXvel = frm2sec(pix2mm(relXvel));
    relYvel = frm2sec(pix2mm(relYvel));
    forwardBack = smooth(relXvel,15)';
    turning = frm2sec(smooth(abs(dataVel(:,3)),15)');
    
    vidPath = fullfile(dataDir,vidName);
    vidObjA = VideoReader(vidPath);
    vidWidth = vidObjA.Width*2;
    vidHeight = vidObjA.Height*2;
    
    offset = 20;%frames to cut on either end
    vis_count = 9;%number of frames
    spacingMode = 2;%1 - frame based   2 - distance based
    mapMode = 1;%1 - turning   2 - speed
    
    
    deltaPos = diff(dataPos(:,1:2)-repmat(dataPos(1,1:2),size(dataPos,1),1),[],1);
    smooth_dist_filt = cumsum(sqrt(deltaPos(offset:end-offset,1).^2+deltaPos(offset:end-offset,2).^2));
    pos_marks = linspace(0,smooth_dist_filt(end),vis_count);
    pos_mark_adj = 1;
    pos_mark_array = repmat(pos_marks.*pos_mark_adj,size(smooth_dist_filt,1),1);
    delta_net_array = repmat(smooth_dist_filt,1,vis_count);
    [~,vis_refs] = min(abs(delta_net_array-pos_mark_array));
    
    if spacingMode == 1
        frmBounds = round(linspace(20,vidObjA.NumberOfFrames-20,vis_count+1));
    else
        frmBounds = vis_refs+offset;
    end
    frmRefs = round(frmBounds(1:end-1)+diff(frmBounds)/2);
    frmCt = numel(frmRefs);
    frmM = uint8(zeros(vidHeight,vidWidth,3,frmCt));
    colrVal = zeros(frmCt,1);
    for iterF = 1:frmCt
        if mapMode == 1
            colorstr = 'hot';
            % colorstr = 'bone';
            %         colorstr = 'jet';
            colrVal(iterF) = mean(turning(frmBounds(iterF:iterF+1)));
        else
            colorstr = 'hot';
            %         colorstr = 'cool';
            %         colrVal(iterF) = mean(forwardBack(frmBounds(iterF:iterF+1)));
            colrVal(iterF) = mean(abs(forwardBack(frmBounds(iterF:iterF+1))));
        end
    end
    colorfun = str2func(colorstr);
    colrVal
    remapper = (uint8(colorfun(256).*255));
    if mapMode == 1
        colrValB = colrVal/max(colrVal);
    else
        %     colrValB = colrVal/max(abs(colrVal))/2+0.5
        colrValB = colrVal/max(colrVal)
        % skewup = 0.2;
        % colrVal = colrVal*(1-skewup)+skewup;
    end
    for iterF = 1:frmCt
        frmR = read(vidObjA,frmRefs(iterF));
        frmS = double(imresize(frmR(:,:,1),2));
        %     frmS = uint8(frmS+abs(frmS-255)*colrVal(iterF));
        frmS = uint8(frmS);
        mapIm = uint8(repmat(zeros(size(frmS)),[1 1 3]));
        for iterR = 1:3
            %         mapIm(:,:,iterR) = intlut(frmS,remapper(:,iterR));
            mapIm(:,:,iterR) = frmS;
        end
        frmM(:,:,:,iterF) = mapIm;
    end
    
    frmF = min(frmM,[],4);
    frmZ = frmF;
    imshow(frmZ)
    hold on
    ptrefs = (frmRefs(1):2:frmRefs(end));
    plot(dataPos(ptrefs,1),dataPos(ptrefs,2),'.','markersize',10,'color',[0 0 0])
    ptrefs = (frmRefs(1):2:frmRefs(end));
    ptrefs(ptrefs < 50 | ptrefs > 100) = [];
    plot(dataPos(ptrefs,1),dataPos(ptrefs,2),'.','markersize',10,'color',[1 0 0])
    for iterF = 1:frmCt
        colr = double(remapper(round(colrValB(iterF)*255)+1,:))/255;
%         plot(dataPos(frmRefs(iterF),1),dataPos(frmRefs(iterF),2),'.','markersize',22,...
%             'color',colr)
        theta = dataPos(frmRefs(iterF),3);
        length = 20;
        u_deprt = cos(theta)*length; v_deprt = -sin(theta)*length;
        quiver(dataPos(frmRefs(iterF),1),dataPos(frmRefs(iterF),2),u_deprt,v_deprt,...
            'MaxHeadSize',12,'LineWidth',1,'AutoScaleFactor',1,'Color',[0 1 0])
        angles = [0 pi-pi/5 pi+pi/5 0];
        scalefacs = [1 .5 .5 1];
        length = 30;
        uA = cos(angles+theta)*length.*scalefacs; vA = -sin(angles+theta)*length.*scalefacs;
%         plot(dataPos(frmRefs(iterF),1)+uA,dataPos(frmRefs(iterF),2)+vA,'color',[0 1 0],'linewidth',2)
    end
    text(30+dataPos(frmRefs(1),1),dataPos(frmRefs(1),2),'start','fontsize',14)
    text(-60+dataPos(frmRefs(end),1),dataPos(frmRefs(end),2),'stop','fontsize',14)
    
%     colormap(colorstr)
%     hcol = colorbar();
%     uiwait(gcf)
end

%%
pdfName = 'LC16_sample_data.pdf';
writePath = fullfile('Z:\CARD LAB SHARE\Ming_temp',pdfName);
export_fig(writePath)

