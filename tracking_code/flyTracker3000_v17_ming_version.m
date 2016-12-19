function flyTracker3000_v17_ming_version
close all
%%%%%%%%%% Display control
showMovie = 1;
showSavedData = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%
op_sys = system_dependent('getos');
if strfind(op_sys,'Microsoft Windows 7')
    dm11Dir = [filesep filesep 'dm11' filesep 'cardlab'];
else
    dm11Dir = [filesep 'Volumes' filesep 'cardlab'];
end
% dataDir = fullfile(dm11Dir,'CARD LAB SHARE','Ming_temp','pBDPGAL4U_CsChrimson');
dataDir = fullfile(dm11Dir,'CARD LAB SHARE','Ming_temp','LC16_CsChrimson');
% dataDir = fullfile(dm11Dir,'CARD LAB SHARE','Ming_temp','unilateral_LC16','all_examples_crop');

resizeFac = 2;
initPosPath = fullfile(dataDir,'initial_positions.mat');
if exist(initPosPath,'file')
    initPos = load(initPosPath);
    initPos = initPos.initPos;
end
%%
vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
rowNames = cellfun(@(x) x(1:end-4),vidList,'uniformoutput',false);
fileCt = numel(rowNames);
varNames = {'headPos','tailPos'};
if ~exist('initPos','var')
    initPos = table(zeros(fileCt,2),zeros(fileCt,2),'VariableNames',varNames,...
        'RowNames',rowNames);
else
    oldNames = initPos.Properties.RowNames;
    newNames = rowNames(cellfun(@(x) ~max(strcmp(oldNames,x)),rowNames));
    fileCt = numel(newNames);
    initPos = [initPos;table(zeros(fileCt,2),zeros(fileCt,2),'VariableNames',varNames,...
        'RowNames',newNames)];
end

rowNames = initPos.Properties.RowNames;
rowNames = rowNames(sum(initPos.headPos,2) == 0);
fileCt = numel(rowNames);

for iterF = 1:fileCt
    fileRef = iterF;
    videoID = rowNames{fileRef};
    vidPath = fullfile(dataDir,[videoID '.avi']);
    try
        vidObjA = VideoReader(vidPath);
    catch
        disp(vidPath)
    end
    video_3D = read(vidObjA,1);
    video_3D = squeeze(video_3D(:,:,1,:));
    frm_B_gray = uint8(video_3D);%avg 1st 10 frames; est. background
    frm_bot = double(frm_B_gray)./255;
    frm_bot = imresize(abs(frm_bot-1),resizeFac,'bilinear');
    
    imshow(frm_bot,'initialmagnification',300)
    hrect = imrect;
    frect = round(getPosition(hrect));
    frm_crop = frm_bot(frect(2):frect(2)+frect(4),frect(1):frect(1)+frect(3));
    close all
    imshow(frm_crop,'initialmagnification',300)
    [x, y] = getpts(gcf);
    x = x+frect(1);
    y = y+frect(2);
    initPos.headPos(videoID,:) = [x(1) y(1)]/resizeFac;
    initPos.tailPos(videoID,:) = [x(2) y(2)]/resizeFac;
    close all
    %         imshow(frm_bot)
    %         hold on
    %         plot(x,y)
    %         return
    save(initPosPath,'initPos')
    disp(['progress: ' num2str(round(fileRef/fileCt*100)) ' %'])
end

load(initPosPath)
rowNames = initPos.Properties.RowNames;
fileCt = numel(rowNames);
%%
% tPosInit = zeros(fileCt,2);
% hPosInit = zeros(fileCt,2);
for fileRef = 1:fileCt
    videoID = rowNames{fileRef};
    disp(videoID)
    vidPath = fullfile(dataDir,[videoID '.avi']);
    try
        vidObjA = VideoReader(vidPath);
    catch
        error(vidPath)
    end
    frameReference = (1:vidObjA.NumberOfFrames);
    
    vidWidth = vidObjA.Width;
    vidHeight = vidObjA.Height;
    frm_count = numel(frameReference);
    
    tracker_data_path = fullfile(dataDir,[videoID '.mat']);
    if showSavedData
        filtered_data_path = fullfile(dataDir,[videoID '_filteredData.mat']);
        load(tracker_data_path)
        ana_varnames = {'filtered_pos','filtered_vel','filtered_accel'};
        analyzer_record = cell2table(cell(1,numel(ana_varnames)),'VariableNames',...
            ana_varnames,'RowNames',{videoID});
    else
        trk_varnames = {'bot_centroid','bot_theta','mvmnt_ref','range_trk','change_trk'};
        tracker_record = cell2table(cell(1,numel(trk_varnames)),'VariableNames',...
            trk_varnames,'RowNames',{videoID});
    end
    
    
    video_3D = read(vidObjA,1);
    video_3D = squeeze(video_3D(:,:,1,:));
    frm_B_gray = uint8(video_3D);%avg 1st 10 frames; est. background
    frm_bot = double(frm_B_gray)./255;
    frm_bot = imresize(abs(frm_bot-1),resizeFac,'bilinear');
    
    hPos = initPos.headPos(videoID,:)*resizeFac;
    tPos = initPos.tailPos(videoID,:)*resizeFac;
    
    vidWidth = vidWidth*resizeFac;
    vidHeight = vidHeight*resizeFac;
    roiPos = [1 vidWidth 1 vidHeight];
    
    
    
    %%%%% Define tracking variables
    trk_hort_ref = 1;
    cntr_pt_init = mean(cat(1,hPos,tPos));
    fly_theta = cart2pol(hPos(1)-tPos(1),-(hPos(2)-tPos(2)));
    fly_length = sqrt(sum((hPos-tPos).^2));
    
    % parent_dir = 'C:\Users\williamsonw\Documents\pezAnalysisRepository';
    % subfun_dir = fullfile(parent_dir,'pezProc_subfunctions');
    % addpath(subfun_dir)
    
    %%%%% Establishing background and re-mapper
    %TEMPL_*** is identical, centered replicates!!! (taken from older, centered frame)
    %SAMPL_*** contains posible new positions and has unknown center!!! (taken from new frame)
    src_fac = 0.5;
    src_leg = round(fly_length*(src_fac)*2);
    [ndxr_struct] = trackingIndexer_spoked4crop_ming_version(src_leg,trk_hort_ref);
    mastr_findr = ndxr_struct.mastr_findr;
    spoke_length = ndxr_struct.spoke_length;
    spoke_leg = ndxr_struct.spoke_leg;
    im_leg = ndxr_struct.im_leg;
    spoke_count = ndxr_struct.spoke_count;
    sampl_ndxr_mastr = (ndxr_struct.sampl_mastr);
    templ_ndxr_mastr = (ndxr_struct.templ_mastr);
    templ_rendxr = (ndxr_struct.templ_rendxr);
    spoke_ndxr_mastr = (ndxr_struct.spoke_mastr);
    spoke_span = spoke_leg*2+1;
    radians_per_spoke = 2*pi/spoke_count;
    layer_ref = round(fly_theta/radians_per_spoke);
    layer_ref(layer_ref > spoke_count) = layer_ref-spoke_count;
    layer_ref(layer_ref <= 0 ) = layer_ref+spoke_count;
    
    %% %%% Initialize bottom frame (under-side view)
    %
%         imshow(frm_bot)
%         hold on
%         plot(cntr_pt_init(1),cntr_pt_init(2),'.','MarkerSize',16,'Color',[0 .7 .7])
%         plot([hPos(1) tPos(1)],[hPos(2) tPos(2)],'.','MarkerSize',14,'Color',[.7 0 0])
%         hold off
%         uiwait(gcf,1)
    %%
    frm_bot = padarray(frm_bot,[im_leg im_leg]);
    se_bot1 = strel('disk',3);
    neg_dim = round(cntr_pt_init);
    pos_dim = round(cntr_pt_init+im_leg*2);
    bot_crop = frm_bot(neg_dim(2):pos_dim(2),neg_dim(1):pos_dim(1));
    bot_open = bot_crop;
    templ_ndxr_layer = templ_ndxr_mastr(spoke_ndxr_mastr(:,layer_ref),:);
    templ_ndxr_layer(templ_ndxr_layer > numel(bot_open)) = numel(bot_open);
    templ_ndxr_layer(templ_ndxr_layer < 1) = 1;
    templ_singl = bot_open(templ_ndxr_layer);
    templ_vis = reshape(templ_singl,spoke_length,spoke_span);
    fly_pos = cntr_pt_init;
    fly_pos_init = fly_pos;
    
%     imshow(cell2mat(templ_vis_cell),'InitialMagnification',300)
    
    %% %%% Video analysis
    % establish tracking variables
    trk_layer = layer_ref;
    trk_theta = trk_layer*radians_per_spoke;
    templ_avg_count = 33;
    templ_array = repmat(templ_singl,[1 templ_avg_count]);
    templ_mean = mean(templ_array,2);
    templ_mean_init = templ_mean;
    templ_trk = templ_mean(templ_rendxr);
    netDeltaTestA = sum((templ_singl(:)-mean(templ_singl(:))).^2);%ss_total
    netDeltaTestB = sum((templ_singl(:)-templ_mean_init(:)).^2);%ss_resid
    changeTest = 1-netDeltaTestB/netDeltaTestA;%sampl_ssd
    
    % preallocate variables, initialize slaves
    botLabels = zeros(frm_count,3);
    mvmnt_ref = zeros(frm_count,1);
    range_trk = zeros(frm_count,1);
    change_trk = zeros(frm_count,1);
    if showSavedData
        botLabels(:,1:2) = tracker_record.bot_centroid{1};
        botLabels(:,3) = tracker_record.bot_theta{1};
        mvmnt_ref = tracker_record.mvmnt_ref{1};
        range_trk = tracker_record.range_trk{1};
        change_trk = tracker_record.change_trk{1};
%         diffWin = 5;
%         diffWin = 19;
        diffWin = 9;
        [Xzero,Xfirst,Xsecond] = golayDifferentiate(botLabels(:,1),diffWin);
        [Yzero,Yfirst,Ysecond] = golayDifferentiate(botLabels(:,2),diffWin);
        [Tzero,Tfirst,Tsecond] = golayDifferentiate(unwrap(botLabels(:,3)),diffWin);
        botLabels = [Xzero(:),Yzero(:),Tzero(:)];
        filtWin = 9;
        botLabels(:,1) = smooth(botLabels(:,1),filtWin)';
        botLabels(:,2) = smooth(botLabels(:,2),filtWin)';
        botLabels(:,3) = smooth(unwrap(botLabels(:,3)),filtWin)';
    end
    % Establish figure for frame-getting
    if showMovie == 1
        scrn_size = get(0, 'ScreenSize');
        fig_position = [scrn_size(3)/10 scrn_size(4)/10 scrn_size(3)*0.6 scrn_size(4)*0.7];
        h_fig = figure('Position',fig_position);
        h_axes1 = axes('Parent',h_fig);
        colormap('gray')
    end
    for frm_ref = 1:frm_count
        
        frameRead = frameReference(frm_ref);
        frameW = read(vidObjA,frameRead);
        frameW = frameW(:,:,1);
        %%%% Process frame for analysis
        frm_raw = frameW;
        frm_gray = double(frm_raw)./255;
        frm_bot = imresize(abs(frm_gray-1),resizeFac);
        frm_bot = padarray(frm_bot,[im_leg im_leg]);
        if ~showSavedData
            %%%% Test whether anything changed from the previous frame
            neg_dim = round(fly_pos);
            pos_dim = round(fly_pos+im_leg*2);
            
            bounds_testA = max([neg_dim(1) pos_dim(1)]) > vidWidth+im_leg*2;
            bounds_testB = max([neg_dim(2) pos_dim(2)]) > vidHeight+im_leg*2;
            bounds_testC = (min([neg_dim pos_dim]) < 1);
            if logical(max([bounds_testA,bounds_testB,bounds_testC]))
                break
            end
            bot_crop = frm_bot(neg_dim(2):pos_dim(2),neg_dim(1):pos_dim(1));
            
            bot_open = imopen(bot_crop,se_bot1);
            %     bot_open = bot_crop;
            
            sampl_ndxr_layer = sampl_ndxr_mastr(spoke_ndxr_mastr(:,trk_layer),:);
            sampl_trk = bot_open(sampl_ndxr_layer);
            ss_total = sum((sampl_trk-repmat(mean(sampl_trk),numel(templ_vis),1)).^2);
            ss_resid = sum((sampl_trk-templ_trk).^2);
            trk_ssd = 1-ss_resid./ss_total;
            [~,max_ref] = max(trk_ssd);
            mastr_refs = mastr_findr(max_ref,:);
            cont_test = sum([mastr_refs frm_ref == 1]);
            
            rangeTestA = sum((sampl_trk(:,1)-mean(sampl_trk(:,1))).^2);%ss_total
            rangeTestB = sum((sampl_trk(:,1)-sampl_trk(:,end)).^2);%ss_resid
            rangeTest = 1-rangeTestB./rangeTestA;
        else
            cont_test = 1;
        end
        if  cont_test ~= 0
            if ~showSavedData
                mvmnt_ref(frm_ref) = frm_ref;
                
                rangeTestA = sum((sampl_trk(:,1)-mean(sampl_trk(:,1))).^2);%ss_total
                rangeTestB = sum((sampl_trk(:,1)-sampl_trk(:,end)).^2);%ss_resid
                rangeTest = 1-rangeTestB./rangeTestA;
                
                netDeltaTestA = sum((templ_singl(:)-mean(templ_singl(:))).^2);%ss_total
                netDeltaTestB = sum((templ_singl(:)-templ_mean_init(:)).^2);%ss_resid
                changeTest = 1-netDeltaTestB/netDeltaTestA;%sampl_ssd
                if changeTest < -1
                    break
                end
                
                
                %%%% Respond to detected changes
                netPointDelta = (fly_pos_init-fly_pos);
                netPosDelta = sqrt(netPointDelta(2)^2+netPointDelta(1)^2);
                if netPosDelta > 30
                    netTest = 1;
                else
                    netTest = 1;
                end
                if netTest
                    templ_singl = sampl_trk(:,max_ref);
                end
                fly_pos = fly_pos+mastr_refs(1:2);
                posTest = [(fly_pos(1)) < roiPos(1)
                    (fly_pos(2)) < roiPos(3)
                    (fly_pos(1)+im_leg*2) > roiPos(2)
                    (fly_pos(2)+im_leg*2) > roiPos(4)];
                if max(posTest)
                    netTest = 0;
                end
                netTest = 1;
                if netTest
                    trk_layer = trk_layer+mastr_refs(3);
                    trk_layer(trk_layer > spoke_count) = trk_layer-spoke_count;
                    trk_layer(trk_layer <= 0 ) = trk_layer+spoke_count;
                    trk_theta = trk_layer*radians_per_spoke;
                    
                    templ_vis = reshape(templ_singl,spoke_length,spoke_span);
                    %             templ_array = circshift(templ_array,[0 1]);
                    %             templ_array(:,1) = templ_singl;
                    %             templ_mean = mean(templ_array(:,end-3:end),2);
                    templ_trk = templ_mean(templ_rendxr);
                end
            else
                fly_pos = botLabels(frm_ref,1:2);
                trk_theta = botLabels(frm_ref,3);
            end
            if showMovie == 1
                if ~showSavedData
                    sampl_vis = reshape(sampl_trk(:,max_ref),spoke_length,spoke_span);
                    templ_vis_mean = reshape(templ_mean,spoke_length,spoke_span);
                    
                    %%%%% Generate new image
                    spcrA = zeros(3,size(bot_open,2));
                    vis_blocks = [spcrA;bot_open;spcrA];
                    spcrB = zeros(3,spoke_span);
                    vis_trks = [sampl_vis;spcrB;templ_vis;spcrB;templ_vis_mean;spcrB];
                    max_dim = max(size(vis_blocks,1),size(vis_trks,1));
                    vis_blocks(max_dim,end) = 0;
                    vis_trks(max_dim,end) = 0;
                    vis_full = [vis_blocks vis_trks];
                    %             vis_full = vis_trks;
                    visual_frame = frm_bot;
                    max_dim = max(size(visual_frame,2),size(vis_full,2));
                    visual_frame(end,max_dim) = 0;
                    vis_full(end,max_dim) = 0;
%                     visual_frame(end-size(vis_full,1)+1:end,end-size(vis_full,2)+1:end) = vis_full;
                    visual_frame(1:size(vis_full,1),end-size(vis_full,2)+1:end) = vis_full;
                else
                    visual_frame = frm_bot;
                    pause(0.1)
                end
                image(visual_frame,'CDataMapping','scaled','Parent',h_axes1);
                axis image, box off, axis off, hold on
                
                
                %%%% Visualize bottom-view labels
                u_com = cos(trk_theta)*fly_length/2;
                v_com = -sin(trk_theta)*fly_length/2;
                
                %%%% Plot data on image
                quiver(fly_pos(1)+im_leg,fly_pos(2)+im_leg,u_com,v_com,'MaxHeadSize',0,'LineWidth',2,...
                    'AutoScaleFactor',1,'Color',[0 .5 .8]);
                plot(fly_pos(1)+im_leg,fly_pos(2)+im_leg,'.','MarkerSize',18,'Color',[.7 0 0])
                if ~showSavedData
                    text(15,75,['Pos Change:  ' num2str(netPosDelta,3)],...
                        'HorizontalAlignment','left','FontWeight','bold',...
                        'Color',[1 1 1],'FontSize',10,'Interpreter','none')
                    text(15,105,['Tmpl Change:  ' num2str(changeTest,3)],...
                        'HorizontalAlignment','left','FontWeight','bold',...
                        'Color',[1 1 1],'FontSize',10,'Interpreter','none')
                    text(15,135,['Frame:  ' int2str(frm_ref)],...
                        'HorizontalAlignment','left','FontWeight','bold',...
                        'Color',[1 1 1],'FontSize',10,'Interpreter','none')
                end
                hold off
                drawnow
            end
            
        end
        if ~showSavedData
            botLabels(frm_ref,:) = [fly_pos,trk_theta];
            range_trk(frm_ref) = rangeTest;
            change_trk(frm_ref) = changeTest;
        end
    end
    if showMovie == 1
        close(h_fig)
    end
    %%
    if ~showSavedData
        tracker_record.bot_centroid = {botLabels(:,1:2)};
        tracker_record.bot_theta = {botLabels(:,3)};
        tracker_record.mvmnt_ref = {mvmnt_ref};
        tracker_record.range_trk = {range_trk};
        tracker_record.change_trk = {change_trk};
        save(tracker_data_path,'tracker_record')
    else
        analyzer_record.filtered_pos = {botLabels};
        analyzer_record.filtered_vel = {[Xfirst(:),Yfirst(:),Tfirst(:)]};
        analyzer_record.filtered_accel = {[Xsecond(:),Ysecond(:),Tsecond(:)]};
        save(filtered_data_path,'analyzer_record')
    end
end