function flyTracker_mingBowl_multiFly(dataDir,showMovie,runMode,skipPrev)
%%
close all

% 1 - initialize fly positions
% 2 - generate indicator vector
% 3 - track flies
% 4 - run after all flies have been tracked
% 5 - diagnostic mode
isInverted = 0;
% 0 - white fly on black background
% 1 - black fly on white background
isFiltered = 0;
makeTemplate4tracking = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%
altStart = 1;
funDir = fileparts(mfilename('fullpath'));
if nargin == 0 || isempty(mfilename)
    funDir = 'Y:\Ming_RubinLab\tracking_code';
    %%%%%%%%%% Display control
    showMovie = 0;
    runMode = 4;
    skipPrev = 0;
    
%     dm11Dir = [filesep filesep 'dm11' filesep 'cardlab'];
    % dataDir = fullfile(dm11Dir,'CARD LAB SHARE','Ming_temp','pBDPGAL4U_CsChrimson');
    % dataDir = fullfile(dm11Dir,'CARD LAB SHARE','Ming_temp','LC16_CsChrimson');
    % dataDir = fullfile(dm11Dir,'CARD LAB SHARE','Ming_temp','unilateral_LC16','all_examples_crop');
    % dataDir = 'Y:\WryanW\Rai2Ryan_all';
%     dataDir = 'D:\rai_tracking\Rai2Ryan_all';
%     dataDir = 'Y:\Ming_RubinLab\LC10_quantification_newAvi';
    dataDir = 'Y:\Ming_RubinLab\newLC10';
    
%     dataDir = 'Y:\Ming_RubinLab\CsChrimson_circular_arena_assay_new\Figure4';
%     dataDir = 'Y:\Ming_RubinLab\CsChrimson_circular_arena_assay_new\Figure5';
%     dataDir = 'Y:\Ming_RubinLab\CsChrimson_circular_arena_assay_new\Figure6\backup_LC16_new';
cd(funDir)
end

if runMode == 1
    initializePositions(dataDir)
    return
end
if runMode == 2
    generateStimVector(dataDir,skipPrev)
    return
end

vidList = dir(fullfile(dataDir,'*.avi'));
vidList = {vidList(:).name};
[~,fileNames] = cellfun(@(x) fileparts(x),vidList,'uniformoutput',false);
fileCt = numel(fileNames);
resizeFac = 1;
if makeTemplate4tracking
    templCell = cell(0,1);
end

if ~showMovie && runMode == 3
    parforArg = 4;
else
    parforArg = 0;
end
% for fileRef = 1:fileCt
parfor (fileRef = 1:fileCt, parforArg)
    videoID = fileNames{fileRef};
%     if ~max(cellfun(@(x) ~isempty(strfind(videoID,x)),{'_v7','_v16','_v31'}))
%         continue
%     end
    disp(videoID)
    metaTablePath = fullfile(dataDir,[videoID 'metaInformation.mat']);
    backgroundPath = fullfile(dataDir,[videoID 'backgroundFrame.tif']);
    bkgfrm = imread(backgroundPath);
    bkgfrm = bkgfrm(:,:,1);
    metaTable = load(metaTablePath);
    metaTable = metaTable.metaTable;
    
    rowNames = metaTable.Properties.RowNames;
    if isempty(metaTable.bot_theta{rowNames{1}}) && runMode > 3
        continue
    end
    if ~isempty(metaTable.bot_theta{rowNames{1}}) && runMode == 3 % skips previously done
        if skipPrev
            continue
        end
    end
    
    flyCt = numel(rowNames);
    for iterF = 1:flyCt
        flyID = rowNames{iterF};
        %%%%% Define tracking variables
        hPos = metaTable.headPos(flyID,:)*resizeFac;
        tPos = metaTable.tailPos(flyID,:)*resizeFac;
        if sum(hPos-tPos) == 0
            metaTable(flyID,:) = [];
        end
    end
    rowNames = metaTable.Properties.RowNames;
    flyCt = numel(rowNames);
    vidPath = fullfile(dataDir,[videoID '.avi']);
    try
        vidObjA = VideoReader(vidPath);
    catch
        error(vidPath)
    end
    frmCt = round(vidObjA.Duration*vidObjA.FrameRate);
    vidWidth = vidObjA.Width;
    vidHeight = vidObjA.Height;
    if runMode == 5
        initFrm = altStart;
    else
        initFrm = 1;
    end
    if runMode == 4
        frameReference = (initFrm:10:frmCt);
    elseif runMode == 5
        frameReference = (initFrm:1:frmCt);
    else
        frameReference = (initFrm:frmCt);
    end
    video_3D = read(vidObjA,1);
    video_3D = squeeze(video_3D(:,:,1,:));
    video_3D = video_3D-bkgfrm;
    frm_B_gray = uint8(video_3D);%avg 1st 10 frames; est. background
    frm_bot = double(frm_B_gray)./255;
    if isInverted
        frm_bot = abs(frm_bot-1);
    end
    frm_bot = imresize(frm_bot,resizeFac,'bilinear');
    
    vidWidth = vidWidth*resizeFac;
    vidHeight = vidHeight*resizeFac;
    roiPos = [1 vidWidth 1 vidHeight];
    flyStruct = struct;
    for iterF = 1:flyCt
        flyID = rowNames{iterF};
        %%%%% Define tracking variables
        hPos = metaTable.headPos(flyID,:)*resizeFac;
        tPos = metaTable.tailPos(flyID,:)*resizeFac;
        fly_length = sqrt(sum((hPos-tPos).^2));
        trk_hort_ref = 1;
        fly_pos = mean(cat(1,hPos,tPos));
        fly_theta = cart2pol(hPos(1)-tPos(1),-(hPos(2)-tPos(2)));
        
        %%%%% Establishing background and re-mapper
        %TEMPL_*** is identical, centered replicates!!! (taken from older, centered frame)
        %SAMPL_*** contains posible new positions and has unknown center!!! (taken from new frame)
        src_fac = 0.6;
        src_leg = round(fly_length*(src_fac)*2);
        posLeg = [1 2 4 7 10 15];
        thetaLeg = [1 3 5];
        theta_ops = [fliplr(thetaLeg*(-1)) 0 thetaLeg];
        heading_factor = 3;
        spoke_count = 360/heading_factor;
        theta_ops = round([theta_ops spoke_count/2 spoke_count/4 -spoke_count/4]);
        pos_ops = [fliplr(posLeg*(-1)) 0 posLeg];
        [ndxr_struct] = trackingIndexer_spoked4crop_ming_version(src_leg,trk_hort_ref,pos_ops,theta_ops,spoke_count);
        
        posLeg = (1:3:40);
        theta_ops = round(linspace(-spoke_count/2,spoke_count/2,37));
        theta_ops = theta_ops(2:end);
        pos_ops = [fliplr(posLeg*(-1)) 0 posLeg];
        [ndxr_struct_highDef] = trackingIndexer_spoked4crop_ming_version(src_leg,trk_hort_ref,pos_ops,theta_ops,spoke_count);
        
        spoke_length = ndxr_struct.spoke_length;
        spoke_leg = ndxr_struct.spoke_leg;
        im_leg = ndxr_struct.im_leg;
        spoke_count = ndxr_struct.spoke_count;
        templ_ndxr_mastr = (ndxr_struct.templ_mastr);
        templ_rendxr = (ndxr_struct.templ_rendxr);
        spoke_ndxr_mastr = (ndxr_struct.spoke_mastr);
        spoke_span = spoke_leg*2+1;
        radians_per_spoke = 2*pi/spoke_count;
        
        
        templPath = fullfile(funDir,'template.mat');
        if makeTemplate4tracking
            templCell = cat(1,templCell,struct('ndxr_struct',ndxr_struct,...
                'frm_bot',frm_bot,'hPos',hPos,'tPos',tPos));
            continue
        else
            %%
            templLoading = load(templPath);
            templUniversal = templLoading.frm_bot;
            ex_pts = templLoading.whole_dests;
            whole_XData = [1 fly_length+im_leg*2];
            whole_YData = whole_XData;
            whole_dests = [im_leg,...
                median(whole_YData);
                fly_length+im_leg,...
                median(whole_YData)];
            Q_tform = cp2tform(ex_pts,whole_dests,'nonreflective similarity');
            frm_bot = imtransform(templUniversal,Q_tform,'XData',whole_XData,...
                'YData',whole_YData,'FillValues',0.01);
            layer_ref = round(pi/radians_per_spoke);
%             imshow(frm_bot)
%             hold on
%             plot(whole_dests(:,1),whole_dests(:,2))
            
            frm_pad = padarray(frm_bot,[im_leg im_leg]);
            se_bot1 = strel('disk',3);
            neg_dim = round(mean(whole_dests));
            pos_dim = round(mean(whole_dests)+im_leg*2);
            bot_crop = frm_pad(neg_dim(2):pos_dim(2),neg_dim(1):pos_dim(1));
            if isFiltered
                bot_open = imopen(bot_crop,se_bot1);
            else
                bot_open = bot_crop;
            end
            templ_ndxr_layer = templ_ndxr_mastr(spoke_ndxr_mastr(:,layer_ref),:);
            templ_singl = bot_open(templ_ndxr_layer);
            templ_vis = reshape(templ_singl,spoke_length,spoke_span);
%             imshow(templ_vis)
        end
        %%%%% Video analysis
        % establish tracking variables
        layer_ref = round(fly_theta/radians_per_spoke);
        layer_ref(layer_ref > spoke_count) = layer_ref-spoke_count;
        layer_ref(layer_ref <= 0 ) = layer_ref+spoke_count;
        trk_layer = layer_ref;
        trk_theta = trk_layer*radians_per_spoke;
        templ_avg_count = 100;% was 33
        templ_array = repmat(templ_singl,[1 templ_avg_count]);
        templ_mean = mean(templ_array,2);
        templ_mean_init = templ_mean;
        
        netDeltaTestA = sum((templ_singl(:)-mean(templ_singl(:))).^2);%ss_total
        netDeltaTestB = sum((templ_singl(:)-templ_mean_init(:)).^2);%ss_resid
        changeTest = 1-netDeltaTestB/netDeltaTestA;%sampl_ssd
        
        % preallocate variables, initialize slaves
        botLabels = zeros(frmCt,3);
        if runMode == 3 % 3
            botLabels(1,1:2) = fly_pos;
            botLabels(1,3) = trk_theta;
            mvmnt_ref = zeros(frmCt,1);
            range_trk = ones(frmCt,1);
            change_trk = zeros(frmCt,1);
            change_trk(1) = changeTest;
        else % 5, 4
            botLabels(:,1:2) = metaTable.bot_centroid{flyID};
            botLabels(:,3) = metaTable.bot_theta{flyID};
            mvmnt_ref = metaTable.mvmnt_ref{flyID};
            range_trk = metaTable.range_trk{flyID};
            change_trk = metaTable.change_trk{flyID};
            if ~max(strcmp(metaTable.Properties.VariableNames,'filtered_pos'))
                metaTable = cat(2,metaTable,cell2table(cell(flyCt,3),'VariableNames',...
                    {'filtered_pos','filtered_vel','filtered_accel'},'RowNames',rowNames));
            end
            try
                diffWin = 9;
                [Xzero,Xfirst,Xsecond] = golayDifferentiate(botLabels(:,1),diffWin);
                [Yzero,Yfirst,Ysecond] = golayDifferentiate(botLabels(:,2),diffWin);
                theta2filt = botLabels(:,3);
                theta2filt = cat(1,0,diff(unwrap(theta2filt(:))));
                theta2filt(theta2filt > pi/12) = 0;
                theta2filt = cumsum(theta2filt);
                [Tzero,Tfirst,Tsecond] = golayDifferentiate(theta2filt,diffWin);
                zeroOrderVectors = [Xzero(:),Yzero(:),Tzero(:)];
                firstOrderVectors = [Xfirst(:),Yfirst(:),Tfirst(:)];
                secondOrderVectors = [Xsecond(:),Ysecond(:),Tsecond(:)];
                metaTable.filtered_pos(flyID) = {zeroOrderVectors};
                metaTable.filtered_vel(flyID) = {firstOrderVectors};
                metaTable.filtered_accel(flyID) = {secondOrderVectors};
            catch ME
                getReport(ME)
                disp('you may need to install the signal processing toolbox')
            end
        end
        flyStruct(iterF).trk_layer = trk_layer;
        flyStruct(iterF).botLabels = botLabels;
        flyStruct(iterF).mvmnt_ref = mvmnt_ref;
        flyStruct(iterF).range_trk = range_trk;
        flyStruct(iterF).change_trk = change_trk;
        flyStruct(iterF).ndxr_struct = ndxr_struct;
        flyStruct(iterF).ndxr_struct_highDef = ndxr_struct_highDef;
        flyStruct(iterF).templ_vis = templ_vis;
        flyStruct(iterF).templ_singl = templ_singl;
        flyStruct(iterF).templ_mean_init = templ_mean_init;
        flyStruct(iterF).templ_mean = templ_mean;
        flyStruct(iterF).templ_array = templ_array;
    end
    if makeTemplate4tracking
        continue
    end
    % Establish figure for frame-getting
    if showMovie == 1
        scrn_size = get(0, 'ScreenSize');
        fig_position = [scrn_size(3)/10 scrn_size(4)/10 scrn_size(3)*0.6 scrn_size(4)*0.7];
        h_fig = figure('Position',fig_position);
        h_axes1 = axes('Parent',h_fig,'nextplot','replace');
        colormap('gray')
        axis image, box off, axis off, hold on
    end
    flyTrkMode = ones(flyCt,1);
    for frm_ref = 2:numel(frameReference)
        if runMode == 4 && showMovie ~= 1
            break
        end
        frameRead = frameReference(frm_ref);
        frameW = read(vidObjA,frameRead);
        frameW = frameW(:,:,1)-bkgfrm;
        %%%% Process frame for analysis
        frm_raw = frameW;
        frm_gray = double(frm_raw)./255;
        if isInverted
            frm_gray = abs(frm_gray-1);
        end
        frm_bot = imresize(frm_gray,resizeFac,'bilinear');
        
        for iterF = 1:flyCt
            %%
            botLabels = flyStruct(iterF).botLabels;
            mvmnt_ref = flyStruct(iterF).mvmnt_ref;
            range_trk = flyStruct(iterF).range_trk;
            change_trk = flyStruct(iterF).change_trk;
            if flyTrkMode(iterF) == 3
                ndxr_struct = flyStruct(iterF).ndxr_struct_highDef;
%                 continue
            elseif flyTrkMode(iterF) == 2
                ndxr_struct = flyStruct(iterF).ndxr_struct_highDef;
            else
                ndxr_struct = flyStruct(iterF).ndxr_struct;
            end
            im_leg = ndxr_struct.im_leg;
            frm_pad = padarray(frm_bot,[im_leg im_leg]);
            if runMode ~= 4
                flyID = rowNames{iterF};
                fly_pos = botLabels(frameReference(frm_ref-1),1:2);
                trk_layer = flyStruct(iterF).trk_layer;
                templ_singl = flyStruct(iterF).templ_singl;
                templ_mean_init = flyStruct(iterF).templ_mean_init;
                templ_mean = flyStruct(iterF).templ_mean;
                templ_array = flyStruct(iterF).templ_array;
                
                mastr_findr = ndxr_struct.mastr_findr;
                spoke_length = ndxr_struct.spoke_length;
                spoke_leg = ndxr_struct.spoke_leg;
                spoke_count = ndxr_struct.spoke_count;
                sampl_ndxr_mastr = (ndxr_struct.sampl_mastr);
                templ_rendxr = (ndxr_struct.templ_rendxr);
                
                templ_trk = templ_mean_init(templ_rendxr);%stays with initial template
                %             templ_trk = templ_mean(templ_rendxr);%uses updating template
                
                spoke_ndxr_mastr = (ndxr_struct.spoke_mastr);
                spoke_span = spoke_leg*2+1;
                radians_per_spoke = 2*pi/spoke_count;
                
                %%%% Test whether anything changed from the previous frame
                neg_dim = round(fly_pos);
                pos_dim = round(fly_pos+im_leg*2);
                bot_crop = frm_pad(neg_dim(2):pos_dim(2),neg_dim(1):pos_dim(1));
                if isFiltered
                    bot_open = imopen(bot_crop,se_bot1);
                else
                    bot_open = bot_crop;
                end
                sampl_ndxr_layer = sampl_ndxr_mastr(spoke_ndxr_mastr(:,trk_layer),:);
                sampl_trk = bot_open(sampl_ndxr_layer);
                ss_total = sum((sampl_trk-repmat(mean(sampl_trk),size(sampl_trk,1),1)).^2);
                ss_resid = sum((sampl_trk-templ_trk).^2);
                trk_ssd = 1-ss_resid./ss_total;
                [~,max_ref] = max(trk_ssd);
                mastr_refs = mastr_findr(max_ref,:);
                
                rangeTestA = sum((templ_singl(:)-mean(templ_singl(:))).^2);%ss_total
                rangeTestB = sum((templ_singl(:)-templ_mean_init(:)).^2);%ss_resid
                rangeTest = 1-rangeTestB./rangeTestA;
                
                netDeltaTestA = sum((templ_singl(:)-mean(templ_singl(:))).^2);%ss_total
                netDeltaTestB = sum((templ_singl(:)-templ_mean(:)).^2);%ss_resid
                changeTest = 1-netDeltaTestB/netDeltaTestA;%sampl_ssd
                if changeTest <= 0.5 && changeTest > -1
                    flyTrkMode(iterF) = 2;
                elseif changeTest <= -1
                    if runMode == 5
                        flyTrkMode(iterF) = 2;
                    else
                        flyTrkMode(iterF) = 3;
                    end
                else
                    flyTrkMode(iterF) = 1;
                end
                %%%% Respond to detected changes
                fly_pos = fly_pos+mastr_refs(1:2);
                templ_singl = sampl_trk(:,max_ref);
                
                trk_layer = trk_layer+mastr_refs(3);
                trk_layer(trk_layer > spoke_count) = trk_layer-spoke_count;
                trk_layer(trk_layer <= 0 ) = trk_layer+spoke_count;
                trk_theta = trk_layer*radians_per_spoke;
                templ_vis = reshape(templ_singl,spoke_length,spoke_span);
                mvmnt_ref(frameRead) = frameRead;
            else
                fly_pos = botLabels(frameRead,1:2);
                trk_theta = botLabels(frameRead,3);
                rangeTest = range_trk(frameRead);
                changeTest = change_trk(frameRead);
            end
            
            if showMovie == 1
                if runMode == 4 % 4
                    visual_frame = frm_bot;
                    changeTest = metaTable.change_trk{flyID}(frameRead);
                else % 5, 3
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
                    %                         vis_full = imresize(vis_full,2);
                    %             vis_full = vis_trks;
                    visual_frame = frm_bot;
                    max_dim = max(size(visual_frame,2),size(vis_full,2));
                    visual_frame(end,max_dim) = 0;
                    vis_full(end,max_dim) = 0;
                    %                     visual_frame(end-size(vis_full,1)+1:end,end-size(vis_full,2)+1:end) = vis_full;
                    visual_frame(1:size(vis_full,1),end-size(vis_full,2)+1:end) = vis_full;
                end
                
%                 set(hIm,'cdata',visual_frame)
                
                if iterF == 1
                    image(visual_frame,'CDataMapping','scaled','Parent',h_axes1);
                    axis image, box off, axis off, hold on
                    set(h_axes1,'nextplot','add')
                end
                
                %%%% Visualize bottom-view labels
                u_com = cos(trk_theta)*fly_length/2;
                v_com = -sin(trk_theta)*fly_length/2;
                
                %%%% Plot data on image
                quiver(fly_pos(1),fly_pos(2),u_com,v_com,'MaxHeadSize',0,'LineWidth',2,...
                    'AutoScaleFactor',1,'Color',[0 .5 .8]);
                plot(fly_pos(1),fly_pos(2),'.','MarkerSize',18,'Color',[.7 0 0]);
                text(1,25,['Frame:  ' int2str(frameRead)],...
                    'HorizontalAlignment','left','FontWeight','bold',...
                    'Color',[1 0 0],'FontSize',20,'Interpreter','none');
                text(-200,iterF*20,['Tmpl Change:  ' num2str(changeTest,3)],...
                    'HorizontalAlignment','left','FontWeight','bold',...
                    'Color',[1 0 0],'FontSize',10,'Interpreter','none');
                text(size(visual_frame,2)+5,iterF*20,['Tmpl Range:  ' num2str(rangeTest,3)],...
                    'HorizontalAlignment','left','FontWeight','bold',...
                    'Color',[1 0 0],'FontSize',10,'Interpreter','none');
                if iterF == flyCt
                    drawnow
                    set(h_axes1,'nextplot','replace')
                end
            end
            if runMode == 3
                botLabels(frameRead,:) = [fly_pos,trk_theta];
            end
            range_trk(frameRead) = rangeTest;
            change_trk(frameRead) = changeTest;
            
            flyStruct(iterF).botLabels = botLabels;
            flyStruct(iterF).mvmnt_ref = mvmnt_ref;
            flyStruct(iterF).range_trk = range_trk;
            flyStruct(iterF).change_trk = change_trk;
            flyStruct(iterF).trk_layer = trk_layer;
            flyStruct(iterF).templ_vis = templ_vis;
            flyStruct(iterF).templ_singl = templ_singl;
            flyStruct(iterF).templ_mean = templ_mean;
            flyStruct(iterF).templ_array = templ_array;
        end
    end
    if showMovie == 1
        close(h_fig)
    end
    %%
    for iterF = 1:flyCt
        flyID = rowNames{iterF};
        botLabels = flyStruct(iterF).botLabels;
        metaTable.bot_centroid(flyID) = {botLabels(:,1:2)};
        metaTable.bot_theta(flyID) = {botLabels(:,3)};
        metaTable.mvmnt_ref(flyID) = {flyStruct(iterF).mvmnt_ref};
        metaTable.range_trk(flyID) = {flyStruct(iterF).range_trk};
        metaTable.change_trk(flyID) = {flyStruct(iterF).change_trk};
        if runMode ~= 5
            parsave_metaTable(metaTablePath,metaTable)
            disp('saved')
        end
    end
end
%%
if makeTemplate4tracking
    makeTemplate(templCell,templPath)
end