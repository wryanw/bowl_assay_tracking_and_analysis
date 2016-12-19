function [graphDataStruct,graphOptions] = getAndPrepTable_ming(graphOptions)
persistent stimStrOld grpListOld plotNameOld graphDataStructOld graphOptionsOld
%%
if isempty(mfilename) || nargin == 0
    stimStrOld = [];
    grpListOld = [];
    plotNameOld = [];
    graphDataStructOld = [];
    graphOptionsOld = [];
    graphOptions = struct('stimStr','all','plotName','escapeFrequency','useManual',false);
end
stimStr = graphOptions.stimStr;
plotName = graphOptions.plotName;
useManual = graphOptions.useManual;
if ~isempty(stimStrOld)
    if isequal(stimStr,stimStrOld) && isequal(plotName,plotNameOld)
        graphDataStruct = graphDataStructOld;
        graphOptions = graphOptionsOld;
        return
    end
end
graphTable = loadDataTable_ming;

grpList = unique(graphTable.exptInfo_stim_type);
grpList = grpList(1);
graphOptions(1).grpList = grpList;

sepPathUse = false;
dirTest = {'leftVSright','forwardVSback','towardVSaway'};
dirChoice = 3;

% if ~isempty(strfind(grpList{1},'chrimson'))
    stimTypeString = 'Photoactivation';
% else
%     stimTypeString = 'Visual_stimulation';%'Visual_stimulation' or 'Photoactivation'
% end

grpCt = numel(grpList);
%     stimTypeLogical = strcmp(graphTable.stim_type,stimTypeString);
stimTypeLogical = true(size(graphTable,1),1);

plotTypeLogical = true(size(graphTable,1),1);
if ~strcmp('Photoactivation',stimTypeString);
    if ~strcmp(stimStr,'all')
        plotTypeLogical(graphTable.exptInfo_elevation ~= 45) = false;
    end
    if strcmp(stimStr,'varied_startstop_size')
        plotTypeLogical(graphTable.exptInfo_elevation ~= 45) = false;
        plotTypeLogical(graphTable.exptInfo_azimuth ~= 0) = false;
    end
    if ~isempty(strfind(stimStr,'lv_tuning_curve'))
        if strcmp(stimStr,'lv_tuning_curve_90')
            plotTypeLogical(graphTable.exptInfo_azimuth ~= 90) = false;
        else
            plotTypeLogical(graphTable.exptInfo_azimuth ~= 0) = false;
        end
        plotTypeLogical(arrayfun(@(x) ~max(x == [5 10 15 20 30 40 80 160]),graphTable.exptInfo_stim_lv)) = false;
    end
    if ~isempty(strfind(stimStr,'lv_tuning_curve')) || ~isempty(strcmp(stimStr,'azimuth_sweep'))
        plotTypeLogical(graphTable.exptInfo_stim_start_degs ~= 10) = false;
        plotTypeLogical(graphTable.exptInfo_stim_stop_degs ~= 180) = false;
    end
    if ~isempty(strfind(stimStr,'azimuth_sweep'))
        plotTypeLogical(arrayfun(@(x) ~max(x == [0 45 90 135 180]),graphTable.exptInfo_azimuth)) = false;
    end
    if ~isempty(strfind(stimStr,'azimuth_sweep')) || strcmp(stimStr,'varied_startstop_size')
        if strcmp(stimStr,'azimuth_sweep_20')
            plotTypeLogical(graphTable.exptInfo_stim_lv ~= 20) = false;
        else
            plotTypeLogical(graphTable.exptInfo_stim_lv ~= 40) = false;
        end
    end
end
for grpNdx = 1:grpCt
    grpRef = grpList{grpNdx};
    groupLogical = strcmp(graphTable.exptInfo_stim_type,grpRef);
    %     groupLogical = true(size(graphTable,1),1);
    
    groupTable = graphTable(groupLogical & stimTypeLogical & plotTypeLogical,:);
    uniqIDsGrp = groupTable.labelB;
    uniqIDlistA = unique(uniqIDsGrp);
    if grpNdx > 1
        uniqNdxA = ismember(uniqIDlistA,uniqIDlistB);
    else
        uniqNdxA = true(size(uniqIDlistA));
    end
    uniqIDlistB = uniqIDlistA(uniqNdxA);
end
uniqIDlist = uniqIDlistB;
if isempty(uniqIDlist)
    disp('no shared uniqIDs')
    return
end

dataCt = size(uniqIDlist,1);

rearrangerPath = 'Z:\CARD LAB SHARE\Ming_temp\LC_activation_with_controls_rearranger.mat';
if exist(rearrangerPath,'file')
    rearranger = load(rearrangerPath);
    rearranger = rearranger.rearranger;
else
    rearranger = (1:dataCt);
end
graphTable.Properties.VariableNames'
%%%%% remove omits
% rearranger(dataOmitLogical) = [];
uniqIDlist = uniqIDlist(rearranger);


jumpChoice = 1;%defaults to 'jumpers'
if ~isempty(strfind(plotName,'nonjumper'))
    jumpChoice = 2;
    sepPathUse = false;
elseif ~isempty(strfind(plotName,'all'))
    jumpChoice = 3;
    sepPathUse = false;
end
jumpOps = {'jumping','nonjumping','all'};
grpCt = numel(grpList);
graphDataStruct = struct;
deg2rad = @(x) x*(pi/180);
rad2deg = @(x) x./(pi/180);
for grpNdx = 1:grpCt
    groupLogical = strcmp(graphTable.exptInfo_stim_type,grpList{grpNdx});
    %     groupLogical = true(size(graphTable,1),1);
    %     stimTypeLogical = strcmp(graphTable.stim_type,stimTypeString);
    stimTypeLogical = true(size(graphTable,1),1);
    groupTable = graphTable(groupLogical & stimTypeLogical,:);
    uniqIDsGrp = groupTable.labelB;
    vidNdcsCell = cellfun(@(x) find(strcmp(uniqIDsGrp,x)),uniqIDlist,'uniformoutput',false);
    graphDataStruct(grpNdx).groupID = grpList(grpNdx);
    for iterE = 1:dataCt
        ndx2plot = vidNdcsCell{iterE};
        exptTable = groupTable(ndx2plot,:);
        if ~strcmp('Photoactivation',stimTypeString);
            %%%%% Was the stimulus within acceptable boundaries??
            stimAzi = exptTable.exptInfo_azimuth(1);
            zeroFlyDept_degs = exptTable.departure_zero_fly;
            zeroStimDept_degs = exptTable.departure_zero_stim;
            zeroFlyStim_degs = zeroFlyDept_degs-zeroStimDept_degs-stimAzi;
            
            while sum(zeroFlyStim_degs < -179) > 0
                zeroFlyStim_degs(zeroFlyStim_degs < -179) = zeroFlyStim_degs(zeroFlyStim_degs < -179)+360;
            end
            while sum(zeroFlyStim_degs > 180)
                zeroFlyStim_degs(zeroFlyStim_degs > 180) = zeroFlyStim_degs(zeroFlyStim_degs > 180)-360;
            end
            zeroFlyStim_degs = zeroFlyStim_degs+90;
            if ~isempty(strfind(plotName,'escapeDirection'))
                stimTol = 30;
            else
                stimTol = 45;
            end
            stimBounds = [-stimTol stimTol]+90;
            stimTest = min([zeroFlyStim_degs > stimBounds(1),zeroFlyStim_degs < stimBounds(2)],[],2);
        else
            stimTest = true(size(exptTable,1),1);%%% removes stim test
        end
        %         plot(zeroFlyStim_degs,'.')
%         disp(exptTable.visStimRef{1})
%         uiwait(gcf)
        if useManual
            fot = exptTable.manFot;
            fowm = exptTable.manFowm;
            shortLogical = (fot-fowm) <= 41;
        else
            fot = exptTable.autoFot;
            fowm = NaN(size(fot));
            shortLogical = NaN(size(fot));
        end
        stimStart = exptTable.stim_start;
        
        stimDur = exptTable.stim_dur;
        graphDataStruct(grpNdx).stimDur(iterE) = stimDur(1);
        if useManual
            jumpTest = exptTable.manJumpTest;
        else
            jumpTest = exptTable.jumpTest;
        end
        notCounting = ~stimTest;
        notCounting(isnan(jumpTest)) = true;
        jumpTest(isnan(jumpTest)) = 0;
        jumpTest = logical(jumpTest);
        nonJumpers = ~jumpTest;
        nonJumpers(notCounting) = false;
        jumpTest(notCounting) = false;
        if strcmp('Photoactivation',stimTypeString);
            threshTime = 200;%in milliseconds
            lateJumpers = fot > round(stimStart+threshTime*6);
            lateJumpers(isnan(lateJumpers)) = false;
            jumpTest(lateJumpers) = false;
            nonJumpers(lateJumpers) = true;
        end
        earlyJumpers = fot < round(stimStart);
        earlyJumpers(isnan(earlyJumpers)) = false;
        jumpTest(earlyJumpers) = false;
        graphDataStruct(grpNdx).jumpCount(iterE) = sum(jumpTest);
        graphDataStruct(grpNdx).dataCount(iterE) = sum(jumpTest)+sum(nonJumpers);
        graphDataStruct(grpNdx).shortUseTotal(iterE) = sum(shortLogical(jumpTest));
%         sum(notCounting)+sum(jumpTest)+sum(nonJumpers)
%         numel(notCounting)
%         if numel(notCounting) ~= sum(notCounting)+sum(jumpTest)+sum(nonJumpers)
%             error('counts problem')
%         end
        grpVarA = exptTable.labelB{1}
        if ~strcmp(grpVarA,uniqIDlist{iterE})
            error('ID mismatch')
        end
        graphDataStruct(grpNdx).dataID(iterE) = {grpVarA};
        graphDataStruct(grpNdx).altID(iterE) = exptTable.labelA(1);
        
        %%%%% scatter plot stuff
        if strcmp(jumpOps{jumpChoice},'jumping')
            data2keep = jumpTest;
        elseif strcmp(jumpOps{jumpChoice},'nonjumping')
            data2keep = nonJumpers;
        else
            data2keep = jumpTest | nonJumpers;
        end
        graphDataStruct(grpNdx).shortUseLogical(iterE) = {shortLogical(data2keep)};
        if ~isempty(strfind(plotName,'escapeDirection'))
            %%%%%%%%%%% Directional stuff
            %             lateNufAccelTest = exptTable.accelMaxFrm > stimStart;
            %             minAccelTest = exptTable.accelMax > 0;
            %             minDistTest = exptTable.departure_XYdist > 0.0;%between 0 and 1 should do
            %             minDistTest = true(size(minDistTest));
            if strcmp(dirTest{dirChoice},'forwardVSback')
                dirResult = zeroFlyDept_degs;% zeros relative to the fly
                dirResult = dirResult+180;%%%%%% Use this for back/forward
            elseif strcmp(dirTest{dirChoice},'leftVSright')
                dirResult = zeroFlyDept_degs;% zeros relative to the fly
                dirResult = dirResult+90;%%%%%% Use this for left/right
            elseif strcmp(dirTest{dirChoice},'towardVSaway')
                dirResult = zeroStimDept_degs;
                dirResult = dirResult+180;%%%%%% Use this for toward/away from stimulus (should zero relative to stimulus)
            end
            while sum(dirResult > 179) > 0
                dirResult(dirResult > 179) = dirResult(dirResult > 179)-360;
            end
            while sum(dirResult < -180) > 0
                dirResult(dirResult < -180) = dirResult(dirResult < -180)+360;
            end
            returnData = dirResult(data2keep);
        elseif ~isempty(strfind(plotName,'frameOfTakeoff'))
            normalFot = fot-stimStart;
            if strcmp(plotName,'frameOfTakeoff_frm')
                if ~strcmp('Photoactivation',stimTypeString);
                    normTo180fot = (normalFot-stimDur)/6;%hard coded for 6000 fps record rate!!!!!!!!!!!
                else
                    normTo180fot = normalFot/6;
                end
            else
                ellovervee = exptTable.exptInfo_stim_lv(1);
                graphDataStruct(grpNdx).ellovervee(iterE) = ellovervee;
                initStimSize = exptTable.exptInfo_stim_start_degs(1);
                graphDataStruct(grpNdx).initStimSize(iterE) = initStimSize;
                finalStimSize = exptTable.exptInfo_stim_stop_degs(1);
                graphDataStruct(grpNdx).finalStimSize(iterE) = finalStimSize;
                if ~strcmp(plotName,'frameOfTakeoff_frm')
                    finalStimSize = 360;
                end
                minTheta = deg2rad(initStimSize);
                maxTheta = deg2rad(finalStimSize);
                stimStartTime = ellovervee/tan(minTheta/2);
                stimEndTime = ellovervee/tan(maxTheta/2);
                fotTime = stimStartTime-normalFot/6;
                fotTheta = 2.*atan(ellovervee./fotTime);
                normTo180fot = rad2deg(fotTheta);
                normTo180fot(fotTime < 0) = normTo180fot(fotTime < 0)+360;
                normTo180fot(fotTime < stimEndTime) = normTo180fot(fotTime < stimEndTime)+360;
            end
            
            returnData = normTo180fot(jumpTest);
        elseif ~isempty(strfind(plotName,'minDeltaZ'))
            returnData = exptTable.zero_minZ_mm(data2keep,1);
        elseif ~isempty(strfind(plotName,'maxDeltaZ'))
            returnData = exptTable.zero_maxZ_mm(data2keep,1);
        elseif ~isempty(strfind(plotName,'deltaHeading'))
            returnData = exptTable.deltaHeadingNet(data2keep,1);
            %             returnData{iterE,grpNdx} = exptTable.deltaHeadingCumSum(scatterTest,1);
        elseif ~isempty(strfind(plotName,'pitch'))
            returnData = exptTable.pitch75frm(data2keep,1);
        elseif ~isempty(strfind(plotName,'accelMax'))
            returnData = exptTable.accelMax(data2keep);
        elseif ~isempty(strfind(plotName,'accelFrm'))
            returnData = (exptTable.accelMaxFrm(jumpTest)-fot(jumpTest))/6;
        elseif ~isempty(strfind(plotName,'velMax'))
            returnData = exptTable.velMax(jumpTest);
        elseif ~isempty(strfind(plotName,'velFrm'))
            returnData = (exptTable.velMaxFrm(jumpTest)-fot(jumpTest))/6;
        elseif strcmp(plotName,'wing2fot')
            sepPathUse = false;
            histVals = (fot(jumpTest)-fowm(jumpTest))/6;
            histVals = log(histVals./6);
            %             x_start = 0;    x_end = 7;  x_step = .075;
            %             histVals = hist(histVals,(x_start:x_step:x_end));
            returnData = histVals+rand(size(histVals));
        elseif ~isempty(strfind(plotName,'spaghettiPlot'))
            if ~isempty(strfind(plotName,'zeroFly'))
                dirResult = zeros(1,sum(data2keep));% zeros relative to the fly
            elseif ~isempty(strfind(plotName,'zeroStim'))
                dirResult = zeroFlyStim_degs(data2keep)-90;% zeros relative to the stimulus
            end
            dirResult = deg2rad(dirResult);
            posArray = exptTable.zeroFlyTheta_XY_mm(data2keep);
            for iterReor = 1:numel(posArray)
                posXY = double(posArray{iterReor}(:,1:2))/1000;
                [posTheta,posRho] = cart2pol(posXY(:,1),posXY(:,2));
                posTheta = posTheta-dirResult(iterReor);
                [posX,posY] = pol2cart(posTheta,posRho);
                posArray{iterReor} = cat(2,[posX;NaN],[posY;NaN]);
            end
            returnData = cat(1,posArray{:});
        elseif ~isempty(strfind(plotName,'wing2fot_histogram'))
            histVals = (fot(jumpTest)-fowm(jumpTest));
            histVals = log(histVals./6);
            returnData = histVals;
        elseif ~isempty(strfind(plotName,'acceleration_overlay')) || ~isempty(strfind(plotName,'speed_overlay'))
            if ~isempty(strfind(plotName,'speed_overlay'))
                vecRef = 1;     scaleFac = 10;
            else
                vecRef = 2;     scaleFac = 100;
            end
            padVal = 60;
            data2keep(fot <= padVal) = false;
            mvmntArray = exptTable.speed_accel_mm(data2keep);
            subfot = fot(data2keep);
            fotNorm = true;
            for iterPad = 1:numel(mvmntArray)
                if fotNorm
                    mvmntVec = [double(mvmntArray{iterPad}(:,vecRef));NaN(padVal,1)];
                    mvmntVec = mvmntVec(subfot(iterPad)-padVal:subfot(iterPad)+padVal);
                else
                    mvmntVec = double(mvmntArray{iterPad}(:,vecRef));
                    mvmntVec = [mvmntVec;NaN(xmax-numel(mvmntVec)+1,1)];
                end
                mvmntArray{iterPad} = mvmntVec/scaleFac;
            end
            returnData = cat(2,mvmntArray{:});
        else % has info for montage
            returnArray = [fowm fot stimStart];
            returnData = returnArray(data2keep,:);
        end
        graphDataStruct(grpNdx).returnData(iterE) = {returnData};
        vidList = exptTable.Properties.RowNames;
        graphDataStruct(grpNdx).videoList(iterE) = {vidList(data2keep)};
    end
end
%%
rearrangerPath = 'Y:\Ming_RubinLab\LC_activation_with_controls_rearranger.mat';
if ~exist(rearrangerPath,'file')
    totalArray = cat(1,graphDataStruct(:).dataCount)';
    jumpArray = cat(1,graphDataStruct(:).jumpCount)';
    if grpCt > 1
        totalArray = sum(totalArray');
        jumpArray = sum(jumpArray');
    end
    jumpPct = jumpArray./totalArray;
    [~,rearranger] = sort(jumpPct);
    save(rearrangerPath,'rearranger')
    [graphDataStruct,graphOptions] = getAndPrepTable_ming(graphOptions);
    graphDataStructOld = graphDataStruct;
    graphOptions.sepPathUse = sepPathUse;
    graphOptions.jumpChoice = jumpChoice;
    graphOptions.dataCt = dataCt;
    graphOptionsOld = graphOptions;
    stimStrOld = stimStr;
    grpListOld = grpList;
    plotNameOld = plotName;
else
    delete(rearrangerPath)
end