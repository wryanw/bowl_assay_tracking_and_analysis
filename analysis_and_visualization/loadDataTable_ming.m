function outputTable = loadDataTable_ming
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
persistent graphTable
%%
if isempty(mfilename)
    graphTable = [];
end
if isempty(graphTable)
    graphTable = load('Y:\Ming_RubinLab\LC_activation_with_controls_graphTable.mat');
    graphTable = graphTable.graphTable;
    
    %%%%% computer and directory variables and information
    op_sys = system_dependent('getos');
    if strfind(op_sys,'Microsoft Windows 7')
        archDir = [filesep filesep 'tier2' filesep 'card'];
        dm11Dir = [filesep filesep 'dm11' filesep 'cardlab'];
    else
        archDir = [filesep 'Volumes' filesep 'card'];
        if ~exist(archDir,'file')
            archDir = [filesep 'Volumes' filesep 'card-1'];
        end
        dm11Dir = [filesep 'Volumes' filesep 'cardlab'];
    end
    if ~exist(archDir,'file')
        error('Archive access failure')
    end
    if ~exist(dm11Dir,'file')
        error('dm11 access failure')
    end
    
    [~,localUserName] = dos('echo %USERNAME%');
    localUserName = localUserName(1:end-1);
    repositoryName = 'pezAnalysisRepository';
    repositoryDir = fullfile('C:','Users',localUserName,'Documents',repositoryName);
    subfun_dir = fullfile(repositoryDir,'pezProc_subfunctions');
    saved_var_dir = fullfile(repositoryDir,'pezProc_saved_variables');
    addpath(subfun_dir,saved_var_dir)
    addpath(fullfile(repositoryDir,'Pez3000_Gui_folder','Matlab_functions','Support_Programs'))
    
    
    
    exptFind = strcmp(graphTable.exptInfo_stim_type,'pulse_General_widthBegin2000_widthEnd2000_cycles1_intensity20');
    stimRep = {'pulse_Williamsonw_widthBegin1000_widthEnd1000_cycles1_intensity20'};
    graphTable.exptInfo_stim_type(exptFind) = repmat(stimRep,sum(exptFind),1);
    exptFind = strcmp(graphTable.exptInfo_stim_type,'pulse_Williamsonw_widthBegin5_widthEnd75_cycles5_intensity30');
    stimRep = {'pulse_Williamsonw_widthBegin1000_widthEnd1000_cycles1_intensity20'};
    graphTable.exptInfo_stim_type(exptFind) = repmat(stimRep,sum(exptFind),1);
    
    graphTable.visStimRef = graphTable.exptInfo_stim_type;
    
    excelDir = 'Y:\WryanW\flyPez_experiment_management';
    excelName = 'bulk_testing_geno.xlsx';
    excelTable = readtable(fullfile(excelDir,excelName),'Sheet','figure_summary_group');
    
    includeVec = logical(excelTable.include);
    colVec = excelTable.Collection_ID(includeVec);
    genoVec = excelTable.Genotype_ID(includeVec);
    protoVec = excelTable.Protocol_ID(includeVec);
    groupNameVec = excelTable.group_name(includeVec);
    lineLabelVec = excelTable.LineLabel(includeVec);
    bulkCt = min([numel(colVec) numel(genoVec) numel(protoVec)]);
    exptIDlist = arrayfun(@(x) sprintf('%016s',num2str(x)),graphTable.exptNum,'uniformoutput',false);
    for iterB = 1:bulkCt
        colStr = sprintf('%04s',num2str(colVec(iterB)));
        genoStr = sprintf('%08s',num2str(genoVec(iterB)));
        protoStr = sprintf('%04s',num2str(protoVec(iterB)));
        exptIDref = [colStr,genoStr,protoStr];
        exptFind = strcmp(exptIDlist,exptIDref);
        graphTable.labelA(exptFind) = repmat(lineLabelVec(iterB),sum(exptFind),1);
        graphTable.labelB(exptFind) = repmat(groupNameVec(iterB),sum(exptFind),1);
    end
    
end
outputTable = graphTable;
