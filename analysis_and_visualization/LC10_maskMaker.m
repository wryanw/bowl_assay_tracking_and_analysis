excelDir = 'Y:\Ming_RubinLab';
excelName = 'Table1_v5_WRW_MW_v2.xlsx';
sheetName = 'graphIDtable_LC10';
excelPath = fullfile(excelDir,excelName);
excelTable = readtable(excelPath,'Sheet',sheetName,'ReadRowNames',true);

vidDir = 'Y:\Ming_RubinLab\montage_vids_500ms';
vidList = dir(fullfile(vidDir,'*.mp4'));
vidList = {vidList(:).name}';
%%
codeList = excelTable.blind_code;
rowNames = excelTable.Properties.RowNames;
for iterR = 1:numel(vidList)
    oldName = vidList{iterR};
    brkpts = strfind(oldName,'_');
    oldLabel = oldName(brkpts(1)+1:brkpts(2)-1);
    labelRef = strcmp(rowNames,oldLabel);
    newName = regexprep(oldName,oldLabel,codeList{labelRef});
    movefile(fullfile(vidDir,oldName),fullfile(vidDir,newName))
end
