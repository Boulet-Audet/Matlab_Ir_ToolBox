function ouput = opus2table(tbl,charCmpNb)
%% Imports the data from mutliple Bruker opus files into a table
if nargin <1 || nargin >2
    disp('Incorect number of arguments');
    disp('ouput = opus2table(tbl,charCmpNb)')
    return
end

if ~istable(tbl)
    disp('Input must be a table')
    return
end

if ~exist('charCmpNb','var') || ~isnumeric(charCmpNb);
    charCmpNb = 11;
end
charCmpNb = abs(round(charCmpNb));

%_________________________
%dat = evalin('base',tableVar);
warning('off','all');

%Check if the field exits
%listField = {};
try idx = iscell(tbl.RunID);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'RunID';
end %try

try idx =iscell(tbl.Bkg);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'Bkg';
end %try
try idx =iscell(tbl.Absorbance) ; 
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'Absorbance';
end %try
try idx = iscell(tbl.TimeIR)  ;
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'TimeIR';
end %try
try idx =iscell(tbl.WaveNb)  ;
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'WaveNb';
end %try
try idx =iscell(tbl.timeRel)  ;
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'timeRel';
end %try


try
    [out] = opusBatchImport; %Imports a batch of Opus files
catch;ouput = []; disp('No valid file selected');return;
end %catch
if  isempty(out)==1
    ouput = tbl;
    disp('No valid blocks in file');
    return; 
end
nbFiles = size(out.ratio,1);

fIdxStart = 1;%Initiates the range
fIdxCur = 1;
cellNB= 0;
for fIdxCur= 1:length(out.ratio.fName) %Runs throught the runID
    [~,nameCur,~] = fileparts(out.ratio.fName{fIdxCur});
    if fIdxCur < length(out.ratio.fName)
    [~,nameNext,~] = fileparts(out.ratio.fName{fIdxCur+1});
    else
    nameNext=['¬',nameCur];
    end %if
    if strcmp(nameCur(1:charCmpNb),nameNext(1:charCmpNb))~=1  %Tests if the next filename is different
        fileTableNb = size(tbl,1);%Evaluates the size of the storage array
        %Checks if the RunID matchs any existing run
        runIDIdx = 0; %Reinitialises the index
        for  idx = 1:fileTableNb 
            try
             fileTableCur = cell2mat(tbl{idx,1});
            catch
                fileTableCur = ' ';
            end %try
            if strncmp(fileTableCur,nameCur,charCmpNb) ==1
                runIDIdx = idx;
            end %if              
        end %for fIdx = 1: filetableNb
        if runIDIdx ==0
            runIDIdx = fileTableNb +1;%Append
            tbl.RunID{runIDIdx} = nameCur;%Appends the name           
        end %if
        tbl.RunID{runIDIdx} = nameCur;
        bkg = cell2mat(out.blocksAll.specRefY(fIdxStart,:));
        tbl.Bkg{runIDIdx} = bkg;
        tbl.Absorbance{runIDIdx} = out.ratio.ratioY(fIdxStart:fIdxCur,:)';
        timeRel = out.ratio.timeRel(fIdxStart:fIdxCur)';
        tbl.TimeIR{runIDIdx} = (timeRel - timeRel(1)); %Saves the range
        tbl.WaveNb{runIDIdx} = out.ratio.ratioX(fIdxStart,:)';

        fIdxStart = fIdxCur+1; %Resets the range
        runIDIdx = runIDIdx+1; %Appends the runIdx
        cellNB = cellNB+1; %Append the number of cells
        fileTableNb = size(tbl,1);%Evaluates the size of the storage array
        disp(['Saved in #',num2str(runIDIdx),'/',num2str(fileTableNb),' | ',nameCur])

    end %if  if strcmp(nameCur,nameNext)~=1
end % for
ouput = tbl;
%Display result
%disp(['Imported ',num2str(nbFiles),' Opus file(s) in ',num2str(cellNB), ' cell(s) after entry : ',num2str(runIDIdx)]);
function [out] = opusBatchImport
%Imports baches of selected Opus files into Matlab with the time vertor 
[fileList,dir] = uigetfile({'*.*','*.*'},'Select *.0000 file(s) for batch Import','MultiSelect', 'on'); 
fclose('all');%Close all open files
if iscell(fileList)==1 %many files selected
    dirFileList = cell(size(fileList));    
    for idxF = 1: size(fileList,2) %Add full path
        currentName = [dir,char(fileList(1,idxF))];
        dirFileList{idxF}= currentName;
    end
    disp([dir,' (',num2str(size(fileList,2)),') files']);
elseif fileList~=0 %Only one file selected
    dirFileList = [dir,fileList];
    fileList= {fileList};%Converts to a single cell
    dirFileList = {dirFileList};%Converts to a single cell
else
    disp('No valid files selected');
    out = [];
    return
end %if
lenL=length(fileList);

%Creates the storage table
empC=cell(lenL,1);
variableList = {'fName','interSmpY','interSmpX','interSmpPara','specRefY','specRefX','specRefPara',...
    'specSmpY','specSmpX','specSmpPara','ratioY','ratioX','ratioPara','ratioTimeAbs','ratioTimeRel'};
out.blocksAll = table(fileList',empC,empC,empC,empC,empC,empC,empC,empC,empC,empC,empC,empC,empC,empC,'VariableNames',variableList);
out.ratio = table(fileList','VariableNames',{'fName'});

%Try to import all the blocks from the selected files
for idxF = 1: lenL;
    if idxF > 49
        if idxF == round(0.2*lenL) || idxF == round(0.4*lenL)||idxF == round(0.6*lenL)|| idxF == round(0.8*lenL)
        disp(['Imported ',num2str(idxF),' of ', num2str(lenL)])
        end %if
    end %if
    fileName = dirFileList{idxF};
    fileID = fopen(fileName);
    importedBlockNames = {};
    
    listName = {'SampleInterferogram','Interferogram','InterferogramChanged','SampleInterferogramChanged'};
    for idx = 1:length(listName)
        try
            [out.blocksAll.interSmpY{idxF},out.blocksAll.interSmpX{idxF},out.blocksAll.interSmpPara{idxF}] = ImportOpus(fileName,listName{idx});
            importedBlockNames{length(importedBlockNames)+1} =listName{idx};
        catch                
        end % try
    end %for

    listName = {'ReferenceSpectrum','ReferenceSpectrumChanged'};
    for idx = 1:length(listName)
        try
            [out.blocksAll.ratioY{idxF},out.blocksAll.ratioX{idxF},out.blocksAll.ratioPara{idxF}] = ImportOpus(fileName,listName{idx});
            importedBlockNames{length(importedBlockNames)+1} =listName{idx};
        catch       
        end %try 
    end
    
    listName = {'SampleSpectrum','SampleSpectrumChanged'};
    for idx = 1:length(listName)
        try %try
            [out.blocksAll.specSmpY{idxF},out.blocksAll.specSmpX{idxF},out.blocksAll.specSmpPara{idxF}] = ImportOpus(fileName,listName{idx});
            importedBlockNames{length(importedBlockNames)+1} =listName{idx};
        catch
        end %try    
    end %For    
   
    listName = {'RatioAbsorption','Ratio','RatioAbsorptionChanged','RatioChanged'};
    for idx = 1:length(listName)       
        try
            [out.blocksAll.ratioY{idxF},out.blocksAll.ratioX{idxF},out.blocksAll.ratioPara{idxF}] = ImportOpus(fileName,listName{idx});
            importedBlockNames{length(importedBlockNames)+1} =listName{idx};
        catch       
        end %try    
    end %for    
   
    if isempty(importedBlockNames) ==1
        out = [];
        disp('No block imported'); return;
    end %if

    if idxF == 1
         out.ratio.Var2 = nan(lenL,1);
         out.ratio.Properties.VariableNames{2} = 'timeAbs';
         out.ratio.Var3 = nan(lenL,1);
         out.ratio.Properties.VariableNames{3} = 'timeRel';
         out.ratio.Var4 = nan(lenL,length(out.blocksAll.ratioY{1}));
         out.ratio.Properties.VariableNames{4} = 'ratioY';
         out.ratio.Var4 = nan(lenL,length(out.blocksAll.ratioX{1}));
         out.ratio.Properties.VariableNames{5} = 'ratioX';
    end
    try
        out.ratio.timeAbs(idxF) = out.blocksAll.ratioPara{idxF}.Instrument.SRT;
        out.ratio.timeRel(idxF) = out.blocksAll.ratioPara{idxF}.Instrument.SRT - out.blocksAll.ratioPara{1}.Instrument.SRT;
    catch
       try
           time = out.blocksAll.ratioPara{idxF}.RatioDataAbsorption.TIM;
           time = str2num(time(1:2))*3600+str2num(time(4:5))*60+str2num(time(7:8));
           day = out.blocksAll.ratioPara{idxF}.RatioDataAbsorption.DAT;
           day = str2num(day(9:10));
           time = time + (day*3600*24);
           out.ratio.timeAbs(idxF) = time;
           out.ratio.timeRel(idxF) = time - out.ratio.timeAbs(1);
       catch;disp('No time stamp found')
       end %try
    end%try
    
    try
    out.ratio.ratioY(idxF,:) = out.blocksAll.ratioY{idxF}';
    out.ratio.ratioX(idxF,:) = out.blocksAll.ratioX{idxF};
    catch
    end %Catch
    fclose('all'); % Closes the current open file
end %for



