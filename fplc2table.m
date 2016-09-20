function output = fplc2table(tbl,charCmpNb)
%% Loads FPLC xls data into the selected table

%UI Variables
if ~exist('charCmpNb','var') || ~isnumeric(charCmpNb);
    charCmpNb = 11;
end
charCmpNb = abs(round(charCmpNb));

if ~istable(tbl)
    disp('Input must be a table')
    return
end

warning ('off','all');%Turn of the warning

%Imports baches of selected Opus files into Matlab with the time vertor 
[fileList,dir] = uigetfile('*.xls','Select *.xls file(s) for batch Chromatogram Import','MultiSelect', 'on'); 
if iscell(fileList)==0
    fileList = {fileList};
else
    fileList= fileList';  
end %if


fileListNb = size(fileList,1);
fileTableNb = size(tbl,1);


%Check if the field exits
try idx = iscell(tbl.FlowRate);  
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'FlowRate';
end %try

try idx =iscell(tbl.Volume);  
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'Volume';
end %try
try idx =iscell(tbl.TimeFPLC) ; 
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'TimeFPLC';
end %try
try idx = iscell(tbl.UV280)  ;
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'UV280';
end %try
try idx =iscell(tbl.Cond)  ;
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'Cond';
end %try
try idx =iscell(tbl.Pressure)  ;
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'Pressure';
end %try
try idx =iscell(tbl.Temp)  ;
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'Temp';
end %try
try idx =iscell(tbl.Mark);  
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'Mark';
end %try
try idx =iscell(tbl.MarkML) ; 
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'MarkML';
end %try
try idx = iscell(tbl.MarkMin);
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'MarkMin';
end %try
try idx = iscell(tbl.RunID);
catch
idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
tbl.Properties.VariableNames{idx} = 'RunID';
end %try


%Runs throught the file list
for fileListIdx= 1:fileListNb
    fileListCur = fileList{fileListIdx};
    fileTableNb = size(tbl,1);
    try
        [~, ~, raw] = xlsread([dir,fileListCur],1,'','');%%Reads the xls file
        importSuccess =1;
    catch error
        disp(['Could not open ',fileListCur])
        disp(error.cause);
        importSuccess = 0;
    end
    
    if importSuccess == 1;
       
        %Check if the fileID matchs one of the database file ID
        fileTableIdx = 0; %Reinitialises the index
        for  idx = 1:fileTableNb 
        fileTableCur = cell2mat(tbl{idx,1});
             if strncmp(fileTableCur,fileListCur,charCmpNb) ==1
            fileTableIdx = idx;
            end %if              
        end %for fIdx = 1: filetableNb
        if fileTableIdx ==0
        fileTableIdx = fileTableNb +1;%Append
        tbl.RunID{fileTableIdx} = fileListCur;
        end %if

        for fieldIdx = 1:size(raw,2)
            fieldCur = cell2mat(raw(3,fieldIdx));         
            if strcmpi(cell2mat(raw(3,fieldIdx)),'ml')==1 && strcmpi(cell2mat(raw(3,fieldIdx+1)),' mAu')==1 %Finds the ml vector 
                tbl.Volume{fileTableIdx} = cell2mat(raw(4:end,fieldIdx));
            end %if
            if strcmpi(fieldCur,'min')==1 && strcmpi(cell2mat(raw(3,fieldIdx+1)),' mAu')==1 %Finds the ml vector 
                tbl.TimeFPLC{fileTableIdx} = cell2mat(raw(4:end,fieldIdx));
            end %if
            if strcmpi(fieldCur,' mAu')==1%Finds the UV vector 
                tbl.UV280{fileTableIdx} = cell2mat(raw(4:end,fieldIdx));
            end %if
            if strcmpi(fieldCur,' mS/cm')==1%Finds the resist vector 
                tbl.Cond{fileTableIdx} = cell2mat(raw(4:end,fieldIdx));
            end %if
            if strcmpi(fieldCur,' MPa')==1%Finds the resist vector 
                tbl.Pressure{fileTableIdx} = cell2mat(raw(4:end,fieldIdx));
            end %if
            if strcmpi(fieldCur,' C')==1%Finds the resist vector 
                tbl.Temp{fileTableIdx} = cell2mat(raw(4:end,fieldIdx));
            end %if
            if strcmpi(fieldCur,' %B')==1%Finds the resist vector 
                tbl.Conc{fileTableIdx} = cell2mat(raw(4:end,fieldIdx));
            end %if
            if strcmpi(fieldCur,'(Set Marks)')==1%Finds the resist vector 
                tbl.Mark{fileTableIdx} = raw(4:end,fieldIdx);
                idx = find(cellfun(@isstr,tbl.Mark{fileTableIdx}),1,'last');%Finds the nan
                tbl.Mark{fileTableIdx} = tbl.Mark{fileTableIdx}(1:idx-1);%Removes the nan
            end %if
            if strcmpi(fieldCur,'min')==1 && strcmpi(cell2mat(raw(3,fieldIdx+1)),'(Set Marks)')==1 %Finds the ml vector 
                vec = cell2mat(raw(4:end,fieldIdx))+1;
                vec(isnan(vec))=0;
                idx = find(vec,1,'last');%Finds the nan
                tbl.MarkMin{fileTableIdx} = vec(1:idx-1)-1;%Removes the nan
            end %if
            if strcmpi(fieldCur,'ml')==1 && strcmpi(cell2mat(raw(3,fieldIdx+1)),'(Set Marks)')==1 %Finds the ml vector 
                vec = cell2mat(raw(4:end,fieldIdx))+1;
                vec(isnan(vec))=0;
                idx = find(vec,1,'last');%Finds the nan
                tbl.MarkML{fileTableIdx} = vec(1:idx-1)-1;%Removes the nan
            end %if
        end  %For fieldIdx = 1:size(raw,2)   

        %Finds the flow rate from the even mark
        if isempty(tbl.Mark{fileTableIdx}) ==0
            idx = strncmp(tbl.Mark{fileTableIdx},'Flow',4);
            idx = find(idx==1);
            idx = idx(1);%Selects the first value
            fileListCurFlow = tbl.Mark{fileTableIdx}{idx};
            idx = strfind(fileListCurFlow,'ml');
            fileListCurFlow = fileListCurFlow(6:idx-2);
            fileListCurFlow = str2double(fileListCurFlow);
            tbl.FlowRate{fileTableIdx} = fileListCurFlow;
        else
            fileListCurFlow = 0.4;%Default rate
            tbl.FlowRate{fileTableIdx} = fileListCurFlow;
        end

        if isempty(tbl.TimeFPLC{fileTableIdx})==1 && isempty(tbl.Volume{fileTableIdx})==0
            tbl.TimeFPLC{fileTableIdx} = tbl.Volume{fileTableIdx}./fileListCurFlow;
        end % If
        
        if isempty(tbl.Volume{fileTableIdx})==1 && isempty(tbl.TimeFPLC{fileTableIdx})==0
            tbl.Volume{fileTableIdx} = tbl.TimeFPLC{fileTableIdx}.*fileListCurFlow;
        end % If

        if isempty(tbl.MarkMin{fileTableIdx})==1 && isempty(tbl.MarkML{fileTableIdx})==0
            tbl.MarkMin{fileTableIdx} = tbl.MarkML{fileTableIdx}./fileListCurFlow;
        end %if
        
        if isempty(tbl.MarkML{fileTableIdx})==1 && isempty(tbl.MarkMin{fileTableIdx})==0
            tbl.MarkML{fileTableIdx} = tbl.MarkMin{fileTableIdx}./fileListCurFlow;
        end %if
        
        
        disp(['Saved ',num2str(fileListIdx),'/',num2str(fileListNb),' in cell ', num2str(fileTableIdx),' | ',fileListCur]);
    
    end %if import sucess

end % for idxList= 1:length(fileListNb)

%Overwrites the array
%assignin('base',tableName, tbl);
output = tbl;
