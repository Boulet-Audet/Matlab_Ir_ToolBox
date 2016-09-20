function [out]=fileList2xl()
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
    %dirFileList = [dir,fileList];
    fileList= {fileList};%Converts to a single cell
    %dirFileList = {dirFileList};%Converts to a single cell
else
    disp('No valid files selected');
    return
end %if
out = fileList;

end %function