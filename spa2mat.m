function output = spa2mat(filename)
%% Import Omnic spa files into Matlab
try  
    % Open the file  
    fid=fopen(filename,'r');  
    % Jump where the values become interesting  
    fseek(fid,hex2dec('11e'),'bof');  
    % Pattern we're looking for  
    pattern = 6228;  
    suspect = 0;  
    while suspect~=pattern  
        oldSuspect = suspect;  
        suspect    = fread(fid,1,'int32');  
    end  
    % The correct address is just before our current suspect  
    address = oldSuspect;  
    % Close the file  
    fclose(fid);  
 catch ex  
    address = 0;  
    disp(ex)  
 end

clc filename='c:\Documents and Settings\User Name\My Documents\Spectral File.SPA';
fid=fopen(filename,'r');% Find the points number
fseek(fid,hex2dec('234'),'bof');
Number_of_DataPoints=fread(fid,1,'int32'); %Find the maximum and minimum of Wavenumber (cm-1)range 
fseek(fid,576,'bof'); Maximum_Wavenumber=fread(fid,1,'single');
Minimum_Wavenumber=fread(fid,1,'single');
Interval=(Maximum_Wavenumber-Minimum_Wavenumber)/(Number_of_DataPoints-1);
Wavenumber=linspace(Minimum_Wavenumber,Maximum_Wavenumber,Number_of_DataPoints);
Wavenumber=flipud(Wavenumber);

%Find the Y-Axis data type:
fseek(fid,hex2dec('360'),'bof');%Transmittance or Absorbance
Y_Label=char(fread(fid,14,'uchar')');% How to define the offset for spectral data still remains unresolved.
fseek(fid,hex2dec('41c'),'bof'); spectrum=fread(fid,Number_of_DataPoints,'single');%'double'); % float64, %real*8
figure(1),plot(Wavenumber,spectrum,'r');
set(gcf,'color','w');
set(gca,'xdir','rev','xcolor','b','ycolor','b','xlim',[round(Minimum_Wavenumber)??,round(Maximum_Wavenumber)]);
xlabel('Wavenumber /cm^{-1}');
ylabel(Y_Label);    

end %function

