function varargout = IRIMGGUI(varargin)
% IRIMGGUI MATLAB code for IRIMGGUI.fig
%      
%%% Gui for Varian IR Imaging

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IRIMGGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @IRIMGGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


function importImage(handles)
global data
global stop

stop =0;
while stop==0

tic % starts the timer

headerFile = get(handles.textHeaderpath,'String');%Gets the string
readHeaderFile(headerFile);%Imports the HDR parameter in the data structure

%% Reads the binary file(s)
%Creates an matrix to contain the all the raw data based on the hdr file
nbFiles = length(data.SmpFileListCur);
try
nbBkg = length(data.bkgListCur);
catch
nbBkg = 0;
imgAllBkg = nan;
displayString = 'No background selelected';
disp(displayString);set(handles.guiEntire,'Name',displayString);
end %try
%Creates an anonymus function to clear variable and free memory
fClear = @(x) clear(inputname(1));

%If data already exist
if isempty(data)~=1
strCondition = questdlg('Do you wish to overwrite previous data of append the dataset?', ...
    'Overwrite or append ?','Overwrite','Append','Overwrite');
end %if


if strcmp(strCondition,'Append')==1
   try
   nbFilesPrevious = size(data.imgAll,1);
   catch
   nbFilesPrevious = 0;
   end %try
else % Overwrite
    nbFilesPrevious = 0;
    try
    fClear(data.imgAll);%Clears the previous data to save memory
    catch
    disp('Data not found to clear memory')
    end %Try
end %if

%Imports the background files
if nbBkg~=nbFiles && nbBkg>1;
nbBkg = 1; 
displayString = ['The number of background is different from the number of samples image,',...
'only the first bkg will be used'];
disp(displayString);set(handles.guiEntire,'Name',displayString);
end %ifnbBkg~= nbFiles


normHigh = str2double(get(handles.editNormHigh,'String'));%Gets the Gui parameters
if isnan(normHigh)==1;normHigh=1900;end;%Tests the GUI parameter
normLow = str2double(get(handles.editNormLow,'String'));%Gets the Gui parameters
if isnan(normLow)==1;normLow=1870;end;%Tests the GUI parameter
if normLow>normHigh;temp=normHigh;normHigh = normLow;normLow=temp;end%Tests the GUI parameter
set(handles.editNormHigh,'String',num2str(normHigh));%Resets the tested GUI value
set(handles.editNormLow,'String',num2str(normLow));%Resets the tested GUI value

normOffsetHigh = str2double(get(handles.editNormOffsetHigh,'String'));%Gets the Gui parameters
if isnan(normOffsetHigh)==1;normOffsetHigh=750;end;%Tests the GUI parameter
normOffsetLow = str2double(get(handles.editNormOffsetLow,'String'));%Gets the Gui parameters
if isnan(normOffsetLow)==1;normOffsetLow=700;end;%Tests the GUI parameter
if normOffsetLow>normOffsetHigh;temp=normOffsetHigh;normOffsetHigh = normOffsetLow;normOffsetLow=temp;end%Tests the GUI parameter
set(handles.editNormOffsetHigh,'String',num2str(normOffsetHigh));%Resets the tested GUI value
set(handles.editNormOffsetLow,'String',num2str(normOffsetLow));%Resets the tested GUI value

if  get(handles.popupmenuRatio,'value')~=1;
    for idxBkgFile = 1:nbBkg % Loop through the files for import/apodisation/phase correction
        displayString = ['(',num2str(idxBkgFile),'/',num2str(nbBkg),') Import backgroud'];
        disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
        try
            currentBkgName=data.bkgListCur{idxBkgFile};
            catch
            displayString = 'No background file selected, import aborded';
            disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            return
        end %try
      
        imgCurrent = readImage(handles,currentBkgName);  %Reads and compute the background file
        
        %Normalises the current image using the selected zero and one range to
        %correct for intensity fluctuations
        checkboxSingleNorm = get(handles.checkboxSingleNorm,'Value');
        if checkboxSingleNorm == 1
            displayString ='Normalising the background single beam'; 
            disp(displayString);%Display operation
            set(handles.guiEntire,'Name',displayString);
            wavenumberVec = data.wavenumberVec;%Gets the wavenumber vector
            [~,normOffsetLowIdx] = min(abs(wavenumberVec-normOffsetLow));%Finds the index
            [~,normOffsetHighIdx] = min(abs(wavenumberVec-normOffsetHigh));%Finds the index
            imgAverageSpectrumOffset = mean(imgCurrent(:,:,normOffsetLowIdx:normOffsetHighIdx),3);%Extracts the offset image
            imgAverageSpectrumOffset = repmat(imgAverageSpectrumOffset,[1 1 size(imgCurrent,3)]);%Replicates the array
            [~,normHighIdx] = min(abs(wavenumberVec-normHigh));%Finds the index
            [~,normLowIdx] = min(abs(wavenumberVec-normLow));%Finds the index
            imgAverageSpectrumNorm = mean(imgCurrent(:,:,normLowIdx:normHighIdx),3);%Extracts the offset image
            imgAverageSpectrumNorm = repmat(imgAverageSpectrumNorm,[1 1 size(imgCurrent,3)]);%Replicates the norm array
            imgAverageSpectrumNorm = imgAverageSpectrumNorm - imgAverageSpectrumOffset;%Substract the offsets
            imgCurrent = (imgCurrent-imgAverageSpectrumOffset)./imgAverageSpectrumNorm;%Normalises the single beam
        end %if checkboxSingleNorm == 1
        
        %Aggregates the image using the selected factor
        %To write

        % Saving computed data
        if idxBkgFile ==1;% if first bkg image
        try
            imgAllBkg= nan(nbBkg,size(imgCurrent,1),size(imgCurrent,2),size(imgCurrent,3));%Creates an array containing all the background spectra images
            catch lasterr
            displayString = ['Not enough memory! ' ,lasterr.message];
            disp(displayString);set(handles.guiEntire,'Name',displayString);
            return
        end

        end % if first image
        %Saves the treated data into the storage array 
        imgAllBkg(idxBkgFile,:,:,:)= imgCurrent;
    end % loops through the bkg files
end %if No ratio is selected

%% Imports the sample files
%Run through the selected files
idxBkgFile =1; %Sets the background index
for idxFile = 1:nbFiles % Loop through the files for import/apodisation/phase correction
    displayString = ['(',num2str(idxFile),'/',num2str(nbFiles),') Import Image'];
    disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
    currentFileName=data.SmpFileListCur{idxFile};
    imgCurrent = readImage(handles,currentFileName); % Read image function
    listing = dir(currentFileName);%Gets the file info
    imgModTimeNum = listing.datenum;%Gets the time number
    imgModTimeStr = listing.date;%Gets the time string
    %Normalises the current image using the selected zero and one range to
    %correct for intensity fluctuations
    checkboxSingleNorm = get(handles.checkboxSingleNorm,'Value');
    if checkboxSingleNorm == 1
        displayString ='Normalising the single beam'; 
        disp(displayString);set(handles.guiEntire,'Name',displayString);%Display operation
        set(handles.guiEntire,'Name',displayString);
            wavenumberVec = data.wavenumberVec;%Gets the wavenumber vector
            [~,normOffsetLowIdx] = min(abs(wavenumberVec-normOffsetLow));%Finds the index
            [~,normOffsetHighIdx] = min(abs(wavenumberVec-normOffsetHigh));%Finds the index
            imgAverageSpectrumOffset = mean(imgCurrent(:,:,normOffsetLowIdx:normOffsetHighIdx),3);%Extracts the offset image
            imgAverageSpectrumOffset = repmat(imgAverageSpectrumOffset,[1 1 size(imgCurrent,3)]);%Replicates the array
            [~,normHighIdx] = min(abs(wavenumberVec-normHigh));%Finds the index
            [~,normLowIdx] = min(abs(wavenumberVec-normLow));%Finds the index
            imgAverageSpectrumNorm = mean(imgCurrent(:,:,normLowIdx:normHighIdx),3);%Extracts the offset image
            imgAverageSpectrumNorm = repmat(imgAverageSpectrumNorm,[1 1 size(imgCurrent,3)]);%Replicates the norm array
            imgAverageSpectrumNorm = imgAverageSpectrumNorm - imgAverageSpectrumOffset;%Substract the offsets
            imgCurrent = (imgCurrent-imgAverageSpectrumOffset)./imgAverageSpectrumNorm;%Normalises the single beam
    end %if checkboxSingleNorm == 1
      
    %% Ratio the spectra
    operationCurrent = get(handles.popupmenuRatio,'String');
    idxValue = get(handles.popupmenuRatio,'Value');
    operationCurrent = operationCurrent{idxValue};

    if strcmp(operationCurrent,'No ratio')==0;
        if size(imgCurrent,1)==size(imgAllBkg,2)&&...
            size(imgCurrent,2)==size(imgAllBkg,3)&&...
                size(imgCurrent,3)==size(imgAllBkg,4)
            displayString = 'Compatible sample and background images for ratio calculation';
            disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            else
            displayString = 'Incompatible sample and background images for ratio calculation';
            disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            operationCurrent = 'No ratio';
        end %  strcmp(handleString(idxValue),'No ratio')==0;
    end %if a ratio mode is selected

    if strcmp(operationCurrent,'Absorbance')==1;
    displayString =['(',num2str(idxFile),'/',num2str(nbFiles),') Absorbance ratio']; 
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    imgCurrent= -1*log10(abs(imgCurrent./squeeze(imgAllBkg(idxBkgFile,:,:,:))));
    data.responseUnit = 'Absorbance';

    elseif strcmp(operationCurrent,'Transmittance')==1;
    displayString =['(',num2str(idxFile),'/',num2str(nbFiles),') Transmittance ratio']; 
    disp(displayString); set(handles.guiEntire,'Name',displayString);
    imgCurrent= (imgCurrent./squeeze(imgAllBkg(idxBkgFile,:,:,:)));
    data.responseUnit = 'Transmittance';

    elseif strcmp(operationCurrent,'Reflectance')==1;
    displayString = ['(',num2str(idxFile),'/',num2str(nbFiles),') Reflectance ratio']; 
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    imgCurrent= -1*log10(imgCurrent./squeeze(imgAllBkg(idxBkgFile,:,:,:)));
    data.responseUnit = 'Reflectance';

    elseif strcmp(operationCurrent,'Substraction')==1;
    displayString = ['(',num2str(idxFile),'/',num2str(nbFiles),') Substration ratio']; 
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    imgCurrent= (imgCurrent-squeeze(imgAllBkg(idxBkgFile,:,:,:)));
    data.responseUnit = 'Substraction';

    elseif strcmp(operationCurrent,'Addition')==1;
    displayString = ['(',num2str(idxFile),'/',num2str(nbFiles),') Addition ratio']; 
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    imgCurrent= (imgCurrent+squeeze(imgAllBkg(idxBkgFile,:,:,:)));
    data.responseUnit = 'Addition';

    else % No ratio
    displayString = ['(',num2str(idxFile),'/',num2str(nbFiles),') No ratio calculation']; 
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    end

    %Increments the background file index
    if nbBkg == nbFiles && nbBkg>1
    idxBkgFile = idxBkgFile +1;
    end

    % Saving computed data
    if idxFile ==1;% if first image
        vecWave = data.wavenumberVec;%Imports the wave vector
        windowLow = str2double(get(handles.editImportImgLow,'string'));%Imports the GUI parameters
        windowHigh = str2double(get(handles.editImportImgHigh,'string'));%Imports the GUI parameters
        if windowLow<min(vecWave)||windowLow>max(vecWave);windowLow=min(vecWave);end;%Test variable
        if windowHigh<min(vecWave)||windowHigh>max(vecWave);windowHigh=max(vecWave);end;%Test variable
        if windowLow>windowHigh;temp=windowLow;windowLow=windowHigh;windowHigh=temp;end;%Swap variable if inverted 
        [~,idxWindowLow]=min(abs(vecWave - windowLow));%Finds the window index
        [~,idxWindowHigh]=min(abs(vecWave - windowHigh));%Finds the window index
        imgCurrent=imgCurrent(:,:,idxWindowLow:idxWindowHigh);%Truncates the spectral range of the first frame  
        if nbFilesPrevious == 0
            try
                data.timeAbs = nan(nbFiles,1);
                data.timeRel = nan(nbFiles,1);
                data.SmpFileList = data.SmpFileListCur;
                data.imgAll= nan(nbFiles,size(imgCurrent,1),size(imgCurrent,2),size(imgCurrent,3));%Preallocates the memory for the storage array
            catch lasterr
                displayString = ['Not enough memory! Use a 64 bit machine! ' ,lasterr.message];
                disp(displayString);set(handles.guiEntire,'Name',displayString);
            end %try
        else
            try
                data.imgAll = cat(1,data.imgAll,nan(nbFiles,size(imgCurrent,1),size(imgCurrent,2),size(imgCurrent,3)));%Appends the array
                data.SmpFileList = [data.SmpFileList,data.SmpFileListCur];
                data.timeAbs = [data.timeAbs,nan(nbFiles,1)];
                data.timeAbs = [data.timeRel,nan(nbFiles,1)];
            catch lasterr
                displayString = ['Not the same size as previous data! ',lasterr.message];
                disp(displayString);set(handles.guiEntire,'Name',displayString);
            end %Catch
         end %if
    else
    imgCurrent=imgCurrent(:,:,idxWindowLow:idxWindowHigh);%Truncates the spectral range  
    end  %if idxFile ==1;% if first image
    %%Saves the treated data into the storage array 
    data.imgAll(idxFile+nbFilesPrevious,:,:,:)= imgCurrent;
    data.timeAbs(idxFile+nbFilesPrevious) = imgModTimeNum;
    data.timeRel(idxFile+nbFilesPrevious) = 24*3600*(imgModTimeNum-data.timeAbs(1));%Saves the time difference and converts in seconds
    
end % Loop through the files for import/apodisation/phase correction

data.imgAllSize = size(data.imgAll);
data.wavenumberVec = vecWave(idxWindowLow:idxWindowHigh);%Saves the new wave vector

timeElapsed = toc; %Ends the timer

displayString = [num2str(length(data.SmpFileList)),' images imported sucessfully in '...
    ,num2str(round(timeElapsed)),' seconds'];
set(handles.guiEntire,'Name',displayString); disp(displayString);drawnow;

enableGui(handles);%Enables the remaining functions
integrateImageRegions(handles) %Integrates the images and refreshed the graphs

stop =1;
end% while stop=='0'

function imgSpectra = readImage(handles,curFullFileName)
global data

displayStringPrevious = get(handles.guiEntire,'Name');

%% Reads the interferogram of the image file selected
try
[~,fileName,ext]= fileparts(curFullFileName);
displayString = ['Reading : ',fileName,ext];
set(handles.guiEntire,'Name',displayString); disp(displayString);drawnow;
img= multibandread(curFullFileName,[data.samples,data.lines,data.bands],...
    'single',data.headerOffest, data.interleave,'ieee-le');
catch lasterr
displayString = ['Not a valid file! ' ,lasterr.message];
disp(displayString);set(handles.guiEntire,'Name',displayString);
return
end % try

waveVector = data.waveVector;%Imports the waveVector
if min(waveVector)>0;%Tests if the file selected is a spectrum
displayString = 'Spectrum selected, saving without computing';
disp(displayString);imgSpectra = img;%Returns the image
data.waveVectorType = 'Wavenumber';
data.wavenumberVec = data.waveVector;%Sets the data vector
return
end % 

% Inverse the interferogram if negative a positive center burst required for FFT
averageInter = squeeze(mean(mean(img,1),2));%Calculates the average Interferogram
averageInterOffset = averageInter -mean(averageInter);%Calculates the average Interferogram
if abs(min(averageInterOffset))>=abs(min(averageInterOffset));
    img=((img-mean(averageInter)).*-1)+mean(averageInter);
end %if

%% Offsets the Interferograms
% The mean of the data is subtracted from the data,
% which brings the baseline of the data to 0. This is needed because
% later the data is zero-filled, and not doing this would introduce a
% discontinuity into the data. A discontinuity, or spike, would result in
% a sine wave moving through the FFT'd data.

if get(handles.checkboxInterZeroOffset,'Value') == 1
    displayString = [displayStringPrevious,'Interferogram offset'];
    disp(displayString);drawnow;
    averageImage = squeeze(mean(img,3));
    averageImage=repmat(averageImage,[1,1,size(img,3)]);
    img = img - averageImage;%Substracts the average intensity of all interferogram in the image    
end %if checkBoxNorma == 1

%% Interferogram base correction fitting
if get(handles.checkboxBaseline,'Value') == 1
displayString = [displayStringPrevious,' | interferogram correction using degree 3 polynomial'];
disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
for idxSample = 1:size(img,1); %Loops through the samples (pixel rows)
for idxLine = 1:size(img,2);%Loops trhough the lines  
cWave = squeeze(img(idxSample,idxLine,:)); 
%Substract a polynome 3 fit
img(idxSample,idxLine,:)=cWave-polyval(polyfit((1:size(img,3)),cWave,3),(1:size(img,3)));
end %Loops trhough the lines    
end %Loops through the samples (pixel rows)
end %if Interferogram baseline correction

%% Centerburst shift correction using spline fitting 
if get(handles.checkboxCburstCorr,'Value')== 1
imgCShift = nan(size(img,1),size(img,2));%Allocates the matrix
[~,cBurst]=max(abs(img),[],3);%Finds the position of the most intense pixel
cBurst=round(mean(mean(cBurst)));%Finds the ZPD/center burst of the most intense pixel
pW1 = 2^(floor(log2(cBurst-1)));%Calculates the largest power of two if the cBurst 
    for idxSample = 1:size(img,1); %Loops through the samples (pixel rows)
        for idxLine = 1:size(img,2);%Loops trhough the lines
            if idxSample== 1 && idxLine ==1
            displayString = [displayStringPrevious,' (0%) Center Burst Shift Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            elseif idxSample== round(size(img,1)/5) && idxLine ==1
            displayString = [displayStringPrevious,' (20%) Center Burst Shift Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            elseif idxSample== round((2*size(img,1))/5) && idxLine ==1
            displayString = [displayStringPrevious,' (40%) Center Burst Shift Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            elseif idxSample== round((3*size(img,1))/5) && idxLine ==1
            displayString = [displayStringPrevious,' (60%) Center Burst Shift Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            elseif idxSample== round((4*size(img,1))/5) && idxLine ==1
            displayString = [displayStringPrevious,' (80%) Center Burst Shift Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
            end %if
            spectrumInterf = squeeze(img(idxSample,idxLine,:));%Extract the interferogram
            %Calculates the width of the interferogram
            cMass = ((cBurst-pW1:cBurst+pW1-1).*spectrumInterf(cBurst-pW1:cBurst+pW1-1)')/mean(spectrumInterf(cBurst-pW1:cBurst+pW1-1)',1);%Finds the exact centre of gravity of the interferogram
            cShift = cMass-cBurst;%Calculates the center Burst shift
            imgCShift(idxSample,idxLine) = cShift; %Saves to the matrix
            spectrumInterfCorr = spline(1:length(spectrumInterf),spectrumInterf,(1:length(spectrumInterf))+cShift)';%Interpolates using a spline the interpolated value
            img(idxSample,idxLine,:) = spectrumInterfCorr; %Saves the interpolate interferogram to the array
        end %Loops trhough the lines   
    end %Loops through the samples (pixel rows)
end % checkboxCburstCorr==1

%% Zero-Fill the Interferogram
% The FFT requires that the number of points in the interferogram
% be a power of 2. If it is not, then the size of the array is increased
% to the next higher power of two and the added points are set to 0
% The array size can then be increased by the zero-fill factor (2, 4 ...).
% For example, if the interferogram contains 16,000 points, the array
% will be padded with zeros to make it 16,384 points. If the zero-fill
% factor is 4 the array will be padded with zeros out to 65,536 points.
% Zero-filling the interferogram before the FFT results in interpolated
% points that lie between the 'actual' points that would result from
% non-zero-filled data, making the resulting data look smoother.
% This technique also improves the photometric accuracy of the data,
% because a computed data point will rarely correspond to the maximum
% absorbtion of a given spectral feature.
% Although zero-filling increases the number of data points in the
% spectrum, it cannot increase the optical resolution of the spectrum.

zeroFillingFactor = get(handles.popupmenuZeroFilling,'Value');
if zeroFillingFactor == 2;
    zeroFillingFactor =3;
    elseif zeroFillingFactor == 3;
    zeroFillingFactor =4;
    elseif zeroFillingFactor == 4;
    zeroFillingFactor =1;
    else %=2x by default
    zeroFillingFactor =2;    
end %zeroFillingFactor == 2;

interTrunc = abs(floor(str2double(get(handles.editForwardTruncation,'string'))));%Get the GUI variable
if isnan(interTrunc)==1;interTrunc=2;end;%Tests the value
if interTrunc>30;interTrunc=2;end%Tests the value
set(handles.editForwardTruncation,'string',num2str(interTrunc));%Sets the tested value

displayString = ['Zero filling using ',num2str(zeroFillingFactor)];
disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;

nbRawPt = size(img,3);
nb0Pt = 2^nextpow2(size(img,3))*zeroFillingFactor; % Next power of 2 from length of y
zeroFilled = zeros(size(img,1),size(img,2),nb0Pt);
zeroFilled(:,:,1:size(img,3)-interTrunc)=img(:,:,1+interTrunc:end);
img = zeroFilled;
%Extrapolate the waveVector
%waveFit=polyfit((1:length(waveVector))',waveVector,1);
%waveVector=(((1:nb0Pt)*(waveFit(1)))+(waveFit(2)))';

boxCheckedCompute = get(handles.checkboxCompute,'Value');
if boxCheckedCompute ==0
    displayString = 'Saving image without computing';
    disp(displayString);
    imgSpectra = img;%Returns the image
    data.waveVectorType = 'Optical retardation';
    data.wavenumberVec = data.waveVector;%Sets the data vector
return ;%Exists the sub
end %if boxCheckedCompute ==0

boxChecked=get(handles.checkboxMertz,'Value');
% Selects the type of apodisation
apodizationType =get(handles.popupmenuApodisation,'Value');

imgSpectra= nan(size(img,1),size(img,2),(size(img,3)/2)); %Creates an array to contain all the spectra images
%imgInterCorr= nan(size(img,1),size(img,2),(size(img,3)));

for idxSample = 1:size(img,1); %Loops through the samples (pixel rows)
    for idxLine = 1:size(img,2);%Loops trhough the lines  
    if idxSample== 1 && idxLine ==1  && boxChecked == 1
        displayString = [displayStringPrevious,' (0%) FFT & Phase Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
        elseif idxSample== round(size(img,1)/5) && idxLine ==1  && boxChecked == 1
        displayString = [displayStringPrevious,' (20%) FFT & Phase Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
        elseif idxSample== round((2*size(img,1))/5) && idxLine ==1 && boxChecked == 1
        displayString = [displayStringPrevious,' (40%) FFT & Phase Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
        elseif idxSample== round((3*size(img,1))/5) && idxLine ==1 && boxChecked == 1
        displayString = [displayStringPrevious,' (60%) FFT & Phase Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
        elseif idxSample== round((4*size(img,1))/5) && idxLine ==1 && boxChecked == 1
        displayString = [displayStringPrevious,' (80%) FFT & Phase Correction'];disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
    end
     
    spectrumInterf = squeeze(img(idxSample,idxLine,:));%Extract the interferogram
    [~,cBurst]=max(abs(spectrumInterf));%Finds the ZPD/ center burst
    if cBurst >(length(spectrumInterf)*0.5); %If the center burst is on the second half of the vector
        spectrumInterf = flipr(spectrumInterf);% flips the vector
        [~,cBurst]=max(abs(spectrumInterf));%Finds the ZPD/ center burst
    end %if
    pW1 = 2^(floor(log2(cBurst-1)));%Calculates the largest power of two of the cBurst 
    pW2 = 2*pW1; %Defines the apodisation window
    nbLen=nbRawPt-cBurst;%Calculates the lenght of the longest branch

        %%Calculates the apodisation function of the interferogram/phase correction vector
    if apodizationType ==4 %Triangular apodisation
        spectrumApodVec = [((abs((1:nbRawPt)'-cBurst)/(nbLen)-1)*-1);zeros(nb0Pt-nbRawPt)];
        phaseApodVec = (abs((1:pW2)'-(pW1))/(pW1)-1)*-1;
        elseif apodizationType ==5 %Box-car apodization
        phaseApodVec = ones(pW2,1);
        spectrumApodVec = [ones(nbRawPt,1);zeros(nb0Pt-nbRawPt,1)];
        elseif apodizationType ==3 %Happ-Genzel apodization
        spectrumApodVec = [(0.54+(0.46*cos(pi*((1:nbRawPt)'-cBurst)/(nbLen))));zeros(nb0Pt-nbRawPt,1)];    
        phaseApodVec=0.54+(0.46*cos(pi*((1:pW2)'-(pW1))/(pW1)));
        elseif apodizationType ==2 %%Blackmann-H3
        spectrumApodVec=0.42323+(0.49755*cos(pi*((1:nbRawPt)'-cBurst)/(nbLen)))+(0.07922*cos(2*pi*((1:nbRawPt)'-cBurst)/(nbLen)));
        spectrumApodVec = [spectrumApodVec;zeros(nb0Pt-nbRawPt,1)];
        phaseApodVec=0.42323+(0.49755*cos(pi*((1:pW2)'-(pW1))/(pW1)))+(0.07922*cos(2*pi*((1:pW2)'-(pW1))/(pW1)));
        else %Blackmann-H4
        spectrumApodVec=0.35875+(0.48829*cos(pi*((1:nbRawPt)'-cBurst)/(nbLen)))+(0.14128*cos(2*pi*((1:nbRawPt)'-cBurst)/(nbLen)))+(0.01168*cos(3*pi*((1:nbRawPt)'-cBurst)/(nbLen)));
        spectrumApodVec = [spectrumApodVec;zeros(nb0Pt-nbRawPt,1)];
        phaseApodVec=0.35875+(0.48829*cos(pi*((1:pW2)'-(pW1))/(pW1)))+(0.14128*cos(2*pi*((1:pW2)'-(pW1))/(pW1)))+(0.01168*cos(3*pi*((1:pW2)'-(pW1))/(pW1)));
    end % Selection
    
    if boxChecked == 1 %If the box is checked perform the Mertz phase correction
    %%Mertz Phase-Correction, apodization and Fourier transform
    % The interferogram contains out-of-phase elements which are introduced
    % by optical path differences in the instrument. The data must be phase
    % corrected or the resulting spectrum will not be photometrically accurate.
    % The phase correction method used here is the Mertz method.
    % An array of the same size as the zero-filled interferogram is prepared
    % and filled with zeros.

    % The points near around ZPD are copied from the interferogram
    % into the new array, but are scaled by the appodization function
    % The interferogram data is rotated so that the right side of the data
    % after the ZPD, including the ZPD point, is moved to the front of the array.    
    phaseCorrAr = spectrumInterf(cBurst-pW1:cBurst+pW1-1);%Creates the vector
    phaseCorrAr= phaseCorrAr.*phaseApodVec;%Applies the apodzation function
    phaseCorrAr=[phaseCorrAr((pW1)+1:end);phaseCorrAr(1:(pW1))];%Rotates the phase correction array
    phaseCorrArComplex=fft(phaseCorrAr);%The data is FFT’d, producing real and imaginary data arrays.
    phaseCorrArReal=real(phaseCorrArComplex);%Extracts the real part 
    phaseCorrArImag=imag(phaseCorrArComplex);%Extracts the imaginary part 
    arPhaseCorr= atan2(phaseCorrArImag,phaseCorrArReal);%Calculates the phase correction for the PW
    spectrumPhase = spline((0:1/((pW2)-1):1)',arPhaseCorr,(0:1/(length(spectrumInterf)-1):1)');%Interpolate the phase curve to full resolution using a spline function
    spectrumCosTerm=cos(spectrumPhase);%Calculate the phase correction COS terms from the partial interferogram
    spectrumSinTerm=sin(spectrumPhase);%Calculate the phase correction SIN terms from the partial interferogram
    spectrumInterf=spectrumInterf.*spectrumApodVec;%Apply the apodisation function on the full interferogram
    spectrumInterf=[spectrumInterf(cBurst:end);spectrumInterf(1:(cBurst-1))];%Rotates the full interferogram
    spectrumComplex=fft(spectrumInterf);%Performs the fourier transform of the rotated interferogram
    spectrumReal=real(spectrumComplex);%Gets the real and imaginary parts of the complex spectrum
    spectrumImag=imag(spectrumComplex);%Gets the real and imaginary parts of the complex spectrum
    spectrumReal = spectrumReal.*spectrumCosTerm;%Corrects for the phase shifts
    spectrumImag = spectrumImag.*spectrumSinTerm;%Corrects for the phase shifts
    curSpectrum = spectrumReal+ spectrumImag;%Combines the real and imaginary parts
    
    if get(handles.checkboxInterNorm,'Value')==1;
        %The imaginary part should become a zero vector and the FFT of the spectrum
        interCorrected = real(fft(curSpectrum)); %Should result in a perfectly symetrical interferogram
        interCorrected = interCorrected./mean(abs(interCorrected(1:length(interCorrected)/8)));%Normalises the interfeogram to 1
        curSpectrum = real(fft(interCorrected));% reperforms the FFT
        %imgInterCorr(idxSample,idxLine,:)= [interCorrected(1+(length(interCorrected)*0.5):end);interCorrected(1:(length(interCorrected)*0.5))]; % Rotates the corrected interferogram
    end %if get(handles.checkboxMertz,'Value')==1;
    %Saves the corrected spectra
    curSpectrum=curSpectrum(1:length(curSpectrum)/2);%Remove the symetric half
 
    else  %No phase correction
    spectrumInterf=spectrumInterf.*spectrumApodVec;%Apply the apodisation function on the full interferogram
    spectrumInterf=[spectrumInterf(cBurst:end);spectrumInterf(1:(cBurst-1))];%Rotates the full interferogram
    spectrumComplex=fft(spectrumInterf);%Performs the fourier transform of the rotated interferogram
    spectrumReal=real(spectrumComplex);%Gets the real and imaginary parts of the complex spectrum
    spectrumImag=imag(spectrumComplex);%Gets the real and imaginary parts of the complex spectrum   
    curSpectrum = spectrumReal+ spectrumImag;    %Combines the real and imaginary parts
    end %if boxChecked == 1
    
    %imgCBurst(idxSample,idxLine)= cBurst; %Saves the Cburst
    %imgSpectra(idxSample,idxLine,:) = interCorrected(1000:length(interCorrected)/2+999); %Saves the single beam into the storage array.
    imgSpectra(idxSample,idxLine,:) = curSpectrum; %Saves the single beam into the storage array.
    end %Loops through the lines  
end %Loops through the samples (pixel rows)

data.waveVectorType = 'Wavenumber';

nbPoints = size(imgSpectra,3);
resolution = round((4*0.224)./max(waveVector));%Estimates the resolution based on the maximum OR
%uDR = round(((2038*16)./resolution)./nbPoints);%Estimates the UDR based on the resolution and number of points
uDR = round(((2038*16)./resolution)./size(waveVector,1));
wavenumberVec = (1:nbPoints)'.*(1/(nbPoints*uDR*6.33*10^(-5))); % Calculates the frequency vector using the HeNe wavelength 633 nm
data.wavenumberVec = wavenumberVec; %Saves the wavenumber vector
data.resolution = resolution; %Saves the resolution
data.uDR = uDR; %Saves the UDR

function readHeaderFile(headerFile)
  global data
   disp('Read header File')
filename = headerFile;
   delimiter = '=';
startRow = 4;
endRow = 13;
% Format string for each line of text:
formatSpec = '%s%s%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);
dataArrayValues = dataArray{1,2};
data.samples = str2double(dataArrayValues{1});
data.lines =  str2double(dataArrayValues{2});
data.bands =  str2double(dataArrayValues{3});
data.headerOffest =  str2double(dataArrayValues{4});
data.fileType = dataArrayValues{5};
data.dataType =  str2double(dataArrayValues{6});
data.interleave = dataArrayValues{7};
data.sensorType = dataArrayValues{8};
data.byteOrder =  str2double(dataArrayValues{9});
data.stretch = dataArrayValues{10};

%Imports the wave vector
delimiter = ',';
startRow = 15;
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
waveVector = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
waveVector=str2double(waveVector{1});
waveVector=waveVector(1:end-1);%removed the last caracter
data.waveVector = waveVector;
if min(waveVector)<0
data.waveVectorType ='Optical retardation';
else
data.waveVectorType ='Wavenumber';
end %if

%% Close the text file.
fclose(fileID); 
disp(['Header file imported successfully']);

function processingSequence(handles)

% vecWave = data.wavenumberVec; %Imports the wavenumber vector
% for idxImg = 1:size(data.imgAll,1)%Loops through the image stack
% arraySize = size(data.imgAll);%Gets the size of the array
% imgCurrent=reshape(data.imgAll(idxImg,:,:,:),arraySize(2:4));%Imports the current image
global data

for idxOp = 1:5 %loopes through the spectra operations

    handleOp = findobj('Tag',['popupmenuOp',num2str(idxOp)]);    %Selects the handles
    handleTopBox = findobj('Tag',['editOp',num2str(idxOp),'Box1']);    %Selects the handles
    handleLowerBox = findobj('Tag',['editOp',num2str(idxOp),'Box2']);    %Selects the handles
    operationCurrent = get(handleOp,'string');%Gets the operation to be performed
    operationCurrent = operationCurrent{get(handleOp,'value')};%Gets the operation to be performed

    switch operationCurrent
    
    % Truncates the spectral range
    case 'Truncate' 
        intervalHigher=str2double(get(handleTopBox,'string'));%Gets the GUI string
        intervalLower=str2double(get(handleLowerBox,'string'));%Gets the GUI string
        if isnan(intervalHigher)==1; intervalHigher= 2000;end;%Tests the input parameters
        if isnan(intervalLower)==1; intervalLower= 900;end;%Tests the input parameters
        if intervalHigher<intervalLower
            intervalTemp = intervalHigher; intervalHigher = intervalLower;
            intervalLower = intervalTemp;%Swap values
        end %if intervalHigher<intervalLower 
        vecWave = data.wavenumberVec;%Imports the wavemumber vector
        if intervalHigher > max(vecWave);intervalHigher=max(vecWave);end
        if intervalLower < min(vecWave);intervalLower=min(vecWave);end
        set(handleTopBox,'string',num2str(intervalHigher));%Resets the values
        set(handleLowerBox,'string',num2str(intervalLower));%Resets the values
        [~,intervalHigherPos]=min(abs(vecWave-intervalHigher));%Identifies the wavenumber interval index
        [~,intervalLowerPos]=min(abs(vecWave-intervalLower));%Identifies the wavenumber interval index
        data.wavenumberVec = vecWave(intervalLowerPos:intervalHigherPos);%Trucates the wavenumber vector
        data.imgAll = data.imgAll(:,:,:,intervalLowerPos:intervalHigherPos);%Truncates the image using the new wavenumber vector
        displayString = ['Data truncated from '...
            ,num2str(intervalHigher),' to ', num2str(intervalLower),' cm-1.']; 
        disp(displayString);%Display operation
        set(handles.guiEntire,'Name',displayString);
        
    %Offset the average absolute intensity of a spectral region to image
    case 'Offset'
        intervalHigher=str2double(get(handleTopBox,'string'));%Gets the GUI string
        intervalLower=str2double(get(handleLowerBox,'string'));%Gets the GUI string
        if isnan(intervalHigher)==1; intervalHigher= 1900;end;%Tests the input parameters
        if isnan(intervalLower)==1; intervalLower= 1850;end;%Tests the input parameters
        if intervalHigher<intervalLower
            intervalTemp = intervalHigher; intervalHigher = intervalLower;
            intervalLower = intervalTemp;
        end %if intervalHigher<intervalLower 
        vecWave = data.wavenumberVec;%Imports the wavemumber vector
        if intervalHigher > max(vecWave);intervalHigher=max(vecWave);end
        if intervalLower < min(vecWave);intervalLower=min(vecWave);end
        set(handleTopBox,'string',num2str(intervalHigher));%Resets the values
        set(handleLowerBox,'string',num2str(intervalLower));%Resets the values
        [~,intervalHigherPos]=min(abs(vecWave-intervalHigher));%Identify the position of the nearest value
        [~,intervalLowerPos]=min(abs(vecWave-intervalLower));%Identify the position of the nearest value
        imgAverage = mean(data.imgAll(:,:,:,intervalLowerPos:intervalHigherPos),4);%Calculates the average intensity of the region selected
        data.imgAll =data.imgAll-repmat(imgAverage,[1,1,1,size(data.imgAll,4)]);%Substract the average intensity from the absolute intensity
        displayString = ['Offset (substract) from '...
            ,num2str(intervalHigher),' to ', num2str(intervalLower),' cm-1 performed.']; 
        disp(displayString);%Display operation
        set(handles.guiEntire,'Name',displayString);
         
    %Normalised the average absolute intensity of a spectral region to image
    case 'Normalise'
    intervalHigher=str2double(get(handleTopBox,'string'));%Gets the GUI string
    intervalLower=str2double(get(handleLowerBox,'string'));%Gets the GUI string
    if isnan(intervalHigher)==1; intervalHigher= 1850;end;%Tests the input parameters
    if isnan(intervalLower)==1; intervalLower= 950;end;%Tests the input parameters
    if intervalHigher<intervalLower
        intervalTemp = intervalHigher; intervalHigher = intervalLower;
        intervalLower = intervalTemp;
    end %if intervalHigher>intervalLower 
    vecWave = data.wavenumberVec;%Imports the wavemumber vector
    if intervalHigher > max(vecWave);intervalHigher=max(vecWave);end;%Tests the input parameters
    if intervalLower < min(vecWave);intervalLower=min(vecWave);end;%Tests the input parameters
    set(handleTopBox,'string',num2str(intervalHigher));%Resets the values
    set(handleLowerBox,'string',num2str(intervalLower));%Resets the values
    [~,intervalHigherPos]=min(abs(vecWave-intervalHigher));%Identify the position of the nearest value
    [~,intervalLowerPos]=min(abs(vecWave-intervalLower)); %Identify the position of the nearest value
    imgAverage = mean(data.imgAll(:,:,:,intervalLowerPos:intervalHigherPos),4);%Calculates the average intensity of the region selected
    imgAverage=reshape(imgAverage,[size(imgAverage,1),size(imgAverage,2),size(imgAverage,3),1]);
    data.imgAll = data.imgAll./repmat(imgAverage,[1,1,1,size(data.imgAll,4)]); %Divides the average intensity from the absolute intensity
    displayString = ['Normalisation from '...
        ,num2str(intervalHigher),' to ', num2str(intervalLower),' cm-1 performed.']; 
    disp(displayString);%Display operation
    set(handles.guiEntire,'Name',displayString);
    
    %Zap a spectral region to image
    case 'Zap'
    intervalHigher=str2double(get(handleTopBox,'string'));%Gets the GUI string
    intervalLower=str2double(get(handleLowerBox,'string'));%Gets the GUI string
    if isnan(intervalHigher)==1; intervalHigher= 2400;end;%Tests the input parameters
    if isnan(intervalLower)==1; intervalLower= 2200;end;%Tests the input parameters
    if intervalHigher<intervalLower
     intervalTemp = intervalHigher; intervalHigher = intervalLower;intervalLower = intervalTemp;%Swap the values
    end %if 
    vecWave = data.wavenumberVec;%Imports the wavemumber vector
    if intervalHigher > max(vecWave);intervalHigher=max(vecWave);end
    if intervalLower < min(vecWave);intervalLower=min(vecWave);end
    set(handleTopBox,'string',num2str(intervalHigher));%Resets the values
    [~,intervalHigherPos]=min(abs(vecWave-intervalHigher));%Identify the idex
    [~,intervalLowerPos]=min(abs(vecWave-intervalLower));%Identify the idex
    for idxImg =1:size(data.imgAll,1);%Loops throught the images to save memory.
        imgLower =data.imgAll(idxImg,:,:,intervalLowerPos);
        imgLower = repmat(imgLower,[1,1,(intervalHigherPos-intervalLowerPos+1)]);
        imgHigher = data.imgAll(idxImg,:,intervalHigherPos);
        imgHigher = repmat(imgHigher,[1,1,(intervalHigherPos-intervalLowerPos+1)]);
        imgSlope = (imgHigher - imgLower)./(intervalHigherPos-intervalLowerPos);%Calculates the linear slope
        imgPos=(0:(intervalHigherPos-intervalLowerPos))';
        imgPos=repmat(imgPos,[1,size(imgCurrent,1),size(imgCurrent,2)]);
        imgPos = permute(imgPos,[2,3,1]);
        imgLine = imgLower+(imgPos.*imgSlope);
        data.imgAll(idxImg,:,:,intervalLowerPos:intervalHigherPos)= 0;%Substract the ROI
        data.imgAll(idxImg,:,:,intervalLowerPos:intervalHigherPos) = imgLine;%Adds the interpolated line
    end %idxImg =1:size(data.imgAll);
    set(handleLowerBox,'string',num2str(intervalLower));%Resets the values
    displayString = ['Zap (replace by straight line) from '...
        ,num2str(intervalHigher),' to ', num2str(intervalLower),' cm-1']; %Display operation
     disp(displayString);;%Display operation
    set(handles.guiEntire,'Name',displayString);
     
    %Interpolates the spectra using a cubic spline function
    case 'Interpolation spline'
    intervalHigher=str2double(get(handleTopBox,'string'));
    vecWave = data.wavenumberVec;%Imports the wavemumber vector
    if isnan(intervalHigher)==1;intervalHigher= 5000;end;%Tests the input parameters
    if round(intervalHigher)>round(length(vecWave)*20);%Tests the input parameters
    intervalHigher = round(length(vecWave)*20);
    elseif intervalHigher<2
    intervalHigher= 2;
    else
    intervalHigher = round(intervalHigher);    
    end %if
    set(handleTopBox,'string',num2str(intervalHigher)); %Sets rounded value
    intervalNew=(vecWave(end)-vecWave(1))./(intervalHigher-1);%Resample the wave vector
    wavenumberVecNew= (vecWave(1)+intervalNew*(0:intervalHigher-1))';
    displayString = ['(',num2str(idxFile),'/',num2str(nbFiles),') Spline interpolation from '...
        ,num2str(length(vecWave)),' to ',num2str(intervalHigher),' points in progress... This might take a while.']; 
    imgTemp=nan(size(data.imgAll));%Doubles the memory allocated... bad
    for idxImg=1:size((data.imgAll),1)
        for idxSample = 1:(size(data.imgAll,2));%Loops through the rows
            for idxLine = 1:(size(data.imgAll,3));%Loops through the columns
                spectrumSpline = squeeze(data.imgAll(idxImg,idxSample,idxLine,:));
                spectrumSpline = spline(vecWave,spectrumSpline,wavenumberVecNew);
                imgTemp(idxImg,idxSample,idxLine,:)=spectrumSpline;    
            end %idxLine = 1:(size(imgCurrent,2))
        end %idxSample = 1:(size(imgCurrent,1))
    end %idxImg=1:size(size(data.imgAll),1)
        data.wavenumberVec=wavenumberVecNew;%Copy the new wavenumber vector
        data.imgAll = imgTemp;%Copy the new data
        fClear(imgTemp);%Clears the temporay array to save
      displayString = ['Spline interpolation from '...
        ,num2str(length(vecWave)),' to ',num2str(intervalHigher),' points performed.']; 
    disp(displayString);%Display operation
    set(handles.guiEntire,'Name',displayString);

    %%% ATR correction
    case 'ATR correction'
    %Calculates the scaling factor from the wavenumber vector
    %The calculation is independant of the refractive index as the value cancel
    %each other out in the calculation of the normalisation factor.
    %This correction does not correct for the abnormal dispersion of the refrative index.
    vecWave = data.wavenumberVec;%Imports the wavemumber vector
    wavenumberDpNorm = max(vecWave)./vecWave;
    wavenumberDpNorm=reshape(wavenumberDpNorm,[1,1,1,length(wavenumberDpNorm)]);
    data.imgAll=data.imgAll./repmat(wavenumberDpNorm,[size(data.imgAll,1),size(data.imgAll,2),size(data.imgAll,3),1]);%Build the 3D array
    displayString = ['Simple ATR correction performed']; 
    disp(displayString);%Display operation
    set(handles.guiEntire,'Name',displayString);
    
    %%% Spc smooth spectra
    case {'Smooth spc S-Golay','Smooth spc linear reg',}
    intervalHigher=str2double(get(handleTopBox,'string'));%Gets the string parameters
    vecWave = data.wavenumberVec;%Imports the wavemumber vector
    if isnan(intervalHigher)==1;intervalHigher=7;end; %Tests the input parameters
    if round(intervalHigher)>round(length(vecWave)/4); %Tests the input parameters
    intervalHigher = round(length(vecWave)/4);
    elseif intervalHigher<2
        intervalHigher= 2;
    else
    intervalHigher = round(intervalHigher);    
    end %if
    set(handleTopBox,'string',num2str(intervalHigher)); %Sets rounded value
    if strcmp(operationCurrent,'Smooth spc S-Golay')==1
    displayString = ['Smooth the spectra using the Savitzky-Golay filter ('...
        ,num2str(intervalHigher),' points)'];
    methodFilter = 'sgolay';
    elseif strcmp(operationCurrent,'Smooth spc linear reg')==1
    displayString = ['Smooth the spectra using a robust linear regression filter ('...
        ,num2str(intervalHigher),' points)'];        
    methodFilter = 'rlowess';
    else %Moving filter
    displayString = ['Smooth the spectra using a moving average filter ('...
        ,num2str(intervalHigher),' points)'];        %Display operation
    methodFilter = 'moving';
    end %if
    disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
    set(handles.guiEntire,'Name',displayString);
    
    for idxImg=1:size(data.imgAll,1);%Loops through the images
        for idxSample = 1:(size(data.imgAll,2));%Loops through the rows
            for idxLine = 1:(size(data.imgAll,3));%Loops through the columns
                spectrumCurrent = squeeze(data.imgAll(idxImg,idxSample,idxLine,:));
                spectrumCurrent = smooth(spectrumCurrent,intervalHigher,methodFilter);
                data.imgAll(idxImg,idxSample,idxLine,:)=spectrumCurrent ;    
            end %idxLine = 1:(size(imgCurrent,2))
        end %idxSample = 1:(size(imgCurrent,1))
    end% for idxImg=1:size((data.imgAll),1)
   
    %%%Smooth 3D
    case 'Smooth 3D'
    intervalHigher=str2double(get(handleTopBox,'string'));
   
    if isnan(intervalHigher)==1;intervalHigher=3;end; %Tests the input parameters
    if round(intervalHigher)>11;
    intervalHigher = 11;
    elseif intervalHigher<3
    intervalHigher= 3;
    else
    intervalHigher = ((round((intervalHigher-1)/2)).*2)+1;%Goes to the next odd number
    end %if
    set(handleTopBox,'string',num2str(intervalHigher)); %Sets rounded value
    %Display operation
    displayString = ['Smooth using a 3D odd box filter ('...
        ,num2str(intervalHigher),' points)'];
    disp(displayString);%Display operation
    
    for idxImg=1:size(data.imgAll,1);%Loops through the images
        data.imgAll(idxImg,:,:,:) = smooth3(data.imgAll(idxImg,:,:,:),'box',intervalHigher);% Smooths the image using a 3D gaussian filter
    end % for idxImg=1:size(data.imgAll,1);%Loops through the images
    
    case 'Pixel aggregate' %%% Aggregate pixels
    % Averages adjacent pixels together and rejects pixels with values 5 times
    % Outside the standard deviation

    %Identidy the pixels with 

    displayString = ['(',num2str(idxFile),'/',num2str(nbFiles),...
        ') Aggregates the pixels('...
        ,num2str(intervalHigher),' points)'];    
    disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
    set(handles.guiEntire,'Name',displayString);
    
    Case 'Water vapour' %%% Water vapour
    %to write
    
    %Gets the First derivative of the spectra
    case 'First derivative'
    vecWave = data.wavenumberVec;%Imports the wavemumber vector
    for idxImg=1:size(data.imgAll,1);%Loops through the images to save RAM
        imgCurrent = squeeze(data.imgAll(idxImg,:,:,:));
        imgCurrent = diff(imgCurrent,1,3);%Diffentiate once along dimension 3
        imgCurrent = imgCurrent./ mean(vecWave); %Divides the values by the data intervals
        imgCurrent = cat(3,imgCurrent,zeros(size(imgCurrent,1),size(imgCurrent,2))); %Add zeros at the end to compensate for the smaller last dimension
        data.imgAll(idxImg,:,:,:) = imgCurrent;
    end % for idxImg=1:size(data.imgAll,1);%Loops through the images
    
    %Gets the Second derivative of the spectra
    case 'Second derivative'
        vecWave = data.wavenumberVec;%Imports the wavemumber vector
    for idxImg=1:size(data.imgAll,1);%Loops through the images to save RAM
        imgCurrent = squeeze(data.imgAll(idxImg,:,:,:));
        imgCurrent = -1.*diff(imgCurrent,2,3);%Diffentiate twice along dimension 3
        imgCurrent = imgCurrent./mean(vecWave); %Divides the values by the data intervals
        imgCurrent = cat(3,imgCurrent,zeros(size(imgCurrent,1),size(imgCurrent,2))); %Add zeros at the end to compensate for the smaller last dimension
        imgCurrent = cat(3,zeros(size(imgCurrent,1),size(imgCurrent,2)),imgCurrent); %Add zeros at the beginningto compensate for the smaller last dimension
        data.imgAll(idxImg,:,:,:) = imgCurrent;
    end % for idxImg=1:size(data.imgAll,1);%Loops through the images
        
    case 'Multiple polygone interpolated offset' %%% Offset the spectra or single beam 
    %%% using an interpolated image based on designated region which should be constant
        
    %Gets the GUI value    
    idxYPos = strsplit(get(handleTopBox,'string'),';');%Imports the x cursor coordinate
    idxXPos = strsplit(get(handleLowerBox,'string'),';');%Imports the y cursor coordinate

    %%Loop through the string
    for idxRegion = 1:length(idxXPos)
        idxXPos{idxRegion}=sscanf(idxXPos{idxRegion},'%d');%Extract the X values
        idxYPos{idxRegion}=sscanf(idxYPos{idxRegion},'%d');%Extract the Y values
        if length(idxYPos{idxRegion})~=length(idxXPos{idxRegion});
                disp('X and Y coordinate do not match'); return;
        end %if
    end %for idxRegion = 1:length(idxXPos)

    %%Tests if any of the coordinates is ouside the bounderies
    for idxRegion = 1:length(idxXPos)
    if isempty(idxXPos{idxRegion})==1;idxXPos{idxRegion}=round(size(data.imgAll,2)/2);end%if
    if isempty(idxYPos{idxRegion})==1;idxYPos{idxRegion}=round(size(data.imgAll,3)/2);end%if
        idxXPos{idxRegion}(idxXPos{idxRegion}<1)=1;idxYPos{idxRegion}(idxYPos{idxRegion}<1)=1;
        idxXPos{idxRegion}(idxXPos{idxRegion}>size(data.imgAll,2))=size(data.imgAll,2);
        idxYPos{idxRegion}(idxYPos{idxRegion}>size(data.imgAll,3))=size(data.imgAll,3);
    end %For
    
    %%Calculates the ROI from given the combined polygones
    ROI = zeros(size(data.imgAll,2),size(data.imgAll,3));%Allocates the Region of Interest (ROI)
    for idxRegion = 1:length(idxXPos)
        ROI = ROI+roipoly(zeros(size(data.imgAll,2),size(data.imgAll,3)),idxXPos{idxRegion},idxYPos{idxRegion});%Ouputs a logical array emcompassing the region of interest
    end % for
    ROI(ROI~=0)=1;%Replaces non zero values by ones
    
    [coorX,coorY]= meshgrid(1:size(data.imgAll,2),1:size(data.imgAll,3));%Generates the position grid
    
    %%Aproximates missing outer frame values using weighted neigboors
    
    filterTopLine=nan((size(data.imgAll,3)),size(data.imgAll,2),size(data.imgAll,3));%Top line
    idxLine = 1;
    for idxColumn = 1:size(data.imgAll,3);
        weightMat = (coorX - idxColumn).^2 + (coorY - idxLine).^2;%Calculates the square of the distance from the pixel
        weightMat(idxLine,idxColumn)=0.0001; %Sets the same pixel weighting
        weightMat = ROI./weightMat;%Mask the non ROI pixels
        weightMat = weightMat./sum(sum(weightMat));%Normalises the valeus
        filterTopLine(idxColumn,:,:) = weightMat;%Stores the weighing factors
    end %for
    
    filterRightLine=nan((size(data.imgAll,2)),size(data.imgAll,2),size(data.imgAll,3));%Right line
    idxColumn = size(data.imgAll,2);
    for idxLine = 1:size(data.imgAll,2);
        weightMat = (coorX - idxColumn).^2 + (coorY - idxLine).^2;%Calculates the square of the distance from the pixel
        weightMat(idxLine,idxColumn)=0.0001;%Sets the same pixel weighting
        weightMat = ROI./weightMat;%Mask the non ROI pixels
        weightMat = weightMat./sum(sum(weightMat));%Normalises the valeus
        filterRightLine(idxLine,:,:) = weightMat;%Stores the weighing factors
    end % for
   
    filterBottomLine=nan((size(data.imgAll,3)),size(data.imgAll,2),size(data.imgAll,3));%Bottom line
    idxLine = size(data.imgAll,3);
    for idxColumn = 1:size(data.imgAll,3);
        weightMat = (coorX - idxColumn).^2 + (coorY - idxLine).^2;%Calculates the square of the distance from the pixel
        weightMat(idxLine,idxColumn)=0.0001; %Sets the same pixel weighting
        weightMat = ROI./weightMat;%Mask the non ROI pixels
        weightMat = weightMat./sum(sum(weightMat));%Normalises the valeus
        filterBottomLine(idxColumn,:,:) = weightMat;%Stores the weighing factors
    end %for
    
    filterLeftLine=nan((size(data.imgAll,2)),size(data.imgAll,2),size(data.imgAll,3));%Right line
    idxColumn = 1;
    for idxLine = 1:size(data.imgAll,2);
        weightMat = (coorX - idxColumn).^2 + (coorY - idxLine).^2;%Calculates the square of the distance from the pixel
        weightMat(idxLine,idxColumn)=0.0001;%Sets the same pixel weighting
        weightMat = ROI./weightMat;%Mask the non ROI pixels
        weightMat = weightMat./sum(sum(weightMat));%Normalises the valeus
        filterLeftLine(idxLine,:,:) = weightMat;%Stores the weighing factors
    end % for

    ROI(ROI==0)=nan; %Replace 0 by nan
    ROIout=ROI;ROIout(1,1:end)=1;ROIout(1:end,end)=1;ROIout(end,1:end)=1;ROIout(1:end,1)=1; %Adds the outline to the ROI
    [imgXcoor,imgYcoor,imgZcoor]= meshgrid(1:size(data.imgAll,2),1:size(data.imgAll,3),1:size(data.imgAll,4));%Generates the position grid
    
    imgXcoorROI = imgXcoor.*repmat(ROIout,[1 1 size(data.imgAll,4)]);%Applies the ROI filter
    imgXcoorROI = reshape(imgXcoorROI,[],1);
    imgXcoorROI(isnan(imgXcoorROI)==1)=[];%Removes the nan values
    imgXcoor = reshape(imgXcoor,[],1);

    imgYcoorROI = imgYcoor.*repmat(ROIout,[1 1 size(data.imgAll,4)]);%Applies the ROI filter
    imgYcoorROI = reshape(imgYcoorROI,[],1);
    imgYcoorROI(isnan(imgYcoorROI)==1)=[];%Removes the nan values
    imgYcoor = reshape(imgYcoor,[],1);

    imgZcoorROI = imgZcoor.*repmat(ROIout,[1 1 size(data.imgAll,4)]);%Applies the ROI filter
    imgZcoorROI = reshape(imgZcoorROI,[],1);
    imgZcoorROI(isnan(imgZcoorROI)==1)=[];%Removes the nan values
    imgZcoor = reshape(imgZcoor,[],1);
     
    for idxImg=1:size(data.imgAll,1);%Loops through the image stack
        displayString =['(',num2str(idxImg),'/',num2str(size(data.imgAll,1)),') Image interpolated offset'];
        disp(displayString);set(handles.guiEntire,'Name',displayString);
        set(handles.guiEntire,'Name',displayString);drawnow;
        
        imgValue =(data.imgAll(idxImg,:,:,:));%Gets the current image
        imgValue = reshape(imgValue,[size(imgValue,2),size(imgValue,3),size(imgValue,4)]);%Removes the singleton dimension
        imgValue = imgValue.*repmat(ROI,[1 1 size(data.imgAll,4)]);%Apply the ROI filter
        
        %Applies the top filter interpolation
        idxLine = 1;
        for idxColumn = 1:size(filterTopLine,1)
            filterCur = squeeze(filterTopLine(idxColumn,:,:));%Selects the current filter for the pixel
            specCur = nanmean(nanmean(imgValue.*repmat(filterCur,[1,1,size(imgValue,3)]),1),2);
            imgValue(idxLine,idxColumn,:)=specCur;%Copies the interpolated spectrum
        end %for

        %Applies the right filter interpolation
        idxColumn = size(data.imgAll,2);
        for idxLine = 1:size(filterRightLine,1)
            filterCur = squeeze(filterRightLine(idxLine,:,:));%Selects the current filter for the pixel
            specCur = nanmean(nanmean(imgValue.*repmat(filterCur,[1,1,size(imgValue,3)]),1),2);
            imgValue(idxLine,idxColumn,:)=specCur;%Copies the interpolated spectrum
        end %for
        
        %Applies the bottom filter interpolation
        idxLine = size(data.imgAll,3);
        for idxColumn = 1:size(filterBottomLine,1)
            filterCur = squeeze(filterBottomLine(idxColumn,:,:));%Selects the current filter for the pixel
            specCur = nanmean(nanmean(imgValue.*repmat(filterCur,[1,1,size(imgValue,3)]),1),2);
            imgValue(idxLine,idxColumn,:)=specCur;%Copies the interpolated spectrum
        end %for

        %Applies the right filter interpolation
        idxColumn = 1;
        for idxLine = 1:size(filterLeftLine,1)
            filterCur = squeeze(filterLeftLine(idxLine,:,:));%Selects the current filter for the pixel
            specCur = nanmean(nanmean(imgValue.*repmat(filterCur,[1,1,size(imgValue,3)]),1),2);
            imgValue(idxLine,idxColumn,:)=specCur;%Copies the interpolated spectrum
        end %for
        
        imgValue = reshape(imgValue,[],1);%Reshapes the 3D array to a vector
        imgValue(isnan(imgValue)==1)=[];%Removes the nan values

        imgInterpolated = griddata(imgXcoorROI,imgYcoorROI,imgZcoorROI,imgValue,imgXcoor,imgYcoor,imgZcoor);%Interpolates the values outside the ROI   
        imgInterpolated=reshape(imgInterpolated,[size(data.imgAll,2),size(data.imgAll,3),size(data.imgAll,4)]);%Reshapes the vector back to a 3D matrix
        imgInterpolated(isnan(imgInterpolated)==1)=0;%Replaces the nan by 0s
        data.imgAll(idxImg,:,:,:)= squeeze(data.imgAll(idxImg,:,:,:)) - imgInterpolated;%Substracts the offset & overwrites the offseted image

    end%for idxImg=1:size(data.imgAll,1);%Loops through the image stack
         
        ROI=reshape(ROI,[size(data.imgAll,2),size(data.imgAll,3),1]);%Reshapes the ROI for expension
        ROI=repmat(ROI,[1,1,size(data.imgAll,4)]);%Extends the ROI along the 4th dimension
        
        for idxImg = 1:size(data.imgAll,1) %Loops throught the images to save memory
            currentImage=reshape(data.imgAll(idxImg,:,:,:),[size(data.imgAll,2),size(data.imgAll,3),size(data.imgAll,4)]); %Averages the region's spectrum by multiplying storage array by the ROI mask (ignores the NaN)
            averageSpectrum = squeeze(nanmean(nanmean((real(currentImage)).*ROI)));
            stackIntegratedRegions(idxImg,idxRegion,:)=averageSpectrum;    
        end %for idxImg = 1:size(data.imgAll)
    
            ROI(ROI==0) = nan;%Replaces all the zero values with NaN
        if sum(sum(ROI))==0%If no pixel within the ROI
            ROI(idxXPos{idxRegion},idxXPos{idxRegion})=1;%Replaces the ROI position index by 1  
        end %sum(sum(ROI))==0
        
        case  'Multiple polygone frame offset'%Offsets the image using the average of the selected areas
            
            %Gets the GUI value    
    idxYPos = strsplit(get(handleTopBox,'string'),';');%Imports the x cursor coordinate
    idxXPos = strsplit(get(handleLowerBox,'string'),';');%Imports the y cursor coordinate

    %%Loop through the string
    for idxRegion = 1:length(idxXPos)
        idxXPos{idxRegion}=sscanf(idxXPos{idxRegion},'%d');%Extract the X values
        idxYPos{idxRegion}=sscanf(idxYPos{idxRegion},'%d');%Extract the Y values
        if length(idxYPos{idxRegion})~=length(idxXPos{idxRegion});
                disp('X and Y coordinate do not match'); return;
        end %if
    end %for idxRegion = 1:length(idxXPos)

    %%Tests if any of the coordinates is ouside the bounderies
    for idxRegion = 1:length(idxXPos)
    if isempty(idxXPos{idxRegion})==1;idxXPos{idxRegion}=round(size(data.imgAll,2)/2);end%if
    if isempty(idxYPos{idxRegion})==1;idxYPos{idxRegion}=round(size(data.imgAll,3)/2);end%if
        idxXPos{idxRegion}(idxXPos{idxRegion}<1)=1;idxYPos{idxRegion}(idxYPos{idxRegion}<1)=1;
        idxXPos{idxRegion}(idxXPos{idxRegion}>size(data.imgAll,2))=size(data.imgAll,2);
        idxYPos{idxRegion}(idxYPos{idxRegion}>size(data.imgAll,3))=size(data.imgAll,3);
    end %For
    
    %%Calculates the ROI from given the combined polygones
    ROI = zeros(size(data.imgAll,2),size(data.imgAll,3));%Allocates the Region of Interest (ROI)
    for idxRegion = 1:length(idxXPos)
        ROI = ROI+roipoly(zeros(size(data.imgAll,2),size(data.imgAll,3)),idxXPos{idxRegion},idxYPos{idxRegion});%Ouputs a logical array emcompassing the region of interest
    end % for
    ROI(ROI~=0)=1;%Replaces non zero values by ones
    ROI(isnan(ROI)==1)=nan;%Replaces all the zero values by nan
     for idxImg=1:size(data.imgAll,1);%Loops through the image stack
        imgCur = data.imgAll(idxImg,:,:,:);%Extract the image
        imgCur=reshape(imgCur,[size(imgCur,2),size(imgCur,3),size(imgCur,4)]);%Removes the first dimension
        meanSpectrum=imgCur.*repmat(ROI,[1,1,size(data.imgAll,4)]);%Applies the ROI mask
        meanSpectrum=squeeze(nanmean(nanmean(meanSpectrum,1),2));%Calculates the average spectrum of the selectect region
        meanSpectrum=reshape(meanSpectrum,[1,1,length(meanSpectrum)]);%Reshapes the array
        imgCur = imgCur-repmat(meanSpectrum,[size(data.imgAll,2),size(data.imgAll,3),1]);%Substract the average spectrum from all the spectra in the image
        data.imgAll(idxImg,:,:,:) = imgCur;%Saves the offseted image to the storage array
     end %for %Loops throught the images
     displayString ='Image frame offset performed on the selected re';
     disp(displayString);set(handles.guiEntire,'Name',displayString);
             
    otherwise
        %Do nothing
end %Switch operationCurrent
end %idxOp = 1:5 

integrateImageRegions(handles) %Integrates the images and refreshes the graphs

function enableGui(handles)

set(handles.menuSave,'enable','on');
set(handles.menuExport,'enable','on');
set(handles.menuProcessing,'enable','on');
set(handles.popupmenuImgSelectionMode,'enable','on');
set(handles.pushbuttonNextSelection,'enable','on');
set(handles.popupmenuSpectrumIntegrationMode,'enable','on');
set(handles.pushbuttonResetRegion,'enable','on');
set(handles.editSpectrumIntervalHigh,'enable','on');
set(handles.popupmenuImageType,'enable','on');
set(handles.uipanelVisualisationLow,'visible','on') %Activates the panel
set(handles.uipanelVisualisation,'visible','on') %Activates the panel
set(handles.uipanelSpectra,'visible','on') %Activates the panel
set(handles.uipanelImage,'visible','on') %Activates the panel
set(handles.uipanelStack,'visible','on') %Activates the panel
set(handles.uipanelHistogram,'visible','on') %Activates the panel


function integrateImageRegions(handles)

global data
tic %Starts the timer
displayString = 'Integrating the images, please wait...';
disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;

waveNb = data.wavenumberVec;%Imports the wavenumber vector

%% Integrates the images using the selected method to reduce the wavenumber dimension to a single variable
try
    spectrumIntegrationType=data.spcInt.type;
    catch
    spectrumIntegrationType = 'Single wavelenght';
    data.spcInt.type = spectrumIntegrationType;
end %Try

switch spectrumIntegrationType

case 'Single peak' % Integrates all the images of the stack using a single peak
    waveHigh=data.spcInt.values(2);%Tests the GUI's paramenters
    waveLow=data.spcInt.values(1);%Tests the GUI's paramenters 
    [~,idxLow]=min(abs(waveNb-waveLow)); %Finds the Low wavenumber index of the interval
    [~,idxHigh]=min(abs(waveNb-waveHigh));%Finds the high wavenumber index of the interval
    imgMean = mean(data.imgAll(:,:,:,idxLow:idxHigh),4);%Gets the means
    imgLowVal=data.imgAll(:,:,:,idxLow);%Gets the low intervals values
    imgHighVal=data.imgAll(:,:,:,idxHigh);%Gets the  high intervals values
    imgHighLowPos = cat(4,imgHighVal,imgLowVal);%Concatenate the fourth dimension
    imgIntervalMin = min(imgHighLowPos,[],4);%Finds the inimum between bouderies
    imgIntervalMax = max(imgHighLowPos,[],4);%Finds the maximum between bouderies
    stackIntImg= imgMean-imgIntervalMin-((imgIntervalMax-imgIntervalMin)/2);%Substracts the min value and half of the diffence between the 2 (linear baseline) 
    data.imgInt = stackIntImg;%Stores the integrated image values

    case 'Peak ratio'

    disp(' case peak ratio Not yet written');
        
    case 'PC1'
    disp(' case PC1 ratio Not yet written');
    
    case 'PC2'
    disp(' case PC2 ratio Not yet written');               
        
    case 'PC3'
    disp(' case PC3 ratio Not yet written');
                
    case 'COG'  %Integrates all the images of the stacks using the COG
    
    waveHigh=data.spcInt.values(2);%Gets the GUI high interval
    waveLow=data.spcInt.values(1);%Gets the GUI low interval
    if waveLow>waveHigh; %Swap values if inverted
     idxTemp = waveHigh; waveHigh = waveLow;waveLow = idxTemp;
    end %if
    [~,idxLow]=min(abs(waveNb-waveLow));waveLow = waveNb(idxLow); %Finds the low wavenumber index of the interval
    [~,idxHigh]=min(abs(waveNb-waveHigh));waveHigh = waveNb(idxHigh);%Finds the high wavenumber index of the interval
    waveArr =reshape(waveNb(idxLow:idxHigh),[1 1 1 size(waveNb(idxLow:idxHigh))]); %Put the wavenumber in the 4rd dimension
    waveArr = repmat(waveArr,size(data.imgAll,1),size(data.imgAll,2),size(data.imgAll,3));... %Repmat the wavevector
    imgCOG = sum(data.imgAll(:,:,:,idxLow:idxHigh).*waveArr,4);%Gets the sum product
    imgCOG = imgCOG./sum(data.imgAll(:,:,:,idxLow:idxHigh),4); %Calculates the COG
    imgCOG(imgCOG<waveLow) = nan;%Removes the values below the minnimum limit
    imgCOG(imgCOG>waveHigh) = nan;%Removes the values above the maximum limit      
    data.imgInt = imgCOG;%Stores the integrated values
           
otherwise %Integrates all the images of the stacks using a single wavelenght 
    idxWaveNb=roundn(str2double(get(handles.editSpectrumSingleSelect,'string')),-1);%Gets and tests the GUI's paramenters
    if isnan(idxWaveNb)==1 || idxWaveNb<min(waveNb)||idxWaveNb>max(waveNb)
        idxWaveNb = waveNb(round(length(waveNb)/2));
    end %if
    data.spcInt.values(1) = idxWaveNb;%Stores the index
    [~,idxWaveValue]=min(abs(waveNb-idxWaveNb));%Finds the index of the selected wavenumber
    stackIntImg = data.imgAll(:,:,:,idxWaveValue);%Extract the single value
    data.imgInt = stackIntImg;%Stores the integrated image values

end %Switch cases

%% Extract the regions of interest (ROI)

%Selects the image integration modes
idxXPos = strsplit(get(handles.editImageXpos,'string'),';');%Imports the x cursor coordinate
idxYPos = strsplit(get(handles.editImageYpos,'string'),';');%Imports the y cursor coordinate

%Loop through the string
for idxRegion = 1:length(idxXPos)
    idxXPos{idxRegion}=sscanf(idxXPos{idxRegion},'%d');%Extract the X values
    idxYPos{idxRegion}=sscanf(idxYPos{idxRegion},'%d');%Extract the Y values
    if length(idxYPos{idxRegion})~=length(idxXPos{idxRegion});
            disp('X and Y coordinate do not match'); return;
    end %if
end %for idxRegion = 1:length(idxXPos)

%Checks if any of the coordinates is ouside the bounderies
for idxRegion = 1:length(idxXPos)
if isempty(idxXPos{idxRegion})==1;idxXPos{idxRegion}=round(size(data.imgAll,2)/2);end%if
if isempty(idxYPos{idxRegion})==1;idxYPos{idxRegion}=round(size(data.imgAll,3)/2);end%if
    idxXPos{idxRegion}(idxXPos{idxRegion}<1)=1;
    idxYPos{idxRegion}(idxYPos{idxRegion}<1)=1;
    idxXPos{idxRegion}(idxXPos{idxRegion}>size(data.imgAll,2))=size(data.imgAll,2);
    idxYPos{idxRegion}(idxYPos{idxRegion}>size(data.imgAll,3))=size(data.imgAll,3);
end %For

data.idxXPos = idxXPos; %Export the position index
data.idxYPos = idxYPos; %Export the position index

nbRegions = length(idxXPos); %Detemines the number of regions
data.nbPixelSelected = nan(nbRegions,1);%Preallocates the number of regions
data.stackRegionAverageSpc=nan(size(data.imgAll,1),nbRegions,size(data.imgAll,4));%Prealocates the region averages spectra
data.stackRegionAverageInt= nan(size(data.imgAll,1),nbRegions);%Prealocates the region averaged integrated value
data.stackRegionsSpectra= cell(size(data.imgAll,1),nbRegions);%Prealocates the all region spectra
data.stackRegionsIntegrated= cell(size(data.imgAll,1),nbRegions);%Prealocates the all region integrated values
for idxRegion = 1:nbRegions %Loops throught the regions
    ROI2D = roipoly(zeros(size(data.imgAll,2),size(data.imgAll,3)),idxXPos{idxRegion},idxYPos{idxRegion});%Ouputs a logical array emcompassing the region of interest
    if sum(sum(ROI2D))==0%If no pixel within the ROI
    ROI2D(idxXPos{idxRegion},idxXPos{idxRegion})=1;%Replaces the ROI position index by 1  
    end %sum(sum(ROI))==0
    ROI2D = ROI2D.*ones(size(data.imgAll,2),size(data.imgAll,3));%Put one within the ROI
    ROI2D(ROI2D==0) = nan;%Replaces all the zero values with NaN
    data.nbPixelSelected(idxRegion) = nansum(nansum(ROI2D));%Calculates the number of selected pixels
    ROI3D=reshape(ROI2D,[size(data.imgAll,2),size(data.imgAll,3),1]);%Reshapes the ROI for repmat
    ROI3D=repmat(ROI3D,[1,1,size(data.imgAll,4)]);%Extends the ROI along the 3rd dimension
    
    for idxImg = 1:size(data.imgAll,1) %Loops throught the images to save memory
        %Calculates the average spectrm for each region of the 3D image
        currentImage=reshape(data.imgAll(idxImg,:,:,:),[size(data.imgAll,2),size(data.imgAll,3),size(data.imgAll,4)]);%Extracts the current image
        curRegionAverageSectrum = squeeze(nanmean(nanmean((real(currentImage)).*ROI3D))); %Averages the region's spectrum by multiplying storage array by the ROI mask (ignores the NaN)
        data.stackRegionAverageSpc(idxImg,idxRegion,:)=curRegionAverageSectrum;%Saves the average spectrum to the array    
        %Selects the region spectra of the 3D image
        currentImage = (real(currentImage)).*ROI3D;%Applies the 3D filter
        currentImage = reshape(currentImage,[],size(currentImage,3)); %Reshapes the array to 2D
        currentImage(any(isnan(currentImage),2),:) = []; %Removes all the Nan values and shrinks the array
        data.stackRegionsSpectra{idxImg,idxRegion} = currentImage;%Saves the Selected ROI to the array     
        %Selects the regions of the integrated image (2D)
        curIntData=reshape(data.imgInt(idxImg,:,:),[size(data.imgAll,2),size(data.imgAll,3)]);%Extracts the current image
        curRegionInt = (real(curIntData)).*ROI2D; %Applies the filter
        curRegionInt = reshape(curRegionInt,[],size(curRegionInt,3)); %Reshapes the array to 2D
        curRegionInt(any(isnan(curRegionInt),2),:) = []; %Removes all the Nan values and shrinks the array
        data.stackRegionsIntegrated{idxImg,idxRegion} = curRegionInt;%Saves the Selected ROI to the array
        data.stackRegionAverageInt(idxImg,idxRegion) = nanmean(curRegionInt);%Saves the average value to the storage array
     end %for idxImg = 1:size(data.imgAll)
 end % for
 
%% Extracts the histograms of the first ROI

% pixelNb=length(data.stackRegionsIntegrated{1,1});%Gets the number of pixels in the first bin
% binNb = abs(round(str2double(get(handles.editBin,'String'))));%Gets the number of bins
% if binNb == 0||isnan(binNb)==1; binNb = round(pixelNb/10);end%Tests the number of bins
% set(handles.editBin,'String',num2str(binNb));%Resets the number of bins
%    
% idxRegion = 1;
% curInt=cell2mat(data.stackRegionsIntegrated(1:end,idxRegion)')';%Extract all the integrated values of the first region
% intMax = max(max(curInt));%Gets the maximum value of the matrix
% intMin = min(min(curInt));%Gets the minimum value of the matrix
% binRange = intMin:(intMax - intMin)/(binNb-1+1):intMax;%Calculates the bin range
% if isempty(binRange);binRange=[1,2];end
% data.stackRegionHistBin = binRange; %Saves the bin range in the global storage array
% binRange(end) = intMax*1.001;%Offest the last bin to include the maximum
% binCount = histc(curInt,binRange,2);%Calculates the histogram for the given interval
% data.stackRegionHist = binCount(1:end,1:end-1); %Saves the data in the global storage array
    
%% Extracts the slice of the first ROI

%Finds the longest vector

%Transform the array into a XYZ coordinate system
%Rotates the arrays along the longest vector
%Loops through the bins
%Avereges Segments the x axis by the number of bins

%Saves the data in the global storage array

%% Fits the image stack kinetic using a sigmoidal function
boxChecked = get(handles.checkboxSigFit,'Value');
if boxChecked==1
    displayString = 'Fitting the curves, please wait...';
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    try
        timeVec = data.timeRel;%Gets the time vector
    catch
        timeVec = (1:size(data.imgAll,1))';   
    end%try
    
    report=cell(17,nbRegions+1);%Generates the emmpty report cell array
    report{1,1}='Filename'; %Writes the header
    report{1,2}='Region number'; %Writes the header
    report{1,3}='A0';report{1,4}='A0(min)';report{1,5}='A0(max)';
    report{1,6}='Aoff';report{1,7}='Aoff(min)';report{1,8}='Aoff(max)';
    report{1,9}='K1';report{1,10}='K1(min)';report{1,11}='K1(max)';
    report{1,12}='K2';report{1,13}='K2(min)';report{1,14}='K2(max)';
    report{1,15}='t half';report{1,16}='t half(min)';report{1,17}='t half(max)';
    report{1,18}='R2';
    
    for idxRegion = 1:nbRegions%Loops throught the regions
        displayString = ['(', num2str(idxRegion),'/',num2str(nbRegions),') region fitting'];
        disp(displayString)
        
        valueVec = data.stackRegionAverageInt(:,idxRegion);

        %Defines the fitting function
        fitModel = fittype('Aoff+(A0-(((K1./K2)+A0)./(1+((K1./(K2*A0))*exp(x.*(K1+(K2*A0)))))))',...
            'coefficients',{'A0','Aoff','K1','K2'},'independent','x');
        guiNbIt = 50;
        %Defines the fitting options
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,min(valueVec),1E-12,1E-12],...
            'Upper',[max(valueVec),max(valueVec),100000.*max(timeVec),100000.*max(timeVec)],...
            'Startpoint',[max(valueVec)-min(valueVec),min(valueVec),10./max(timeVec),10000./max(timeVec)],...
            'Robust', 'LAR',... %or 'On' select 'LAR'
            'Algorithm','Trust-Region',... %
            'MaxIter', guiNbIt,...
            'DiffMinChange',1e-15.*max(valueVec),...
            'DiffMaxChange',1e15.*max(valueVec),...
            'Display','Off'); %or 'Notify'; 
        
        try %Fits the x and y
            [fittedValues,gof,fitPar] = fit(timeVec,valueVec,fitModel,fitopt);
            arrCoeffs = coeffvalues(fittedValues);%Gets the FIT coefficient
            arrCoeffsInterval = confint(fittedValues,0.95);%95 Confidence interval
            [~,name,~] = fileparts(data.SmpFileList{1});%Gets thje filename
            report{idxRegion+1,1}= name;
            report{idxRegion+1,2}=idxRegion;
            A0= arrCoeffs(1);report{idxRegion+1,3}=A0;
            A0_min= arrCoeffsInterval(1,1);report{idxRegion+1,4}=A0_min;
            A0_max= arrCoeffsInterval(2,1);report{idxRegion+1,5}=A0_max;
            Aoff= arrCoeffs(2);report{idxRegion+1,6}=Aoff;
            Aoff_min= arrCoeffsInterval(1,2);report{idxRegion+1,7}=Aoff_min;
            Aoff_max= arrCoeffsInterval(2,2);report{idxRegion+1,8}=Aoff_max;
            K1= arrCoeffs(3);report{idxRegion+1,9}=K1;
            K1_min= arrCoeffsInterval(1,3);report{idxRegion+1,10}=K1_min;
            K1_max= arrCoeffsInterval(2,3);report{idxRegion+1,11}=K1_max;
            K2= arrCoeffs(4);report{idxRegion+1,12}=K2;
            K2_min= arrCoeffsInterval(1,4);report{idxRegion+1,13}=K2_min;
            K2_max= arrCoeffsInterval(2,4);report{idxRegion+1,14}=K2_max;
            tHalf= (log((K2*A0)/K1))/(K1+(K2*A0));%Calculates the t1/2
            
            report{idxRegion+1,15}= tHalf;%Saves the t half
            tHalf_min= (log((K2_min*A0_min)/K1_max))/(K1_max+(K2_max*A0_max));%Calculates the t1/2
            report{idxRegion+1,16}= tHalf_min;%Saves the t half
            tHalf_max= (log((K2_max*A0_max)/K1_min))/(K1_min+(K2_min*A0_min));%Calculates the t1/2
            report{idxRegion+1,17}= tHalf_max;%Saves the t half
            
            timeRelInter = (min(timeVec):((max(timeVec)-min(timeVec))/499):max(timeVec))';%Interpolates the time vector
            data.timeRelInter = timeRelInter;
            data.fittedInterValues{idxRegion} = Aoff+(A0-(((K1./K2)+A0)./(1+((K1./(K2*A0))*exp(timeRelInter.*(K1+(K2*A0)))))));
            data.fittedValues{idxRegion} = Aoff+(A0-(((K1./K2)+A0)./(1+((K1./(K2*A0))*exp(timeVec.*(K1+(K2*A0)))))));

            report{idxRegion+1,18}=gof.rsquare;%Saves the R square coefficient
           
        catch lasterr
            displayString = ['Fit error! ' ,lasterr.message];
            disp(displayString);
            set(handles.checkboxSigFit,'Value',0);   
        end %Catch

    end % For through the files

data.fitReport = report;
    
end %if boxChecked==1

timeElapsed = toc; %Ends the timer

displayString = ['Integration completed in '...
    ,num2str(round(timeElapsed)),' seconds'];
disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
refreshGraph(handles);%Refreshes the graphs


function [firstVariable,scondVariable]=swapValues(firstVariable,scondVariable) %Swap the variable values, used for min max variables
temp = firstVariable;
firstVariable = scondVariable;
scondVariable = temp;

function refreshGraph(handles)
%Refreshes the graphs

global data
global hAll

%Imports the data from the global storage array
try
    stackRegionAverageSpc =data.stackRegionAverageSpc;
catch
    stackRegionAverageSpc = data.stackIntegratedRegions; 
end %try
try
    stackRegionAverageInt = data.stackRegionAverageInt;
catch
    stackRegionAverageInt = data.stackIntRegionsSpectra;
end
try
stackIntImg = data.imgInt;
catch
stackIntImg = data.stackIntImg;  
end
spectrumIntegrationType = data.spcInt.type;

%Calculates the colour map
   rCMap(1,:)= 0.01+(0.99.*(1-abs(cos((1:100-1)./1.5))));
   rCMap(2,:)= 0.01+(0.99.*(1-abs(cos((1:100-1)./2.8))));
   rCMap(3,:)= 0.01+(0.99.*(1-abs(cos((1:100-1)./3.5))));

%Creates an anonymus function to clear variable and free memory
fClear = @(x) clear(inputname(1));

hAll = handles;

set(handles.uipanelVisualisation, 'Visible', 'off');%Disables the GUI whilst it refreshes the graphs
set(handles.uipanelVisualisationLow, 'Visible', 'off');%Disables the GUI whilst it refreshes the graphs
drawnow;

waveNb = data.wavenumberVec;%Imports the wavenumner vector
nbFiles =  size(data.imgAll,1);%Gets the number of images

%Gets and tests the parameters from the GUI
if size(data.imgAll,1) > 1
set(handles.sliderFrame,'Min', 1);
set(handles.sliderFrame,'Max', nbFiles);
set(handles.sliderFrame,'SliderStep',[(1/((nbFiles)-1)),(1/(nbFiles-1))]);
set(handles.sliderFrame,'Enable', 'on'); %Activates the slider
else  
set(handles.sliderFrame,'Enable', 'off'); %Deactivates the slider since only one frame
end

idxImg = floor(get(handles.sliderFrame,'value'));% Gets the selected image
if idxImg > size(data.imgAll,1) || idxImg <1%Tests the frame index
idxImg =1;
end
set(handles.editFrameSelect,'string',num2str(idxImg));
set(handles.sliderFrame,'value',idxImg);


%% Plots the region average spectra (First axis)

spectrumResponseMax=str2double(get(handles.editSpectrumResponseMax,'string'));%Gets the GUI parameter
if isnan(spectrumResponseMax)==1;spectrumResponseMax = max(max(stackRegionAverageSpc));end%Tests GUI parameter
spectrumResponseMin=str2double(get(handles.editSpectrumResponseMin,'string'));%Gets the GUI parameter
if isnan(spectrumResponseMin)==1;spectrumResponseMin = min(min(stackRegionAverageSpc));end%Tests GUI parameter
if spectrumResponseMin>spectrumResponseMax;%Tests GUI parameter
[spectrumResponseMin,spectrumResponseMax]=swapValues(spectrumResponseMin,spectrumResponseMax);    
elseif spectrumResponseMin==spectrumResponseMax
spectrumResponseMin = min(min(min(min(data.imgAll))));
spectrumResponseMax = max(max(max(max(data.imgAll))));       
end %if
set(handles.editSpectrumResponseMax,'string',num2str(spectrumResponseMax));%Sets the tested GUI paremeter
set(handles.editSpectrumResponseMin,'string',num2str(spectrumResponseMin));%Sets the tested GUI paremeter

spectrumWaveMax = str2double(get(handles.editSpectrumWaveMax,'string'));
if isnan(spectrumWaveMax)==1;spectrumWaveMax = max(waveNb);end
spectrumWaveMin = str2double(get(handles.editSpectrumWaveMin,'string'));
if isnan(spectrumWaveMin)==1;spectrumWaveMin = min(waveNb);end
if spectrumWaveMin>spectrumWaveMax;
[spectrumWaveMin,spectrumWaveMax]=swapValues(spectrumWaveMin,spectrumWaveMax);    
elseif spectrumWaveMin==spectrumWaveMax
spectrumWaveMax = max(waveNb);spectrumWaveMin = min(waveNb);  
end %if
set(handles.editSpectrumWaveMax,'string',num2str(spectrumWaveMax));%Sets the tested GUI paremeter
set(handles.editSpectrumWaveMin,'string',num2str(spectrumWaveMin));%Sets the tested GUI paremeter

hold(handles.axesSpectrum,'off');
for idxR = 1:size(stackRegionAverageSpc,2);
    vecCur = squeeze(stackRegionAverageSpc(idxImg,idxR,:));
    plot(handles.axesSpectrum,waveNb,vecCur,'LineWidth',2,'LineStyle','-',...
    'Color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);
    set(handles.axesSpectrum,'XDir','reverse','xlim',[spectrumWaveMin spectrumWaveMax])
    set(handles.axesSpectrum,'ylim',[spectrumResponseMin spectrumResponseMax])
    xlabel(handles.axesSpectrum,'Wavenumber (cm^{-1})');
    try
    reponseUnit= data.responseUnit;
    catch
    reponseUnit = 'Absorbance';
    end %try
    ylabel(handles.axesSpectrum,reponseUnit);
    %'YTickLabel',{'-1.0','-0.5','.0','0.5','1.0','1.50'})

    hold(handles.axesSpectrum,'on')
    
    switch spectrumIntegrationType
    case 'Single peak'
    waveLow=data.spcInt.values(1);
    waveHigh=data.spcInt.values(2);
    [~,idxLow] =min(abs(waveNb-waveLow));
    [~,idxHigh] =min(abs(waveNb-waveHigh));
    
    plot(handles.axesSpectrum,[waveLow,waveHigh],...
        [stackRegionAverageSpc(idxImg,idxR,idxLow),stackRegionAverageSpc(idxImg,idxR,idxHigh)],...
    'LineWidth',2,'Color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);        

    case 'Peak ratio'
    case 'PC1'
    case 'PC2'
    case 'PC3'
        
    case 'COG' % Plot COG cursor
    waveLow=data.spcInt.values(1);
    waveHigh=data.spcInt.values(2);
    [~,idxLow] =min(abs(waveNb-waveLow));
    [~,idxHigh] =min(abs(waveNb-waveHigh));
    
    plot(handles.axesSpectrum,[waveLow,waveHigh],...
        [stackRegionAverageSpc(idxImg,idxR,idxLow),stackRegionAverageSpc(idxImg,idxR,idxHigh)],...
    'LineWidth',2,'Color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);      
    
    wavenumberValue=stackRegionAverageSpc(idxImg,idxR);
    [~,idxWaveNb] =min(abs(waveNb-wavenumberValue));
    valY = stackRegionAverageSpc(idxImg,idxR,idxWaveNb);
    plot(handles.axesSpectrum,wavenumberValue,valY,'MarkerSize',15,'Marker','s','LineWidth',2,'Color',[0.99 0.99 0.99]);
    plot(handles.axesSpectrum,wavenumberValue,valY,'MarkerSize',13,'Marker','s','LineWidth',2,'Color',[0.51 0.51 0.51]);
    plot(handles.axesSpectrum,wavenumberValue,valY,'MarkerSize',11,'Marker','s','LineWidth',2,'Color',[0.01 0.01 0.01]);       
  
    otherwise %Single point integration
    %Plot cursor
    wavenumberValue = data.spcInt.values(1);
    [~,idxWaveNb] =min(abs(waveNb-wavenumberValue));
    plot(handles.axesSpectrum,wavenumberValue,stackRegionAverageSpc(idxImg,idxR,idxWaveNb),...
    'MarkerSize',15,'Marker','s','LineWidth',2,'Color',[0.99 0.99 0.99]);
    plot(handles.axesSpectrum,wavenumberValue,stackRegionAverageSpc(idxImg,idxR,idxWaveNb),...
    'MarkerSize',13,'Marker','s','LineWidth',2,'Color',[0.51 0.51 0.51]);
    plot(handles.axesSpectrum,wavenumberValue,stackRegionAverageSpc(idxImg,idxR,idxWaveNb),...
    'MarkerSize',11,'Marker','s','LineWidth',2,'Color',[0.01 0.01 0.01]);

    end %Switch between case

end %for idxRegions = 1:size(spectrumRegionAverage,2);

%% Plots the selected integrated image (Second axis)

imgXMax = str2double(get(handles.editImageXMax,'string'));%Gets the GUI parameter
imgYMax = str2double(get(handles.editImageYMax,'string'));%Gets the GUI parameter
imgYMin = str2double(get(handles.editImageYMin,'string'));%Gets the GUI parameter
imgXMin = str2double(get(handles.editImageXMin,'string'));%Gets the GUI parameter
imgZMax = str2double(get(handles.editImageZMax,'string'));%Gets the GUI parameter
imgZMin = str2double(get(handles.editImageZMin,'string'));%Gets the GUI parameter

if isnan(imgXMax)==1;imgXMax = size(data.imgAll,2);end%Tests the GUI parameter
if isnan(imgXMin)==1;imgXMin = 1;end 
if imgXMin>imgXMax;%Tests the GUI parameter
[imgXMin,imgXMax]=swapValues(imgXMin,imgXMax);
elseif imgXMin==imgXMax
imgXMax = size(data.imgAll,2);imgXMin = 1;   
end %if

if isnan(imgYMax)==1;imgYMax = size(data.imgAll,3);end%Tests the GUI parameter
if isnan(imgYMin)==1;imgYMin = 1;end
if imgYMin>imgYMax;%Tests the GUI parameter
[imgYMin,imgYMax]=swapValues(imgYMin,imgYMax);
elseif imgYMin==imgYMax
imgYMax = size(data.imgAll,3);imgYMin = 1;   
end %if

if isnan(imgZMax)==1;imgZMax =max(max(max(stackIntImg(idxImg,:,:,:))));end%Tests the GUI parameter
if isnan(imgZMin)==1;imgZMin = min(min(min(stackIntImg(idxImg,:,:,:))));end
if imgZMin>imgZMax;%Tests the GUI parameter
[imgZMin,imgZMax]=swapValues(imgZMin,imgZMax);
elseif imgZMin==imgZMax
imgZMax = max(max(max(max(stackIntImg(idxImg,:,:,:)))));imgZMin = min(min(min(min(stackIntImg(idxImg,:,:,:)))));   
end %if

set(handles.editImageYMax,'string',num2str(imgYMax));%Resets the GUI parameter
set(handles.editImageYMin,'string',num2str(imgYMin));%Resets the GUI parameter
set(handles.editImageXMax,'string',num2str(imgXMax));%Resets the GUI parameter
set(handles.editImageXMin,'string',num2str(imgXMin));%Resets the GUI parameter

%Loads and sets the color map
cMap= [0 0 0;0.036363 0 0.036363635212183;0.072727270424366 0 0.072727270424366;0.109090909361839 0 0.109090909361839;0.145454540848732 0 0.145454540848732;0.181818187236786 0 0.181818187236786;0.218181818723679 0 0.218181818723679;0.254545450210571 0 0.254545450210571;0.290909081697464 0 0.290909081697464;0.327272742986679 0 0.327272742986679;0.363636374473572 0 0.363636374473572;0.400000005960464 0 0.400000005960464;0.363636374473572 0 0.454545468091965;0.327272742986679 0 0.509090900421143;0.290909081697464 0 0.563636362552643;0.254545450210571 0 0.618181824684143;0.218181818723679 0 0.672727286815643;0.181818187236786 0 0.727272748947144;0.145454540848732 0 0.781818211078644;0.109090909361839 0 0.836363613605499;0.072727270424366 0 0.890909075737;0.036363635212183 0 0.9454545378685;0 0 1;0 0.0666666701436043 1;0 0.133333340287209 1;0 0.200000002980232 1;0 0.266666680574417 1;0 0.333333343267441 1;0 0.400000005960464 1;0 0.466666668653488 1;0 0.533333361148834 1;0 0.600000023841858 1;0 0.666666686534882 1;0 0.733333349227905 1;0 0.800000011920929 1;0 0.866666674613953 1;0 0.933333337306976 1;0 1 1;0.0769230797886848 1 0.923076927661896;0.15384615957737 1 0.846153855323792;0.230769231915474 1 0.769230782985687;0.307692319154739 1 0.692307710647583;0.384615391492844 1 0.615384638309479;0.461538463830948 1 0.538461565971375;0.538461565971375 1 0.461538463830948;0.615384638309479 1 0.384615391492844;0.692307710647583 1 0.307692319154739;0.769230782985687 1 0.230769231915474;0.846153855323792 1 0.15384615957737;0.923076927661896 1 0.0769230797886848;1 1 0;1 0.923076927661896 0;1 0.846153855323792 0;1 0.769230782985687 0;1 0.692307710647583 0;1 0.615384638309479 0;1 0.538461565971375 0;1 0.461538463830948 0;1 0.384615391492844 0;1 0.307692319154739 0;1 0.230769231915474 0;1 0.15384615957737 0;1 0.0769230797886848 0;1 0 0];
data.cMap=cMap;%Exports the colour map used
set(IRIMGGUI,'colormap',cMap);
imgCurrent = stackIntImg(idxImg,:,:,:);
imgCurrent = reshape(imgCurrent,size(imgCurrent,2),size(imgCurrent,3));%Removes the first dimension
nbMapLevels = 64;%Sets the number of levels

vecMapLevelList = imgZMin:((imgZMax-imgZMin)/(nbMapLevels-1)):imgZMax;%Generates the colour vector map
vecTickList = cell(6,1);
for i = 1:6
    vecTickList{i}= num2str(((i*(10/nbMapLevels)*(imgZMax-imgZMin))+imgZMin),'%f');
end %for
data.vecTickList = vecTickList; %Saves the tick list
[imgXPos,imgYpos]=meshgrid(1:size(imgCurrent,1),1:size(imgCurrent,2));%Creates a mesh grid of the image dimension

imgType = get(handles.popupmenuImageType,'string');
imgType =imgType{get(handles.popupmenuImageType,'value')};

drawnow;
hold(handles.axesImage,'off')
if strcmp(imgType,'Contour map')==1
contourf(handles.axesImage,imgXPos,imgYpos,imgCurrent,'LevelList',vecMapLevelList,'Fill','on','lineStyle','none');
handles.colorbar = colorbar('peer',handles.axesImage); %This line adds a colour scale
else % Pixel image
image(1+((imgCurrent-imgZMin).*((nbMapLevels-1)/(imgZMax-imgZMin))),'Parent',handles.axesImage); % Plots the image
axes(handles.axesImage); %Set current axes
%handles.colorbar = colorbar('peer',handles.axesImage,'TickDir','out','YTickLabel',vecTickList); %Adds a colour scale
%handles.colorbar = colorbar('TickDir','out','YTickLabel',vecTickList); %Adds a colour scale
end %if strcmp(imgType,'')==1
hold(handles.axesImage,'on')

idxXPos=data.idxXPos;%Imports the integrated region indexes
idxYPos=data.idxYPos;%Imports the integrated region indexes

for idxR = 1:length(idxXPos)
    plot(handles.axesImage,[idxXPos{idxR};idxXPos{idxR}(1)],[idxYPos{idxR};idxYPos{idxR}(1)],...
        'LineWidth',4,'Color',[0.99 0.99 0.99]);%Plots the wide white lines.
    plot(handles.axesImage,[idxXPos{idxR};idxXPos{idxR}(1)],[idxYPos{idxR};idxYPos{idxR}(1)],...
        'LineWidth',2,'Color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);%Plots teh narrow black lines.
%     text(handles.axesImage,mean(idxXPos{idxR}(:)),mean(idxYPos{idxR}(:)),num2str(idxR),...
%         'FontSize',18,'color',[0.99,0.99,0.99],'FontWeight','Bold');%Adds the polygone index white number
%     text(handles.axesImage,mean(idxXPos{idxR}(:)),mean(idxYPos{idxR}(:)),num2str(idxR),...
%         'FontSize',14,'color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);%Adds the polygone index white number
end %for idxRegion = 1:length(idxXPos)

set(handles.axesImage,'XDir','normal','xlim',[imgXMin imgXMax]);%Scalew the X axis
set(handles.axesImage,'YDir','normal','ylim',[imgYMin imgYMax])%Scalew the Y axis
xlabel(handles.axesImage,'X pixel'); %Adds the axis label
ylabel(handles.axesImage,'Y pixel'); %Adds the axis label
drawnow;

%% Plots the integrated regions of all images (Third axis)

selectedFrameUnit=get(handles.popupmenuStackAxisUnit,'String');%Get the handle's string
selectedFrameUnit=selectedFrameUnit{get(handles.popupmenuStackAxisUnit,'Value')};%Gets the frame axis unit value
if strcmp(selectedFrameUnit,'Relative time (s)')==1
try
stackXvec = data.timeRel;%Sets the relative time as x Axis
catch
stackXvec = (1:nbFiles);%Frame number only   
set(handles.popupmenuStackAxisUnit,'Value',1);%Resets the value
end %Try    
else
stackXvec = (1:nbFiles);%Frame number only
end %if
    
stackXMax = str2double(get(handles.editStackXMax,'string'));
if isnan(stackXMax)==1||stackXMax>max(stackXvec);stackXMax = max(stackXvec);end
stackXMin = str2double(get(handles.editStackXMin,'string'));
if isnan(stackXMin)==1||stackXMin<min(stackXvec);stackXMin=min(stackXvec);end
if stackXMin>stackXMax;
    [stackXMin,stackXMax]=swapValues(stackXMin,stackXMax);
    elseif stackXMin==stackXMax
     stackXMin=min(stackXvec);
     stackXMax = max(stackXvec);
     if  stackXMax ==1;  stackXMax = 2;end;
end
stackYMax = str2double(get(handles.editStackYMax,'string'));
if isnan(stackYMax)==1;stackYMax = max(max(stackRegionAverageInt));end
if isnan(stackYMax)==1;stackYMax = 1;end
stackYMin = str2double(get(handles.editStackYMin,'string'));
if isnan(stackYMin)==1;stackYMin = min(min(stackRegionAverageInt));end
if isnan(stackYMin)==1;stackYMin = 0;end
if stackYMin>stackYMax;%Tests that min is less than max
    [stackYMin,stackYMax]=swapValues(stackYMin,stackYMax);
    elseif stackYMin==stackYMax
    stackYMin = min(min(stackRegionAverageInt))-(0.2*(min(min(stackRegionAverageInt))));
    stackYMax = max(max(stackRegionAverageInt))+(0.1*(max(max(stackRegionAverageInt))));
end %if
 
set(handles.editStackXMax,'string',num2str(stackXMax));%Resets the GUI parameter
set(handles.editStackXMin,'string',num2str(stackXMin));%Resets the GUI parameter
set(handles.editStackYMax,'string',num2str(stackYMax));%Resets the GUI parameter
set(handles.editStackYMin,'string',num2str(stackYMin));%Resets the GUI parameter

hold(handles.axesStack,'off')
boxChecked = get(handles.checkboxSigFit,'Value');

for idxR = 1:size(stackRegionAverageInt,2)
    plot(handles.axesStack,stackXvec,stackRegionAverageInt(:,idxR),...
    'LineWidth',1,'LineStyle','none','MarkerSize',5,'Marker','o',...
    'Color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)],'MarkerFaceColor',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);
    hold(handles.axesStack,'on')
    if boxChecked==1
        if strcmp(selectedFrameUnit,'Relative time (s)')==1
        plot(handles.axesStack,data.timeRelInter,data.fittedInterValues{idxR},'LineWidth',2,'Color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);
        else
        plot(handles.axesStack,(1:length(data.fittedValues{idxR})),data.fittedValues{idxR},'LineWidth',2,'Color',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);
        end %if
 
    end %if boxChecked==1
end

plot(handles.axesStack,[stackXvec(idxImg),stackXvec(idxImg)],[stackYMin,stackYMax],...
'LineWidth',2,'LineStyle','--','Color',[0.01,0.01,0.01]);%Plot cursor line
hold(handles.axesStack,'on')

xlabel(handles.axesStack,selectedFrameUnit);%Write the x axis label
%set(handles.axesStack,'XDir','normal','xlim',[stackXMin,stackXMax]);%Scalew the X axis
%set(handles.axesStack,'YDir','normal','ylim',[stackYMin,stackYMax])%Scalew the Y axis

%Gets the handles of all the axes
% Make ticks outside
set(handles.axesSpectrum,'XMinorTick','on','TickDir','out','YMinorTick','on','TickDir','out');
set(handles.axesImage,'XMinorTick','off','TickDir','out','YMinorTick','off','TickDir','out');
set(handles.axesStack,'XMinorTick','on','TickDir','out','YMinorTick','on','TickDir','out');

%% Plots the histogram (Fourth axis)

% %Get the offset ratio from the GUI
% 
% histXMax = str2double(get(handles.editHistXMax,'string'));%Gets the GUI parameter
% histYMax = str2double(get(handles.editHistYMax,'string'));%Gets the GUI parameter
% histYMin = str2double(get(handles.editHistYMin,'string'));%Gets the GUI parameter
% histXMin = str2double(get(handles.editHistXMin,'string'));%Gets the GUI parameter
% offsetRatio = str2double(get(handles.editOffsetRatio,'string'));%Gets the GUI parameter
% 
% if offsetRatio < 0 || offsetRatio > 1; offsetRatio = 0.5;end %Tests the offset ratio
% 
% if isnan(histXMax)==1;histXMax = max(max(data.stackRegionHist))*offsetRatio*size(data.imgAll,1);end%Tests the GUI parameter
% if isnan(histXMin)==1;histXMin = 0;end 
% if histXMin>histXMax; %Tests the GUI parameter
% [histXMin,histXMax]=swapValues(histXMin,histXMax);
% elseif histXMin==histXMax
% histXMax = max(max(data.stackRegionHist))*offsetRatio*size(data.imgAll,1);histXMin = 0;   
% end %if
% 
% if isnan(histYMax)==1;histYMax = max(data.stackRegionHistBin);end%Tests the GUI parameter
% if isnan(histYMin)==1;histYMin = min(data.stackRegionHistBin);end
% if histYMin>histYMax;%Tests the GUI parameter
% [histYMin,histYMax]=swapValues(histYMin,histYMax);
% elseif histYMin==histYMax;histYMax = max(data.stackRegionHistBin);histYMin = min(data.stackRegionHistBin);
% end %if
% 
% histBinNb = round(str2double(get(handles.editBin,'string')));%Gets the number of bins
% if histBinNb<2 || isnan(histBinNb)==1; histBinNb = 2;end %Test the number of bins
% 
% set(handles.editHistYMax,'string',num2str(histYMax));%Resets the GUI parameter
% set(handles.editHistYMin,'string',num2str(histYMin));%Resets the GUI parameter
% set(handles.editHistXMax,'string',num2str(histXMax));%Resets the GUI parameter
% set(handles.editHistXMin,'string',num2str(histXMin));%Resets the GUI parameter
% set(handles.editBin,'string',num2str(histBinNb));%Resets the GUI parameter
% set(handles.editOffsetRatio,'string',num2str(offsetRatio));%Resets the GUI parameter
% 
% vecX = data.stackRegionHistBin(1:end-1);%Gets the bin vector
% stackRegionHist = data.stackRegionHist;%Gets the histogram data to plot
% 
% hold(handles.axesHisto,'off')
% for idx = 1:size(data.imgAll,1)
%     offset = (idx-1) .* offsetRatio .* max(max(data.stackRegionHist)); %Calculates the offset
%     plot(handles.axesHisto,stackRegionHist(idx,:)+ offset,vecX,...
%     'LineWidth',1,'LineStyle','-','MarkerSize',5,'Marker','.',...
%     'Color',[rCMap(1,idx) rCMap(2,idx) rCMap(3,idx)],'MarkerFaceColor',[rCMap(1,idxR) rCMap(2,idxR) rCMap(3,idxR)]);
%     hold(handles.axesHisto,'on')
% end
% 
% offset = (idxImg-1).* offsetRatio .* max(max(data.stackRegionHist));
% plot(handles.axesHisto,[offset,offset],[min(data.stackRegionHistBin),max(data.stackRegionHistBin)],...
%     'LineWidth',2,'LineStyle','--','Color',[0.01,0.01,0.01]);%Plot cursor line
%     hold(handles.axesHisto,'on')
% 
% xlabel(handles.axesHisto,'Pixel Count');%Write the x axis label
% ylabel(handles.axesHisto,'Integrated Value');%Write the y axis label
% 
% set(handles.axesHisto,'XDir','normal','xlim',[histXMin,histXMax]);%Scalew the X axis
% set(handles.axesHisto,'YDir','normal','ylim',[histYMin,histYMax])%Scalew the Y axis

%Offest each histogram by 0.5 max 

%% Resets the button down functions which were reseted by the plot functions
set(handles.axesSpectrum,'ButtonDownFcn','IRIMGGUI(''clickSpectrum'')');
set(handles.axesImage,'ButtonDownFcn','IRIMGGUI(''clickImage'')');
set(handles.axesStack,'ButtonDownFcn','IRIMGGUI(''clickStack'')');

displayString =[data.SmpFileList{idxImg},' | Refreshed'];  
disp(displayString);set(handles.guiEntire,'Name',displayString);

set(handles.uipanelCompute,'Visible','off');%Hides uipanel
set(handles.uipanelProcessing,'Visible','off');%Hides uipanel

set(handles.uipanelVisualisation,'Visible', 'on');%Shows uipanel
set(handles.uipanelVisualisationLow,'Visible', 'on');%Shows uipanel

set(handles.uipanelSpectra,'Visible','on');%Shows uipanel
set(handles.uipanelImage,'Visible','on');%Shows uipanel
set(handles.uipanelStack,'Visible','on');%Shows uipanel
set(handles.uipanelHistogram,'Visible','on');%Shows uipanel
set(handles.uipanelSlice,'Visible','on');%Shows uipanel


%Saves the index in case the mat file is loaded later 

function clickSpectrum

global data
global hAll

handles = hAll;

h = handles.popupmenuSpectrumIntegrationMode;
spectrumIntegrationType = get(h,'string');
spectrumIntegrationType = spectrumIntegrationType{get(h,'value')};  

set(handles.uipanelSpectra,'visible','on') %Shows the graph
set(handles.uipanelImage,'visible','off') %Hides the other graphs
set(handles.uipanelStack,'visible','off') %Hides the other graphs
set(handles.uipanelHistogram,'Visible','off');%Hides the other graphs

switch spectrumIntegrationType
    case 'Single peak'
        
set(handles.editSpectrumIntervalHigh,'visible','on');
set(handles.editSpectrumIntervalLow,'visible','on');
set(handles.editSpectrumSingleSelect,'visible','off');


displayString =('Please select a two extremity points on the spectrum');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
set(handles.guiEntire,'Name',displayString);
[cursorX,~]=ginput(2);
cursorX=roundn(cursorX,-1);% Rounds to the first decimal digit
set(handles.editSpectrumIntervalHigh,'string',num2str(max(cursorX)));
set(handles.editSpectrumIntervalLow,'string',num2str(min(cursorX)));

data.spcInt.type=spectrumIntegrationType;%Saves the integartion type in the storage structure
data.spcInt.values(1) = min(cursorX);
data.spcInt.values(2) = max(cursorX);%Saves the range to the storage 

displayString =[num2str(max(cursorX)),' to ',num2str(min(cursorX)),' cm-1 selected'];  
disp(displayString);set(handles.guiEntire,'Name',displayString);
    
    case 'Peak with baseline'
        %to write

    case 'Peak ratio'
        %to write

    case 'PC1'
                %to write
    case 'PC2'
                %to write
    case 'PC3'
                %to write
                
    case 'COG'
set(handles.editSpectrumIntervalHigh,'visible','on'); %Hides the other panes
set(handles.editSpectrumIntervalLow,'visible','on'); %Hides the other panes
set(handles.editSpectrumSingleSelect,'visible','off'); %Hides the other panes

displayString =('To calculate the COG, please select a two extremity points on the spectrum');  
disp(displayString); set(handles.guiEntire,'Name',displayString);
[cursorX,~]=ginput(2);
cursorX=roundn(cursorX,-1);% Rounds to the first decimal digit
set(handles.editSpectrumIntervalHigh,'string',num2str(max(cursorX)));
set(handles.editSpectrumIntervalLow,'string',num2str(min(cursorX)));

data.spcInt.type=spectrumIntegrationType;%Saves the integartion type in the storage structure
data.spcInt.values(1) = min(cursorX);
data.spcInt.values(2) = max(cursorX);%Saves the range to the storage 

displayString =[num2str(max(cursorX)),' to ',num2str(min(cursorX)),' cm-1 selected'];  
disp(displayString); set(handles.guiEntire,'Name',displayString);      
                
    otherwise %Single wavelength
displayString =('Please select a single pixel on the spectrum');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
[cursorX,~]=ginput(1);
cursorX=roundn(cursorX,-1);
set(handles.editSpectrumSingleSelect,'string',num2str(cursorX))
data.spcInt.type=spectrumIntegrationType;
data.values = cursorX;
displayString =[num2str(cursorX),' cm-1 selected'];  
disp(displayString); set(handles.guiEntire,'Name',displayString);

end % Switch case

set(handles.uipanelSpectra,'visible','on') %Shows the panel
set(handles.uipanelImage,'visible','on')  %Shows the panel
set(handles.uipanelHistogram,'visible','on')  %Shows the panel
set(handles.uipanelStack,'visible','on')  %Shows the panel

integrateImageRegions(handles) % Integrates the stack using the GUI parameters
refreshGraph(handles) %Refresh graphs 

function clickImage

global data
global hAll
handles = hAll;

h = handles.popupmenuImgSelectionMode;
imgIntegrationType = get(h,'string');
imgIntegrationType = imgIntegrationType{get(h,'value')};  

%Hides the two other graphs
set(handles.uipanelSpectra,'visible','off') %Shows the panel
set(handles.uipanelImage,'visible','on') %Shows  the panel
set(handles.uipanelHistogram,'visible','off')  %Hides the panel
set(handles.uipanelStack,'visible','off') %Hides the panel

switch imgIntegrationType
    
    case 'Single polygone'
    displayString =('Please select polygones by clicking on the image before pressing ENTER');  
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    waitfor(msgbox(displayString,'Message')),

    [cursorX,cursorY] = ginput;%Selects a single polygone
    cursorX = round(cursorX);%Rounds the X cursor position
    cursorY=round(cursorY);%Rounds the Y cursor position

    data.idxXPos{1} = cursorX;    %Saves the X cursor position in the storage array
    data.idxYPos{1} = cursorY;    %Saves the Y cursor position in the storage array

    cursorXstring=''; %Initiates the cursor string
    cursorYstring=''; %Initiates the cursor string
    for idxPoint=1:length(cursorX)
    cursorXstring = [cursorXstring,' ',num2str(cursorX(idxPoint))];%Creates the X cursor string
    cursorYstring = [cursorYstring,' ',num2str(cursorY(idxPoint))];%Creates the X cursor string
    end %for
    cursorXstring = cursorXstring(2:end);
    cursorYstring = cursorYstring(2:end);

    displayString =['Single polygone selected'];  
    disp(displayString);set(handles.guiEntire,'Name',displayString);

    set(handles.editImageXpos,'string',num2str(cursorXstring));
    set(handles.editImageYpos,'string',num2str(cursorYstring));

        case 'Single square'
            %To write
        case 'Multiple pixels'
           %To write
        case 'Multiple polygones'

    displayString =('Please select polygone(s) by clicking on the image before pressing ENTER');  
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    waitfor(msgbox(displayString,'Message')),

    strCondition = 'Yes'; idxRegion = 0;polygone=[]; %Sets the initial parameters
    while strcmp(strCondition,'Yes')==1 && idxRegion<50

    [cursorX,cursorY] = ginput; %Selects a polygone using the dialogue
    idxRegion =  idxRegion+1; %Increments the number of polygone
    cursorX = round(cursorX); %Rounds the pixel position
    cursorY=round(cursorY); %Rounds the pixel position
    data.idxXPos{idxRegion} = cursorX;%Saves the data in the storage structure
    data.idxYPos{idxRegion} = cursorY;%Saves the data in the storage structure
    cursorX = cat(1,cursorX,cursorX(1));cursorY = cat(1,cursorY,cursorY(1));%Concatenate the array to add the first value
    polygoneCur=line(cursorX,cursorY,'LineWidth',3,'Color',[1 1 1]);%Plots the selected polygone
    polygone=cat(1,polygone,polygoneCur);%Appends the array
    polygoneCur=line(cursorX,cursorY,'LineWidth',2,'Color',[0 0 0]);%Plots the selected polygone
    polygone=cat(1,polygone,polygoneCur);%Appends the array
    
    strCondition = questdlg([num2str(idxRegion),' polygone(s) selected. Do you wish to select another one ?'], ...
    'Select another polygone ?','Yes','No','Yes');
    end % While loop
    nbRegion = idxRegion;%Saves the number of regions

    %Creates the GUI string
    cursorXstring='';cursorYstring=''; %Initiates the cursor string
    for idxRegion = 1:nbRegion
        for idxPoint=1:length(data.idxXPos{idxRegion})
            cursorXstring = [cursorXstring,' ',num2str(data.idxXPos{idxRegion}(idxPoint))];%Adds the next point    
            cursorYstring = [cursorYstring,' ',num2str(data.idxYPos{idxRegion}(idxPoint))];%Adds the next point  
        end %for
        cursorXstring = [cursorXstring,';'];% Adds the delimeter of the polygone
        cursorYstring = [cursorYstring,';'];% Adds the delimeter of the polygone
    end %for
    cursorXstring = cursorXstring(2:end-1);%Removes the first and last caracter
    cursorYstring = cursorYstring(2:end-1);%Removes the first and last caracter

    displayString =[idxRegion,'Polygone(s) selected'];  
    disp(displayString);set(handles.guiEntire,'Name',displayString);


    set(handles.editImageXpos,'string',num2str(cursorXstring));%Writes the cursor string in the GUI
    set(handles.editImageYpos,'string',num2str(cursorYstring));%Writes the cursor string in the GUI

        case 'Multiple squares'

        otherwise %Single pixel
    displayString =('Please select a single pixel on the image');  
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    [cursorX,cursorY]=ginput(1);
    cursorX=round(cursorX);
    cursorY=round(cursorY);
    set(handles.editImageXpos,'string',num2str(cursorX));
    set(handles.editImageYpos,'string',num2str(cursorY));
    data.idxXPos{1} = cursorX;
    data.idxYPos{1} = cursorY;
    displayString =['X = ',num2str(cursorX),' , Y = ',num2str(cursorY),' pixel selected'];  
    disp(displayString);set(handles.guiEntire,'Name',displayString);

end % switch imgIntegrationType
        
set(handles.uipanelSpectra,'visible','on') %Shows the panel
set(handles.uipanelImage,'visible','on')  %Shows the panel
set(handles.uipanelHistogram,'visible','on')  %Shows the panel
set(handles.uipanelStack,'visible','on')  %Shows the panel

data.imgIntegrationType = imgIntegrationType;%Export the integration type
integrateImageRegions(handles) % Integrates the stack using the GUI parameters


function clickStack

global data
global hAll

handles = hAll;

%Hides the two other graphs
set(handles.uipanelSpectra,'visible','off') %Hides the panel
set(handles.uipanelImage,'visible','off') %Hides the panel
set(handles.uipanelHistogram,'visible','off') %Hides the panel
set(handles.uipanelStack,'visible','on') %Shows the panel

displayString =('Please select an image from the stack');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
[cursorX,~]=ginput(1);
cursorX=round(cursorX);
set(handles.editFrameSelect,'string',num2str(cursorX));
data.idxImg = cursorX;
displayString =['Image ',num2str(cursorX),' selected'];  
disp(displayString);set(handles.guiEntire,'Name',displayString);

set(handles.uipanelSpectra,'visible','on') %Shows the panel
set(handles.uipanelImage,'visible','on')  %Shows the panel
set(handles.uipanelHistogram,'visible','on')  %Shows the panel
set(handles.uipanelStack,'visible','on')  %Shows the panel

refreshGraph(handles)%Refresh graphs 

function computeEnable(handles,onOrOff)

handleGroup=[handles.checkboxBaseline,...
    handles.popupmenuApodisation,...
    handles.popupmenuZeroFilling,...
    handles.checkboxMertz,...
    handles.checkboxInterZeroOffset,...
    handles.checkboxSingleNorm,...
    handles.checkboxSingleNorm,...
    handles.checkboxCburstCorr,...
    handles.editNormHigh,...
    handles.editNormLow,...
    handles.editForwardTruncation,...
    handles.editImportImgLow,...
    handles.editImportImgHigh,...
    handles.checkboxInterNorm,...
    ];
if strcmp(onOrOff,'off') == 0;onOrOff='on';end
set(handleGroup,'enable',onOrOff);%Sets the property of the handle group

checkboxSingleNorm_Callback([],[],handles);

function refreshGui(handles)%Refreshes the GUI
%  List of spectral operation
for idxOp = 1:5

handleOp = findobj('Tag',['popupmenuOp',num2str(idxOp)]);
handleTopBox = findobj('Tag',['editOp',num2str(idxOp),'Box1']);
handleLowerBox = findobj('Tag',['editOp',num2str(idxOp),'Box2']);
handleClickButton = findobj('Tag',['pushbuttonClickSelectOp',num2str(idxOp)]);
handleNextOp  = findobj('Tag',['popupmenuOp',num2str(idxOp+1)]);

operationCurrent = get(handleOp,'string');
operationCurrent = operationCurrent{get(handleOp,'value')};

if strcmp(operationCurrent,'Truncate')==1 || strcmp(operationCurrent,'Normalise')==1 ||...
strcmp(operationCurrent,'Offset')==1 || strcmp(operationCurrent,'Zap')==1
set(handleOp,'enable','on');
set(handleTopBox,'enable','on');
set(handleLowerBox,'enable','on');
set(handleNextOp,'enable','on');
set(handleClickButton,'enable','on');

elseif strcmp(operationCurrent,'Interpolation spline')==1
set(handleOp,'enable','on');
set(handleTopBox,'enable','on');
set(handleLowerBox,'enable','off');
set(handleNextOp,'enable','on');   
set(handleLowerBox,'string','');
set(handleClickButton,'enable','off');

elseif strcmp(operationCurrent,'Smooth spc S-Golay')==1||...
        strcmp(operationCurrent,'Smooth spc moving')==1||...
         strcmp(operationCurrent,'Smooth spc linear reg')==1||...
         strcmp(operationCurrent,'Smooth 3D')
set(handleOp,'enable','on');
set(handleTopBox,'enable','on');
set(handleLowerBox,'enable','off');
set(handleNextOp,'enable','on');
set(handleLowerBox,'string','');
set(handleClickButton,'enable','off');

elseif strcmp(operationCurrent,'Smooth image')==1
set(handleOp,'enable','on');
set(handleTopBox,'enable','on');
set(handleLowerBox,'enable','off');
set(handleNextOp,'enable','on');
set(handleLowerBox,'string','');
set(handleClickButton,'enable','off');

elseif strcmp(operationCurrent,'ATR correction')==1
set(handleOp,'enable','on');set(handleTopBox,'enable','off');
set(handleLowerBox,'enable','off');set(handleNextOp,'enable','on');

elseif strcmp(operationCurrent,'Water vapour')==1
set(handleOp,'enable','on');set(handleTopBox,'enable','off');
set(handleLowerBox,'enable','off');set(handleNextOp,'enable','on');       

elseif strcmp(operationCurrent,'Multiple polygone frame offset')==1||...
        strcmp(operationCurrent,'Multiple polygone interpolated offset')
set(handleOp,'enable','on');set(handleTopBox,'enable','on');
set(handleLowerBox,'enable','on');set(handleNextOp,'enable','on');
set(handleClickButton,'enable','on');

else %No operation
set(handleTopBox,'enable','off');
set(handleLowerBox,'enable','off');
set(handleNextOp,'enable','off');
set(handleNextOp,'value',find(strcmp(get(handleNextOp,'string'),'No operation')));
set(handleClickButton,'enable','off');
end %if

end %idxOp = 1:5 

% --- Executes just before IRIMGGUI is made visible.
function IRIMGGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = IRIMGGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in buttonImport.
function buttonImport_Callback(hObject, eventdata, handles)
importImage(handles)


% --- Executes on button press in buttonHeaderFile.
function buttonHeaderFile_Callback(hObject, eventdata, handles)
%Selects the header file and loads its path
[headerFile,directory] = uigetfile('*.hdr','Select hdr header file','Select hdr ENVI header file');
%Tests if a header file is selected
if headerFile == 0
displayString = 'No header file selected';
set(handles.guiEntire,'Name',displayString);
set(handles.buttonSelectImage,'enable','off');
set(handles.textImgFilePath,'enable','off');
 return
end %aData.infoFile == []
headerFileFull = [directory,headerFile];%Concatenates the directory string to filename
set(handles.textHeaderpath,'String',headerFileFull);%Sets the selected file as string
displayString = ['Header file selected : ',headerFile];
set(handles.guiEntire,'Name',displayString); disp(displayString)
set(handles.buttonSelectImage,'enable','on');
set(handles.textImgFilePath,'enable','on');


% --- Executes on button press in buttonSelectImage.
function buttonSelectImage_Callback(hObject, eventdata, handles)
global data

[data.SmpFileListCur,smpDirectory] = uigetfile({'*.seq';'*.dat'},'Select *.seq or *.dat file(s)','MultiSelect', 'on'); 
% Now you can set the "string" property of the static textbox.
% Puts only the filename if only one file was selected
if iscell(data.SmpFileListCur)==1 %many files sselected
for fileNumber = 1: size(data.SmpFileListCur,2) %Add full path
currentName = [smpDirectory,char(data.SmpFileListCur(1,fileNumber))];
data.SmpFileListCur{fileNumber}= currentName;
end
displayString =[smpDirectory,' (',num2str(size(data.SmpFileListCur,2)),') files']; 
set(handles.textImgFilePath,'String',displayString)
set(handles.guiEntire,'Name',displayString)
elseif data.SmpFileListCur~=0 %Only one file selected
data.SmpFileListCur = [smpDirectory,data.SmpFileListCur];
displayString = data.SmpFileListCur;
set(handles.textImgFilePath,'String',displayString);
set(handles.guiEntire,'Name',displayString);
data.SmpFileListCur = {data.SmpFileListCur};%Converts to a cell
else
displayString =('No valid files selected');
disp(displayString);set(handles.guiEntire,'Name',displayString);
set(handles.textImgFilePath,'String',displayString)
set(handles.guiEntire,'Name',displayString)

set(handles.checkboxCompute,'enable','off');
set(handles.buttonImport,'enable','off');
set(handles.popupmenuRatio,'enable','off');
computeEnable(handles,'off');
return
end

set(handles.checkboxCompute,'enable','on');
set(handles.buttonImport,'enable','on');
set(handles.popupmenuRatio,'enable','on');
computeEnable(handles,'on');

refreshGui(handles)

% --- Executes on button press in pushbuttonSelectBkg.
function pushbuttonSelectBkg_Callback(hObject, eventdata, handles)
global data
[data.bkgListCur,bkgDirectory] = uigetfile({'*.seq';'*.dat'},'Select background *.seq or *.dat file(s)','Select bkg file(s)','MultiSelect', 'on'); 
% Now you can set the "string" property of the static textbox.
% Puts only the filename if only one file was selected
if iscell(data.bkgListCur)==1 %many files sselected
for fileNumber = 1: size(data.bkgList,2) %Add full path
currentName = [bkgDirectory,char(data.bkgListCur(1,fileNumber))];
data.bkgList{fileNumber}= currentName;
end
displayString =[bkgDirectory,' (',num2str(size(data.bkgListCur,2)),') files']; 
set(handles.textSelectBkg,'String',displayString);
set(handles.display,'String',displayString)
elseif data.bkgListCur~=0 %Only one file selected
data.bkgListCur = [bkgDirectory,data.bkgListCur];
displayString = data.bkgListCur;
data.bkgListCur = {data.bkgListCur}; % Converts to cells
set(handles.textSelectBkg,'String',displayString);
set(handles.guiEntire,'Name',displayString);
else
displayString =('No valid files selected');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
set(handles.textSelectBkg,'String',displayString)
return

end
% --- Executes on selection change in popupmenuRatio.
function popupmenuRatio_Callback(hObject, eventdata, handles)
handleString = get(handles.popupmenuRatio,'String');
idxValue = get(handles.popupmenuRatio,'Value');
if strcmp(handleString(idxValue),'No ratio')==1;
set(handles.pushbuttonSelectBkg,'enable','off');
set(handles.textSelectBkg,'enable','off');
displayString =('No ratio will be applied to the images');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
else
set(handles.pushbuttonSelectBkg,'enable','on');
set(handles.textSelectBkg,'enable','on');
displayString =('A ratio will be applied to the images');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
end

% --- Executes during object creation, after setting all properties.
function popupmenuRatio_CreateFcn(hObject, eventdata, handles)

% --- Executes on selection change in popupmenuZeroFilling.
function popupmenuZeroFilling_Callback(hObject, eventdata, handles)
global data
nbPoints = data.bands ;
nbPoints = 2^nextpow2(nbPoints);
handleString = get(handles.popupmenuZeroFilling,'String');
idxValue = get(handles.popupmenuZeroFilling,'Value');
if strcmp(handleString(idxValue),'Zero filling 1x')==1;
displayString =['Next power of 2 multiplied by 1 ( ', num2str(nbPoints*0.5),' data points )']; 
disp(displayString);set(handles.guiEntire,'Name',displayString);
elseif strcmp(handleString(idxValue),'Zero filling 3x')==1;
displayString =['Next power of 2 multiplied by 3 ( ', num2str(nbPoints*1.5),' data points )'];
disp(displayString);set(handles.guiEntire,'Name',displayString);
elseif strcmp(handleString(idxValue),'Zero filling 4x')==1;
displayString =['Next power of 2 multiplied by 4 ( ', num2str(nbPoints*2),' data points )'];
disp(displayString);set(handles.guiEntire,'Name',displayString);
else
displayString =['Next power of 2 multiplied by 2 ( ', num2str(nbPoints*1),' data points )'];  
disp(displayString);set(handles.guiEntire,'Name',displayString);
end

% --- Executes during object creation, after setting all properties.
function popupmenuZeroFilling_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuApodisation.
function popupmenuApodisation_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupmenuApodisation_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxMertz.
function checkboxMertz_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkboxBaseline.
function checkboxBaseline_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkboxCompute.
function checkboxCompute_Callback(hObject, eventdata, handles)
ticked = get(handles.checkboxCompute,'Value');
if ticked==0
computeEnable(handles,'off');
displayString =('No computing will be performed on image(s)');  
else
computeEnable(handles,'on');
displayString =('Computing will be performed on image(s)');  
end %if
disp(displayString);set(handles.guiEntire,'Name',displayString);


function edit1_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
refreshGraph(handles)


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuOp1.
function popupmenuOp1_Callback(hObject, eventdata, handles)
refreshGui(handles)


% --- Executes during object creation, after setting all properties.
function popupmenuOp1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuOp2.
function popupmenuOp2_Callback(hObject, eventdata, handles)
refreshGui(handles)


% --- Executes during object creation, after setting all properties.
function popupmenuOp2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuOp3.
function popupmenuOp3_Callback(hObject, eventdata, handles)
refreshGui(handles)


% --- Executes during object creation, after setting all properties.
function popupmenuOp3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOp1Box1_Callback(hObject, eventdata, handles)
refreshGui(handles)


% --- Executes during object creation, after setting all properties.
function editOp1Box1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editOp1Box2_Callback(hObject, eventdata, handles)
refreshGui(handles)


% --- Executes during object creation, after setting all properties.
function editOp1Box2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonReprocess.
function pushbuttonReprocess_Callback(hObject, eventdata, handles)
processingSequence(handles)
set(handles.uipanelCompute,'Visible','off');
set(handles.uipanelProcessing,'Visible','off');
set(handles.uipanelVisualisation,'Visible','on');


% --- Executes on button press in pushbuttonRefreshImport.
function pushbuttonRefreshImport_Callback(hObject, eventdata, handles)
% To write



% --- Executes on selection change in popupmenuOp4.
function popupmenuOp4_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuOp4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuOp5.
function popupmenuOp5_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuOp5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOp2Box1_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp2Box1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOp2Box2_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp2Box2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOp3Box1_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp3Box1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOp3Box2_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp3Box2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editOp4Box1_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp4Box1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOp4Box2_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp4Box2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOp5Box1_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp5Box1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editOp5Box2_Callback(hObject, eventdata, handles)
refreshGui(handles)

% --- Executes during object creation, after setting all properties.
function editOp5Box2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderFrame_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function sliderFrame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on selection change in popupmenuImgSelectionMode.
function popupmenuImgSelectionMode_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuImgSelectionMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonNextSelection.
function pushbuttonNextSelection_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes on button press in pushbuttonResetRegion.
function pushbuttonResetRegion_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes on selection change in popupmenuSpectrumIntegrationMode.
function popupmenuSpectrumIntegrationMode_Callback(hObject, eventdata, handles)
h = handles.popupmenuSpectrumIntegrationMode;
spectrumIntegrationType = get(h,'string');
spectrumIntegrationType = spectrumIntegrationType{get(h,'value')};  

if strcmp(spectrumIntegrationType,'Single peak')==1
        
set(handles.editSpectrumIntervalHigh,'visible','on');
set(handles.editSpectrumIntervalLow,'visible','on');
set(handles.editSpectrumSingleSelect,'visible','off');

elseif strcmp(spectrumIntegrationType,'COG')==1
       
set(handles.editSpectrumIntervalHigh,'visible','on');
set(handles.editSpectrumIntervalLow,'visible','on');
set(handles.editSpectrumSingleSelect,'visible','off');

else %Single wavelength
set(handles.editSpectrumIntervalHigh,'visible','off');
set(handles.editSpectrumIntervalLow,'visible','off');
set(handles.editSpectrumSingleSelect,'visible','on');
end %ife
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuSpectrumIntegrationMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on mouse press over axes background.
function axesImage_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes on mouse press over axes background.
function axesStack_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes on selection change in popupmenuImageType.
function popupmenuImageType_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuImageType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSpectrumIntervalHigh_Callback(hObject, eventdata, handles)
integrateImageRegions(handles) %Integrates the images and refreshed the graphs


% --- Executes during object creation, after setting all properties.
function editSpectrumIntervalHigh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSpectrumIntervalLow_Callback(hObject, eventdata, handles)
integrateImageRegions(handles) %Integrates the images and refreshed the graphs

function editSpectrumIntervalLow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editImageXpos_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editImageXpos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editImageYpos_Callback(hObject, eventdata, handles)
refreshGraph(handles)


% --- Executes during object creation, after setting all properties.
function editImageYpos_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editFrameSelect_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editFrameSelect_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSpectrumWaveMax_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editSpectrumWaveMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSpectrumWaveMin_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editSpectrumWaveMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSpectrumResponseMax_Callback(hObject, eventdata, handles)
refreshGraph(handles)


% --- Executes during object creation, after setting all properties.
function editSpectrumResponseMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSpectrumResponseMin_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editSpectrumResponseMin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editImageXMin_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editImageXMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editImageXMax_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editImageXMax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editImageYMax_Callback(hObject, eventdata, handles)
refreshGraph(handles)


% --- Executes during object creation, after setting all properties.
function editImageYMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editImageYMin_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editImageYMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editStackXMin_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editStackXMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editStackXMax_Callback(hObject, eventdata, handles)
refreshGraph(handles)

function editStackXMax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
refreshGraph(handles)

function editSpectrumSingleSelect_Callback(hObject, eventdata, handles)
integrateImageRegions(handles) %Integrates the images and refreshed the graphs

% --- Executes during object creation, after setting all properties.
function editSpectrumSingleSelect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editStackYMax_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editStackYMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editStackYMin_Callback(hObject, eventdata, handles)
refreshGraph(handles)

% --- Executes during object creation, after setting all properties.
function editStackYMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes6_ButtonDownFcn(hObject, eventdata, handles)
disp('haha')


% --- Executes on button press in pushbuttonSpectrumClick.
function pushbuttonSpectrumClick_Callback(hObject, eventdata, handles)
clickSpectrum

% --- Executes on button press in pushbuttonSpectrumCrop.
function pushbuttonSpectrumCrop_Callback(hObject, eventdata, handles)

%Hides the two other graphs
set(handles.uipanelSpectra,'visible','on') %Activates the panel
set(handles.uipanelImage,'visible','off') %Activates the panel
set(handles.uipanelStack,'visible','off') %Activates the panel
set(handles.uipanelHistogram,'visible','off') %Activates the panel

displayString =('Please select a spectrum window');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
[cursorX,cursorY]=ginput(2);
rangeXMax = max(cursorX);
rangeXMin = min(cursorX);
rangeYMax = max(cursorY);
rangeYMin = min(cursorY);
h=handles.axesSpectrum;
line1=line([rangeXMin,rangeXMax],[rangeYMax,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
line2=line([rangeXMin,rangeXMax],[rangeYMin,rangeYMin],'LineWidth',2,'color',[0.1 0.1 0.1]);
line3=line([rangeXMin,rangeXMin],[rangeYMin,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
line4=line([rangeXMax,rangeXMax],[rangeYMin,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
pause(0.75)
delete([line1,line2,line3,line4]);

set(handles.editSpectrumResponseMax,'string',num2str(rangeYMax))
set(handles.editSpectrumResponseMin,'string',num2str(rangeYMin))
set(handles.editSpectrumWaveMax,'string',num2str(rangeXMax))
set(handles.editSpectrumWaveMin,'string',num2str(rangeXMin))

displayString ='Spectrum cropped';  
disp(displayString);set(handles.guiEntire,'Name',displayString);
       
%Shows the two other graphs
set(handles.uipanelSpectra,'visible','on') %Activates the panel
set(handles.uipanelImage,'visible','on') %Activates the panel
set(handles.uipanelStack,'visible','on') %Activates the panel
set(handles.uipanelHistogram,'visible','on') %Activates the panel

refreshGraph(handles);%Refresh graphs 

% --- Executes on button press in pushbuttonImageCrop.
function pushbuttonImageCrop_Callback(hObject, eventdata, handles)
%Hides the two other graphs
set(handles.uipanelSpectra,'visible','off') %Activates the panel
set(handles.uipanelImage,'visible','on') %Activates the panel
set(handles.uipanelStack,'visible','off') %Activates the panel
set(handles.uipanelHistogram,'visible','off') %Activates the panel

displayString =('Please select an image window');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
[cursorX,cursorY]=ginput(2);
rangeXMax = max(cursorX);
rangeXMin = min(cursorX);
rangeYMax = max(cursorY);
rangeYMin = min(cursorY);
h=handles.axesImage;
line1=line([rangeXMin,rangeXMax],[rangeYMax,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
line2=line([rangeXMin,rangeXMax],[rangeYMin,rangeYMin],'LineWidth',2,'color',[0.1 0.1 0.1]);
line3=line([rangeXMin,rangeXMin],[rangeYMin,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
line4=line([rangeXMax,rangeXMax],[rangeYMin,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
pause(0.75)
delete([line1,line2,line3,line4]);

set(handles.editImageXMax,'string',num2str(rangeXMax))
set(handles.editImageXMin,'string',num2str(rangeXMin))
set(handles.editImageYMax,'string',num2str(rangeYMax))
set(handles.editImageYMin,'string',num2str(rangeYMin))

displayString ='Image cropped';  
disp(displayString);set(handles.guiEntire,'Name',displayString);
       
%Shows the two other graphs
set(handles.uipanelSpectra,'visible','on') %Activates the panel
set(handles.uipanelImage,'visible','on') %Activates the panel
set(handles.uipanelStack,'visible','on') %Activates the panel
set(handles.uipanelHistogram,'visible','on') %Activates the panel

%Refresh graphs 
refreshGraph(handles)

% --- Executes on button press in pushbuttonImageClick.
function pushbuttonImageClick_Callback(hObject, eventdata, handles)
clickImage %Launch

% --- Executes on button press in pushbuttonStackClick.
function pushbuttonStackClick_Callback(hObject, eventdata, handles)
clickStack %Launch

% --- Executes on button press in pushbuttonStackCrop.
function pushbuttonStackCrop_Callback(hObject, eventdata, handles)
%Hides the two other graphs
set(handles.uipanelSpectra,'visible','off') %Activates the panel
set(handles.uipanelImage,'visible','off') %Activates the panel
set(handles.uipanelStack,'visible','on') %Activates the panel
set(handles.uipanelHistogram,'visible','off') %Activates the panel

displayString =('Please select an stack window');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
[cursorX,cursorY]=ginput(2);
rangeXMax = max(cursorX);
rangeXMin = min(cursorX);
rangeYMax = max(cursorY);
rangeYMin = min(cursorY);
h=handles.axesStack;
line1=line([rangeXMin,rangeXMax],[rangeYMax,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
line2=line([rangeXMin,rangeXMax],[rangeYMin,rangeYMin],'LineWidth',2,'color',[0.1 0.1 0.1]);
line3=line([rangeXMin,rangeXMin],[rangeYMin,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
line4=line([rangeXMax,rangeXMax],[rangeYMin,rangeYMax],'LineWidth',2,'color',[0.1 0.1 0.1]);
pause(0.75)
delete([line1,line2,line3,line4]);

set(handles.editStackXMax,'string',num2str(rangeXMax))
set(handles.editStackXMin,'string',num2str(rangeXMin))
set(handles.editStackYMax,'string',num2str(rangeYMax))
set(handles.editStackYMin,'string',num2str(rangeYMin))

displayString ='Image cropped';  
disp(displayString);set(handles.guiEntire,'Name',displayString);
       
%Shows the two other graphs
set(handles.uipanelSpectra,'visible','on') %Activates the panel
set(handles.uipanelImage,'visible','on') %Activates the panel
set(handles.uipanelStack,'visible','on') %Activates the panel
set(handles.uipanelHistogram,'visible','on') %Activates the panel

%Refresh graphs 
refreshGraph(handles)

% --------------------------------------------------------------------
function uipanelImage_ButtonDownFcn(hObject, eventdata, handles)
clickImage


function editRegionSelected_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editRegionSelected_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menuImport_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menuLoad_Callback(hObject, eventdata, handles)

userChoice = questdlg('Unsaved data will be lost. Are you sure you wish to load a *.mat file?',...
	'Load data','Load data','Cancel','Cancel');
switch userChoice
    case  'Cancel'
displayString =('Data not loaded');  
disp(displayString);set(handles.guiEntire,'Name',displayString); drawnow;
case  'Load data'
[fileName,pathName,~] = uigetfile('*.mat','Select a *.mat file');     
if fileName==0
displayString =('No *.mat selected, data not loaded');  
disp(displayString);set(handles.guiEntire,'Name',displayString);  
else
clear data %Clears the old data from the work space
displayString =('Loading *.mat file, please wait...');  
disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
fullname = [pathName,fileName];
load(fullname);
displayString =(['Data loaded ',fullname]);  
disp(displayString);set(handles.guiEntire,'Name',displayString); 

enableGui(handles); %Enables the menus of the GUI
set(handles.guiEntire,'Name',fileName)

idxXPosStr = '';
idxYPosStr = '';
for idxRegion=1:size(data.idxXPos,2)%Loops through the regions
    idxXPosStr = [idxXPosStr,';',mat2str(data.idxXPos{idxRegion})];
    idxYPosStr = [idxYPosStr,';',mat2str(data.idxYPos{idxRegion})];         
end
set(handles.editImageXpos,'string',idxXPosStr(2:end));%Imports the cursor coordinate
set(handles.editImageYpos,'string',idxYPosStr(2:end));%Imports the cursor coordinate

set(handles.guiEntire,'Name',fileName);%Gives the filename to the gui
set(handles.uipanelVisualisationLow,'visible','on') %Activates the panel
set(handles.uipanelVisualisation,'visible','on') %Activates the panel
refreshGraph(handles); %Refresh the graph
end
end %Switch case


% --------------------------------------------------------------------
function menuSave_Callback(hObject, eventdata, handles)
global data
fullName = data.SmpFileList{1};
[currentDirectory,name,ext] = fileparts(fullName);
defaultName = [currentDirectory,'\',name,'.mat'];
[fileName,pathName,~] =  uiputfile(defaultName,'Save data into *.mat file');
if fileName==0
displayString =('No *.mat selected, data not saved');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
else
displayString ='Saving data, please wait...';  
disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
fileName = [pathName,fileName];
save(fileName,'data') %saves only the variables or fields of a structure array specified by variables.
set(handles.guiEntire,'Name',fileName);%Gives the filename to the gui
displayString =['Data saved as ',name];  
disp(displayString);set(handles.guiEntire,'Name',displayString);
set(handles.guiEntire,'Name',displayString)
end


% --------------------------------------------------------------------
function menuExportPDF_Callback(hObject, eventdata, handles)
global data

set(handles.uipanelSpectra,'visible','on') %Activates the panel
set(handles.uipanelImage,'visible','on') %Activates the panel
set(handles.uipanelHistogram,'visible','on') %Activates the panel
set(handles.uipanelStack,'visible','on') %Activates the panel

fullName = data.SmpFileList{1};
[currentDirectory,name,ext] = fileparts(fullName);
defaultName = [currentDirectory,'\',name,'.pdf'];
[fileName,pathName,ext] =  uiputfile(defaultName,'Save Graph as *.pdf');
fullName = [pathName,fileName];
if fileName==0
    displayString =('No *.pdf selected, data not saved');  
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    set(handles.guiEntire,'Name',displayString)
    else
    displayString =['Saving ',fileName,' Please wait'];  
    set(handles.guiEntire,'Name',displayString)
    disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;

    gcf = figure; %Creates a new figure
    set(gcf,'colormap',data.cMap);
    set(gcf,'Position',[300 300 1700 400]);% Sets the size of the figure      
    set(gcf,'PaperOrientation','landscape')%Sets the paper orientation
    set(gcf,'paperunits','inch');%Sets the paper units
    set(gcf,'papersize',[8,22]);%Sets the paper size large enough to fit the GUI
    set(gcf,'paperposition',[0.5,15,22,22]);
    
    h = copyobj(handles.axesSpectrum,gcf);
    positionVal = get(h,'Position');
    positionVal(1) = positionVal(1) + 0;
        set(h,'Position',positionVal);
    
    h = copyobj(handles.axesImage,gcf);
    positionVal = get(h,'Position');
    positionVal(1) = positionVal(1) + 100;
    positionVal(2) = positionVal(2) -0.5;
    set(h,'Position',positionVal);
    colorbar('peer',h);% Creates the color bar
    format shortE
    colorbar('peer',h,'TickDir','out','YTickLabel',data.vecTickList); %Adds a colour scale

    h = copyobj(handles.axesStack,gcf);
    positionVal = get(h,'Position');
    positionVal(1) = positionVal(1) + 180;
    set(h,'Position',positionVal);
    
    h = copyobj(handles.axesHisto,gcf);
    positionVal = get(h,'Position');
    positionVal(1) = positionVal(1) + 270;
    positionVal(2) = positionVal(2) + 2;
    set(h,'Position',positionVal);
    
	saveas(gcf,fullName,'pdf');
    close(gcf);%Closes the figure
    
    pause(0.1)
    displayString =['Graphs saved as ',fileName];  
    disp(displayString);set(handles.guiEntire,'Name',displayString)
end %If


% --------------------------------------------------------------------
function menuExportVideo_Callback(hObject, eventdata, handles)
global data

%GUI handle parameters for GUI building
frameInterpolation = 1;%1 for yes
interFactorFrame = 10;
interFactorXY = 3;
frameRate = 20;
colourImage = 1;%1 for yes
scaleImage =1;%1 for yes

fileName = data.SmpFileList{1};
[currentDirectory,name,~] = fileparts(fileName);
defaultName = [currentDirectory,'\',name,'.mp4'];
[fileName,pathName,~] =  uiputfile({'*.mp4';'*.avi'},'Save images as',defaultName);
[~,~,ext] = fileparts(fileName);
if strcmp(ext,'.mp4')==1
videoFormat = 'MPEG-4';
elseif strcmp(ext,'.avi')==1   
videoFormat = 'Motion JPEG AVI';
else
videoFormat = 'MPEG-4';
end%if
fileName = [pathName,fileName];
if fileName(1)==0
    displayString =('No file selected, video not saved');  
    disp(displayString);
set(handles.guiEntire,'Name',displayString)
else

displayString ='Exporting video, please wait...';  
disp(displayString);set(handles.guiEntire,'Name',displayString);drawnow;
    
images = data.imgInt;%Imports the images data
images=permute(images,[2 3 1]);%Permutes the dimensions

%Scales the data
if scaleImage == 1;
minumumValue=min(min(min(images)));%Finds the minimum of the array
maximumValue=max(max(max(images)));%Finds the maximum of the array
images = (images-minumumValue)./(maximumValue-minumumValue);%Scales the arrays
end %if scaleImage == 1;

%Interpolate frames
if frameInterpolation == 1;
[dim1,dim2,dim3]=meshgrid(1:(1/interFactorXY):size(images,1),...
    1:(1/interFactorXY):size(images,2),1:(1/interFactorFrame):size(images,3));%Generates the interpolated dridimensional space
images = interp3(images,dim1,dim2,dim3);%Imperpolates in 3d the image stack
end % interpolation == 1

%Loads the color map
try
cMap = data.cMap;
catch
cMap = jet(64);
end %Try
if colourImage == 1;
[xq,yq]=meshgrid(1:3,1:((length(cMap)-1)/(255-1)):64);%Interpolates the colour map to 255
cMap=interp2(cMap,xq,yq);%Interpolates the colour map to 255
imgAllRGB = nan(size(images,1),size(images,2),3,size(images,3));
for idxImg =1:size(images,3)
imgRGB = ind2rgb(gray2ind(images(:,:,idxImg),255),cMap);%Converts the from grey scale to colormap
imgAllRGB(:,:,:,idxImg) = imgRGB(:,:,:);
end % for
images = imgAllRGB;
else %No colour added
images = reshape(images,[size(images,1),size(images,2),1,size(images,3)]);%Adds an additional dimension
end    %colourImage == 1;

video = VideoWriter(fileName,videoFormat);
video.FrameRate = frameRate;
open(video);%Opens the video
writeVideo(video,images);%Writes the video data
close(video);%Closes the video

displayString =[name,' ',videoFormat,' video saved'];  
disp(displayString);set(handles.guiEntire,'Name',displayString);

end %if isempty(fileName)==1

% --------------------------------------------------------------------
function menuProcessing_Callback(hObject, eventdata, handles)
set(handles.uipanelCompute,'Visible','off');
set(handles.uipanelProcessing,'Visible','off');

% --------------------------------------------------------------------
function test1_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menuImportImage_Callback(hObject, eventdata, handles)
set(handles.uipanelCompute,'Visible','on');
set(handles.uipanelProcessing,'Visible','off');
set(handles.uipanelVisualisation,'Visible','off');
set(handles.uipanelVisualisationLow,'Visible','off');
set(handles.uipanelSpectra,'Visible','off');
set(handles.uipanelImage,'Visible','off');
set(handles.uipanelStack,'Visible','off');
set(handles.uipanelHistogram,'Visible','off');

% --------------------------------------------------------------------
function menuSequenceProcessing_Callback(hObject, eventdata, handles)
set(handles.uipanelCompute,'Visible','off');
set(handles.uipanelVisualisation,'Visible','off');
set(handles.uipanelVisualisationLow,'Visible','off');
set(handles.uipanelSpectra,'Visible','off');
set(handles.uipanelImage,'Visible','off');
set(handles.uipanelStack,'Visible','off');
set(handles.uipanelHistogram,'Visible','off');
set(handles.uipanelProcessing,'Visible','on');

% --------------------------------------------------------------------
function menuTruncate_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menuOffset_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menuAtrCorrection_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menuSmooth_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menuExport_Callback(hObject, eventdata, handles)


% --- Executes when guiEntire is resized.
function guiEntire_ResizeFcn(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonCancelProcessing.
function pushbuttonCancelProcessing_Callback(hObject, eventdata, handles)
global stop
stop =1;%Send the value to stop the processing
set(handles.uipanelCompute,'Visible','off');
set(handles.uipanelProcessing,'Visible','off');
set(handles.uipanelVisualisation,'Visible','on');
set(handles.uipanelVisualisationLow,'Visible','on');
set(handles.uipanelSpectra,'Visible','on');
set(handles.uipanelImage,'Visible','on');
set(handles.uipanelStack,'Visible','on');
set(handles.uipanelHistogram,'Visible','on');



% --- Executes on button press in pushbuttonCancelImportFile.
function pushbuttonCancelImportFile_Callback(hObject, eventdata, handles)
global stop
stop =1;%Send the value to stop the import
set(handles.uipanelCompute,'Visible','off');
set(handles.uipanelProcessing,'Visible','off');
set(handles.uipanelVisualisation,'Visible','on');
set(handles.uipanelVisualisationLow,'Visible','on');
set(handles.uipanelSpectra,'Visible','on');
set(handles.uipanelImage,'Visible','on');
set(handles.uipanelStack,'Visible','on');
set(handles.uipanelHistogram,'Visible','on');



function editImportImgHigh_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function editImportImgHigh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editImportImgLow_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editImportImgLow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxInterZeroOffset.
function checkboxInterZeroOffset_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkboxSingleNorm.
function checkboxSingleNorm_Callback(hObject, eventdata, handles)
handleGroup=[handles.editNormOffsetHigh,handles.editNormOffsetLow...
    handles.editNormHigh,handles.editNormLow];
if get(handles.checkboxSingleNorm,'value')==1&&get(handles.checkboxCompute,'value')==1
set(handleGroup,'enable','on');
else
set(handleGroup,'enable','off');
end %if

% --- Executes on button press in checkboxCburstCorr.
function checkboxCburstCorr_Callback(hObject, eventdata, handles)


function editNormHigh_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editNormHigh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editNormLow_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function editNormLow_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNormOffsetHigh_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editNormOffsetHigh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNormOffsetLow_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function editNormOffsetLow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editForwardTruncation_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editForwardTruncation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonClickSelectOp1.
function pushbuttonClickSelectOp1_Callback(hObject, eventdata, handles)
idxOp =1;
clickSelect(idxOp,handles);


% --- Executes on button press in pushbuttonClickSelectOp2.
function pushbuttonClickSelectOp2_Callback(hObject, eventdata, handles)
idxOp =2;
clickSelect(idxOp,handles);


% --- Executes on button press in pushbuttonClickSelectOp3.
function pushbuttonClickSelectOp3_Callback(hObject, eventdata, handles)
idxOp =3;
clickSelect(idxOp,handles);


% --- Executes on button press in pushbuttonClickSelectOp4.
function pushbuttonClickSelectOp4_Callback(hObject, eventdata, handles)
idxOp =4;
clickSelect(idxOp,handles);


% --- Executes on button press in pushbuttonClickSelectOp5.
function pushbuttonClickSelectOp5_Callback(hObject, eventdata, handles)
idxOp =5;
clickSelect(idxOp,handles);

function clickSelect(idxOp,handles)

%Selects the current operation handles
handleOp = findobj('Tag',['popupmenuOp',num2str(idxOp)]);
handleTopBox = findobj('Tag',['editOp',num2str(idxOp),'Box1']);
handleLowerBox = findobj('Tag',['editOp',num2str(idxOp),'Box2']);

operationCurrent = get(handleOp,'string');%Gets theoperation string array
operationCurrent = operationCurrent{get(handleOp,'value')};%Selected operation


if strcmp(operationCurrent,'Truncate')==1 || strcmp(operationCurrent,'Normalise')==1 ||...
    strcmp(operationCurrent,'Offset')==1 || strcmp(operationCurrent,'Zap')==1
%Two points click cases    
set(handles.uipanelSpectra,'visible','on') %Activates the panel
set(handles.uipanelImage,'visible','off') %Deactivates the panel
set(handles.uipanelStack,'visible','off') %Deactivates the panel
set(handles.uipanelHistogram,'visible','off') %Deactivates the panel
displayString =('Please select a two extremity points on the region of interest');  
disp(displayString);set(handles.guiEntire,'Name',displayString);
waitfor(msgbox(displayString,'Message'));%Message box

[cursorX,~]=ginput(2);%Two point GUI input
%cursorX=roundn(cursorX,-1);% Rounds to the first decimal digit
set(handleTopBox,'string',num2str(max(cursorX)));%Sets the GUI value
set(handleLowerBox,'string',num2str(min(cursorX)));%Sets the GUI value
displayString =[num2str(max(cursorX)),' to ',num2str(min(cursorX)),' cm-1 selected'];  
disp(displayString);set(handles.guiEntire,'Name',displayString);
set(handles.uipanelSpectra,'visible','off') %Deactivates the panel

    elseif strcmp(operationCurrent,'Multiple polygone frame offset')==1||...
        strcmp(operationCurrent,'Multiple polygone interpolated offset')
    set(handleOp,'enable','on');set(handleTopBox,'enable','on');
    set(handleLowerBox,'enable','on')

    set(handles.uipanelSpectra,'visible','off') %Activates the panel
    set(handles.uipanelImage,'visible','on') %Deactivates the panel
    set(handles.uipanelStack,'visible','off') %Deactivates the panel
    
    displayString =('Please select polygone(s) by clicking on the image before pressing ENTER');  
    disp(displayString);set(handles.guiEntire,'Name',displayString);
    waitfor(msgbox(displayString,'Message'));%Message box

    strCondition = 'Yes'; idxRegion = 0;polygone=[]; %Sets the initial parameters
    while strcmp(strCondition,'Yes')==1 && idxRegion<50
    [cursorX,cursorY] = ginput; %Selects a polygone using the dialogue
    idxRegion =  idxRegion+1; %Increments the number of polygone
    cursorX = round(cursorX); %Rounds the pixel position
    cursorY = round(cursorY); %Rounds the pixel position
    idxXPos{idxRegion}=cursorX;%Saves the data
    idxYPos{idxRegion}=cursorY;%Saves the data
    cursorX=cat(1,cursorX,cursorX(1));cursorY=cat(1,cursorY,cursorY(1));%Add the first point at the end
    cgf = handles.axesImage;%Selects the current figure
    polygoneCur=line(cursorX,cursorY);%Plots the selected polygone
    polygone=cat(1,polygone,polygoneCur);%Appends the array
    polygoneCur=line(cursorX,cursorY,'LineWidth',3,'Color',[1 1 1]);%Plots the selected polygone
    polygone=cat(1,polygone,polygoneCur);%Appends the array
    polygoneCur=line(cursorX,cursorY,'LineWidth',2,'Color',[0 0 0]);%Plots the selected polygone
    polygone=cat(1,polygone,polygoneCur);%Appends the array
    strCondition = questdlg([num2str(idxRegion),' polygone(s) selected. Do you wish to select another one ?'], ...
    'Select another polygone ?','Yes','No','Yes');
    end % While loop
    nbRegion = idxRegion;%Saves the number of regions
    set(polygone(:),'Visible','off');    %Clears the selected polygones
    %Creates the GUI string
    cursorXstring='';cursorYstring=''; %Initiates the cursor string
    for idxRegion = 1:nbRegion
        for idxPoint=1:length(idxXPos{idxRegion})
            cursorXstring = [cursorXstring,' ',num2str(idxXPos{idxRegion}(idxPoint))];%Adds the next point    
            cursorYstring = [cursorYstring,' ',num2str(idxYPos{idxRegion}(idxPoint))];%Adds the next point  
        end %for
        cursorXstring = strcat(cursorXstring,';');% Adds the delimeter of the polygone
        cursorYstring = strcat(cursorYstring,';');% Adds the delimeter of the polygone
    end %for
    cursorXstring = cursorXstring(2:end-1);%Removes the first and last caracter
    cursorYstring = cursorYstring(2:end-1);%Removes the first and last caracter

    displayString =[idxRegion,'Polygone(s) selected'];  
    disp(displayString);set(handles.guiEntire,'Name',displayString);

    set(handleLowerBox,'string',num2str(cursorXstring));%Writes the cursor string in the GUI
    set(handleTopBox,'string',num2str(cursorYstring));%Writes the cursor string in the GUI
    set(handles.uipanelImage,'visible','off') %Deactivates the panel
    
    else %No operation
    disp('No click selection needed');
end %if


% --- Executes on selection change in popupmenuStackAxisUnit.
function popupmenuStackAxisUnit_Callback(hObject, eventdata, handles)
stackXMax = 1;
stackXMin = 1;
stackYMax = 1;
stackYMin = 1;
set(handles.editStackXMax,'string',num2str(stackXMax));%Resets the GUI parameter
set(handles.editStackXMin,'string',num2str(stackXMin));%Resets the GUI parameter
set(handles.editStackYMax,'string',num2str(stackYMax));%Resets the GUI parameter
set(handles.editStackYMin,'string',num2str(stackYMin));%Resets the GUI parameter
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function popupmenuStackAxisUnit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxSigFit.
function checkboxSigFit_Callback(hObject, eventdata, handles)
integrateImageRegions(handles) %Integrates the images and refreshed the graphs

% --------------------------------------------------------------------
function menuExportFit_Callback(hObject, eventdata, handles)
global data

%Exports the FITed data as a .xls file
fileName = data.SmpFileList{1};
[currentDirectory,name,~] = fileparts(fileName);
defaultName = [currentDirectory,'\',name,'.xls'];
[fileName,pathName,~] =  uiputfile('*.xls','Save report as',defaultName);
xlswrite([pathName,fileName],data.fitReport);
displayString =[name,'.xls report saved'];  
disp(displayString);set(handles.guiEntire,'Name',displayString);



function editHistXMin_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editHistXMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editHistXMax_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editHistXMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editHistYMax_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 

% --- Executes during object creation, after setting all properties.
function editHistYMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editHistYMin_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editHistYMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuHistoSelect.
function popupmenuHistoSelect_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuHistoSelect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editBin_Callback(hObject, eventdata, handles)
integrateImageRegions(handles) %Integrates the images and refreshed the graphs


% --- Executes during object creation, after setting all properties.
function editBin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonFullScaleSpectra.
function pushbuttonFullScaleSpectra_Callback(hObject, eventdata, handles)
responseMax = 1;
responseMin = 1;
spectrumWaveMax = 1;
spectrumWaveMin = 1;
set(handles.editSpectrumResponseMax,'String',num2str(responseMax))
set(handles.editSpectrumResponseMin,'String',num2str(responseMin))
set(handles.editSpectrumWaveMax,'string',num2str(spectrumWaveMax));%Sets the tested GUI paremeter
set(handles.editSpectrumWaveMin,'string',num2str(spectrumWaveMin));%Sets the tested GUI paremeter
refreshGraph(handles);%Refresh graphs 


% --- Executes on button press in pushbuttonFullScaleImage.
function pushbuttonFullScaleImage_Callback(hObject, eventdata, handles)
imgYMax = 1;
imgYMin = 1;
imgXMax = 1;
imgXMin = 1;
set(handles.editImageYMax,'string',num2str(imgYMax));%Resets the GUI parameter
set(handles.editImageYMin,'string',num2str(imgYMin));%Resets the GUI parameter
set(handles.editImageXMax,'string',num2str(imgXMax));%Resets the GUI parameter
set(handles.editImageXMin,'string',num2str(imgXMin));%Resets the GUI parameter
refreshGraph(handles);%Refresh graphs 


% --- Executes on button press in pushbuttonFullScaleHisto.
function pushbuttonFullScaleHisto_Callback(hObject, eventdata, handles)
histYMax = 1;
histYMin = 1;
histXMax = 1;
histXMin = 1;
set(handles.editHistYMax,'string',num2str(histYMax));%Resets the GUI parameter
set(handles.editHistYMin,'string',num2str(histYMin));%Resets the GUI parameter
set(handles.editHistXMax,'string',num2str(histXMax));%Resets the GUI parameter
set(handles.editHistXMin,'string',num2str(histXMin));%Resets the GUI parameter
refreshGraph(handles);%Refresh graphs 


% --- Executes on button press in pushbuttonFullScaleStack.
function pushbuttonFullScaleStack_Callback(hObject, eventdata, handles)
stackXMax = 1;
stackXMin = 1;
stackYMax = 1;
stackYMin = 1;
set(handles.editStackXMax,'string',num2str(stackXMax));%Resets the GUI parameter
set(handles.editStackXMin,'string',num2str(stackXMin));%Resets the GUI parameter
set(handles.editStackYMax,'string',num2str(stackYMax));%Resets the GUI parameter
set(handles.editStackYMin,'string',num2str(stackYMin));%Resets the GUI parameter
refreshGraph(handles);%Refresh graphs 



function editSliceXMin_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editSliceXMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSliceXMax_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editSliceXMax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSliceYMax_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editSliceYMax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSliceYMin_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editSliceYMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBinSlice_Callback(hObject, eventdata, handles)
integrateImageRegions(handles) %Integrates the images and refreshed the graphs


% --- Executes during object creation, after setting all properties.
function editBinSlice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonFullScaleSlice.
function pushbuttonFullScaleSlice_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes on button press in checkboxInterNorm.
function checkboxInterNorm_Callback(hObject, eventdata, handles)



function editOffsetRatio_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editOffsetRatio_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editImageZMin_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 


% --- Executes during object creation, after setting all properties.
function editImageZMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editImageZMax_Callback(hObject, eventdata, handles)
refreshGraph(handles);%Refresh graphs 

% --- Executes during object creation, after setting all properties.
function editImageZMax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
