function allKin=csvreader; %
%MATLAB code to import the varian kinetic files into matlab
%Written by Maxime Boulet-Audet

%% User selected independant variables
guiGeaderOffset= 5;
guiDelimiter = ',';
guiStartRow = 6;
guiSubFirst = 0 ; %Substract the first spectrum from each series, 1=yes 0=no
guiIntegrateRange = 0 ; %Integrates the selected range, 1=yes 0=no
guiPlotFit = 1 ; 
guiReductionType = 'Range'; %Range = Integrate a range, PLS = partial least square method
guiIntLow = 1250; %Lower boundery of the integration range
guiIntHigh = 1350; % Higher bboundery of the integration range
guiFit = 0 ; %Perform the sigmoidal fit of the region 1=yes 0=no
guiNbIt = 50;
guiDefaultTime = 26.6;
dfaultFilename = 'C:\Users\mert2619\Wk\XP_FTIR_Projects\IgG solution kin\14-01-20_1524 2055  igg 1i2i2 ph4i0 70c.csv';

%% Tests the user input variables

if guiIntLow>guiIntHigh
    guiIntTemp=guiIntLow;
    guiIntLow = guiIntHigh;
    guiIntHigh = guiIntTemp; clear guiIntTemp;
end
guiNbIt = floor(guiNbIt);
if guiNbIt<1||guiNbIt>500;guiNbIt=25;end

%% Import the kinetic file names 

smpFileList = {};
[smpFileNameList,smpDirectory] = uigetfile({'*.csv'},'Select *.csv file(s)','MultiSelect', 'on'); 
% Puts only the filename if only one file was selected
if iscell(smpFileNameList)==1 %many files sselected
    for idxFile = 1:length(smpFileNameList) %Add full path
        smpFileList{idxFile} = [smpDirectory,char(smpFileNameList(1,idxFile))];
    end
    displayString =['Selected ',smpDirectory,' (',num2str(size(smpFileNameList,2)),') files'];
    disp(displayString);
elseif smpFileNameList~=0 %Only one file selected
    smpFileList = [smpDirectory,smpFileNameList];
    displayString = ['Selected ',smpFileList];
    disp(displayString);
    smpFileList = {smpFileList};%Converts to a cell
    smpFileNameList = {smpFileNameList};%Converts to a cell
else
    displayString =('No valid files selected');
    disp(displayString);
    return %Exits the function without import
end %if

%Stores the file list in the storage array
allKin.fileList=smpFileList;
allKin.smpFileNameList = smpFileNameList;

%Creates the storage array
allKin.spectraKin = cell(1,length(smpFileList));

%% Imports the kinetic data

%Loops throught the file list
for idxFile = 1:length(smpFileList);
    curfilename = smpFileList{idxFile};
    
    % Tries to determine the meta data from the filename
    strPrefix = ' ph';    strSulfix = ' ';
    curStr = findstrmidle(curfilename,strPrefix,strSulfix);
    %Subsitutes a character for the period
    if strcmp(curStr,'')==1;curStr='0';end
    curStr = strrep(curStr, 'i','.');curStr = strrep(curStr, '|','.');
    allKin.pH{idxFile}=str2double(curStr);
         
    strPrefix = ' tz';strSulfix = ' ';
    curStr = findstrmidle(curfilename,strPrefix,strSulfix);    
    %Subsitutes a character for the period
    if strcmp(curStr,'')==1;curStr=num2str(guiDefaultTime);end
    curStr = strrep(curStr, 'i','.');curStr = strrep(curStr, '|','.');
    allKin.timeInterval{idxFile}=str2double(curStr);
    
    strPrefix = ' nacl';strSulfix = ' ';
    curStr = findstrmidle(curfilename,strPrefix,strSulfix);    
    %Subsitutes a character for the period
    curStr = strrep(curStr, 'i','.');curStr = strrep(curStr, '|','.');
    allKin.naCl{idxFile}=str2double(curStr);
        
    %
    fid=fopen(curfilename,'r');

    %Determines the lenght and width of the file and imports
    nbLines = 0;
    tline = fgetl(fid);
    while ischar(tline)
      tline = fgetl(fid);
      nbLines = nbLines+1;
          if nbLines == guiGeaderOffset
            nbCol = sum(tline==guiDelimiter)+1;
          end %if
    end
    
    fclose(fid); %Closes the file
    
    %Creates an array of the file size
    allSpectra=nan(nbLines-guiGeaderOffset,nbCol);

    %Creates the format spec
    formatSpec = [];
    for idxCol = 1:nbCol
    formatSpec = [formatSpec,'%f'];
    end

    %Re-opens the file
    fid=fopen(curfilename,'r');
    dataArray = textscan(fid, formatSpec, 'Delimiter', guiDelimiter,...
        'EmptyValue' ,NaN,'HeaderLines' ,guiStartRow-1, 'ReturnOnError', false);
    dataArray=cell2mat(dataArray(1,1:length(dataArray)));
    
    fclose(fid); %Closes the file

     waveVector = dataArray(1:end,1);
     allSpectra =  dataArray(:,2:end);
     allSpectra=flipdim(allSpectra,2); % Flips the array
     %Stores the data in the storage array
     allKin.spectraKin{idxFile} = allSpectra;
     allKin.waveVec{idxFile} = waveVector;
     %Stores the rime interval as a vector
     timeInterval = allKin.timeInterval{idxFile};
     if isnan(timeInterval)==1; timeInterval=guiDefaultTime;end
     timeVec = ((1:size(allSpectra,2))-1)';
     timeVec = timeVec.*timeInterval;
     allKin.timeVec{idxFile} = timeVec;
     
end % Looping throught the files

%% Substracts the first spectrum from the kinetics

if strcmpi(guiReductionType,'Range')==1

%% Integrates the selected range

    %Loops throught the files
    allKin.reduced = cell(1,length(smpFileList));
    for idxFile = 1:length(smpFileList)
        % Locates the index of the bounderies
        [~,idxIntLow]=min(abs(allKin.waveVec{idxFile}-guiIntLow));
        [~,idxIntHigh]=min(abs(allKin.waveVec{idxFile}-guiIntHigh));
        %Gets the average of the spectral range
        allKin.reduced{idxFile} = mean(allKin.spectraKin{idxFile}(idxIntLow:idxIntHigh,:),1);
    end %Loops throught the files

elseif strcmpi(guiReductionType,'PLS')==1
%% Performs a Partial Least Square (PLS) on the data set to determine the concentration

    %Loops throught the files
    allKin.reduced = cell(1,length(smpFileList));
    for idxFile = 1:length(smpFileList)
        
        % Locates the index of the bounderies
        [~,idxIntLow]=min(abs(allKin.waveVec{idxFile}-guiIntLow));
        [~,idxIntHigh]=min(abs(allKin.waveVec{idxFile}-guiIntHigh));
        spectraToEstimate = allKin.spectraKin{idxFile}(idxIntLow:idxIntHigh,:)';
        
        ncomp = 1 ; %The number of component used n - number of training observations
        %save('plsTrainingIgGConcentration','trainSpectra','trainConc','trainWave');
        trainSpectra =nan; % contains the training spectra, one obsevation per row
        trainConc = nan; %  contains the training concetrations, one obsevation per row
        trainWave = nan; % constains the training spectra wavenumber vactor
        load('plsTrainingIgGConcentration') %Loads the PLS training spectra and concentrations.
        
        % Locates the index of the bounderies
        [~,idxIntLow]=min(abs(trainWave-guiIntLow));
        [~,idxIntHigh]=min(abs(trainWave-guiIntHigh));
        trainSpecRange = trainSpectra(:,idxIntLow:idxIntHigh);
        
        if size(trainSpecRange,2)~=size(spectraToEstimate,2)
        disp('The training matrix and the estimated spectra do not have the same step, interpolation would be needed');
        return %exits the script
        end %if condition
        
        % Returns the PLS regression coefficients BETA.
        % BETA isa (p+1)-by-m matrix, containing intercept terms in the first row:
        % MSE The first row of MSE contains mean-squared errorsfor the predictor variables in X, and the secondrow contains mean-squared errors for the response variable(s) in Y.
        [predictorLoad,responseLoad,predictorScores,responseScores,beta,PCTVAR,MSE] = plsregress(trainSpecRange,trainConc,ncomp);
        clear predictorLoad responseLoad predictorScores responseScores PCTVAR
        %Loads the kinetic data

        %Estimates the concentration based on the beta matrix
         responseEstimated = [ones(size(spectraToEstimate,1),1),spectraToEstimate]*beta;
         %Saves the PLS data in the storage array
         allKin.reduced{idxFile} = responseEstimated';
    end %Loops throught the files

end %If guiReductionType
 
%% Fits the kinetics with a custom function

if guiFit ==1
    
    %Calculate the colour scheme for the data points
    cPoint = 0.01:(0.3/(length(smpFileList)-1)):0.31;
    %Calculates the colour scheme for the fitted curves
    cCurve = 0.51:(0.3/(length(smpFileList)-1)):0.81;
    
    for idxFile = 1:length(smpFileList)
        displayString = ['Fitting file ', num2str(idxFile),'/',num2str(length(smpFileList))];
        disp(displayString)
        
        x = allKin.timeVec{idxFile};
        y = allKin.reduced{idxFile}';

        %Defines the fitting function
        fitModel = fittype('Aoff+(A0-(((K1./K2)+A0)./(1+((K1./(K2*A0))*exp(x.*(K1+(K2*A0)))))))',...
            'coefficients',{'A0','Aoff','K1','K2'},'independent','x');

        %Defines the fitting options
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0,1E-12,1E-12],...
            'Upper',[2.*max(y),2.*max(y),100000.*max(x),100000.*max(x)],...
            'Startpoint',[0,max(y)-min(y),1./max(x),1000./max(x)],...
            'Robust', 'LAR',... %or 'On' select 'LAR'
            'Algorithm','Trust-Region',... %
            'MaxIter', guiNbIt,...
            'DiffMinChange',1e-15.*max(y),...
            'DiffMaxChange',1e15.*max(y),...
            'Display','Off'); %or 'Notify'; 
        %Fits the data of the radial bin
        
        try
        %Fits the x and y
            [fitedValues,gof,fitPar] = fit(x,y,fitModel,fitopt);
            arrCoeffs = coeffvalues(fitedValues);
            A0= arrCoeffs(1);allKin.a0{idxFile}=A0;
            Aoff= arrCoeffs(2);allKin.Aoff{idxFile}=Aoff;
            K1= arrCoeffs(3);allKin.k1{idxFile}=K1;
            K2= arrCoeffs(4);allKin.k2{idxFile}=K2;
            yFitted = Aoff+(A0-(((K1./K2)+A0)./(1+((K1./(K2*A0))*exp(x.*(K1+(K2*A0)))))));
            arrCoeffsInterval = confint(fitedValues,0.95);
            allKin.r2{idxFile}=gof.rsquare;

            %Calculates the t1/2
            tHalf= (log((K2*A0)/K1))/(K1+(K2*A0));
            allKin. tHalf{idxFile}=tHalf;
                        
            if guiPlotFit == 1
            %Plots the result
            if idxFile==1;
            %    hold off;
            hFig = figure('Color',[1 1 1]);
            end %if
            
            hdle= line(x,y);
            box('on');
            set(hdle,'Marker','.','MarkerSize',10,'LineStyle','none','Color',[cPoint(idxFile),cPoint(idxFile),cPoint(idxFile)]);
             %Draws the fitting lines
            hdle(1) = line(x,yFitted);
            hdle(2) = line([tHalf,tHalf],[Aoff+A0,Aoff-(K1/K2)]);
            hdle(3) = line([tHalf,max(x)],[A0+Aoff,A0+Aoff]);
            hdle(4) = line([tHalf,min(x)],[Aoff-(K1/K2),Aoff-(K1/K2)]);
            hdle(5) = line([tHalf-(2/(K2*A0)),tHalf+(2/(K2*A0))],[Aoff-(K1/K2),A0+Aoff]);
            set(hdle,'LineWidth',3,'LineStyle','-','Color',[cCurve(idxFile),cCurve(idxFile),cCurve(idxFile)]);
            
            allAxes = findall(hFig,'Type','axes');
            set(allAxes,'XMinorTick','on','TickDir','out','YMinorTick','on','TickDir','out');
            
            title([smpFileNameList{idxFile},...
                ' | R2 = ',num2str(roundn(gof.rsquare,-3)),...
                ' | tHalf = ',num2str(roundn(tHalf,0)),...
                ' | K2 = ',num2str(roundn((K2),-5))]);
            xlabel('Time (s)');% Creates ylabel
            ylabel('Concentration fraction');%Creates xLabel
            pause(1);
            hold on
            %close(hFig);
            end % If
            
        catch lasterr
            displayString = ['Fit error! ' ,lasterr.message];
            disp(displayString);
        end %Catch

    end % For through the files
end % if entry condition

 disp('Import succesful');
 
    function strFound = findstrmidle(strSearch,strPrefix,strSulfix)
    
    %Returns the string betweeen the last instance of the prefix and a
    %sulfix
    
    idxPrefix = strfind(strSearch,strPrefix);
    if isempty(idxPrefix)==1
        strFound= '';return
    else
        idxPrefix = idxPrefix(end);
    end
        
    idxSulfix = strfind(strSearch,strSulfix);
    if isempty(idxSulfix)==1
    idxSulfix = length(strSearch);
    end
    
    idxSulfix(idxSulfix<=(idxPrefix+length(strPrefix))) = nan;
    idxSulfix = min(idxSulfix);
    if isnan(idxSulfix)==1
    idxSulfix = length(strSearch);
    end    
    
    strFound  = strSearch((idxPrefix+length(strPrefix)):(idxSulfix-length(strSulfix)));

 
            
            