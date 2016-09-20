function out = tableIrAnalysis(tbl,idList)
%% Inputs the user parameters
plsModelIdx = 2;
offsetRange = [1780,1820];
baselineRange = [1780,1820;840,880];%Leave empty to avoid
derDegree = 1;
absVal = 0; 
rangePLS = [1490,1610;900,1100]; %[1610,1670];%; Empty for skip PLS analysis
nbPC = 2;
cycleAn = 1; %Analysis of the cycles
cycleNbFrame = 5;%Number of frames to average
markTimeOffset = 1;
waterVapourSubIdx = 3; % 194; % Use 0 for no sub
waterVapourSubRg =[1900,1800];
subRange = [];%[40,45];%[60,62]; % Use [] for no subtraction
spcRefIdx = []; % Index of the reference spectrum    
rangeCog = [1670,1615;1580,1402;1312,1285];
rangeInt = [1670,1615;1580,1402];  
movingAverage = 1;
    
%% Stored the user parameters in a structure variable
par.plsModelIdx = plsModelIdx;
par.offsetRange = offsetRange;
par.baselineRange = baselineRange;
par.derDegree = derDegree;
par.absVal = absVal ; 
par.rangePLS = rangePLS;
par.nbPC = nbPC;
par.waterVapourSubIdx = waterVapourSubIdx;
par.waterVapourSubRg = waterVapourSubRg;
par.subRange = subRange;
par.rangeCOG = rangeCog;
par.rangeInt = rangeInt;  
par.movingAverage = movingAverage;
par.markTimeOffset = markTimeOffset;
%% Test inputs
tic %Starts the timer

tblName = inputname(1);
if ~istable(tbl)
    disp([tblName,' is not a table'])
    return
end     

if nargin < 2 
    idList = [1,size(tbl,1)]; %Analyse the entire table
    disp(['Analysing table ',tblName,'...Please wait'])
end
idList =abs(round(min(idList))):1:abs(round(max(idList)));
if min(abs(idList-plsModelIdx))~=0
    idList = [plsModelIdx,idList];
end % if
    
    %% Check if the table fields exists

try idx = iscell(tbl.RefAutoSubkSol);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'RefAutoSubkSol';
end %try

try idx = iscell(tbl.AbsorbanceRaw);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'AbsorbanceRaw';
end %try

try idx = iscell(tbl.AbsorbanceFinal);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'AbsorbanceFinal';
end %try

try idx = iscell(tbl.RsquaredPLS);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'RsquaredPLS';
end %try

try idx = iscell(tbl.rangePLSInt);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'rangePLSInt';
end %try

try idx = iscell(tbl.AbsorbanceDiff);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'AbsorbanceDiff';
end %try

try idx = iscell(tbl.FTIR1int);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'FTIR1int';
end %try

try idx = iscell(tbl.ftirIntRg);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'ftirIntRg';
end %try

try idx = iscell(tbl.par);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'par';
end %try

try idx = iscell(tbl.cog);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'cog';
end %try 

try idx = iscell(tbl.int);  
catch
    idx = size(tbl,2)+1; tbl.(idx){1} = [];%Adds a field
    tbl.Properties.VariableNames{idx} = 'int';
end %try 

%% Spectra pre-treatment          
for idListIdx = 1:length(idList) %Loops througth the file id
    idIdx = idList(idListIdx);
    % Loads parameters
    xWave = tbl.WaveNb{idIdx};
    absorbance = tbl.Absorbance{idIdx};        
    xTimeIR = tbl.TimeIR{idIdx};
    if isempty(xTimeIR)
        xTimeIR = zeros(1,size(absorbance,2));
    end %if

    % Offset
    try
        offsetRangeIdx=nan(size(offsetRange));
        [~,offsetRangeIdx(1)] = min(abs(xWave-min(offsetRange)));
        [~,offsetRangeIdx(2)] =  min(abs(xWave-max(offsetRange)));
        absorbanceOffset = repmat(mean(absorbance(min(offsetRangeIdx):max(offsetRangeIdx),:),1),[size(absorbance,1) 1]);
        absorbance = absorbance-absorbanceOffset;
    catch error
        disp(['RunID ',num2str(idIdx),' RunID | Offset failed! ',error.message,' at line ',num2str(error.stack.line)]);
    end % try

    % Baseline subtraction
    if ~isempty(par.baselineRange);
        try
            temp = repmat(mean(par.baselineRange,2),1,size(absorbance,2));
            baselinePoint1X = temp(1,:);
            baselinePoint2X = temp(2,:);
            [~,rangeIdx(1)]=min(abs(xWave-min(baselineRange(1,:))));
            [~,rangeIdx(2)]=min(abs(xWave-max(baselineRange(1,:))));
            baselinePoint1Y = mean(absorbance(min(rangeIdx):max(rangeIdx),:),1);
            [~,rangeIdx(1)]=min(abs(xWave-min(baselineRange(2,:))));
            [~,rangeIdx(2)]=min(abs(xWave-max(baselineRange(2,:))));
            baselinePoint2Y = mean(absorbance(min(rangeIdx):max(rangeIdx),:),1);
            baselineSlope = (baselinePoint2Y - baselinePoint1Y)./(baselinePoint2X-baselinePoint1X);%Caluates a linear slope a
            baselineAsympt = baselinePoint1Y - (baselineSlope.*baselinePoint1X);%Calulate the asymptote b
            baselineMat = (repmat(xWave,1,size(baselineSlope,2))).*(repmat(baselineSlope,size(xWave,1),1)) + repmat(baselineAsympt,size(xWave),1);
            absorbance = absorbance - baselineMat; %Subtracts the baseline
          catch error
            disp(['RunID ',num2str(idIdx),' Baseline failed! ',error.message,' at line ',num2str(error.stack.line)]);
        end %catch
    end

    % Frame subtration
    if ~isempty(subRange); % Use 0 for no subtraction
        try
            [~,subRangeIdx(1)]=min(abs(xTimeIR-min(subRange(1,:))));%Finds the index of the selected range
            [~,subRangeIdx(2)]=min(abs(xTimeIR-max(subRange(1,:))));%Finds the index of the selected range
            absorbanceOffset = mean(absorbance(:,min(subRangeIdx):max(subRangeIdx)),2);%Calculates the average spectrum of the X range
            absorbanceOffset = repmat(absorbanceOffset,[1 size(absorbance,2)]);
            absorbance = absorbance-absorbanceOffset;%Subtract the average spectrum
        catch
            disp(['RunID ',num2str(idIdx),' Frame substraction failed! ',error.message,' at line ',num2str(error.stack.line)]);    
        end %catch
    end %if 

    % Water vapour auto substraction
    if ~isempty(waterVapourSubIdx);
        try
            sRef = tbl.Absorbance{waterVapourSubIdx};
            [~,rangeIdx(1)]=min(abs(xWave-min(waterVapourSubRg)));
            [~,rangeIdx(2)]=min(abs(xWave-max(waterVapourSubRg)));
            rangeIdx = sort(rangeIdx);
            sRefCrop = sRef(min(rangeIdx):max(rangeIdx));
            for vecIdx = 1:size(absorbance,2) %Loop through the spectra
                sInCrop = absorbance(min(rangeIdx):max(rangeIdx),vecIdx);
                [~,kSol] = specAutoSub(sInCrop,sRefCrop);%Autosubtract
                sIn  = absorbance(:,vecIdx) - (kSol.*sRef);
                absorbance(:,vecIdx) = sIn;
            end %for        
        catch error
            disp(['RunID ',num2str(idIdx),' Water subtraction failed! ',error.message,' at line ',num2str(error.stack.line)]);
        end %try
    end %if

    % Reference spectrum autosubtraction
    if ~isempty(spcRefIdx);
        tbl.refAutoSubkSol = nan(size(absorbance,2),1);
        try
            sRef = tbl.Absorbance{spcRefIdx};
            sRef = squeeze(mean(sRef,2));%Averages the first dimension
            for idx = 1:size(absorbance,2) %Loop through the spectra
            sIn = absorbance(:,idx);  
            [sIn,kSol] = specAutoSub(sIn,sRef);%Autosubtract
            refAutoSubkSol(idx)= kSol;%Saves the subtarction factor
            absorbance(:,idx) = sIn;
            end %for
            tbl.refAutoSubkSol = refAutoSubkSol;
            tbl.RefAutoSubkSol{idIdx}=refAutoSubkSol;
        catch error
            disp(['RunID ',num2str(idIdx),' Ref subtraction failed! ',error.message,' at line ',num2str(error.stack.line)]);
        end %try            
    end %if

    % Moving average smoothing
    if movingAverage > 1 && ~isempty(movingAverage);
        try
            absorbanceSmooth = nan(size(absorbance));
            for waveIdx=1:size(absorbance,1)
                absorbanceSmooth(waveIdx,:) = smooth(absorbance(waveIdx,:),movingAverage,'moving');
            end
            absorbance = absorbanceSmooth; %Copy the matrix
        catch error
            disp(['RunID ',num2str(idIdx),' Moving avg failed! ',error.message,' at line ',num2str(error.stack.line)]);  
        end           
    end %if par.movingAverage > 1;

    % Derivative
    if derDegree ~=0 && ~isempty(derDegree)
        try
            absorbance=diff(absorbance,derDegree,1);
            absorbance= [absorbance;zeros(par.derDegree,size(absorbance,2))];%Patches missing end values
        catch error
            disp(['RunID ',num2str(idIdx),' Derivative failed ! ',error.message,' at line ',num2str(error.stack.line)]);  
        end            
    end %if

    % Make spectra absolute
    if par.absVal == 1;
        absorbance = abs(absorbance);            
    end % if par.absVal = 1
    tbl.AbsorbanceFinal{idIdx}=absorbance;% Saves the pre-treated spectra
    tbl.par{idIdx} = par; %Saves the parameters

end %For loop through the ID

    %% Multivariate and monovariate analysis
for idListIdx = 1:length(idList) %Loops througth the run id
    idIdx = idList(idListIdx);
    xWave = tbl.WaveNb{idIdx};
    absorbance = tbl.AbsorbanceFinal{idIdx};
    if isempty(absorbance)
        absorbance = tbl.Absorbance{idIdx};
    end %if
    
    % PLS Multivariate analysis
    if ~isempty(par.rangePLS) && ~isempty(tbl.Response{par.plsModelIdx})
        rangePLSIdx=nan(2,size(rangePLS,2));
        plsModelIdx = par.plsModelIdx;
        for idx = 1:size(rangePLS,2)
            try
                [~,rangePLSIdx(idx,1)]=min(abs(xWave-min(rangePLS(idx,:))));%Finds the lower bound index
                [~,rangePLSIdx(idx,2)]=min(abs(xWave-max(rangePLS(idx,:))));%Finds the upper bound index
            catch error
                disp(['RunID ',num2str(idIdx),' Invalid PLS Range! ',error.message,' at line ',num2str(error.stack.line)])
            end %try

        end %for

        try
            for responseIdx= 1:size(tbl.Response{plsModelIdx},1)
                absorbance = tbl.AbsorbanceFinal{idIdx};
                if isempty(absorbance)
                    absorbance = tbl.Absorbance{idIdx};
                end %if
                xPredictor = tbl.AbsorbanceFinal{plsModelIdx};
                xPredictor = xPredictor(min(rangePLSIdx(responseIdx,:)):max(rangePLSIdx(responseIdx,:)),:)';
                yResponse = tbl.Response{plsModelIdx}';
                yResponse  = yResponse(:,responseIdx);
                xQuant = tbl.AbsorbanceFinal{idIdx};
                xQuant = xQuant(min(rangePLSIdx(responseIdx,:)):max(rangePLSIdx(responseIdx,:)),:)';
                [xLoadings,yLoadings,xScores,yScores,bPLS,PLSPctVar] = plsregress(xPredictor,yResponse,nbPC);% Models the dataset
                yfitPLSPredictor = [ones(size(xPredictor,1),1), xPredictor]*bPLS;%Fits the model data sets
                %figure;scatter(response,yfitPLS);line([0 max(response)],[0 max(response)]);
                TSS = sum((yResponse-mean(yResponse)).^2);
                RSS_PLS = sum((yResponse-yfitPLSPredictor).^2);
                if responseIdx == 1
                    tbl.RsquaredPLS{idIdx} = nan(1,size(tbl.Response{plsModelIdx},1)); % Allocation
                end %if
                tbl.RsquaredPLS{idIdx}(responseIdx)= 1 - RSS_PLS/TSS;% Calculates the R2
                yfitPLSQuant = [ones(size(xQuant,1),1),xQuant]*bPLS;%Fits the model data sets
                if responseIdx == 1
                    tbl.rangePLSInt{idIdx} = nan(size(tbl.Response{plsModelIdx},1),size(absorbance,2)); % Allocation
                end %if                
                tbl.rangePLSInt{idIdx}(responseIdx,:)=yfitPLSQuant;
            end%idx= 1:size(responsePLSModel,1)                
        catch error
            disp(['RunID ',num2str(idIdx),' PLS modeling failed! ',error.message,' at line ',num2str(error.stack.line)])
        end %try            
    end %if

    %Calculated the Centre of Gravity (COG)
    if ~isempty(rangeCog);
        try
            rangeCogIdx = nan(size(rangeCog));%Create the index array
            cogVal = nan(size(rangeCog,1),size(absorbance,2));%Create the index array
            for rangeIdx = 1:size(rangeCog,1) % Loops throught the ranges
                [~,rangeCogIdx(rangeIdx,1)]=min(abs(xWave-min(rangeCog(rangeIdx,:))));%Finds the lower bound index
                [~,rangeCogIdx(rangeIdx,2)]=min(abs(xWave-max(rangeCog(rangeIdx,:))));%Finds the upper bound index
                for vecIdx = 1:size(absorbance,2) % Loops throught the vectors
                    cogAbs = absorbance(min(rangeCogIdx(rangeIdx,:)):max(rangeCogIdx(rangeIdx,:)),:);%Gets the absorbance range
                    cogWave = xWave(min(rangeCogIdx(rangeIdx,:)):max(rangeCogIdx(rangeIdx,:)),:);% Gets the xRange
                    if size(cogWave,2) ~= size(cogAbs,2)
                        cogWave  = repmat(cogWave(:,1),[1 size(cogAbs,2)]);
                    end %if
                    cogVal(rangeIdx,:) = sum(cogAbs.*cogWave,1)./sum(cogAbs,1);%Calculates the COG
                end %for
            end %for
                tbl.cog{idIdx} = cogVal;%Saves the analysis data
        catch error
            disp(['RunID ',num2str(idIdx),' Failed COG | ',error.message])
        end %try
    end %if

    if ~isempty(rangeInt);
        try 
            rangeIntIdx = nan(size(rangeInt));%Create the index array
            intVal = nan(size(rangeInt,1),size(absorbance,2));%Create the index array
            for rangeIdx = 1:size(rangeInt,1) % Loops throught the ranges
                [~,rangeIntIdx(rangeIdx,1)]=min(abs(xWave-min(rangeInt(rangeIdx,:))));%Finds the lower bound index
                [~,rangeIntIdx(rangeIdx,2)]=min(abs(xWave-max(rangeInt(rangeIdx,:))));%Finds the upper bound index
                for vecIdx = 1:size(absorbance,2) % Loops throught the vectors
                    intAbs = absorbance(min(rangeIntIdx(rangeIdx,:)):max(rangeIntIdx(rangeIdx,:)),:);%Gets the absorbance range
                    intWave = xWave(min(rangeIntIdx(rangeIdx,:)):max(rangeIntIdx(rangeIdx,:)),:);% Gets the xRange
                    if size(intWave,2) ~= size(intAbs,2)
                        intWave  = repmat(intWave(:,1),[1 size(intAbs,2)]);
                    end %if
                    intVal(rangeIdx,:) = mean(intAbs);%Calculates the Integral
                end %for
            end %for
            tbl.int{idIdx} = cogVal;%Saves the analysis data
        catch error
            disp(['RunID ',num2str(idIdx),' Failed Integrate | ',error.message])
        end %try
    end %if
    
end %for idIdx = min(idList):max(idList)

%% Cycle analysis
if cycleAn == 1 %Cycle analysis condition
    try
        for idListIdx = 1:length(idList) %Loops througth the run id
            idIdx = idList(idListIdx);
            markArr = tbl.Mark{idIdx};
            tFTIR = tbl.TimeIR{idIdx};
            if isempty(markArr) || isempty(tFTIR)
                continue %skips the run id
            end
            tFTIR = tFTIR./60;
            markMinutes = tbl.MarkMin{idIdx};
            markMinutes = markMinutes - markTimeOffset;
            markInj = strcmp(markArr,'Injection Valve Inj');
            markPos4 = strcmp(markArr,'Buffer Valve Pos 4');
            markComp = markInj + markPos4;
            markArrSel = cell(sum(markComp)+1,4); 
            markArrSel{1,1} = '0';
            markComp = markMinutes .* markComp;
            cellFilIdx = 1;
            for cellIdx = 1:length(markComp)
                if markComp(cellIdx)>0
                    markArrSel{cellFilIdx,2} = markComp(cellIdx);
                    markArrSel{cellFilIdx+1,1} = markArr(cellIdx);
                    cellFilIdx = cellFilIdx + 1;
                end %if
            end
            markArrSel{size(markArrSel,1),2} = max(tFTIR);%Add a cycles at the end
            markPlsDat = tbl.rangePLSInt{idIdx};
            markInjIdx = 1;
            markCipIdx = 1;
            for markIdx = 1:size(markArrSel,1)%Runs through the marks
                 [~,markRangePLSIdx] = min(abs(tFTIR-markArrSel{markIdx,2}));%Finds the lower bound index
                 markRangeIdxMin = markRangePLSIdx - cycleNbFrame +1;
                 markRangeIdxMax = markRangePLSIdx + 0;
                 if (markRangePLSIdx - cycleNbFrame) < 1 
                    markArrSel{markIdx,3} = markPlsDat(1,markRangeIdxMax);%Integrates the FTIR reponse
                    markArrSel{markIdx,4} = 0;%Integrates the FTIR reponse
                 else
                    markArrSel{markIdx,3} = mean(markPlsDat(1,markRangeIdxMin:markRangeIdxMax),2);%Integrates the FTIR reponse
                    markArrSel{markIdx,4} = std(markPlsDat(1,markRangeIdxMin:markRangeIdxMax),0,2);%Integrates the FTIR reponse
                 end %if
                 if strcmp(markArrSel{markIdx,1},'Injection Valve Inj') ==1;
                    markArrSel{markIdx,1} = ['Cyl',num2str(markInjIdx)];
                    markInjIdx =  markInjIdx +1;
                 elseif strcmp(markArrSel{markIdx,1},'Buffer Valve Pos 4')==1;
                    markArrSel{markIdx,1} = ['CIP',num2str(markCipIdx)];
                    markCipIdx = markCipIdx +1;
                 end %if             
            end %for
            tbl.FTIR1int{idIdx} = markArrSel;%Saves the result back to the table        
        end %for
    catch error
            disp(['RunID ',num2str(idIdx),'Failed cycle analysis | ',error.message])
    end %try
end %if

%% Exports the results
out = tbl; %Exports the table    
disp(['Table ',tblName,' saved after ', num2str(toc * 1000,'%.1f'),' ms']);
   
end %function

function [out1,out2] = v1lowerthan2(var1,var2)
   if var1 > var2
       tmp = var2;var2 = var1;var1 = tmp;
   end
   out1 = var1;out2 = var2;   
end %fct


