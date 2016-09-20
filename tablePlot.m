function figureHandles = tablePlot(tbl,idIdx)
%% Fixed parameters
closesPreviousFigures = 1 ;%Close previous graphs
plotSurfBod = 0; % 1 plots a surface plot, 0 does not
plotSpcBod = 1; % 1 plots the spectrum plot, 0 does not
xAxisRangeIntRange = [400,500];
xAxisRangeIntStepsNB = 2;
xAxisRangeIntWidth = 1;
yAxisRangeSpcLim = [-0.001,0.002];
xAxisRangeSpcLim = [1900,900];
yAxisSpcOffset = 0.0000;
plotProfileBod = 1 ;
xAxisRange = [];
yAxisRangeUV = [];
yAxisRangeCond = [];
yAxisRangeFTIR = [-1.0,6.0];
plotDifBod = 0; % 1 plots the spectrum plot, 0 does not
yAxisDiffOffset = 0.0002;
baseType = 'Minutes'; %baseType = 'Minutes'; %baseType = 'Seconds';
figureHorPosition = 1900;
cropRange = [900,1900];
plotCyles = 1 ; %Plots the cycles int values

%% Test inputs
if ~istable(tbl)
    disp('Input must be a table')
    return
end

if ~exist('idIdx','var') || ~isnumeric(idIdx);
    idIdx = size(tbl,1);
    disp('No ID selected, plotting the last entry')
end
idIdx = abs(round(idIdx));
tblName = inputname(1);
tic %Start timer

%% Loads parameters

id = tbl.RunID{idIdx};
tFPLC = tbl.TimeFPLC{idIdx};
FlowRate = tbl.FlowRate{idIdx};
uv280 = tbl.UV280{idIdx};
cond = tbl.Cond{idIdx};
pFTIR = tbl.rangePLSInt{idIdx};
tFTIR = tbl.TimeIR{idIdx};
xWave = tbl.WaveNb{idIdx};
absorbance = tbl.AbsorbanceFinal{idIdx};
bkg = tbl.Bkg{idIdx};

%% Converts time axes

try
    %Converts IR time from absolute seconds into relative minutes  
    %timeFPLC = (timeFPLC - repmat(timeFPLC(1),[length(timeFPLC) 1]));
catch; disp('no Time FPLC')
end %if

if length(tFTIR)~= size(absorbance,2)
    tFTIR= tFTIR(1):((tFTIR(end)-tFTIR(1))/(size(absorbance,2)-1)):tFTIR(end);    
end %if lenght(timeIR)==1

if size(tFTIR,2)>size(tFTIR,1)
    tFTIR = tFTIR';%Transpose
end
%Converts the x-axis unit
if strcmp(baseType,'Seconds') == 1;
    tFPLC = tFPLC.*60;   
elseif strcmp(baseType,'mL') == 1;
    tFTIR = (tFTIR - tFTIR(1))./60;%Converts IR time from absolute seconds into relative minutes  
    tFTIR = tFTIR.*FlowRate;
    tFPLC = tFPLC.*FlowRate;
else %Minutes
    tFTIR = (tFTIR - tFTIR(1))./60;%Converts IR time from absolute seconds into relative minutes  
end %if

if isempty(xAxisRange)
    xAxisRange = [min([min(tFPLC),min(tFTIR)]),max([max(tFPLC),max(tFTIR)])];
end

if isempty(yAxisRangeCond)
    yAxisRangeCond = [min(cond)-(0.05*(max(cond)-min(cond))),max(cond)+(0.05*(max(cond)-min(cond)))];
end %if

if isempty(yAxisRangeUV)
    yAxisRangeUV = [min(uv280)-(0.05*(max(uv280)-min(uv280))),max(uv280)+(0.05*(max(uv280)-min(uv280)))];
end %if

if closesPreviousFigures == 1
    close all % Close all previous figures
end %if

%% Crop FTIR spectra
try
    cropRangeIdx = nan(size(cropRange));
    [~,cropRangeIdx(1)]=min(abs(xWave-min(cropRange)));
    [~,cropRangeIdx(2)]=min(abs(xWave-max(cropRange)));
    absorbance = absorbance(min(cropRangeIdx):max(cropRangeIdx),:);
    xWave = xWave(min(cropRangeIdx):max(cropRangeIdx));
catch error
    disp(['No FTIR !',error.message]);
end %Try

%% Generates the surface plot
if plotSurfBod ==1
    try
    figureSurf = figure('Position',[100 20 800 800],'Color',[1 1 1],'Name',id,'NumberTitle','off');
    colormap('hot');
    axes1 = axes('Parent',figureSurf);
    hold(axes1,'on');
    surf(tFTIR,xWave,absorbance,'Parent',axes1,'LineStyle','none');
    view(axes1,[-75 35]); % view(axes1,[-77.5 35.6]);
    grid(axes1,'on');
    axis(axes1,'ij');
    title(axes1,id);
    xlabel('Time (min)');% Create xlabel
    ylabel('Wavenumber (1/cm)');% Create ylabel
    zlabel('Absorbance (O.D.)');% Create zlabel    
    catch error
        disp(['No valid data to plot',error.message])
    end
end %if

%% Generates a spectrum plot
if plotSpcBod == 1
    try
        fSpc = figure('Position',[figureHorPosition+1400 500 500 500],'Color',[1 1 1],'Name',id,'NumberTitle','off');
        figureHandles.fSpc = fSpc;
        axesSpc = axes('Parent',fSpc);
        hold(axesSpc,'on');
        title(axesSpc,id);
        xlabel('Wavenumber (1/cm)');% Create xlabel
        ylabel('Absorbance (O.D.)');% Create ylabel
        set(axesSpc,'FontSize',12,'LineWidth',2,'TickDir','out','XDir','reverse');
        ylim(axesSpc,[min(yAxisRangeSpcLim) max(yAxisRangeSpcLim)]); 
        xlim(axesSpc,[min(xAxisRangeSpcLim) max(xAxisRangeSpcLim)]);
        absorbanceInt = nan(size(xWave,1),xAxisRangeIntStepsNB);
        sizeIdx = xAxisRangeIntStepsNB;
        xAxisRangeIntStep = (max(xAxisRangeIntRange) - min(xAxisRangeIntRange))/(xAxisRangeIntStepsNB-1);
        xAxisRangeInt = [(((1:xAxisRangeIntStepsNB)-1)*xAxisRangeIntStep)+min(xAxisRangeIntRange)-(0.5*xAxisRangeIntWidth);(((1:xAxisRangeIntStepsNB)-1)*xAxisRangeIntStep)+min(xAxisRangeIntRange)+(0.5*xAxisRangeIntWidth)]';
        for rgIdx  = 1: xAxisRangeIntStepsNB %Loops throught the range
            [~,xAxisRangeIntIdx(1)]=min(abs(tFTIR-min(xAxisRangeInt(rgIdx,:))));
            [~,xAxisRangeIntIdx(2)]=min(abs(tFTIR-max(xAxisRangeInt(rgIdx,:))));
            absorbanceInt(:,rgIdx) = mean(absorbance(:,min(xAxisRangeIntIdx):max(xAxisRangeIntIdx)),2);%Averages the x range
            if sizeIdx < 2
                colourVal = [0,0,1];
            else
                colourVal = [((rgIdx-1)/(sizeIdx-1)),0,((sizeIdx-rgIdx)/(sizeIdx-1))];
            end %if   
            line(xWave,absorbanceInt(:,rgIdx)+(yAxisSpcOffset*(rgIdx-1)),'Parent',axesSpc,'LineWidth',2,'Color',colourVal);
        end %for
    catch error
        disp(['No spectra to plot ! ',error.message])
    end   
 end % if plotSpcBod == 1
 
%% Creates the spectra difference plot
 
if plotDifBod == 1 && plotSpcBod == 1; % 1 plots the spectrum plot, 0 does not
    try
        fDif = figure('Position',[250 30 800 600],'Color',[1 1 1],'Name',id,'NumberTitle','off');
        figureHandles.fDif = fDif;
        axesSpc = axes('Parent',fDif);
        hold(axesSpc,'on');
        title(axesSpc,id);
        xlabel('Wavenumber (1/cm)');% Create xlabel
        ylabel('Diff Abs. (O.D.)');% Create ylabel
        set(axesSpc,'FontSize',12,'LineWidth',2,'TickDir','out','XDir','reverse');
        ylim(axesSpc,[min(yAxisRangeSpcLim) max(yAxisRangeSpcLim)]); 
        xlim(axesSpc,[min(xAxisRangeSpcLim) max(xAxisRangeSpcLim)]);
        absorbanceDiff = diff(absorbanceInt,1,2);
        sizeIdx = size(absorbanceDiff,2);
            for rgIdx  = 1:sizeIdx %Loops throught the range
                if sizeIdx < 2
                    colourVal = [0,0,1];
                else
                    colourVal = [((rgIdx-1)/(sizeIdx-1)),0,((sizeIdx-rgIdx)/(sizeIdx-1))];
                end %if    
            line(xWave,absorbanceDiff(:,rgIdx)+(yAxisDiffOffset*(rgIdx-1)),'Parent',axesSpc,'LineWidth',2,'Color',colourVal);
            end %for
        tbl.AbsorbanceDiff{idIdx}= absorbanceDiff;% Output
    catch error
            disp(['No spectra to plot ! ',error.message])
    end  %try 
    
end %if
 
%% Creates the profile figure
if plotProfileBod == 1
    fProfile = figure('Position',[figureHorPosition 50 1200 900],'Color',[1 1 1],'Name',id,'NumberTitle','off');
    figureHandles.fProfile = fProfile;
    % Create axes
    axes1 = axes('Parent',fProfile);
    % Set the remaining axes properties
    set(axes1,'OuterPosition',[0 0.75 1 0.25],'LineWidth',2);
    try
        line(tFPLC,uv280,'Parent',axes1,'LineWidth',2,'Color',[0 0 0.85]);
    catch error
        disp(['No UV 280 data ! ',error.message]);
    end %try
%     sizeIdx = size(xAxisRangeInt,1);
%     if sizeIdx <= 1; sizeIdx = 2;end;
%     for intIdx = 1:sizeIdx
%         colourVal = [((intIdx-1)/(sizeIdx-1)),0,((sizeIdx-intIdx)/(sizeIdx-1))];
%         try
%         xAxisRangeUV = repmat(mean(xAxisRangeInt(intIdx,:)),1,2);
%         line(xAxisRangeUV ,yAxisRangeUV,'Parent',axes1,'LineWidth',3,'Color',colourVal);
%         end
%     end %for
    xlim(axes1,[min(xAxisRange) max(xAxisRange)]);
    ylim(axes1,[min(yAxisRangeUV) max(yAxisRangeUV)]);
    ylabel({'UV280','(mA)'});
    box(axes1,'off');
    set(axes1,'TickDir','out');
    % Create axes
    axes2 = axes('Parent',fProfile);
    % Set the remaining axes properties
    set(axes2,'OuterPosition',[0 0.50 1 0.25],'LineWidth',2);
    try
        line(tFPLC,cond,'Parent',axes2,'LineWidth',2,'Color',[(255/255),(165/255),0]);
    catch error
        disp(['No Cond data ! ',error.message]);
    end %try
%     for intIdx = 1:sizeIdx
%         colourVal = [((intIdx-1)/(sizeIdx-1)),0,((sizeIdx-intIdx)/(sizeIdx-1))];
%         try
%         line(repmat(mean(xAxisRangeInt(intIdx,:)),1,2),yAxisRangeCond,'Parent',axes2,'LineWidth',3,'Color',colourVal);
%         end
%     end %for
    xlim(axes2,[min(xAxisRange) max(xAxisRange)]);
    ylim(axes2,[min(yAxisRangeCond) max(yAxisRangeCond)]);
    ylabel({'Conductivity','(mS/cm)'});
    box(axes2,'off');
    set(axes2,'TickDir','out');
    % Create axes
    axes3 = axes('Parent',fProfile);
    % Set the remaining axes properties
    set(axes3,'OuterPosition',[0 0.00 1 0.50],'LineWidth',2);
    try
        rgIdx = 1;
        line(tFTIR,pFTIR(rgIdx,:),'Parent',axes3,'LineWidth',2,'Color',[0.85 0 0]);
    catch error
        disp(['No FTIR data ! ',error.message]);
    end %try
    xlim(axes3,[min(xAxisRange) max(xAxisRange)]);
    ylim(axes3,[min(yAxisRangeFTIR) max(yAxisRangeFTIR)]);
    ylabel({'In-column concentration','mg/mL'});
    xlabel(baseType);
    if plotSpcBod == 1
        for intIdx = 1:sizeIdx
            %xAxisRangeInt(intIdx,:)
            colourVal = [((intIdx-1)/(sizeIdx-1)),0,((sizeIdx-intIdx)/(sizeIdx-1))];
            try
                line(repmat(mean(xAxisRangeInt(intIdx,:)),1,2),yAxisRangeFTIR,'Parent',axes3,'LineWidth',3,'Color',colourVal);
            catch
            end
        end %for
    end %if
    box(axes3,'off')
    set(axes3,'TickDir','out');

    title(axes1,id);
end %if plotProfileBod == 1

%% Plot column chart
if plotCyles == 1
    try
    cycleDat = cell2mat(tbl.FTIR1int{idIdx}(:,3));
    cycleDatErr = cell2mat(tbl.FTIR1int{idIdx}(:,4));
    cycleX = 0:length(cycleDat)-1;
    cycleXLabels = tbl.FTIR1int{idIdx}(:,1);
    cycleXLabels = cycleXLabels';%Transpose
    fCycle = figure('Position',[figureHorPosition+1100 50 900 300],'Color',[1 1 1],'Name',id,'NumberTitle','off');
    figureHandles.fCycle = fCycle;
    axesCycle = axes('Parent',fCycle);    
    hBar = bar(cycleX,cycleDat,0.6,'FaceColor',[0.24 0.24 0.24]);
    hold on
    hBar = errorbar(cycleX,cycleDat,cycleDatErr,'.','LineWidth',2,'Color',[0 0 0]);
    xlabel('Cycle number');
    ylabel('In-Column protein concentration');
    set(axesCycle,'TickDir','out','XTick',cycleX);
    set(axesCycle,'XTick',cycleX,'XTickLabel',cycleXLabels);
    title(axesCycle,id);
    catch error
        disp(error.message)
    end    
end %if

%% Closing
vName=@(x) inputname(1);
if ~exist(vName(figureHandles),'var')
    disp('No figure plotted');
    return
end


disp([tblName,' Run ID ',num2str(idIdx),' plotted after ', num2str(toc * 1000,'%.1f'),' ms']);
disp(id);
end %function

function [out1,out2] = v1lowerthan2(var1,var2)
   if var1 > var2
       tmp = var2;var2 = var1;var1 = tmp;
   end
   out1 = var1;out2 = var2;   
end %fct


