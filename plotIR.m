function figId = plotIR(wavenumber,absorbance,xRange,yRange)
%% Test input
if nargin < 2
    disp('Not enough input arguments');
    disp('plotIR(wavenumber,absorbance,[xRange],[yRange])');
    return
end


%% Ploting parameters
if exist('xRange','var') ==0
    xRange = [min(min(wavenumber)),max(max(wavenumber))];
else
   xRange = [min(xRange), max(xRange)];
end

if exist('yRange','var') ==0
    ymin = min(min(absorbance));
    ymax = max(max(absorbance));
    yRange = [(ymin-0.1*(ymax-ymin)),(ymax+0.1*(ymax-ymin))];
else
    yRange = [min(yRange), max(yRange)];
end

%% Create figure
tVal  = clock;
tStr  = [num2str(tVal(1),'%04d'),num2str(tVal(2),'%02d'),...
    num2str(tVal(3),'%02d'),' ',num2str(tVal(4),'%02d'),num2str(tVal(5),'%02d'),num2str(round(tVal(6)),'%02d')];
figId = figure('Color',[1 1 1],'NumberTitle','off','Name',tStr);

% Create axes
axes1 = axes('Parent',figId);
hold(axes1,'on');
hold on;
colourVec = [0.8*rand,0.8*rand,0.8*rand];
plot(wavenumber,absorbance,'LineWidth',2,'Color',colourVec);% Create multiple lines using matrix input to plot
xlabel('Wavenumber (1/cm)');% Create xlabel
ylabel('Absorbance');% Create ylabel
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Arial','FontSize',14.0,'FontSmoothing','off',...
    'FontUnits','centimeters','TickDir','out','XColor',[0 0 0],'XDir','reverse',...
    'XMinorTick','on','YColor',[0 0 0],'LineWidth',2,'YMinorTick','on','ZColor',[0 0 0]);

ylim(axes1,yRange); 
xlim(axes1,xRange);

function figId = plotUV(valX,valY,xRange,yRange)
%% Test input
if nargin < 2
    disp('Not enough input arguments...');
    disp('figId = plotUv(time,absorbance,[xRange],[yRange])');
    return
end

%% Ploting parameters
if ~exist('xRange','var')
    xRange = [min(min(valX)),max(max(valX))];
else
   xRange = [min(xRange), max(xRange)];
end

if ~exist('yRange','var')
    ymin = min(min(valY));
    ymax = max(max(valY));
    yRange = [(ymin-0.1*(ymax-ymin)),(ymax+0.1*(ymax-ymin))];
else
    yRange = [min(yRange), max(yRange)];
end

%% Create figure
tVal  = clock;
tStr  = [num2str(tVal(1),'%04d'),num2str(tVal(2),'%02d'),...
    num2str(tVal(3),'%02d'),' ',num2str(tVal(4),'%02d'),num2str(tVal(5),'%02d'),num2str(round(tVal(6)),'%02d')];
figId = figure('Color',[1 1 1],'NumberTitle','off','Name',tStr);

% Create axes
axes1 = axes('Parent',figId);
hold(axes1,'on');
hold on;
colourVec = [0.8*rand,0.8*rand,0.8*rand];
plot(valX,valY,'LineWidth',2,'Color',colourVec);% Create multiple lines using matrix input to plot
xlabel('Time (s)');% Create xlabel
ylabel('UV 280 nm (mA)');% Create ylabel
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Arial','FontSize',14.0,'FontSmoothing','off',...
    'FontUnits','centimeters','TickDir','out','XColor',[0 0 0],...
    'XMinorTick','on','YColor',[0 0 0],'LineWidth',2,'YMinorTick','on','ZColor',[0 0 0]);

ylim(axes1,yRange); 
xlim(axes1,xRange);

