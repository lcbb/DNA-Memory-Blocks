%% Load my settings
set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 0.5, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultLineLineWidth', 0.5, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 12, ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.025 0.020]);

% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

%% Load data
clear all

[dogORbldg,dogOrbldgHdr] = fca_readfcs('export_Specimen_001_DOG_OR_BLDG_005_Comp-PE-Texas Red-A+.fcs');
[NotdogORbldg,NotdogOrbldgHdr] = fca_readfcs('export_Specimen_001_DOG_OR_BLDG_005_Comp-PE-Texas Red-A-.fcs');

%% Plot data
close all

figure(1)
set(gcf,'Position',[1333 471 200 200])

x = NotdogORbldg(:,19);
y = dogORbldg(:,19);

% We ignore negative values
for i = 1:length(x)
    if x(i) < 0
       x(i) = NaN;
    end
end

for i = 1:length(y)
    if y(i) < 0
       y(i) = NaN;
    end
end

[~,edges_x] = histcounts(log10(x));
[~,edges_y] = histcounts(log10(y));

hold on
h1 = histogram(x,10.^edges_x);
h1.FaceColor = [255 0 0]./255;
h1.LineStyle = 'None';
h2 = histogram(y,10.^edges_y);
h2.FaceColor = [0 0 255]./255;
h2.LineStyle = 'None';
ax = gca;
ax.XScale = 'log'; 
ax.XLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5]; 
xlabel('"dog","building"-TAMRA [a.u.]')
ylabel('Particles [counts]');
box off

tightfig;