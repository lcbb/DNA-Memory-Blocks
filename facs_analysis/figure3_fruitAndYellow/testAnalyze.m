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

[fruit,fruitHdr] = fca_readfcs('export_Specimen_001_FRUIT_AND_YELLOW_006_Q1 Comp-PE-Texas Red-A- , Comp-APC-A+.fcs');
[yellow,yellowHdr] = fca_readfcs('export_Specimen_001_FRUIT_AND_YELLOW_006_Q3 Comp-PE-Texas Red-A+ , Comp-APC-A-.fcs');
[others,othersHdr] = fca_readfcs('export_Specimen_001_FRUIT_AND_YELLOW_006_Q4 Comp-PE-Texas Red-A- , Comp-APC-A-.fcs');
[banana,bananaHdr] = fca_readfcs('export_Specimen_001_FRUIT_AND_YELLOW_006_Q2 Comp-PE-Texas Red-A+ , Comp-APC-A+.fcs');

%% Plot data
close all

figure(1)
set(gcf,'Position',[1333 471 200 200])

scatter(banana(:,19),banana(:,7),'.','MarkerEdgeColor',[0.78 0.082 0.522])
hold on
scatter(yellow(:,19),yellow(:,7),'.','MarkerEdgeColor',[0.6 0.6 0.6])
scatter(fruit(:,19),fruit(:,7),'.','MarkerEdgeColor',[0.4 0.4 0.4])
scatter(others(:,19),others(:,7),'.','MarkerEdgeColor',[0.2 0.2 0.2])
ax = gca;
ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [5e2 3e4];
ax.YLim = [1e2 3e4];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.YAxis.TickValues = [1e2 1e3 1e4 1e5];
xlabel('"yellow"-TAMRA [a.u.]'); 
ylabel('"fruit"-AF647 [a.u.]');
box on

tightfig;
