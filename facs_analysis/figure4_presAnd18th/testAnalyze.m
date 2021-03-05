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

[wash,washHdr] = fca_readfcs('export_Specimen_001_LINCOLN_007_Comp-PE-Texas Red-A, Comp-APC-A subset.fcs');
[lin,linHdr] = fca_readfcs('export_Specimen_001_LINCOLN_007_Comp-PE-Texas Red-A, Comp-APC-A subset-1.fcs');
[others,othersHdr] = fca_readfcs('export_Specimen_001_LINCOLN_007_Comp-PE-Texas Red-A, Comp-APC-A subset-2.fcs');

%% Plot data
close all

figure(1)
set(gcf,'Position',[1333 471 200 200])

scatter(wash(:,19),wash(:,7),'.','MarkerEdgeColor',[0.4 0.4 0.4])
hold on
scatter(lin(:,19),lin(:,7),'.','MarkerEdgeColor',[43 182 115]./255)
scatter(others(:,19),others(:,7),'.','MarkerEdgeColor',[0.2 0.2 0.2])
ax = gca;
ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [1e2 1e5];
ax.YLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.YAxis.TickValues = [1e2 1e3 1e4 1e5];
xlabel('"president"-TAMRA [a.u.]'); 
ylabel('"18th century"-AF647 [a.u.]');
box on

tightfig;
