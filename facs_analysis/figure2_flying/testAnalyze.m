%% Load my settings
set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 14, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 14, ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.025 0.020]);

% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

%% Load data
[plane1,plane1Hdr] = fca_readfcs('export_Specimen_001_1to1flying_006_Comp-FITC-A, Comp-APC-A subset.fcs');
[Nplane1,Nplane1Hdr] = fca_readfcs('export_Specimen_001_1to1flying_006_Comp-FITC-A, Comp-APC-A subset-1.fcs');

[plane100,plane100Hdr] = fca_readfcs('export_Specimen_001_1to100flying_007_Comp-FITC-A, Comp-APC-A subset.fcs');
[Nplane100,Nplane100Hdr] = fca_readfcs('export_Specimen_001_1to100flying_007_Comp-FITC-A, Comp-APC-A subset-1.fcs');

[plane10k,plane10kHdr] = fca_readfcs('export_Specimen_001_1to10000flying_008_Comp-FITC-A, Comp-APC-A subset.fcs');
[Nplane10k,Nplane10kHdr] = fca_readfcs('export_Specimen_001_1to10000flying_008_Comp-FITC-A, Comp-APC-A subset-1.fcs');

[plane1M,plane1MHdr] = fca_readfcs('export_Specimen_001_1to1Mflying_009_Comp-FITC-A, Comp-APC-A subset.fcs');
[Nplane1M,Nplane1MHdr] = fca_readfcs('export_Specimen_001_1to1Mflying_009_Comp-FITC-A, Comp-APC-A subset-1.fcs');
%% Plot data
close all

figure(1)
set(gcf,'Position',[1333 471 800 400])

subplot(1,4,1)
scatter(plane1(:,11),plane1(:,7),'.','MarkerEdgeColor',[0.65 0.17 0.64])
hold on
scatter(Nplane1(:,11),Nplane1(:,7),'.','MarkerEdgeColor',[0.2 0.2 0.2])
ax = gca;
ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [1e3 1e5];
ax.YLim = [1e1 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e3 1e4 1e5];
ax.YAxis.TickValues = [1e1 1e2 1e3 1e4 1e5];
xlabel('FITC [a.u.]'); 
ylabel('{\it flying}-AF647 [a.u.]');
box off

subplot(1,4,2)
scatter(plane100(:,11),plane100(:,7),'.','MarkerEdgeColor',[0.65 0.17 0.64])
hold on
scatter(Nplane100(:,11),Nplane100(:,7),'.','MarkerEdgeColor',[0.2 0.2 0.2])
ax = gca;
ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [1e3 1e5];
ax.YLim = [1e1 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e3 1e4 1e5];
ax.YAxis.TickValues = [1e1 1e2 1e3 1e4 1e5];
xlabel('FITC [a.u.]'); 
ylabel('{\it flying}-AF647 [a.u.]');
box off

subplot(1,4,3)
scatter(plane10k(:,11),plane10k(:,7),'.','MarkerEdgeColor',[0.65 0.17 0.64])
hold on
scatter(Nplane10k(:,11),Nplane10k(:,7),'.','MarkerEdgeColor',[0.2 0.2 0.2])
ax = gca;
ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [1e3 1e5];
ax.YLim = [1e1 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e3 1e4 1e5];
ax.YAxis.TickValues = [1e1 1e2 1e3 1e4 1e5];
box off

subplot(1,4,4)
scatter(plane1M(:,11),plane1M(:,7),'.','MarkerEdgeColor',[0.65 0.17 0.64])
hold on
scatter(Nplane1M(:,11),Nplane1M(:,7),'.','MarkerEdgeColor',[0.2 0.2 0.2])
ax = gca;
ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [1e3 1e5];
ax.YLim = [1e1 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e3 1e4 1e5];
ax.YAxis.TickValues = [1e1 1e2 1e3 1e4 1e5];
xlabel('FITC [a.u.]'); 
ylabel('{\it flying}-AF647 [a.u.]');
box off

tightfig;

