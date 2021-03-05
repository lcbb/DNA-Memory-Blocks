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
clear all

[dog,dogHdr] = fca_readfcs('export_Specimen_001_Wolf (3)_007_Q3 Comp-PE-Texas Red-A+ , Comp-APC-A-.fcs');
[wild,wildHdr] = fca_readfcs('export_Specimen_001_Wolf (3)_007_Q1 Comp-PE-Texas Red-A- , Comp-APC-A+.fcs');
[wolf,wolfHdr] = fca_readfcs('export_Specimen_001_Wolf (3)_007_Q2 Comp-PE-Texas Red-A+ , Comp-APC-A+.fcs');
[others,othersHdr] = fca_readfcs('export_Specimen_001_Wolf (3)_007_Q4 Comp-PE-Texas Red-A- , Comp-APC-A-.fcs');

%% Plot data
close all

figure(1)
set(gcf,'Position',[1333 471 200 200])

scatter(dog(:,19),dog(:,7),'.','MarkerEdgeColor',[0.4 0.4 0.4])
hold on
scatter(wild(:,19),wild(:,7),'.','MarkerEdgeColor',[0.6 0.6 0.6])
scatter(wolf(:,19),wolf(:,7),'.','MarkerEdgeColor',[0.78 0.082 0.522])
scatter(others(:,19),others(:,7),'.','MarkerEdgeColor',[0.2 0.2 0.2])
ax = gca;
ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [1e2 1e5];
ax.YLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.YAxis.TickValues = [1e2 1e3 1e4 1e5];
xlabel('"dog"-TAMRA [a.u.]'); 
ylabel('"wild"-AF647 [a.u.]');
box on

tightfig;

%%
figure(2)
set(gcf,'Position',[1333 471 200 200])

x = wolf(:,9);

% We ignore negative values
for i = 1:length(x)
    if x(i) < 0
       x(i) = NaN;
    end
end

[~,edges_x] = histcounts(log10(x));

h1 = histogram(x,10.^edges_x);
h1.FaceColor = [0.78 0.082 0.522];
h1.LineStyle = 'None';
ax = gca;
ax.XScale = 'log'; 
ax.XLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5]; 
xlabel('"black & white"-TYE705 [a.u.]')
ylabel('Particles [counts]');
box off

%%
figure(3)
set(gcf,'Position',[1333 471 472 372])

scatter3(dog(:,19),dog(:,7),dog(:,9),'.','MarkerEdgeColor',[0.4 0.4 0.4])
hold on
scatter3(wild(:,19),wild(:,7),wild(:,9),'.','MarkerEdgeColor',[0.6 0.6 0.6])
scatter3(wolf(:,19),wolf(:,7),wolf(:,9),'.','MarkerEdgeColor',[0.78 0.082 0.522])
scatter3(others(:,19),others(:,7),others(:,9),'.','MarkerEdgeColor',[0.2 0.2 0.2])

plot3(dog(:,19),dog(:,7),ones(size(dog(:,9))).*100,'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(wild(:,19),wild(:,7),ones(size(wild(:,9))).*100,'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(wolf(:,19),wolf(:,7),ones(size(wolf(:,9))).*100,'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(others(:,19),others(:,7),ones(size(others(:,9))).*100,'.','MarkerEdgeColor',[0.9 0.9 0.9])

plot3(dog(:,19),ones(size(dog(:,7))).*100,dog(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(wild(:,19),ones(size(wild(:,7))).*100,wild(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(wolf(:,19),ones(size(wolf(:,7))).*100,wolf(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(others(:,19),ones(size(others(:,7))).*100,others(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])

plot3(ones(size(dog(:,19))).*100,dog(:,7),dog(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(ones(size(wild(:,19))).*100,wild(:,7).*100,wild(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(ones(size(wolf(:,19))).*100,wolf(:,7),wolf(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])
plot3(ones(size(others(:,19))).*100,others(:,7),others(:,9),'.','MarkerEdgeColor',[0.9 0.9 0.9])

ax = gca;
ax.XScale = 'log'; ax.YScale = 'log'; ax.ZScale = 'log';
ax.XLim = [1e2 1e5]; ax.YLim = [1e2 1e5]; ax.ZLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.YAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.ZAxis.TickValues = [1e2 1e3 1e4 1e5];
xlabel('"dog"-TAMRA [a.u.]'); 
ylabel('"wild"-AF647 [a.u.]');
zlabel('"black & white"-TYE705 [a.u.]')
box on
view(60,35)
set(gca,'Ydir','reverse', ...
        'GridLineStyle', 'none', ...
        'MinorGridLineStyle', 'none');


tightfig;

%% 
figure(4)
set(gcf,'Position',[1333 471 1200 360])

subplot(1,3,1)
scatter3(dog(:,19),dog(:,7),dog(:,9),'.','MarkerEdgeColor',[0.4 0.4 0.4])
hold on
scatter3(wild(:,19),wild(:,7),wild(:,9),'.','MarkerEdgeColor',[0.6 0.6 0.6])
scatter3(wolf(:,19),wolf(:,7),wolf(:,9),'.','MarkerEdgeColor',[0.78 0.082 0.522])
scatter3(others(:,19),others(:,7),others(:,9),'.','MarkerEdgeColor',[0.2 0.2 0.2])

ax = gca;
ax.XScale = 'log'; ax.YScale = 'log'; ax.ZScale = 'log';
ax.XLim = [1e2 1e5]; ax.YLim = [1e2 1e5]; ax.ZLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.YAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.ZAxis.TickValues = [1e2 1e3 1e4 1e5];
xlabel('"dog"-TAMRA [a.u.]'); 
ylabel('"wild"-AF647 [a.u.]');
zlabel('"black & white"-TYE705 [a.u.]')
box on
set(gca, ...
        'GridLineStyle', 'none', ...
        'MinorGridLineStyle', 'none');
view(90,0)

subplot(1,3,2)
scatter3(dog(:,19),dog(:,7),dog(:,9),'.','MarkerEdgeColor',[0.4 0.4 0.4])
hold on
scatter3(wild(:,19),wild(:,7),wild(:,9),'.','MarkerEdgeColor',[0.6 0.6 0.6])
scatter3(wolf(:,19),wolf(:,7),wolf(:,9),'.','MarkerEdgeColor',[0.78 0.082 0.522])
scatter3(others(:,19),others(:,7),others(:,9),'.','MarkerEdgeColor',[0.2 0.2 0.2])

ax = gca;
ax.XScale = 'log'; ax.YScale = 'log'; ax.ZScale = 'log';
ax.XLim = [1e2 1e5]; ax.YLim = [1e2 1e5]; ax.ZLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.YAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.ZAxis.TickValues = [1e2 1e3 1e4 1e5];
xlabel('"dog"-TAMRA [a.u.]'); 
ylabel('"wild"-AF647 [a.u.]');
zlabel('"black & white"-TYE705 [a.u.]')
box on
set(gca,'Xdir','reverse', ...
        'GridLineStyle', 'none', ...
        'MinorGridLineStyle', 'none');
view(180,0)

subplot(1,3,3)
scatter3(dog(:,19),dog(:,7),dog(:,9),'.','MarkerEdgeColor',[0.4 0.4 0.4])
hold on
scatter3(wild(:,19),wild(:,7),wild(:,9),'.','MarkerEdgeColor',[0.6 0.6 0.6])
scatter3(wolf(:,19),wolf(:,7),wolf(:,9),'.','MarkerEdgeColor',[0.78 0.082 0.522])
scatter3(others(:,19),others(:,7),others(:,9),'.','MarkerEdgeColor',[0.2 0.2 0.2])

ax = gca;
ax.XScale = 'log'; ax.YScale = 'log'; ax.ZScale = 'log';
ax.XLim = [1e2 1e5]; ax.YLim = [1e2 1e5]; ax.ZLim = [1e2 1e5];
ax.Color = 'none';
ax.XAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.YAxis.TickValues = [1e2 1e3 1e4 1e5];
ax.ZAxis.TickValues = [1e2 1e3 1e4 1e5];
xlabel('"dog"-TAMRA [a.u.]'); 
ylabel('"wild"-AF647 [a.u.]');
zlabel('"black & white"-TYE705 [a.u.]')
box on
set(gca, ...
        'GridLineStyle', 'none', ...
        'MinorGridLineStyle', 'none');
view(0,90)

tightfig;
