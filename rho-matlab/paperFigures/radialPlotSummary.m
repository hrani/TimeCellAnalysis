%% Radial Plot Summary - run after running paperFiguresSynth.m
% Requires edits before usable as standalone

disp('Plotting Radial Plot Summary ...')

fig7 = figure(7);
clf
set(fig7, 'Position', [1500, 300, 900, 1200])
subplot(4, 3, 1)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(1, :)/min(dependenceAll(1, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('tiMean', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 2)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(2, :)/min(dependenceAll(2, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('tiBoot', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 3)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(3, :)/min(dependenceAll(3, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('tiBoth', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 4)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(4, :)/min(dependenceAll(4, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('r2bMean', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 5)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(5, :)/min(dependenceAll(5, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('r2bBoot', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 6)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(6, :)/min(dependenceAll(6, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('peq', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 7)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(7, :)/min(dependenceAll(7, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('tiMean-O', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 8)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(5, :)/min(dependenceAll(1, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('tiBase-O', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 9)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(5, :)/min(dependenceAll(1, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('r2bMean-O', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 10)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(5, :)/min(dependenceAll(1, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('r2bBase-O', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

subplot(4, 3, 11)
theta = linspace(0, 2*pi, input.nParams);
rho = transpose(dependenceAll(5, :)/min(dependenceAll(1, :))); %Remember to transpose
pol = secplot(theta, rho);
ax = gca;
ax.RTickLabel = num2cell(sort(dependence2));
ax.ThetaGrid = 'off';
ax.ThetaTick = rad2deg(theta);
ax.ThetaTickLabel = paramLabels;
title('peqBase-O', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
set(gca, 'FontSize', figureDetails.fontSize)

txt = {'Parameter Ranges -', ...
    'Noise: 10% to 70%', ...
    'EW: 30th to 90th %ile', ...
    'Imp.: 0 to 66 frames', ...
    'HTR: 33% to 99%', ...
    'Back. \lambda = 0.2255 to 1.1740'};
text(0, 3, txt, 'FontSize', figureDetails.fontSize, 'FontWeight', 'bold')

print(sprintf('%s/figs/RadialSummary-%i-%i-%i-%i-%i_%i', ...
    HOME_DIR2, ...
    input.gDate, ...
    input.gRun, ...
    input.nDatasets, ...
    workingOnServer), ...
    '-dpng')

disp('... Radial Plot Summary Plotted.')