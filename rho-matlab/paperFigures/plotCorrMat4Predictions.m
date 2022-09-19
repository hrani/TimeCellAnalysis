clear
close all
tic
figureDetails = compileFigureDetails(11, 2, 5, 0.2, 'inferno'); %(fontSize, lineWidth, markerSize, transparency, colorMap)
%algoLabels = {'tiMean', 'tiBoot', 'tiBoth', 'r2bMean', 'r2bBoot', 'peq'};
algoLabels = {'tiMean', 'tiBoot', 'tiBoth', 'r2bMean', 'r2bBoot', 'peq', 'tiMean-O', 'tiBase-O', 'r2bBase-O', 'peqBase-O'};
predictionTable = readtable('prediction.csv', 'Range', 'A2:J76546', 'ReadRowNames', false);
prediction = table2array(predictionTable);
correlationMatrix4Predictions = corrcoef(prediction(56296:76545, :), 'Rows', 'pairwise'); %Restricted to the physiology datasets (w/Background Activity)
fig1 = figure(1);
clf
set(fig1, 'Position', [100, 100, 800, 800])
h = heatmap(correlationMatrix4Predictions, ...
    'Colormap', linspecer, ...
    'Title', 'Correlation Coefficients - Predictions', ...
    'CellLabelColor','none');
h.XDisplayLabels = algoLabels;
h.YDisplayLabels = algoLabels;
%h.YDisplayLabels = {'', '', '', '', '', ''};
set(gca, 'FontSize', figureDetails.fontSize)
print('corrMat4Predictions', ...
    '-dpng')
toc