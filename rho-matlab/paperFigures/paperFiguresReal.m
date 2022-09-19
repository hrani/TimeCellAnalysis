% Paper Figures
% AUTHOR: Kambadur Ananthamurthy
% DETAILS: 13 real datasets were analysed on the basis
% of a variety of numerical procedures. Cells in each dataset were given
% analog scores on the same basis.
% Load the Consolidated Analysis and look for patterns in the plots.

close all
%clear

tic

%%
input.nMethods = 4; % 2x ti, r2b, peq
input.nAlgos = 10; %detection algorithms; 3x ti, 2x r2b, 4x Otsu, peq
nConcordanceAlgos = input.nAlgos;
peakSortOtherCells = 1;

workingOnServer = 2; %Current
% Directory config
if workingOnServer == 1
    HOME_DIR = '/home/bhalla/ananthamurthy/';
    saveDirec = strcat(HOME_DIR, 'Work/Analysis/Imaging/');
elseif workingOnServer == 2
    HOME_DIR = '/home/ananth/Documents/';
    HOME_DIR2 = '/media/ananth/Storage/';
    saveDirec = strcat(HOME_DIR2, 'Work/Analysis/RHO/');
else
    HOME_DIR = '/home/ananth/Documents/';
    HOME_DIR2 = '/home/ananth/Desktop/';
    saveDirec = strcat(HOME_DIR2, 'Work/Analysis/RHO/');
end
%Additinal search paths
addpath(genpath(strcat(HOME_DIR, 'rho-matlab/CustomFunctions')))
%addpath(genpath(strcat(HOME_DIR, 'paperFigures'))) %or just 'cd' to this
%addpath(genpath(strcat(HOME_DIR, 'TimeCellAnalysis/TcPy'))) %python/C++ analysis; not required if the .csv files are copied to HOME_DIR/AnalysisOutput/
addpath(genpath(strcat(HOME_DIR, 'AnalysisOutput'))) % location for .csv (python/C++ analysis outputs)
addpath(genpath(strcat(HOME_DIR, 'rho-matlab/localCopies')))
make_db_realBatch

diaryOn = 0;
if diaryOn
    if workingOnServer == 1
        diary (strcat(HOME_DIR, '/logs/benchmarksDiary'))
    else
        diary (strcat(HOME_DIR2, '/logs/benchmarksDiary_', num2str(gDate), '_', num2str(gRun)))
    end
    diary on
end

%Load .csv from Analysis Output dir
r2b = load(sprintf('%s/AnalysisOutput/r2b_realData.csv', HOME_DIR));
ti = load(sprintf('%s/AnalysisOutput/ti_realData.csv', HOME_DIR));
peq = load(sprintf('%s/AnalysisOutput/peq_realData.csv', HOME_DIR));

%nObservations = input.nDatasets*input.nCells;
nObservations = size(ti, 1);
predictor = zeros(nObservations, input.nMethods);
prediction = zeros(nObservations, input.nAlgos);

%Predictors
predictor(:, 1) = ti(:, 3);
predictor(:, 2) = ti(:, 4);
predictor(:, 3) = r2b(:, 4)./r2b(:, 3); %r2b
predictor(:, 4) = peq(:, 4); %peqScore

%Predictions (X)
prediction(:, 1) = ti(:, 6);
prediction(:, 2) = ti(:, 7);
prediction(:, 3) = and(ti(:, 6), ti(:, 7)); % same as and(predictions(:, 1), predictions(:, 2))
prediction(:, 4) = r2b(:, 6);
prediction(:, 5) = r2b(:, 7);
prediction(:, 6) = peq(:, 3);

%Now to add the Otsu based classifications
[T1, EM1] = graythresh(predictor(:, 1));
[T2, EM2] = graythresh(predictor(:, 2));
[T3, EM3] = graythresh(predictor(:, 3));
[T4, EM4] = graythresh(predictor(:, 4));

prediction(:, 7) = predictor(:, 1) > T1; %meanScore_ti_otsu
prediction(:, 8) = predictor(:, 2) > T2; %baseScore_ti_otsu
prediction(:, 9) = predictor(:, 3) > T3; %baseScore_r2b_otsu
prediction(:, 10) = predictor(:, 4) > T4; %peqScore_otsu

X = prediction;

% Search for uniquely identified Time Cells with some threshold
unique1 = zeros(size(X, 1), input.nAlgos);
unique2 = zeros(size(X, 1), input.nAlgos);
unique3 = zeros(size(X, 1), input.nAlgos);
for obs = 1:size(X, 1)
    if sum(squeeze(X(obs, :))) == 1
        unique1(obs, :) = X(obs, :);
    elseif sum(squeeze(X(obs, :))) == 2
        unique2(obs, :) = X(obs, :);
    elseif sum(squeeze(X(obs, :))) == 3
        unique3(obs, :) = X(obs, :);
    end
end

figureDetails = compileFigureDetails(11, 2, 5, 0.2, 'magma'); %(fontSize, lineWidth, markerSize, transparency, colorMap)
%Extra colormap options: inferno/plasma/viridis/magma
%C = distinguishable_colors(input.nMethods);
C = magma(input.nMethods*2);
%Normalize
maximaN = zeros(input.nMethods, 1);
minimaN = zeros(input.nMethods, 1);
normPredictor = zeros(size(predictor));
for method = 1:input.nMethods
    maximaN(method) = max(predictor(:, method));
    %minimaN(method) = min(predictor(:, method));

    normPredictor(:, method) = predictor(:, method)/maximaN(method);
end

%Z-Score
zPredictor = zscore(predictor);

% Normalize Z-Score
maximaZ = zeros(input.nMethods, 1);
minimaZ = zeros(input.nMethods, 1);
nzPredictor = zeros(size(predictor));
for method = 1:input.nMethods
    maximaZ(method) = max(zPredictor(:, method));
    minimaZ(method) = min(zPredictor(:, method));

    posi = find(zPredictor(:, method) > 0);
    negi = find(zPredictor(:, method) < 0);
    nzPredictor(posi, method) = zPredictor(posi, method)/maximaN(method);
    nzPredictor(negi, method) = -1 * zPredictor(negi, method)/minimaN(method);
end

correlationMatrix = corrcoef(normPredictor, 'Rows', 'pairwise');
correlationMatrix4Predictions = corrcoef(prediction, 'Rows', 'pairwise');

%Labels
algoLabels = {'tiMean', 'tiBoot', 'tiBoth', 'r2bMean', 'r2bBoot', 'peq', 'tiMean-O', 'tiBase-O', 'r2bBase-O', 'peqBase-O'};
concordanceAlgoLabels = {'>=1', '>=2', '>=3', '>=4', '>=5', '>=6', '>=7', '>=8', '>=9', '=10'};
equivalenceLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
methodLabels = {'tiMean', 'tiBase', 'r2bBase', 'peqBase'};
procedureLabels = {'Synth.', 'TI', 'R2B', 'PEQ'};

paramLabels = {'Noise (%)', 'EW (%ile)', 'Imp. (frames)', 'HTR (%)', 'Back. (\lambda)'};
metricLabels = {'TNR', 'FNR', 'FPR', 'TPR'};

metricLabels_pos = {'TP', 'FP'};
metricLabels_neg = {'TN', 'FN'};

metricLabels_true = {'TP', 'TN'};
metricLabels_false = {'FP', 'FN'};
metricLabels_timeCells = {'True Pos', 'False Neg'};
metricLabels_otherCells = {'True Neg', 'False Pos'};
metricLabels2 = {'Recall', 'Precision', 'F1 Score'};

%% Scores and Comparison - All

disp('Plotting Scores and Comparisions ...')
fig7 = figure(7);
clf
set(fig7, 'Position', [1000, 300, 900, 1200])

for method = 1:input.nMethods

    if method == 1
        subplotNum = [1, 6];
    elseif method == 2
        subplotNum = [2, 7];
    elseif method == 3
        subplotNum = [21, 26];
    elseif method == 4
        subplotNum = [22, 27];
    end

    subplot(36, 5, subplotNum)
    h1 = histogram(predictor(:, method), ...
        'EdgeColor', C(5, :), ...
        'FaceColor', C(5, :));
    if method == 1
        h1.NumBins = 10000;
    elseif method == 2
        h1.NumBins = 6000;
    elseif method == 3
        h1.NumBins = 400;
    elseif method == 4
        h1.NumBins = 10000;
        %     elseif method == 5
        %         h.NumBins = 10000;
    end
    binWidth = h1.BinWidth;
    %     ylabel('Counts', ...
    %         'FontSize', figureDetails.fontSize, ...
    %         'FontWeight', 'bold')
    hold off
    title(sprintf('%s', char(methodLabels(method))), ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')

    if method == 3
        xlabel('Score', ...
            'FontSize', figureDetails.fontSize, ...
            'FontWeight', 'bold')
        ylabel('Counts', ...
            'FontSize', figureDetails.fontSize, ...
            'FontWeight', 'bold');
    else
        ylabel('')
        xlabel('')
    end

    set(gca, 'FontSize', figureDetails.fontSize)
end
hold off

subplot(36, 5, [4:5, 9:10, 14:15, 19:20, 29:30])
h1 = heatmap(correlationMatrix, ...
    'Colormap', linspecer, ...
    'Title', 'Correlation Coefficients - Scores', ...
    'CellLabelColor','none');
h1.XDisplayLabels = methodLabels;
h1.YDisplayLabels = methodLabels;
%h.YDisplayLabels = {'', '', '', '', '', ''};

set(gca, 'FontSize', figureDetails.fontSize)

subplot(36, 5, [(51:52), (56:57), (61:62), (66:67), (71:72), (76:77)])
h2 = heatmap(correlationMatrix4Predictions, ...
    'Colormap', linspecer, ...
    'Title', 'Correlation Coefficients - Predictions', ...
    'CellLabelColor', 'none');
h2.XDisplayLabels = algoLabels;
h2.YDisplayLabels = algoLabels;
set(gca, 'FontSize', figureDetails.fontSize)

timeCellPredictions = (sum(prediction, 1)/size(prediction, 1))*100;
subplot(36, 5, [(54:55), (59:60), (64:65), (69:70), (74:75), (79:80)])
b = barh(timeCellPredictions);
set(gca, 'YDir','reverse')
%b.FaceColor = C(1, :);
b.FaceColor = C(5, :);
axis tight
title(sprintf('Time Cell Predictions (N=%i)', size(prediction, 1)), ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xlabel('Counts (%)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
%yticks([])
yticklabels(algoLabels)
set(gca, 'FontSize', figureDetails.fontSize)

load(sprintf('%s/Work/Analysis/Imaging/realData.mat', HOME_DIR2))
input.nDatasets = length(realData);
selectedAlgo = 3; %Using TI-Bootstrap w/ Activity Filter
normalizeEachCell = 1;
allCellsTrialAvg = [];
for dataseti = 1:input.nDatasets
    currentDatasetTrialAvg = squeeze(mean(realData(dataseti).dfbf, 2));
    allCellsTrialAvg = [allCellsTrialAvg; currentDatasetTrialAvg];
end
if normalizeEachCell
    for celli = 1:size(allCellsTrialAvg, 1)
        maxVal = max(squeeze(allCellsTrialAvg(celli, :)));
        allCellsTrialAvg(celli, :) = allCellsTrialAvg(celli, :)/maxVal;
    end
end
allTimeCellsTrialAvg = allCellsTrialAvg(find(prediction(:, selectedAlgo)), :);
allOtherCellsTrialAvg = allCellsTrialAvg(find(~prediction(:, selectedAlgo)), :);

%Sort based on Peak - Time Cells
[~, I1] = max(allTimeCellsTrialAvg, [], 2); %find peak for each cell, along frames; [M, I]
[~, I2] = sort(I1);
sortedAllTimeCellsTrialAvg = allTimeCellsTrialAvg(I2, :);

%Sort based on Peak - Other Cells
[~, I3] = max(allOtherCellsTrialAvg, [], 2); %find peak for each cell, along frames; [M, I]
[~, I4] = sort(I3);
sortedAllOtherCellsTrialAvg = allOtherCellsTrialAvg(I4, :);

subplot(36, 5, (106:140))
imagesc(sortedAllTimeCellsTrialAvg*100)
colormap(figureDetails.colorMap)
title(sprintf('Real Physiology | %s | Peak Sorted', algoLabels{selectedAlgo}), ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xlabel('', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
ylabel('All Time Cells', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
z = colorbar;
if normalizeEachCell
    ylabel(z,'Norm. dF/F (%)', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
else
    ylabel(z,'dF/F (%)', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
end
set(gca, 'FontSize', figureDetails.fontSize)

randomI = sort(randi([1 1478], 1, size(sortedAllTimeCellsTrialAvg, 1)));
subplot(36, 5, (146:180))
imagesc(sortedAllOtherCellsTrialAvg(randomI, :)*100)
colormap(figureDetails.colorMap)
xlabel('Frames', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
ylabel('Other Cells (e.g.)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
z = colorbar;
if normalizeEachCell
    ylabel(z,'Norm. dF/F (%)', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
else
    ylabel(z,'dF/F (%)', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
end
minVal = min(min(sortedAllTimeCellsTrialAvg))*100;
maxVal = max(max(sortedAllTimeCellsTrialAvg))*100;
caxis([minVal maxVal])
set(gca, 'FontSize', figureDetails.fontSize)

print(sprintf('%s/figs/7-RealDataScoresNPredictions-%i-%i', ...
    HOME_DIR2, ...
    nObservations, ...
    workingOnServer), ...
    '-dpng')

disp('... Scores and Comparisons Plotted')
%%
elapsedTime2 = toc;
fprintf('Elapsed Time: %.4f seconds\n', elapsedTime2)
disp('... All Done!')