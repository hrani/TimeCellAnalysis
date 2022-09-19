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

workingOnServer = 0; %Current
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

figureDetails = compileFigureDetails(11, 2, 5, 0.2, 'inferno'); %(fontSize, lineWidth, markerSize, transparency, colorMap)
%Extra colormap options: inferno/plasma/viridis/magma
%C = distinguishable_colors(input.nMethods);
C = linspecer(input.nMethods);
%C2 = linspecer(4);

%Normalize
maximaN = zeros(input.nMethods, 1);
minimaN = zeros(input.nMethods, 1);
normPredictor = zeros(size(predictor));
for method = 1:input.nMethods
    maximaN(method) = max(predictor(:, method));
    %minimaN(method) = min(predictor(:, method));

    %     posi = find(predictor(:, method) > 0);
    %     negi = find(predictor(:, method) < 0);
    %     normPredictor(posi, method) = predictor(posi, method)/maximaN(method);
    %     normPredictor(negi, method) = -1 * predictor(negi, method)/minimaN(method);

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

disp('Plotting Trial-Avg Activity ...')
load(sprintf('%s/Work/Analysis/Imaging/realData.mat', HOME_DIR2))
input.nDatasets = length(realData);
allCellsTrialAvg = [];

for algo = 1:input.nAlgos
    for dataseti = 1:input.nDatasets
        currentDatasetTrialAvg = squeeze(mean(realData(dataseti).dfbf, 2));
        allCellsTrialAvg = [allCellsTrialAvg; currentDatasetTrialAvg];
    end
    allTimeCellsTrialAvg = allCellsTrialAvg(find(prediction(:, algo)), :);
    allOtherCellsTrialAvg = allCellsTrialAvg(find(~prediction(:, algo)), :);

    %Sort based on Peak - Time Cells
    [~, I1] = max(allTimeCellsTrialAvg, [], 2); %find peak for each cell, along frames; [M, I]
    [~, I2] = sort(I1);
    sortedAllTimeCellsTrialAvg = allTimeCellsTrialAvg(I2, :);

    %Sort based on Peak - Other Cells
    [~, I3] = max(allOtherCellsTrialAvg, [], 2); %find peak for each cell, along frames; [M, I]
    [~, I4] = sort(I3);
    sortedAllOtherCellsTrialAvg = allOtherCellsTrialAvg(I4, :);

    subplot(16, 5, (1:40))
    % imagesc(allTimeCellsTrialAvg*100)
    % colormap(linspecer)
    % h3 = heatmap(sortedAllTimeCellsTrialAvg*100, ...
    %     'Title', 'tiBoot | Peak Sorted', ...
    %     'Colormap', linspecer, ...
    %     'CellLabelColor','none');
    % h3.GridVisible = 'off';
    % %h3.XLabel(customXLabels2);
    % h3.YLabel('All Time Cells', ...
    %     'FontSize', figureDetails.fontSize, ...
    %     'FontWeight', 'bold')
    % % hAx=h3.NodeChildren(3); % return the heatmap underlying axes handle
    % % hAx.FontWeight='bold';
    % %set(gca, 'FontSize', figureDetails.fontSize)
    %
    % cYLabels = 1:size(allTimeCellsTrialAvg, 1);
    % customYLabels = string(cYLabels);
    % customYLabels(mod(cYLabels, 50) ~= 0) = " ";
    % h3.YDisplayLabels = customYLabels;
    %
    % cXLabels = 1:size(allTimeCellsTrialAvg, 2);
    % customXLabels = string(cXLabels);
    % customXLabels(mod(cXLabels, 50) ~= 0) = " ";
    %
    % customXLabels2 = zeros(1, size(customXLabels,2));
    % customXLabels2(1, :) = " ";
    % %h3.XDisplayLabels = customXLabels;
    % h3.XDisplayLabels = customXLabels2;
    % s = struct(h3); %weird hack on a weird hack
    % s.XAxis.TickLabelRotation = 0;  %
    if length(sortedAllTimeCellsTrialAvg) < (size(prediction, 1))/2
        imagesc(sortedAllTimeCellsTrialAvg*100)
        ylabel('All Time Cells', ...
            'FontSize', figureDetails.fontSize, ...
            'FontWeight', 'bold')
    else
        imagesc(sortedAllTimeCellsTrialAvg*100)
        ylabel('All Time Cells', ...
            'FontSize', figureDetails.fontSize, ...
            'FontWeight', 'bold')
    end
    colormap(linspecer)
    title(sprintf('%s | Peak Sorted', algoLabels{algo}), ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
    xlabel('', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
    z = colorbar;
    ylabel(z,'dF/F (%)', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
    set(gca, 'FontSize', figureDetails.fontSize)

    subplot(16, 5, (41:80))
    % h4 = heatmap(sortedAllOtherCellsTrialAvg(1:size(sortedAllTimeCellsTrialAvg,1), :)*100, ...
    %     'ColorLimits',[min(min(sortedAllTimeCellsTrialAvg))*100 max(max(sortedAllTimeCellsTrialAvg))*100], ...
    %     'Colormap', linspecer, ...
    %     'CellLabelColor','none');
    % h4.GridVisible = 'off';
    % h4.XLabel('Frames', ...
    %     'FontSize', figureDetails.fontSize, ...
    %     'FontWeight', 'bold');
    % h4.YLabel('Other Cells (e.g.)', ...
    %     'FontSize', figureDetails.fontSize, ...
    %     'FontWeight', 'bold')
    % % hAx=h4.NodeChildren(3); % return the heatmap underlying axes handle
    % % hAx.FontWeight='bold';
    % %set(gca, 'FontSize', figureDetails.fontSize)
    % h4.YDisplayLabels = customYLabels;
    % h4.XDisplayLabels = customXLabels;
    % s = struct(h4); %weird hack on a weird hack
    % s.XAxis.TickLabelRotation = 0; %vertical
    if length(sortedAllOtherCellsTrialAvg) >= length(sortedAllTimeCellsTrialAvg)
        imagesc(sortedAllOtherCellsTrialAvg(1:size(sortedAllTimeCellsTrialAvg, 1), :)*100)
        ylabel('Other Cells (e.g.)', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
    else
        imagesc(sortedAllOtherCellsTrialAvg*100)
        ylabel('All Other Cells', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
    end
    colormap(linspecer)
    xlabel('Frames', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
    z = colorbar;
    ylabel(z,'dF/F (%)', ...
        'FontSize', figureDetails.fontSize, ...
        'FontWeight', 'bold')
    minVal = min(min(sortedAllTimeCellsTrialAvg))*100;
    maxVal = max(max(sortedAllTimeCellsTrialAvg))*100;
    caxis([minVal maxVal])
    set(gca, 'FontSize', figureDetails.fontSize)

    print(sprintf('%s/figs/%i-%s-RealDataScoresNPredictions-%i-%i', ...
        HOME_DIR2, ...
        algo, ...
        algoLabels{algo}, ...
        nObservations, ...
        workingOnServer), ...   
        '-dpng')
end

%%
elapsedTime2 = toc;
fprintf('Elapsed Time: %.4f seconds\n', elapsedTime2)
disp('... All Done!')