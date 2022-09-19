% Paper Figures - For Synthetic Data (Fig1, Fig4, Fig5, Fig6, supFig6, Fig8, supFig8)
% AUTHOR: Kambadur Ananthamurthy
% DETAILS: 537 uniquely tagged synthetic datasets were analysed on the basis
% of a variety of numerical procedures. Cells in each dataset were given
% analog scores on the same basis.
% Load the Consolidated Analysis and look for patterns in the plots.
% RGB Triplet list for MATLAB:https://in.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html

close all
%clear

tic
%%
input.nCells = 135;
input.nMethods = 4; % 2x ti, r2b, peq
input.nAlgos = 10; %detection algorithms; 3x ti, 2x r2b, 4x Otsu
nConcordanceAlgos = input.nAlgos;
input.nParams = 5; %Noise, EW, Imp., HTR, and Background

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
make_db

saveFolder = strcat(saveDirec, db.mouseName, '/', db.date, '/');

datasetCatalog = 0; %Only to select the batch for datasets
if datasetCatalog == 0
    %Synthetic Dataset Details
    input.gDate = 20220618; %generation date; int
    input.gRun = 1; %generation run number; int
    input.nDatasets = 567; %int
elseif datasetCatalog == 1
    % Synthetic Dataset Details
    input.gDate = 20220308; %generation date
    input.gRun = 1; %generation run number
    input.nDatasets = 537;
end

%Local variables
gDate = input.gDate; %int
gRun = input.gRun; %int
nDatasets = input.nDatasets; %int

diaryOn = 0;
if diaryOn
    if workingOnServer == 1
        diary (strcat(HOME_DIR, '/logs/benchmarksDiary'))
    else
        diary (strcat(HOME_DIR2, '/logs/benchmarksDiary_', num2str(gDate), '_', num2str(gRun)))
    end
    diary on
end

% Operations
plotRefQ = 1;
plotAnalysedQs = 1;
plotDatasetCheck = 0;
bandClasses = 1; %To bundle (TP, FP) and (TN, FN)

%Load .csv from Analysis Dir
r2b = load(sprintf('%sAnalysisOutput/r2b_%i_gRun%i_batch_%i.csv', HOME_DIR, gDate, gRun, nDatasets));
ti = load(sprintf('%sAnalysisOutput/ti_%i_gRun%i_batch_%i.csv', HOME_DIR, gDate, gRun, nDatasets));
peq = load(sprintf('%sAnalysisOutput/peq_%i_gRun%i_batch_%i.csv', HOME_DIR, gDate, gRun, nDatasets));

%nObservations = input.nDatasets*input.nCells;
nObservations = size(ti, 1);
predictor = zeros(nObservations, input.nMethods);
prediction = zeros(nObservations, input.nAlgos);

%Get Response and Predictors
%response = r2b(:, 3); %same as ti(:, 3)
response = ti(:, 3); %same as r2b(:, 3) or peq(:, 3)
predictor(:, 1) = ti(:, 4);
predictor(:, 2) = ti(:, 5);
predictor(:, 3) = r2b(:, 5)./r2b(:, 4); %r2b
predictor(:, 4) = peq(:, 5); %peqScore

%Get Response (Y) and Predictions (X)
%[Y, X] = consolidatePredictions(input, sdo_batch, cData);
Y = response;
prediction(:, 1) = ti(:, 7);
prediction(:, 2) = ti(:, 8);
prediction(:, 3) = and(ti(:, 7), ti(:, 8)); % same as and(predictions(:, 1), predictions(:, 2))
prediction(:, 4) = r2b(:, 7);
prediction(:, 5) = r2b(:, 8);
prediction(:, 6) = peq(:, 4);

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

input.removeNaNs = 0;
input.saveFolder = saveFolder;

% TP, FN, FP, TN
%nCases = length(Y);
%confusionMatrix = zeros(input.nAlgos, 2, 2);

clear results1
clear results2
clear results3
[results1, ~, results3] = compareAgainstTruth(X, Y, input);

% Search for uniquely identified Time Cells with some threshold
unique1 = zeros(length(Y), input.nAlgos);
unique2 = zeros(length(Y), input.nAlgos);
unique3 = zeros(length(Y), input.nAlgos);
for obs = 1:length(Y)
    if sum(squeeze(X(obs, :))) == 1
        unique1(obs, :) = X(obs, :);
    elseif sum(squeeze(X(obs, :))) == 2
        unique2(obs, :) = X(obs, :);
    elseif sum(squeeze(X(obs, :))) == 3
        unique3(obs, :) = X(obs, :);
    end
end

figureDetails = compileFigureDetails(11, 2, 8, 10, 0.2, 'magma'); %(fontSize, lineWidth, markerSize, capSize, transparency, colorMap)
%Extra colormap options: inferno/plasma/viridis/magma
%C = distinguishable_colors(input.nAlgos);
C = magma(input.nMethods*2);

selectedRGBs = [4, 1, 3, 5];
C1 = zeros(length(selectedRGBs), 3);
for rgbi = 1:length(selectedRGBs)
    C1(rgbi, :) = C(selectedRGBs(rgbi), :);
end
C1rev = C1(length(selectedRGBs):-1:1, :);

%disp('Generating subplot 1 ...')
input.removeNaNs = 1;

% %Prepare the Look Up Table (LUT)
% [response, predictor] = consolidatePredictors(sdo_batch, cData, input); %NOTE: "Y" may be identical to "response".

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
%correlationMatrix = corrcoef(predictor);

%Labels
algoLabels = {'tiMean', 'tiBoot', 'tiBoth', 'r2bMean', 'r2bBoot', 'peq', 'tiMean-O', 'tiBase-O', 'r2bBase-O', 'peqBase-O'};
algoLabels2 = {'tme', 'tbo', 'tbh', 'rme', 'rbo', 'peq', 'tmeO', 'tbaO', 'rbaO', 'peqO'};
concordanceAlgoLabels = {'>=1', '>=2', '>=3', '>=4', '>=5', '>=6', '>=7', '>=8', '>=9', '=10'};
equivalenceLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
methodLabels = {'tiMean', 'tiBase', 'r2bBase', 'peqBase'};
procedureLabels = {'Synth.', 'TI', 'R2B', 'PEQ'};

impX = 0; %0 for off, or non-zero multiplier
paramLabels = {'Noise (%)', 'EW (%ile)', 'Imp. (frames)', 'HTR (%)', 'Back. (\lambda)'};
paramLabels_mod = {'Noise (%)', 'EW (%ile)', [num2str(impX),'x Imp.'], 'HTR (%)', 'Back. (\lambda)'};

paramLabels2 = {'Noise', 'EW', 'Imp.', 'HTR', 'Back.'};
paramLabels3 = {'Noise (%)', 'EW (%ile)', 'EW (SDX)', 'Imp. (frames)', 'Imp. Type', 'HTR (%)', 'HT Assign.', 'Trial Order', 'Back. (\lambda)', 'TC (%)', 'Cell Order' };
paramLabels3_mod = {'Noise (%)', 'EW (%ile)', 'EW (SDX)', [num2str(impX), ' x Imp.'], 'Imp. Type', 'HTR (%)', 'HT Assign.', 'Trial Order', 'Back. (\lambda)', 'TC (%)', 'Cell Order' };

metricLabels = {'TNR', 'FNR', 'FPR', 'TPR'};
metricLabels_pos = {'TP', 'FP'};
metricLabels_neg = {'TN', 'FN'};

metricLabels_true = {'TP', 'TN'};
metricLabels_false = {'FP', 'FN'};
metricLabels_timeCells = {'True Pos', 'False Neg'};
metricLabels_otherCells = {'True Neg', 'False Pos'};
metricLabels2 = {'Recall', 'Precision', 'F1 Score'};

selectedAlgo = 3; %Using TI-Bootstrap w/ Activity Filter
%% Example Schematic
addpath(genpath(sprintf('%srho-matlab/src', HOME_DIR)))
%generateSyntheticDataExample(workingOnServer, figureDetails)

% %% Scores and Comparison - All
% plotOptimalPoint = 0;
% disp('Plotting Scores and Comparisons ...')
% fig4 = figure(4);
% clf
% set(fig4, 'Position', [1, 1, 900, 1200])
% for method = 1:input.nMethods
%     ptcScores = normPredictor(response == 1, method);
%     ocScores = normPredictor(response == 0, method);
% 
%     if method == 1
%         subplotNum = [(1:2), (10:11)];
%     elseif method == 2
%         subplotNum = [(4:5), (13:14)];
%     elseif method == 3
%         subplotNum = [(55:56), (64:65)];
%     elseif method == 4
%         subplotNum = [(58:59), (67:68)];
%     end
% 
%     subplot(21, 9, subplotNum)
%     h = histogram(ptcScores, ...
%         'EdgeColor', C(5, :), ...
%         'FaceColor', C(5, :));
%     ylim([0 108])
% 
%     set(gca,'xticklabel',{[]})
%     if method == 1
%         h.NumBins = 10000;
%     elseif method == 2
%         h.NumBins = 6000;
%     elseif method == 3
%         h.NumBins = 400;
%     elseif method == 4
%         h.NumBins = 10000;
%     end
%     binWidth = h.BinWidth;
% 
%     set(gca, 'yticklabel', {'0', '100', '200'})
% 
%     hold off
%     title(sprintf('%s', char(methodLabels(method))), ...
%         'FontSize', figureDetails.fontSize, ...
%         'FontWeight', 'bold')
%     ax = gca;
%     ax.YAxis.TickLabelFormat = '%,g';
% 
%     if method == 1
%         xlim([-0.1, 1])
%     elseif method == 2
%         xlim([-0.2, 0.6])
%     elseif method == 3
%         xlim([-0.001, 0.04])
%     elseif method == 4
%         xlim([-0.05, 1])
%     end
%     if method == 1
%         lgd = legend({'Time Cells'}, ...
%             'Location', 'northeast');
%         lgd.FontSize = figureDetails.fontSize-3;
%     end
%     set(gca, 'FontSize', figureDetails.fontSize)
% 
%     if method == 1
%         subplotNum = [(19:20), (28:29)];
%     elseif method == 2
%         subplotNum = [(22:23), (31:32)];
%     elseif method == 3
%         subplotNum = [(73:74), (82:83)];
%     elseif method == 4
%         subplotNum = [(76:77), (85:86)];
%     end
% 
%     subplot(21, 9, subplotNum)
%     histogram(ocScores, ...
%         'BinWidth', binWidth, ...
%         'EdgeColor', C(1, :), ...
%         'FaceColor', C(1, :))
%     ylim([0 108])
% 
%     if method == 3
%         ylabel('Counts', ...
%             'FontSize', figureDetails.fontSize, ...
%             'FontWeight', 'bold')
%     end
%     hold off
%     if method == 3
%         xlabel('Binned Norm. Scores', ...
%             'FontSize', figureDetails.fontSize, ...
%             'FontWeight', 'bold')
%     end
%     set(gca, 'yticklabel', {'0', '100', '200'})
% 
%     if method == 1
%         xlim([-0.1, 1])
%     elseif method == 2
%         xlim([-0.2, 0.6])
%     elseif method == 3
%         xlim([-0.001, 0.04])
%     elseif method == 4
%         xlim([-0.05, 1])
%     end
% 
%     if method == 1
%         lgd = legend({'Other Cells'}, ...
%             'Location', 'northeast');
%         lgd.FontSize = figureDetails.fontSize-3;
%     end
%     set(gca, 'FontSize', figureDetails.fontSize)
% 
% end
% hold off
% 
% subplot(21, 9, [(61:63), (70:72), (79:81), (88:90)])
% %for method = input.nMethods:-1:1
% for method = 1:input.nMethods
%     %disp(method)
%     %Linear Regression
%     clear Mdl
%     Mdl = fitglm(predictor(:, method), response, ...
%         'Distribution', 'binomial', ...
%         'Link', 'logit');
% 
%     %Get the scores
%     score_log = Mdl.Fitted.Probability;
% 
%     %ROC Curve coordinates
%     [Xlog, Ylog, Tlog, AUClog, optOP] = perfcurve(response, score_log, 1);
%     %disp(Xlog)
%     %disp(Ylog)
%     %disp(optOP)
% 
%     %     if isnan(Xlog) || isnan(Ylog)
%     %         warning('Probably skipping ...')
%     %     end
% 
%     %Plots
%     plot(Xlog, Ylog, '-.', ...
%         'Color', C1(method, :), ...
%         'LineWidth', figureDetails.lineWidth)
%     hold on
% 
%     if plotOptimalPoint
%         try
%             plot(optOP(1), optOP(2), 'Color', C1(method, :), 'o')
%         catch
%             warning('Unable to plot optimal operational point')
%         end
%         hold on
%     end
% end
% hold off
% % title(sprintf('ROC Curves (N=%i)', input.nDatasets), ...
% %     'FontSize', figureDetails.fontSize, ...
% %     'FontWeight', 'bold')
% title('ROC Curves', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xlabel('False Positive Rate', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% ylabel('True Positive Rate', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% lgd1 = legend(methodLabels, ...
%     'Location', 'southeast');
% lgd1.FontSize = figureDetails.fontSize-3;
% 
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(21, 9, [(7:9), (16:18), (25:27), (34:36)])
% % h = heatmap(correlationMatrix, ...
% %     'Colormap', linspecer, ...
% %     'Title', 'Corr. Coeff. - Scores', ...
% %     'CellLabelColor','none');
% % h.XDisplayLabels = methodLabels;
% % h.YDisplayLabels = methodLabels;
% imagesc(correlationMatrix)
% colormap(figureDetails.colorMap)
% %colormap(C)
% %colormap(jet)
% colorbar
% title('Corr. Coeff. - Scores', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xticklabels(methodLabels)
% yticklabels(methodLabels)
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% if workingOnServer == 2
%     normalizeEachCell = 1;
%     allCellsTrialAvg = [];
%     if ~exist('sdo_batch', 'var')
%         disp('[INFO] Loading a large batch of datasets ...')
%         load(sprintf('%sWork/Analysis/RHO/M26/20180514/synthDATA_%i_gRun%i_batch_%i.mat', HOME_DIR2, gDate, gRun, nDatasets))
%         disp('[INFO] ... done!')
%     end
%     for dataseti = 1:input.nDatasets
%         currentDatasetTrialAvg = squeeze(mean(sdo_batch(dataseti).syntheticDATA, 2));
%         allCellsTrialAvg = [allCellsTrialAvg; currentDatasetTrialAvg]; %Not properly preallocated because of unknown no. of classified TCs
%     end
%     if normalizeEachCell
%         for celli = 1:size(allCellsTrialAvg, 1)
%             maxVal = max(squeeze(allCellsTrialAvg(celli, :)));
%             allCellsTrialAvg(celli, :) = allCellsTrialAvg(celli, :)/maxVal;
%         end
%     end
% 
%     allTimeCellsTrialAvg = allCellsTrialAvg(find(prediction(:, selectedAlgo)), :);
%     allOtherCellsTrialAvg = allCellsTrialAvg(find(~prediction(:, selectedAlgo)), :);
% 
%     %Sort based on Peak - Time Cells
%     [~, I1] = max(allTimeCellsTrialAvg, [], 2); %find peak for each cell, along frames; [M, I]
%     [~, I2] = sort(I1);
%     sortedAllTimeCellsTrialAvg = allTimeCellsTrialAvg(I2, :);
% 
%     %Sort based on Peak - Other Cells
%     [~, I3] = max(allOtherCellsTrialAvg, [], 2); %find peak for each cell, along frames; [M, I]
%     [~, I4] = sort(I3);
%     sortedAllOtherCellsTrialAvg = allOtherCellsTrialAvg(I4, :);
% 
%     subplot(21, 9, (109:144))
%     imagesc(sortedAllTimeCellsTrialAvg*100)
%     colormap(figureDetails.colorMap)
%     %colormap(C)
%     title(sprintf('Synthetic Data | %s | Peak Sorted', algoLabels{selectedAlgo}), ...
%         'FontSize', figureDetails.fontSize, ...
%         'FontWeight', 'bold')
%     xlabel('', ...
%         'FontSize', figureDetails.fontSize, ...
%         'FontWeight', 'bold')
%     ylabel('All Time Cells', ...
%         'FontSize', figureDetails.fontSize, ...
%         'FontWeight', 'bold')
%     z = colorbar;
%     if normalizeEachCell
%         ylabel(z, 'Norm. dF/F (%)', ...
%             'FontSize', figureDetails.fontSize, ...
%             'FontWeight', 'bold')
%     else
%         ylabel(z,'dF/F (%)', ...
%             'FontSize', figureDetails.fontSize, ...
%             'FontWeight', 'bold')
%     end
%     set(gca, 'FontSize', figureDetails.fontSize)
% 
%     randomI = sort(randi([1 size(sortedAllOtherCellsTrialAvg, 1)], 1, size(sortedAllTimeCellsTrialAvg, 1)));
%     subplot(21, 9, (155:189))
%     imagesc(sortedAllOtherCellsTrialAvg(randomI, :)*100)
%     colormap(figureDetails.colorMap)
%     %colormap(C)
%     xlabel('Frames', ...
%         'FontSize', figureDetails.fontSize, ...
%         'FontWeight', 'bold')
%     ylabel('Other Cells (e.g.)', ...
%         'FontSize', figureDetails.fontSize, ...
%         'FontWeight', 'bold')
%     z = colorbar;
%     if normalizeEachCell
%         ylabel(z, 'Norm. dF/F (%)', ...
%             'FontSize', figureDetails.fontSize, ...
%             'FontWeight', 'bold')
%     else
%         ylabel(z,'dF/F (%)', ...
%             'FontSize', figureDetails.fontSize, ...
%             'FontWeight', 'bold')
%     end
%     minVal = min(min(sortedAllTimeCellsTrialAvg))*100;
%     maxVal = max(max(sortedAllTimeCellsTrialAvg))*100;
%     %caxis([minVal maxVal])
%     set(gca, 'FontSize', figureDetails.fontSize)
% end
% print(sprintf('%s/figs/4-ComparingScores-%i-%i-%i-%i-%i_%i', ...
%     HOME_DIR2, ...
%     input.gDate, ...
%     input.gRun, ...
%     input.nDatasets, ...
%     workingOnServer), ...
%     '-dpng')
% 
% disp('... Scores and Comparisons Plotted')

% %% Performance Metrics
% 
% disp('Plotting Performance Metrics ...')
% %results1 is ordered as TP, FP, TN, FN
% results1_mod(:, (1:2))= results1(:, [1, 4]);
% results1_mod(:, (3:4)) = results1(:, [3, 2]);
% results1_mod(:, [2, 4]) = results1_mod(:, [2, 4])*(-1);
% 
% results1_mod2 = results1;
% results1_mod2(:, [2, 4]) = results1(:, [2, 4])*(-1);
% 
% fig5 = figure(5);
% clf
% set(fig5, 'Position', [500, 500, 900, 1200])
% subplot(10, 8, (1:16))
% if bandClasses
%     d1 = bar(results1_mod2(:, (1:2)), 'stacked');
% else
%     d1 = bar(results1_mod(:, (1:2)), 'stacked');
% end
% d1(1).FaceColor = C1(3, :);
% d1(2).FaceColor = C1(4, :);
% %xlim([0 5])
% axis tight
% title(sprintf('Ground Truth: Time Cells (N=%i)', sum(response)), ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% % xlabel('All Algorithms', ...
% % 'FontSize', figureDetails.fontSize, ...
% %     'FontWeight', 'bold')
% ylabel('Rate', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% ylim([-1, 1])
% yticks(-1:0.5:1)
% yticklabels({'1.0', '0.5', '0.0', '0.5', '1.0'})
% xticklabels({'', '', '', '', '', '', '', '', })
% xticklabels(algoLabels)
% %xtickangle(45)
% if bandClasses
%     lgd = legend(metricLabels_pos, ...
%         'Location', 'southeastoutside');
% else
%     lgd = legend(metricLabels_timeCells, ...
%         'Location', 'southeastoutside');
% end
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(11, 8, (25:40))
% if bandClasses
%     d2 = bar(results1_mod2(:, (3:4)), 'stacked');
% else
%     d2 = bar(results1_mod(:, (3:4)), 'stacked');
% end
% d2(1).FaceColor = C1(2, :);
% d2(2).FaceColor = C1(1, :);
% %xlim([0 5])
% axis tight
% title(sprintf('Ground Truth: Other Cells (N=%i)', length(response)-sum(response)), ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% % xlabel('All Algorithms', ...
% % 'FontSize', figureDetails.fontSize, ...
% %     'FontWeight', 'bold')
% ylabel('Rate', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% ylim([-1, 1])
% yticks(-1:0.5:1)
% yticklabels({'1.0', '0.5', '0.0', '0.5', '1.0'})
% xticklabels({'', '', '', '', '', '', '', '', })
% xticklabels(algoLabels)
% xtickangle(45)
% if bandClasses
%     lgd = legend(metricLabels_neg, ...
%         'Location', 'southeastoutside');
% else
%     lgd = legend(metricLabels_otherCells, ...
%         'Location', 'southeastoutside');
% end
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% %subplot(11, 8, (25:40))
% subplot(10, 8, (49:56))
% d = bar(results3);
% modICol = [3 4 2 1]; %[1:4] for default
% colororder(C1(modICol, :))
% %xlim([0 5])
% axis tight
% title(sprintf('Predictive Performance Metrics (N=%i)', length(response)), ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% % xlabel('All Algorithms', ...
% % 'FontSize', figureDetails.fontSize, ...
% %     'FontWeight', 'bold')
% ylabel('Rate', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% ylim([0, 1])
% yticks(0:0.5:1)
% yticklabels({'0.0', '0.5', '1.0'})
% xticklabels({'', '', '', '', '', '', '', '', })
% xticklabels(algoLabels)
% xtickangle(45)
% lgd = legend(metricLabels2, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(10, 8, [65:66, 73:74])
% corrMatPredictionsAll = corrcoef(prediction, 'Rows', 'pairwise');
% h = heatmap(corrMatPredictionsAll, ...
%     'Title', 'Corr. Coeff. - Predictions', ...
%     'XDisplayLabels', algoLabels, ...
%     'YDisplayLabels', algoLabels);
% s = struct(h);
% s.XAxis.TickLabelRotation = 45;
% colormap(figureDetails.colorMap)
% colorbar
% % title('Corr. Coeff. - Predictions', ...
% %     'FontSize', figureDetails.fontSize, ...
% %     'FontWeight', 'bold')
% % xlabel(algoLabels, ...
% %     'FontSize', figureDetails.fontSize, ...
% %     'FontWeight', 'bold')
% % ylabel(algoLabels, ...
% %     'FontSize', figureDetails.fontSize, ...
% %      'FontWeight', 'bold')
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% meanInUse = [0.0000000, 8002165, 8002031, 8002165; 0, 16035957, 16035894, 16035957]; %bytes
% runTime = [0.0000000, 1.013099479675293, 1.9082653522491455, 0.03733084201812744; 0, 2.0923821687698365, 3.8474811792373655, 0.06629612445831298]; %secs
% 
% %Populate Synthesis Numbers
% if exist('totalMem', 'var')
%     meanInUse(2, 1) = (totalMem/input.nDatasets)*(1024^2); %in bytes now; dirty hack
%     meanInUse(1, 1) = meanInUse(2, 1)/2;
% 
% else
%     synthMemUsage = load(sprintf('%sWork/Analysis/RHO/M26/20180514/synthDATA_20220618_gRun1_batch_567.mat', HOME_DIR2), 'totalMem');
%     meanInUse(2, 1) = (synthMemUsage.totalMem/input.nDatasets)*(1024^2); %in bytes now; dirty hack
%     meanInUse(1, 1) = meanInUse(2, 1)/2;
% end
% 
% if exist('elapsedTime', 'var')
%     runTime(2, 1) = elapsedTime/input.nDatasets;
%     runTime(1, 1) = runTime(2, 1)/2;
% else
%     synthRuntime = load(sprintf('%sWork/Analysis/RHO/M26/20180514/synthDATA_20220618_gRun1_batch_567.mat', HOME_DIR2), 'elapsedTime');
%     runTime(2, 1) = synthRuntime.elapsedTime/input.nDatasets;
%     runTime(1, 1) = runTime(2, 1)/2;
% end
% 
% %Usage
% subplot(10, 8, [69, 77])
% b1 = bar(meanInUse'/(1024^2));
% %b1.FaceColor = [0.8500 0.3250 0.0980];
% %b1.FaceColor = C(2,:);
% modICol = [3 4 2 1]; %[1:4] for default
% colororder(C1(modICol, :))
% axis tight
% title('Memory Usage', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% ylabel('Workspace (MB)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% %xticks([1, 2, 3, 4])
% ylim([0, 42])
% %yticks([0, 200, 400, 600])
% % lgd = legend({'135 cells', '540 cells', '1350 cells'}, ...
% %     'Location', 'southeastoutside');
% % lgd.FontSize = figureDetails.fontSize-4;
% % lgd = legend({'', ''}, ...
% %     'Location', 'southwestoutside');
% legend('hide')
% xticklabels(procedureLabels)
% xtickangle(45)
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(10, 8, [71:72, 79:80])
% %boxplot(runTime, procedureLabels)
% %errorbar(meanRunTime', stdRunTime', 'r*', 'CapSize', 10)
% b2 = bar(runTime');
% %b2.FaceColor = [0.8500 0.3250 0.0980];
% %b2.FaceColor = C(2,:);
% modICol = [3 4 2 1]; %[1:4] for default
% colororder(C1(modICol, :))
% axis tight
% title('Runtimes', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% ylabel('Time (sec)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% ylim([0, 5])
% lgd = legend({'67 cells', '135 cells'}, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% xticklabels(procedureLabels)
% xtickangle(45)
% 
% yticks([0, 1, 2, 3, 4])
% yticklabels({'0', '1', '2', '3', '4'})
% 
% set(gca, 'FontSize', figureDetails.fontSize)
% print(sprintf('%s/figs/5-PerformanceEvaluation-%i-%i-%i_%i', ...
%     HOME_DIR2, ...
%     input.gDate, ...
%     input.gRun, ...
%     input.nDatasets, ...
%     workingOnServer), ...
%     '-dpng')
% 
% disp('... Performance Metrics Plotted')

%% Sensitivity and Resource

disp('Plotting Sensitivity and Resource ...')
xTicks = 1:input.nAlgos;
joinBundledResults3_base = [];
joinBundledResults3_mod = [];
joinAllResults3 = [];
nSets = 3;

if datasetCatalog == 0
    %Cell wise - Multiple dataset by nCells (=135) to get the following
    iNoise1 = (56296:57645); %low
    iNoise2 = (57646:58995); %medium
    iNoise3 = (58996:60345); %high

    iEW1 = (60346:61695); %small
    iEW2 = (61696:63045); %medium
    iEW3 = (63046:64395); %large

    iImp1 = (64396:65745); %low
    iImp2 = (65746:67095); %medium
    iImp3 = (67096:68445); %high

    iHTR1 = (68446:69795); %low
    iHTR2 = (69796:71145); %medium
    iHTR3 = (71146:72495); %high

    iBackgnd1 = (72496:73845); %low
    iBackgnd2 = (73846:75195); %medium
    iBackgnd3 = (75196:76545); %high

    %Parameter Ranges - see configSynth.m or setupSynthDataParams8.m
    noiseVal1 = 10; % in %
    noiseVal2 = 40; % in %
    noiseVal3 = 70; % in %

    ewVal1 = 30; % in %ile
    ewVal2 = 60; % in %ile
    ewVal3 = 90; % in %ile

    impVal1 = 0; % in frames
    impVal2 = 33; % in frames
    impVal3 = 66; % in frames

    htrVal1 = 33; % in %
    htrVal2 = 66; % in %
    htrVal3 = 99; % in %

    backVal1 = 0.15; %average lambda value
    backVal2 = 0.6; %average lambda value
    backVal3 = 1.05; %lambda value

    caseResponse = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3];
end
if datasetCatalog == 0
    nShuffles = 10; %same throughout
else
    nShuffles = 5; %same throughout
end

dependenceAll = zeros(input.nAlgos, input.nParams);
groupAlgos = 1:input.nAlgos;
groupParams = 1:input.nParams;
bundledResults3 = zeros(nShuffles, input.nAlgos, nSets);
allResults3 = zeros(nShuffles, input.nAlgos, nSets, 3);

X1 = prediction(iNoise1, :);
X2 = prediction(iNoise2, :);
X3 = prediction(iNoise3, :);

Y1 = Y(iNoise1);
Y2 = Y(iNoise2);
Y3 = Y(iNoise3);

for iSet = 1:nSets
    if iSet == 1
        Xeffect = X1;
        Yeffect = Y1;
    elseif iSet == 2
        Xeffect = X2;
        Yeffect = Y2;
    elseif iSet == 3
        Xeffect = X3;
        Yeffect = Y3;
    end

    for shuffle = 1:nShuffles
        myCells = ((input.nCells*(shuffle-1))+1):(input.nCells*shuffle);
        clear results1
        clear results2
        clear results3
        [results1, results2, results3] = compareAgainstTruth(Xeffect(myCells, :), Yeffect(myCells), input);
        bundledResults3(shuffle, :, iSet) = results3(:, 3);
        allResults3(shuffle, :, iSet, :) = results3;
    end
end
joinAllResults3 = [joinAllResults3; squeeze(allResults3(:, :, 1, :))]; %For baseline - Noise expts.
joinBundledResults3_base = [joinBundledResults3_base; squeeze(bundledResults3(:, :, 1))]; %For baseline - Noise expts.
joinBundledResults3_mod = [joinBundledResults3_mod; bundledResults3];
meanBundledResults3 = squeeze(mean(bundledResults3, 1, 'omitnan'));
stdBundledResults3 = squeeze(std(bundledResults3, 1, 'omitnan'));

mBR_Noise = meanBundledResults3; %For Supplementary
sBR_Noise = stdBundledResults3; %For Supplementary

X1 = prediction(iEW1, :);
X2 = prediction(iEW2, :);
X3 = prediction(iEW3, :);

Y1 = Y(iEW1);
Y2 = Y(iEW2);
Y3 = Y(iEW3);

for iSet = 1:nSets
    if iSet == 1
        Xeffect = X1;
        Yeffect = Y1;
    elseif iSet == 2
        Xeffect = X2;
        Yeffect = Y2;
    elseif iSet == 3
        Xeffect = X3;
        Yeffect = Y3;
    end

    for shuffle = 1:nShuffles
        myCells = ((input.nCells*(shuffle-1))+1):(input.nCells*shuffle);
        clear results1
        clear results2
        clear results3
        [results1, results2, results3] = compareAgainstTruth(Xeffect(myCells, :), Yeffect(myCells), input);
        bundledResults3(shuffle, :, iSet) = results3(:, 3);
        allResults3(shuffle, :, iSet, :) = results3;
    end
end
joinAllResults3 = [joinAllResults3; squeeze(allResults3(:, :, 2, :))]; %For baseline - EW expts.
joinBundledResults3_base = [joinBundledResults3_base; squeeze(bundledResults3(:, :, 2))]; %For baseline - EW expts.
joinBundledResults3_mod = [joinBundledResults3_mod; bundledResults3];
meanBundledResults3 = squeeze(mean(bundledResults3, 1, 'omitnan'));
stdBundledResults3 = squeeze(std(bundledResults3, 1, 'omitnan'));

mBR_EW = meanBundledResults3; %For Supplementary
sBR_EW = stdBundledResults3; %For Supplementary

X1 = prediction(iImp1, :);
X2 = prediction(iImp2, :);
X3 = prediction(iImp3, :);

Y1 = Y(iImp1);
Y2 = Y(iImp2);
Y3 = Y(iImp3);

for iSet = 1:nSets
    if iSet == 1
        Xeffect = X1;
        Yeffect = Y1;
    elseif iSet == 2
        Xeffect = X2;
        Yeffect = Y2;
    elseif iSet == 3
        Xeffect = X3;
        Yeffect = Y3;
    end

    for shuffle = 1:nShuffles
        myCells = ((input.nCells*(shuffle-1))+1):(input.nCells*shuffle);
        clear results1
        clear results2
        clear results3
        [results1, results2, results3] = compareAgainstTruth(Xeffect(myCells, :), Yeffect(myCells), input);
        bundledResults3(shuffle, :, iSet) = results3(:, 3);
        allResults3(shuffle, :, iSet, :) = results3;
    end
end
joinAllResults3 = [joinAllResults3; squeeze(allResults3(:, :, 1, :))]; %For baseline - Imprecision expts.
joinBundledResults3_base = [joinBundledResults3_base; squeeze(bundledResults3(:, :, 1))]; %For baseline - Imprecision expts.
joinBundledResults3_mod = [joinBundledResults3_mod; bundledResults3];
meanBundledResults3 = squeeze(mean(bundledResults3, 1, 'omitnan'));
stdBundledResults3 = squeeze(std(bundledResults3, 1, 'omitnan'));

mBR_Imp = meanBundledResults3; %For Supplementary
sBR_Imp = stdBundledResults3; %For Supplementary

X1 = prediction(iHTR1, :);
X2 = prediction(iHTR2, :);
X3 = prediction(iHTR3, :);

Y1 = Y(iHTR1);
Y2 = Y(iHTR2);
Y3 = Y(iHTR3);

for iSet = 1:nSets
    if iSet == 1
        Xeffect = X1;
        Yeffect = Y1;
    elseif iSet == 2
        Xeffect = X2;
        Yeffect = Y2;
    elseif iSet == 3
        Xeffect = X3;
        Yeffect = Y3;
    end

    for shuffle = 1:nShuffles
        myCells = ((input.nCells*(shuffle-1))+1):(input.nCells*shuffle);
        clear results1
        clear results2
        clear results3
        [results1, results2, results3] = compareAgainstTruth(Xeffect(myCells, :), Yeffect(myCells), input);
        bundledResults3(shuffle, :, iSet) = results3(:, 3);
        allResults3(shuffle, :, iSet, :) = results3;
    end
end
joinAllResults3 = [joinAllResults3; squeeze(allResults3(:, :, 2, :))]; %For baseline - HTR expts.
joinBundledResults3_base = [joinBundledResults3_base; squeeze(bundledResults3(:, :, 2))]; %For baseline - HTR expts.
joinBundledResults3_mod = [joinBundledResults3_mod; bundledResults3];
meanBundledResults3 = squeeze(mean(bundledResults3, 1, 'omitnan'));
stdBundledResults3 = squeeze(std(bundledResults3, 1, 'omitnan'));

mBR_HTR = meanBundledResults3; %For Supplementary
sBR_HTR = stdBundledResults3; %For Supplementary

X1 = prediction(iBackgnd1, :);
X2 = prediction(iBackgnd2, :);
X3 = prediction(iBackgnd3, :);

Y1 = Y(iBackgnd1);
Y2 = Y(iBackgnd2);
Y3 = Y(iBackgnd3);

for iSet = 1:nSets
    if iSet == 1
        Xeffect = X1;
        Yeffect = Y1;
    elseif iSet == 2
        Xeffect = X2;
        Yeffect = Y2;
    elseif iSet == 3
        Xeffect = X3;
        Yeffect = Y3;
    end

    for shuffle = 1:nShuffles
        myCells = ((input.nCells*(shuffle-1))+1):(input.nCells*shuffle);
        clear results1
        clear results2
        clear results3
        [results1, results2, results3] = compareAgainstTruth(Xeffect(myCells, :), Yeffect(myCells), input);
        bundledResults3(shuffle, :, iSet) = results3(:, 3);
        allResults3(shuffle, :, iSet, :) = results3;
    end
end
joinAllResults3 = [joinAllResults3; squeeze(allResults3(:, :, 2, :))]; %For baseline - Back. Act. expts.
joinBundledResults3_base = [joinBundledResults3_base; squeeze(bundledResults3(:, :, 2))]; %For baseline - Back. Act. expts.
joinBundledResults3_mod = [joinBundledResults3_mod; bundledResults3];
meanBundledResults3 = squeeze(mean(bundledResults3, 1, 'omitnan'));
stdBundledResults3 = squeeze(std(bundledResults3, 1, 'omitnan'));

mBR_Back = meanBundledResults3; %For Supplementary
sBR_Back = stdBundledResults3; %For Supplementary

meanJoinAllResults3 = squeeze(mean(joinAllResults3, 1, 'omitnan'));
stderrJoinAllResults3 = squeeze(std(joinAllResults3, 1, 'omitnan')/sqrt(size(joinAllResults3, 1)));
meanJoinBundledResults3 = squeeze(mean(joinBundledResults3_base, 1, 'omitnan'));
stderrJoinBundledResults3 = squeeze(std(joinBundledResults3_base, 1, 'omitnan')/sqrt(size(joinBundledResults3_base, 1)));

depSig %Simple Linear Regression using fitlm() - Custom

fig6 = figure(6);
clf
set(fig6, 'Position', [1050, 1000, 900, 1200])

subplot(21, 2, (1:4))
b1 = bar(meanJoinAllResults3);
modICol = [3 4 2 1]; %[1:4] for default
colororder(C1(modICol, :))
%colororder(C1)
axis tight
title('Baseline - Physiological Regime', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
ylabel('Rate', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels(algoLabels);
xtickangle(45)
lgd = legend(metricLabels2, ...
    'Location', 'southeastoutside');
lgd.FontSize = figureDetails.fontSize-3;
set(gca, 'FontSize', figureDetails.fontSize)

subplot(21, 2, [9, 11, 13])
myYR = reshape(squeeze(joinBundledResults3_mod(1:10, selectedAlgo, :)), [], 1); %For Noise (%) using tiBoth
mdl = fitlm(caseResponse, myYR, 'linear');
h = plot(mdl);
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = C(5, :);
dataHandle.Marker = 'o';
dataHandle.MarkerSize = figureDetails.markerSize;
fitHandle.Color = C(1, :);
set(cbHandles, 'Color', C(2, :), 'LineWidth', figureDetails.lineWidth)
% hold on
% %plot([1, 0.5, 0])
% plot([mBR_Noise(selectedAlgo, 1), (mBR_Noise(selectedAlgo, 3)+mBR_Noise(selectedAlgo, 1))/2, mBR_Noise(selectedAlgo, 3)], '--k')
% hold off
title('Dependence Schematic', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xlabel('Noise (%)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels({num2str(noiseVal1), '', num2str(noiseVal2), '', num2str(noiseVal3)})
ylabel('Adj. F1 Score', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
lgd = legend({algoLabels{selectedAlgo}, 'Fit', '95% PI'}, ...
    'Location', 'southeastoutside');
lgd.FontSize = figureDetails.fontSize-3;
set(gca, 'FontSize', figureDetails.fontSize)

dependence1 = slopes(1, :);
dependence1_err_neg = allRMSE(1, :)/sqrt(9);
dependence1_err_pos = dependence1_err_neg;
 
subplot(21, 2, [12, 14])
%yyaxis right
b2 = bar(dependence1);
%b2.FaceColor = [0.8500 0.3250 0.0980];
b2.FaceColor = C1(2, :);
ylabel('\Delta F1 Score/\Delta Noise', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
hold on
er = errorbar(xTicks, dependence1, dependence1_err_neg, dependence1_err_pos, ...
    'MarkerSize', figureDetails.markerSize, ...
    'CapSize', figureDetails.capSize);
er.Color = C(5, :);
er.LineStyle = 'none';
hold off
axis tight
title('Noise (%)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels(algoLabels);
xtickangle(45)
set(gca, 'FontSize', figureDetails.fontSize)
ylim([-0.7, 0.1])

dependence2 = slopes(2, :);
dependence2_err_neg = allRMSE(2, :)/sqrt(9);
dependence2_err_pos = dependence2_err_neg;

subplot(21, 2, [21, 23])
%yyaxis right
b2 = bar(dependence2);
%b2.FaceColor = [0.8500 0.3250 0.0980];
b2.FaceColor = C1(2, :);
ylabel('\Delta F1 Score/\Delta EW', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
%set(b2,'FaceAlpha', 0.5)
hold on
er = errorbar(xTicks, dependence2, dependence2_err_neg, dependence2_err_pos, ...
    'MarkerSize', figureDetails.markerSize, ...
    'CapSize', figureDetails.capSize);
er.Color = C(5, :);
er.LineStyle = 'none';
hold off
axis tight
title('EW (%ile)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels(algoLabels);
xtickangle(45)
set(gca, 'FontSize', figureDetails.fontSize)
ylim([-0.022, 0.18])

dependence3 = slopes(3, :);
dependence3_stderr_neg = allRMSE(3, :)/sqrt(9);
dependence3_stderr_pos = dependence3_stderr_neg;

subplot(21, 2, [22, 24])
b2 = bar(dependence3);
%b2.FaceColor = [0.8500 0.3250 0.0980];
b2.FaceColor = C1(2, :);
ylabel('\Delta F1 Score/\Delta Imp.', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
hold on
er = errorbar(xTicks, dependence3, dependence3_stderr_neg, dependence3_stderr_pos, ...
    'MarkerSize', figureDetails.markerSize, ...
    'CapSize', figureDetails.capSize);
er.Color = C(5, :);
er.LineStyle = 'none';
hold off
axis tight
title('Imp. (frames)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels(algoLabels);
xtickangle(45)
set(gca, 'FontSize', figureDetails.fontSize)
ylim([-0.72, 0.1])

dependence4 = slopes(4, :);
dependence4_err_neg = allRMSE(4, :)/sqrt(9);
dependence4_err_pos = dependence4_err_neg;

subplot(21, 2, [31, 33])
b2 = bar(dependence4);
%b2.FaceColor = [0.8500 0.3250 0.0980];
b2.FaceColor = C1(2, :);
ylabel('\Delta F1 Score/\Delta HTR', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
hold on
er = errorbar(xTicks, dependence4, dependence4_err_neg, dependence4_err_pos, ...
    'MarkerSize', figureDetails.markerSize, ...
    'CapSize', figureDetails.capSize);
er.Color = C(5, :);
er.LineStyle = 'none';
hold off
axis tight
title('HTR (%)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels(algoLabels);
xtickangle(45)
set(gca, 'FontSize', figureDetails.fontSize)
ylim([-0.03, 0.26])

dependence5 = slopes(5, :);
dependence5_err_neg = allRMSE(5, :)/sqrt(9);
dependence5_err_pos = dependence5_err_neg;

subplot(21, 2, [32, 34])
b2 = bar(dependence5);
%b2.FaceColor = [0.8500 0.3250 0.0980];
b2.FaceColor = C1(2, :);
ylabel('\Delta F1 Score/\Delta Back.', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
hold on
er = errorbar(xTicks, dependence5, dependence5_err_neg, dependence5_err_pos, ...
    'MarkerSize', figureDetails.markerSize, ...
    'CapSize', figureDetails.capSize);
er.Color = C(5, :);
er.LineStyle = 'none';
hold off
axis tight
title('Background (\lambda)', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels(algoLabels);
xtickangle(45)
set(gca, 'FontSize', figureDetails.fontSize)
ylim([-0.22, 0.03])

% Concordance - slightly different algorithm from just using compareAgainstTruth()
X4 = prediction; %All cells
sumX = sum(X4, 2);
nCases = length(Y);

% Now, for any cell, if sumX >= some threshold, we suggest calling it a Time Cell
X4 = zeros(nCases, nConcordanceAlgos);

%Preallocation
allTP = zeros(nCases, nConcordanceAlgos);
allTN = zeros(nCases, nConcordanceAlgos);
allFP = zeros(nCases, nConcordanceAlgos);
allFN = zeros(nCases, nConcordanceAlgos);
TPR = zeros(nConcordanceAlgos, 1); %Sensitivity, Recall, or True Positive Rate
FPR = zeros(nConcordanceAlgos, 1); %Fall-out or False Positive Rate
TNR = zeros(nConcordanceAlgos, 1); %Specificity or True Negative Rate
FNR = zeros(nConcordanceAlgos, 1); %Miss Rate or False Negative Rate
PPV = zeros(input.nAlgos, 1); %Precision or Positive Predictive Value

results4 = zeros(nConcordanceAlgos, 3); %3 because we'll look at Recall, Precision, and F1 Score
startCase = 56296; %Onwards down the row
for algo = 1:nConcordanceAlgos

    X4(:, algo) = (sumX >= algo);

    for myCase = startCase:nCases
        if Y(myCase) %Positive Cases
            if X4(myCase, algo)
                allTP(myCase, algo) = 1;
            else
                allFN(myCase, algo) = 1;
            end
        else %Negative Cases
            if X4(myCase, algo)
                allFP(myCase, algo) = 1;
            else
                allTN(myCase, algo) = 1;
            end
        end
    end

    TP = sum(allTP(:, algo));
    FN = sum(allFN(:, algo));
    FP = sum(allFP(:, algo));
    TN = sum(allTN(:, algo));

    if TP ~= 0
        TPR(algo) = TP/(TP + FN); %Recall
        PPV(algo) = TP/(TP + FP); %Precision
    else
        TPR(algo) = 0; %Recall
        PPV(algo) = 0; %Precision
    end

    if FP ~= 0
        FPR(algo) = FP/(FP + TN);
    else
        FPR(algo) = 0;
    end

    if TN ~= 0
        TNR(algo) = TN/(TN + FP);
    else
        TNR(algo) = 0;
    end

    if FN ~= 0
        FNR(algo) = FN/(FN + TP);
    else
        FNR(algo) = 0;
    end

    results4(algo, 1) = TPR(algo); %Recall
    results4(algo, 2) = PPV(algo); %Precision

    if TP == 0
        results4(algo, 3) = 0;
    else
        results4(algo, 3) = 2 * results4(algo, 1) * results4(algo, 2)/(results4(algo, 1) + results4(algo, 2)); %F1 Score
    end
end

subplot(21, 2, (39:42))
d = bar(results4);
modICol = [3 4 2 1]; %[1:4] for default
colororder(C1(modICol, :))
%colororder(C1)
xlim([0 5])
axis tight
xlabel('Classification Concordance Threshold', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
ylabel('Rate', ...
    'FontSize', figureDetails.fontSize, ...
    'FontWeight', 'bold')
xticklabels(concordanceAlgoLabels)
xtickangle(45)
lgd = legend(metricLabels2, ...
    'Location', 'southeastoutside');
lgd.FontSize = figureDetails.fontSize-3;
set(gca, 'FontSize', figureDetails.fontSize)

print(sprintf('%s/figs/6-Resource-%i-%i-%i-%i-%i_%i', ...
    HOME_DIR2, ...
    input.gDate, ...
    input.gRun, ...
    input.nDatasets, ...
    workingOnServer), ...
    '-dpng')

disp('... Sensitivity and Resource Plotted.')

%Save Dependence values for python/MATLAB
dependenceTable = array2table(dependenceAll);
writetable(dependenceTable, 'physiologyDependence.csv');

% %% Sup. Fig 1
% disp('Plotting Sup. Fig 1 ...')
% %Supplementary Figure - Fig16
% fig16 = figure(16);
% clf
% set(fig16, 'Position', [100, 100, 900, 1200])
% %colororder();
% subplot(10, 2, [1, 3])
% %Equivalence
% myCellIds = 56296:76545;
% [equivalence, equivalenceSplits] = getEquivalenceMatches(prediction(myCellIds, 1:6), response(myCellIds, 1));
% %h = histogram(equivalenceSplits, 'FaceColor', C([6, 1], :));
% b = bar(equivalenceSplits, '');
% %colororder(C1(modICol, :))
% axis tight
% title('Equivalence (XNOR)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% %xticks = ([1, 2]);
% xticklabels({'', 'No Match', '', 'Match'})
% ylabel('Counts', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% lgd = legend({'Time Cells', 'Other Cells'}, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(10, 2, [2, 4])
% errorbar(mBR_Noise', sBR_Noise', ...
%     'LineWidth', figureDetails.lineWidth, ...
%     'CapSize', 10)
% title('Dependence', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xlabel('Noise (%)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xticklabels({num2str(noiseVal1), '', num2str(noiseVal2), '', num2str(noiseVal3)})
% ylabel('F1 Score', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% lgd = legend(algoLabels, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(10, 2, [9, 11])
% errorbar(mBR_EW', sBR_EW', ...
%     'LineWidth', figureDetails.lineWidth, ...
%     'CapSize', 10)
% title('Dependence', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xlabel('Event Width (%ile)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xticklabels({num2str(ewVal1), '', num2str(ewVal2), '', num2str(ewVal3)})
% ylabel('F1 Score', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% lgd = legend(algoLabels, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(10, 2, [10, 12])
% errorbar(mBR_Imp', sBR_Imp', ...
%     'LineWidth', figureDetails.lineWidth, ...
%     'CapSize', 10)
% title('Dependence', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xlabel('Imp (frames)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xticklabels({num2str(impVal1), '', num2str(impVal2), '', num2str(impVal3)})
% ylabel('F1 Score', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% lgd = legend(algoLabels, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(10, 2, [17, 19])
% errorbar(mBR_HTR', sBR_HTR', ...
%     'LineWidth', figureDetails.lineWidth, ...
%     'CapSize', 10)
% title('Dependence', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xlabel('Hit Trial Ratio (%)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xticklabels({num2str(htrVal1), '', num2str(htrVal2), '', num2str(htrVal3)})
% ylabel('F1 Score', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% lgd = legend(algoLabels, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% subplot(10, 2, [18, 20])
% errorbar(mBR_Back', sBR_Back', ...
%     'LineWidth', figureDetails.lineWidth, ...
%     'CapSize', 10)
% title('Dependence', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xlabel('Background (\lambda)', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% xticklabels({num2str(backVal1), '', num2str(backVal2), '', num2str(backVal3)})
% xtickangle(0)
% ylabel('F1 Score', ...
%     'FontSize', figureDetails.fontSize, ...
%     'FontWeight', 'bold')
% lgd = legend(algoLabels, ...
%     'Location', 'southeastoutside');
% lgd.FontSize = figureDetails.fontSize-3;
% set(gca, 'FontSize', figureDetails.fontSize)
% 
% print(sprintf('%s/figs/sup6-Resource-%i-%i-%i-%i-%i_%i', ...
%     HOME_DIR2, ...
%     input.gDate, ...
%     input.gRun, ...
%     input.nDatasets, ...
%     workingOnServer), ...
%     '-dpng')
% 
% disp('... Sup. Fig1 plotted.')

%Normalize Dependence Estimates
for parami = 1:size(dependenceAll, 2)
    maxVal = max(squeeze(dependenceAll(:, parami)));
    normDependenceAll(:, parami) = dependenceAll(:, parami)/maxVal;
    fprintf('Param: %i Max: %d\n', parami, maxVal)
end

% %% Final Summary
% disp('Plotting summary spider ...')
% %for algo = 1:size(dependenceAll, 1)
% multipliers = zeros(size(dependenceAll, 2), 1);
% for parami = 1:size(dependenceAll, 2)
%     if parami == 1 || parami == 3 || parami == 5
%         minVal = min(squeeze(dependenceAll(:, parami)));
%         normDependenceAll(:, parami) = dependenceAll(:, parami)/minVal;
%         multipliers(parami) = minVal;
%     else
%         maxVal = max(squeeze(dependenceAll(:, parami)));
%         normDependenceAll(:, parami) = dependenceAll(:, parami)/maxVal;
%         multipliers(parami) = maxVal;
%     end
%     paramLabels2{parami} = [sprintf(' %0.4f', multipliers(parami)) ' / ' paramLabels{parami}];
% end
% 
% fig8 = figure(8);
% % clf
% set(fig8, 'Position', [1050, 1000, 900, 800])
% spider_plot(normDependenceAll(1:6, :), ...
%     'FillOption', 'on', ...
%     'AxesLabels', paramLabels2)
% title('Physiological Regime Dependencies for Time Cell Detection', ...
%     'FontSize', figureDetails.fontSize+1, ...
%     'FontWeight', 'bold')
% lgd = legend(algoLabels);
% lgd.FontSize = figureDetails.fontSize+1;
% print(sprintf('%s/figs/8-Summary-%i-%i-%i-%i-%i_%i', ...
%     HOME_DIR2, ...
%     input.gDate, ...
%     input.gRun, ...
%     input.nDatasets, ...
%     workingOnServer), ...
%     '-dpng')
% 
% disp('... Summary spider plotted.')

%% Supplementary to Summary
if ~exist('sdcp', 'var')
    load(sprintf('%sWork/Analysis/RHO/M26/20180514/synthDATA_%i_gRun%i_batch_%i.mat', HOME_DIR2, gDate, gRun, nDatasets), 'sdcp')
end
synthBatchProfile = getSynthBatchProfile(sdcp, 1, input.nDatasets, impX);

unphysi = 1:297;
canoni = 298:417;
physi = 418:567;

% disp('Plotting Sup. Fig 2 ...')
% fig18 = figure(18);
% clf
% set(fig18, 'Position', [100, 100, 900, 1200])
% subplot(12, 12, [1, 2, 3, 13, 14, 15, 25, 26, 27, 37, 38, 39])
% imagesc(synthBatchProfile.dProfileM2Norm(unphysi, :))
% colorbar
% colormap(figureDetails.colorMap)
% title('Modulation', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% ylabel('"Unphys." Datasets', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10 11])
% if impX ~= 0
%     xticklabels(paramLabels3_mod)
% else
%     xticklabels(paramLabels3)
% end
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% subplot(12, 12, [73, 74, 75, 85, 86, 87])
% imagesc(synthBatchProfile.dProfileM2Norm(canoni, :))
% colorbar
% colormap(figureDetails.colorMap)
% title('Modulation', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% ylabel('"Canon." Datasets', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10 11])
% if impX ~= 0
%     xticklabels(paramLabels3_mod)
% else
%     xticklabels(paramLabels3)
% end
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% subplot(12, 12, [121, 122, 123, 133, 134, 135])
% imagesc(synthBatchProfile.dProfileM2Norm(physi, :))
% colorbar
% colormap(figureDetails.colorMap)
% title('Modulation', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% ylabel('"Phys." Datasets', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10 11])
% if impX ~= 0
%     xticklabels(paramLabels3_mod)
% else
%     xticklabels(paramLabels3)
% end
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% %Confusion Matrices
% [dTP, dFP, dFN, dTN] = getConfusionMatrixPerDataset(prediction, response);
% 
% %Unphys
% %subplot(12, 12, [8, 9, 20, 21, 32, 33, 44, 45])
% subplot(12, 12, [5, 6, 7, 17, 18, 19, 29, 30, 31, 41, 42, 43])
% imagesc(dFP(unphysi, :)/input.nCells)
% colorbar
% colormap(figureDetails.colorMap)
% caxis([0 1])
% title('False Positives', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10 11])
% xticklabels(algoLabels)
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% subplot(12, 12, [9, 10, 11, 21, 22, 23, 33, 34, 35, 45, 46, 47])
% imagesc(dFN(unphysi, :)/input.nCells)
% colorbar
% colormap(figureDetails.colorMap)
% caxis([0 1])
% title('False Negatives', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10])
% xticklabels(algoLabels)
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% %Canon
% subplot(12, 12, [77, 78, 79, 89, 90, 91])
% imagesc(dFP(canoni, :)/input.nCells)
% colorbar
% colormap(figureDetails.colorMap)
% caxis([0 1])
% title('False Positives', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10])
% xticklabels(algoLabels)
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% subplot(12, 12, [81, 82, 83, 93, 94, 95])
% imagesc(dFN(canoni, :)/input.nCells)
% colorbar
% colormap(figureDetails.colorMap)
% caxis([0 1])
% title('False Negatives', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10])
% xticklabels(algoLabels)
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% %Phys
% subplot(12, 12, [125, 126, 127, 137, 138, 139])
% imagesc(dFP(physi, :)/input.nCells)
% colorbar
% colormap(figureDetails.colorMap)
% caxis([0 1])
% title('False Positives', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10])
% xticklabels(algoLabels)
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% subplot(12, 12, [129, 130, 131, 141, 142, 143])
% imagesc(dFN(physi, :)/input.nCells)
% colorbar
% colormap(figureDetails.colorMap)
% caxis([0 1])
% title('False Negatives', ...
%     'FontSize', figureDetails.fontSize-1, ...
%     'FontWeight', 'bold')
% clear xticks
% xticks([1 2 3 4 5 6 7 8 9 10])
% xticklabels(algoLabels)
% xtickangle(90)
% set(gca, 'FontSize', figureDetails.fontSize-1)
% 
% print(sprintf('%s/figs/sup8Full-Summary-%i-%i-%i-%i-%i_%i', ...
%     HOME_DIR2, ...
%     input.gDate, ...
%     input.gRun, ...
%     input.nDatasets, ...
%     workingOnServer), ...
%     '-dpng')
% disp('... Sup Fig. 2 plotted!')

%% Finish
elapsedTime2 = toc;
fprintf('Elapsed Time: %.4f seconds\n', elapsedTime2)
disp('... All Done!')