%Paper Figures 2 - Kambadur Ananthamurthy
%Loads analysis output from .csv
%Plots the method scores and predictive performance
%The plots here are prepare the ground truth per dataset instead of per
%cell.
%Here we look at regime splits -
% 1. "Unphysiological",
% 2. "Canonical",
% 3. "Physiological"

clear
close all

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

%Additional search paths
%addpath(genpath(strcat(HOME_DIR, 'rho-matlab/CustomFunctions')))
%addpath(genpath(strcat(HOME_DIR, 'paperFigures'))) %
%addpath(genpath(strcat(HOME_DIR, 'TimeCellAnalysis/TcPy'))) %python/C++ analysis
addpath(genpath(strcat(HOME_DIR, 'AnalysisOutput'))) % location for .csv (python/C++ analysis outputs)
addpath(genpath(strcat(HOME_DIR, 'rho-matlab/localCopies')))

gDate = '20220618'; %string
gRun = 1; %int
nDatasets = 567; %int

%Load .csv
r2b = load(sprintf('r2b_%s_gRun%i_batch_%i.csv', gDate, gRun, nDatasets));
ti = load(sprintf('ti_%s_gRun%i_batch_%i.csv', gDate, gRun, nDatasets));
peq = load(sprintf('peq_%s_gRun%i_batch_%i.csv', gDate, gRun, nDatasets));

%Predictions
sigMean_ti = ti(:, 6);
sigBootstrap_ti = ti(:, 7);
%theAND_ti = sigMean_ti & sigBootstrap_ti;
sigMean_r2b = r2b(:, 6);
sigBootstrap_r2b = r2b(:, 7);
sigPeq = peq(:, 4);

%Ground Truth - In this case TN, FN, FP, TP - per Dataset
groundTruth = load(sprintf('groundTruth_%s_gRun%i_batch_%i.csv', gDate, gRun, nDatasets));

%%
nCells = 135;
nMethods = 6;
methodLabels = {'Mau-Pk', 'Mau-TI', 'Mau-Pk&&TI', 'R2B-thr', 'R2B-boot', 'PEQ'};
confLabels = {'TNR', 'FNR', 'FPR', 'TPR'};
perfLabels = {'Recall', 'Precision', 'F1 Score'};
%% Unphysiological - Datasets 1 to 297
disp('Unphysiological Regime ...')
datasetStart = 1;
datasetEnd = 297;
nDats = datasetEnd - datasetStart + 1;
totalCells = nDats*nCells;
[results1, results2] = getResults(datasetStart, datasetEnd, nMethods, groundTruth);

fig1 = figure(1);
clf
set(fig1, 'Position', [100 100 1200 1000])
subplot(2, 1, 1)
bar(results1)
title(sprintf('Unphysiological Regime - %s gRun%i nCells: %i', gDate, gRun, totalCells), 'FontSize', 16, 'FontWeight', 'bold')
xticklabels(methodLabels)
ylabel('Rate', 'FontSize', 16, 'FontWeight', 'bold')
legend(confLabels, 'Location','best')
set(gca, 'FontSize', 14)
%print(sprintf('confMat_Unphys_%s_gRun%i_batch_%i', gDate, gRun, (datasetEnd-datasetStart+1)), '-dpng')

subplot(2, 1, 2)
bar(results2)
title(sprintf('Unphysiological Regime - %s gRun%i nCells: %i', gDate, gRun, totalCells), 'FontSize', 16, 'FontWeight', 'bold')
xticklabels(methodLabels)
ylabel('Rate', 'FontSize', 16, 'FontWeight', 'bold')
legend(perfLabels, 'Location', 'best')
set(gca, 'FontSize', 14)
print(sprintf('unphys_%s_gRun%i_batch_%i', gDate, gRun, (datasetEnd-datasetStart+1)), '-dpng')

%% Canonical - Datasets 298 to 417
disp('Canonical Regime ...')
datasetStart = 298;
datasetEnd = 417;
nDats = datasetEnd - datasetStart + 1;
totalCells = nDats*nCells;
[results1, results2] = getResults(datasetStart, datasetEnd, nMethods, groundTruth);

fig2 = figure(2);
clf
set(fig2, 'Position', [100 100 1200 1000])
subplot(2, 1, 1)
bar(results1)
title(sprintf('Canonical Regime - %s gRun%i nCells: %i', gDate, gRun, totalCells), 'FontSize', 16, 'FontWeight', 'bold')
xticklabels(methodLabels)
ylabel('Rate', 'FontSize', 16, 'FontWeight', 'bold')
legend(confLabels, 'Location', 'best')
set(gca, 'FontSize', 14)
%print(sprintf('confMat_Canon_%s_gRun%i_batch_%i', gDate, gRun, (datasetEnd-datasetStart+1)), '-dpng')

subplot(2, 1, 2)
bar(results2)
title(sprintf('Canonical Regime - %s gRun%i nCells: %i', gDate, gRun, totalCells), 'FontSize', 16, 'FontWeight', 'bold')
xticklabels(methodLabels)
ylabel('Rate', 'FontSize', 16, 'FontWeight', 'bold')
legend(perfLabels, 'Location', 'best')
set(gca, 'FontSize', 14)
print(sprintf('canon_%s_gRun%i_batch_%i', gDate, gRun, (datasetEnd-datasetStart+1)), '-dpng')

%% Physiological Regime - Dataset 418 to 537
disp('Physiological Regime ...')
datasetStart = 418;
datasetEnd = 567;
nDats = datasetEnd - datasetStart + 1;
totalCells = nDats*nCells;
[results1, results2] = getResults(datasetStart, datasetEnd, nMethods, groundTruth);

fig3 = figure(3);
clf
set(fig3, 'Position', [100 100 1200 1000])
subplot(2, 1, 1)
bar(results1)
title(sprintf('Physiological Regime - %s gRun%i nCells: %i', gDate, gRun, totalCells), 'FontSize', 16, 'FontWeight', 'bold')
xticklabels(methodLabels)
ylabel('Rate', 'FontSize', 16, 'FontWeight', 'bold')
legend(confLabels, 'Location', 'best')
set(gca, 'FontSize', 14)
%print(sprintf('confMat_Phys_%s_gRun%i_batch_%i', gDate, gRun, (datasetEnd-datasetStart+1)), '-dpng')

subplot(2, 1, 2)
bar(results2)
title(sprintf('Physiological Regime - %s gRun%i nCells: %i', gDate, gRun, totalCells), 'FontSize', 16, 'FontWeight', 'bold')
xticklabels(methodLabels)
ylabel('Rate', 'FontSize', 16, 'FontWeight', 'bold')
legend(perfLabels, 'Location', 'best')
set(gca, 'FontSize', 14)
print(sprintf('phys_%s_gRun%i_batch_%i', gDate, gRun, (datasetEnd-datasetStart+1)), '-dpng')