function [results1, results2] = getResults(datasetStart, datasetEnd, nMethods, groundTruth)

results1 = zeros(nMethods, 4); %TN, FN, FP, TP
results2 = zeros(nMethods, 3); %recall, precision, f1Score
for methodi = 1:nMethods
    allTN = 0;
    allFN = 0;
    allFP = 0;
    allTP = 0;
    for dataseti = datasetStart:datasetEnd
        fprintf('[INFO] Method: %i, Dataset: %i\n', methodi, dataseti);
        if methodi == 1
            %sigMean_ti - Vars2-5
            TN = squeeze(groundTruth(dataseti, 2));
            FN = squeeze(groundTruth(dataseti, 3));
            FP = squeeze(groundTruth(dataseti, 4));
            TP = squeeze(groundTruth(dataseti, 5));

        elseif methodi == 2
            %sigBootstrap_ti - Vars6-9
            TN = squeeze(groundTruth(dataseti, 6));
            FN = squeeze(groundTruth(dataseti, 7));
            FP = squeeze(groundTruth(dataseti, 8));
            TP = squeeze(groundTruth(dataseti, 9));

        elseif methodi == 3
            %theAND_ti - Vars10-13
            TN = squeeze(groundTruth(dataseti, 10));
            FN = squeeze(groundTruth(dataseti, 11));
            FP = squeeze(groundTruth(dataseti, 12));
            TP = squeeze(groundTruth(dataseti, 13));

        elseif methodi == 4
            %sigMean_r2b - Vars14-18
            TN = squeeze(groundTruth(dataseti, 14));
            FN = squeeze(groundTruth(dataseti, 15));
            FP = squeeze(groundTruth(dataseti, 16));
            TP = squeeze(groundTruth(dataseti, 17));

        elseif methodi == 5
            %sigBootstrap_r2b - Vars19-21
            TN = squeeze(groundTruth(dataseti, 18));
            FN = squeeze(groundTruth(dataseti, 19));
            FP = squeeze(groundTruth(dataseti, 20));
            TP = squeeze(groundTruth(dataseti, 21));
        
        elseif methodi == 6
            %sigPeq - Vars22-25
            TN = squeeze(groundTruth(dataseti, 22));
            FN = squeeze(groundTruth(dataseti, 23));
            FP = squeeze(groundTruth(dataseti, 24));
            TP = squeeze(groundTruth(dataseti, 25));
        end
        
        fprintf('TN:%i, FN:%i, FP:%i, TP:%i\n', TN, FN, FP, TP)
        allTN = allTN + TN;
        allFN = allFN + FN;
        allFP = allFP + FP;
        allTP = allTP + TP;
    end
    fprintf('allTN:%i, allFN:%i, allFP:%i, allTP:%i\n', allTN, allFN, allFP, allTP)

    TNR = allTN/(allTN + allFP);
    FNR = allFN/(allFN + allTP);
    FPR = allFP/(allFP + allTN);
    TPR = allTP/(allTP + allFN); %Recall
    PPV = allTP/(allTP + allFP); %Precision
    f1Score = 2*TPR*PPV/(TPR+PPV);

    fprintf('TNR:%0.4f, FNR:%0.4f, FPR:%0.4f, TPR:%0.4f\n', TNR, FNR, FPR, TPR)

    %Adding through datasets in the batch
    results1(methodi, 1) = TNR;
    results1(methodi, 2) = FNR;
    results1(methodi, 3) = FPR;
    results1(methodi, 4) = TPR;

    results2(methodi, 1) = TPR; %Recall
    results2(methodi, 2) = PPV; %Precision
    results2(methodi, 3) = f1Score;
end
end