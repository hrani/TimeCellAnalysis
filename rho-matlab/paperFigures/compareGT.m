function [TP, FN, FP, TN] = compareGT(predictions, gTruth, nCells)
TP = 0;
FN = 0;
FP = 0;
TN = 0;

for cell = 1:nCells
    if gTruth(cell) == 1 %GT - True Cases
        if predictions(cell) == 1 %
            TP = TP+1;
        else
            FN = FN+1;
        end
    else %GT - False Cases
        if predictions(cell) == 1
            FP = FP+1;
        else
            TN = TN+1;
        end
    end
end
fprintf('------> TN: %i FN: %i FP: %i TP: %i\n', TN, FN, FP, TP)
end