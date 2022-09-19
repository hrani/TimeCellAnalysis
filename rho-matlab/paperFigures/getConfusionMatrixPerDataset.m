function [dTP, dFP, dFN, dTN] = getConfusionMatrixPerDataset(prediction, response)
nCases = size(prediction, 1); %length(response)
nAlgos = size(prediction, 2);

nDatasets = nCases/135; %nCells/dataset currently set to 135
dTP = zeros(nDatasets, nAlgos);
dFP = zeros(nDatasets, nAlgos);
dFN = zeros(nDatasets, nAlgos);
dTN = zeros(nDatasets, nAlgos);

for algo = 1: nAlgos
    dataset = 1;
    tp = 0;
    fp = 0;
    fn = 0;
    tn = 0;

    for myCase = 1:nCases
        if response(myCase) %True Cases
            if prediction(myCase, algo)
                tp = tp+1;
            else
                fn = fn+1;
            end
        else %False Cases
            if prediction(myCase, algo)
                fp = fp+1;
            else
                tn = tn+1;
            end
        end

        dTP(dataset, algo) = tp;
        dFP(dataset, algo) = fp;
        dFN(dataset, algo) = fn;
        dTN(dataset, algo) = tn;

        if ~mod(myCase, 135) %nCells/dataset set to 135
            dataset = dataset+1;
            tp = 0;
            fp = 0;
            fn = 0;
            tn = 0;
        end
    end
end
end