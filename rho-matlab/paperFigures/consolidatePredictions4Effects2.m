
% close all
% clear

% tic
function [Y, X] = consolidatePredictions4Effects2(input, ptcList, pred1, pred2, pred3, pred4, pred5, indices)

nCells = input.nCells;
nAlgos = input.nAlgos;
nDatasets = length(indices);
%figureDetails = compileFigureDetails(16, 2, 5, 0.5, 'inferno'); %(fontSize, lineWidth, markerSize, transparency, colorMap)
%Extra colormap options: inferno/plasma/viridis/magm
%C = distinguishable_colors(nAlgos);

% Prepare Look Up Table (lut)
disp('Creating Confusion Matrix ...')
preLUT = zeros(nCells*nDatasets, nAlgos+1);

%reality = zeros(nDatasets, nCells);

count = 0;
for dnum = indices(1):indices(nDatasets)
    count = count + 1;
    start = ((count-1)*nCells + 1);
    finish = count*nCells;
    %fprintf('Start = %i, Finish = %i\n', start, finish)
    
    reality = zeros(nCells, 1); %Preallocate
    reality(ptcList(dnum).ptcList) = 1; %Ground Truth
    
    for algo = 1:nAlgos+1
        if algo == 1
            preLUT(start:finish, algo) = reality;
        elseif algo == 2
            preLUT(start:finish, algo) = pred1;
        elseif algo == 3
            preLUT(start:finish, algo) = squeeze(cData.methodA.mAOutput_batch(dnum).timeCells2);
        elseif algo == 4
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells1);
        elseif algo == 5
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells2);
    end
end
if input.removeNaNs
    %Edit out NaNs as O
    reshapedLUT = reshape(preLUT, [], 1);
    reshapedLUT(isnan(reshapedLUT)) = 0;
    LUT = reshape(reshapedLUT, size(preLUT));
else
    LUT = preLUT;
end

Y = LUT(:, 1); %Truth
X = LUT(:, 2:input.nAlgos+1);

disp('... done!')

end
% toc
