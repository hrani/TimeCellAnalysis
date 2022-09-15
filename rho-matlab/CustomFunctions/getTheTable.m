function [response, predictor] = getTheTable(sdo_batch, cData, input)

% Prepare Look Up Table (LUT)
disp('Creating Look Up Table ...')

nCells = input.nCells;
nMethods = input.nMethods;
nDatasets = input.nDatasets;

preLUT = zeros(nCells*nDatasets, nMethods+2);

count = 0;
for dnum = 1:nDatasets
    count = count + 1;
    start = ((count-1)*nCells + 1);
    finish = count*nCells;
    %fprintf('Start = %i, Finish = %i\n', start, finish)
    
    reality = zeros(nCells, 1); %Preallocate
    reality(sdo_batch(dnum).ptcList) = 1; %Ground Truth
    
    for column = 1:nMethods+2
        if column == 1
            preLUT(start:finish, column) = reality; %Ground Truth - True Class Labels
        elseif column == 2
            preLUT(start:finish, column) = squeeze(sdo_batch(dnum).Q);
        elseif column == 3
            preLUT(start:finish, column) = squeeze(cData.methodA.mAOutput_batch(dnum).Q); %scores
        elseif column == 4
            preLUT(start:finish, column) = squeeze(cData.methodB.mBOutput_batch(dnum).Q); %scores
        elseif column == 5
            preLUT(start:finish, column) = squeeze(cData.methodC.mCOutput_batch(dnum).Q2); %scores
        elseif column == 6
            preLUT(start:finish, column) = squeeze(cData.methodD.mDOutput_batch(dnum).Q); %scores
        elseif column == 7
            preLUT(start:finish, column) = squeeze(cData.methodE.mEOutput_batch(dnum).Q); %scores
        elseif column == 8
            preLUT(start:finish, column) = squeeze(cData.methodF.mFOutput_batch(dnum).Q2); %scores
        end

        if algo == 1
            preLUT(start:finish, algo) = reality;
        elseif algo == 2
            preLUT(start:finish, algo) = squeeze(cData.methodA.mAOutput_batch(dnum).timeCells1);
        elseif algo == 3
            preLUT(start:finish, algo) = squeeze(cData.methodA.mAOutput_batch(dnum).timeCells2);
        elseif algo == 4
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells1);
        elseif algo == 5
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells2);
        elseif algo == 6
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells3);
        elseif algo == 7
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells4);
        elseif algo == 8
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells5);
        elseif algo == 9
            preLUT(start:finish, algo) = squeeze(cData.methodB.mBOutput_batch(dnum).timeCells6);
        elseif algo == 10
            preLUT(start:finish, algo) = squeeze(cData.methodC.mCOutput_batch(dnum).timeCells1);
        elseif algo == 11
            preLUT(start:finish, algo) = squeeze(cData.methodC.mCOutput_batch(dnum).timeCells2);
        elseif algo == 12
            preLUT(start:finish, algo) = squeeze(cData.methodD.mDOutput_batch(dnum).timeCells1);
        elseif algo == 13
            preLUT(start:finish, algo) = squeeze(cData.methodE.mEOutput_batch(dnum).timeCells1);
        elseif algo == 14
            preLUT(start:finish, algo) = squeeze(cData.methodF.mFOutput_batch(dnum).timeCells1);
        elseif algo == 15
            preLUT(start:finish, algo) = squeeze(cData.methodF.mFOutput_batch(dnum).timeCells2);
        end
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
disp('... done!')

response = LUT(:, 1); % Ground Truth - True Class Labels
predictor = LUT(:, 3:input.nMethods+2);
end
% toc
