%{
Here are the various fields within the MATLAB structure 'sdcp' -->
1. timeCellPercent
2. cellOrder
3. maxHitTrialPercent
4. hitTrialPercentAssignment
5. trialOrder
6. eventWidth
7. eventAmplificationFactor
8. eventTiming
9. startFrame
10. endFrame
11. imprecisionFWHM
12. imprecisionType
13. noise
14. noisePercent
15. randomseed
16. addBackgroundSpikes4ptc
17. addBackgroundSpikes4oc
18. backDistLambda
19. comment
%}

function synthBatchProfile = getSynthBatchProfile(sdcp, dStart, dEnd, impX)
nDatasets = dEnd-dStart+1;
nParams = 5; %Noise, EW, Imp., HTR, and Background
nParams2 = 11;

s.noisePercent = [];
s.eventWidth = [];
s.imprecisionFWHM = [];
s.maxHitTrialPercent = [];
s.backDistLambda = [];
dProfile1 = repmat(s, 1, nDatasets);
clear s

s.noisePercent = [];
s.eventWidth = [];
s.eventWidthStdDevX = [];
s.imprecisionFWHM = [];
s.imprecisionType = [];
s.maxHitTrialPercent = [];
s.hitTrialPercentAssignent = [];
s.trialOrder = [];
s.backDistLambda = [];
s.timeCellPercent = [];
s.cellOrder = [];
dProfile2 = repmat(s, 1, nDatasets);
clear s


%Preallocation
dProfileM1 = zeros(nDatasets, nParams);
dProfileM2 = zeros(nDatasets, nParams2);
dProfileM1Norm = zeros(nDatasets, nParams);
dProfileM2Norm = zeros(nDatasets, nParams2);

for dataseti = dStart:dEnd
    %1. Noise (%)
    dProfile1(dataseti).noisePercent = sdcp(dataseti).noisePercent;
    dProfile2(dataseti).noisePercent = sdcp(dataseti).noisePercent;

    dProfileM1(dataseti, 1) = sdcp(dataseti).noisePercent;
    
    dProfileM2(dataseti, 1) = sdcp(dataseti).noisePercent;

    %2. EW (%ile)
    dProfile1(dataseti).eventWidth = sdcp(dataseti).eventWidth{1};
    dProfile2(dataseti).eventWidth = sdcp(dataseti).eventWidth{1};
    dProfile2(dataseti).eventWidthStdDevX = sdcp(dataseti).eventWidth{2};
    
    dProfileM1(dataseti, 2) = sdcp(dataseti).eventWidth{1};

    dProfileM2(dataseti, 2) = sdcp(dataseti).eventWidth{1};
    dProfileM2(dataseti, 3) = sdcp(dataseti).eventWidth{2};

    %3. Imp. (frames)
    if impX == 0
        dProfile1(dataseti).imprecisionFWHM = sdcp(dataseti).imprecisionFWHM;
        dProfile2(dataseti).imprecisionFWHM = sdcp(dataseti).imprecisionFWHM;

        dProfileM1(dataseti, 3) = sdcp(dataseti).imprecisionFWHM;

        dProfileM2(dataseti, 4) = sdcp(dataseti).imprecisionFWHM;
    else
        dProfile1(dataseti).imprecisionFWHM = sdcp(dataseti).imprecisionFWHM * impX;
        dProfile2(dataseti).imprecisionFWHM = sdcp(dataseti).imprecisionFWHM * impX;

        dProfileM1(dataseti, 3) = sdcp(dataseti).imprecisionFWHM * impX;

        dProfileM2(dataseti, 4) = sdcp(dataseti).imprecisionFWHM * impX;
    end

    if strcmpi(sdcp(dataseti).imprecisionType, 'uniform')
        dProfile2(dataseti).imprecisionType = 1;
        
        dProfileM2(dataseti, 5) = 1;

    elseif strcmpi(sdcp(dataseti).imprecisionType, 'normal')
        dProfile2(dataseti).imprecisionType = 0.5;
        
        dProfileM2(dataseti, 5) = 0.5;

    elseif strcmpi(sdcp(dataseti).imprecisionType, 'none')
        dProfile2(dataseti).imprecisionType = 0;
        
        dProfileM2(dataseti, 5) = 0;

    end

    %4. HTR (%)
    dProfile1(dataseti).maxHitTrialPercent = sdcp(dataseti).maxHitTrialPercent;
    dProfile2(dataseti).maxHitTrialPercent = sdcp(dataseti).maxHitTrialPercent;
    

    dProfileM1(dataseti, 4) = sdcp(dataseti).maxHitTrialPercent;

    dProfileM2(dataseti, 6) = sdcp(dataseti).maxHitTrialPercent;

    if strcmpi(sdcp(dataseti).hitTrialPercentAssignment, 'fixed')
        dProfile2(dataseti).hitTrialPercentAssignment = 0;
        dProfileM2(dataseti, 7) = 0;
    elseif strcmpi(sdcp(dataseti).hitTrialPercentAssignment, 'random')
        dProfile2(dataseti).hitTrialPercentAssignment = 1;
        dProfileM2(dataseti, 7) = 1;
    end
    
    if strcmpi(sdcp(dataseti).trialOrder, 'random')
        dProfile2(dataseti).trialOrder = 1;
        
        dProfileM2(dataseti, 8) = 1;
    elseif strcmpi(sdcp(dataseti).trialOrder, 'basic')
        dProfile2(dataseti).trialOrder = 0;

        dProfileM2(dataseti, 8) = 0;
    end

    %5. Background (\lambda)
    if sdcp(dataseti).addBackgroundSpikes4ptc
        dProfile1(dataseti).backDistLambda = sdcp(dataseti).backDistLambda;
        dProfile2(dataseti).backDistLambda = sdcp(dataseti).backDistLambda;
        
        dProfileM1(dataseti, 5) = sdcp(dataseti).backDistLambda;

        dProfileM2(dataseti, 9) = sdcp(dataseti).backDistLambda;
    else
        dProfile1(dataseti).backDistLambda = 0;
        dProfile2(dataseti).backDistLambda = 0;
    
        dProfileM2(dataseti, 9) = 0;
    end

    % Time Cell Percentages
    dProfile2(dataseti).timeCellPercent = sdcp(dataseti).timeCellPercent;
    
    dProfileM2(dataseti, 10) = sdcp(dataseti).timeCellPercent;

    if strcmpi(sdcp(dataseti).cellOrder, 'basic')
        dProfile2(dataseti).cellOrder = 0;

        dProfileM2(dataseti, 11) = 0;
    elseif strcmpi(sdcp(dataseti).cellOrder, 'random')
        dProfile2(dataseti).cellOrder = 1;
        
        dProfileM2(dataseti, 11) = 0;
    end

    for parami = 1:nParams
        dProfileM1Norm(:, parami) = dProfileM1(:, parami)/max(squeeze(dProfileM1(:, parami)));
    end

    for parami = 1:nParams2
        dProfileM2Norm(:, parami) = dProfileM2(:, parami)/max(squeeze(dProfileM2(:, parami)));
    end
end
synthBatchProfile.dProfile1 = dProfile1;
synthBatchProfile.dProfile2 = dProfile2;
synthBatchProfile.dProfileM1 = dProfileM1;
synthBatchProfile.dProfileM2 = dProfileM2;
synthBatchProfile.dProfileM1Norm = dProfileM1Norm;
synthBatchProfile.dProfileM2Norm = dProfileM2Norm;

end