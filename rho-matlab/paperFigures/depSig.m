%{
BASELINE F1 SCORE: joinBundledResults3_base has 50 datasets (Ns) by 10 algorithms
MODULATION F1 SCORE: joinBundledResults3_mod has 50 datasets by 10 algorithms, each in 3 sets (low, med, high)
%}

clear mdl
doLinearFit = 2;
doANOVA = 0;
doTTest = 0;
g1 = 1:size(joinBundledResults3_mod, 2); %Effect of Algorithm
nTests = 5;
allC = zeros(input.nParams, nTests, 45, 6);
% fig2 = figure(2);
% clf
% set(fig2, 'Position', [100, 100, 1200, 900])
alpha = 0.01;
fprintf('[INFO] Level of Significance set to %.4f\n', alpha)
count = 0;

linFitLabels = {'Data', 'Fit', '95% PI'};
linFitLabels2 = {'Data', 'Fit'};
if doLinearFit == 1
    x = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3];
    slopes = zeros(input.nParams, input.nAlgos);
    intercepts = zeros(input.nParams, input.nAlgos);
    fig162 = figure(162);
    clf
    set(fig162, 'Position', [1, 1, 1200, 900])

    for algoi = 1: input.nAlgos
        for parami = 1:input.nParams
            count = count + 1;
            fprintf('[INFO] Param: %i ...\n', parami)
            datasets = (parami-1)*10 + (1:10);
            y = squeeze(joinBundledResults3_mod(datasets, :, :));

            %yR = reshape(squeeze(mean(squeeze(y(:, algoi, :)), 1)), [], 1); %reshaped
            yR = reshape(sort(squeeze(y(:, algoi, :))), [], 1);
            [p,S,mu] = polyfit(x, yR, 1);
            slopes(parami, algoi) = p(1);
            intercepts(parami, algoi) = p(2);
            [y_fit1, delta1] = polyval(p, x, S, mu);
            deltas(parami, algoi) = delta1;

            subplot(input.nAlgos, input.nParams, count)
            plot(x,yR,'bo')
            hold on
            plot(x,y_fit1,'k-')
            plot(x,y_fit1+2*delta1,'m--',x,y_fit1-2*delta1,'m--')
            title(sprintf('%s | %s', algoLabels{algoi}, paramLabels{parami}))
            lgd = legend('Data','Lin. Fit','95% PI', ...
                'Location', 'northwest');
            lgd.FontSize = figureDetails.fontSize-3;
        end
    end
    set(gca, 'FontSize', figureDetails.fontSize)
    print(sprintf('%s/figs/SupFig162-LinearFits-%i-%i-%i-%i-%i_%i', ...
        HOME_DIR2, ...
        input.gDate, ...
        input.gRun, ...
        input.nDatasets, ...
        workingOnServer), ...
        '-dpng')

    fig163 = figure(2);
    clf
    set(fig163, 'Position', [100, 100, 1200, 600])
    subplot(1, 2, 1)
    imagesc(flip(slopes))
    title('Slopes', 'FontSize', figureDetails.fontSize, 'FontWeight', 'bold')
    xticks(1:input.nAlgos)
    xticklabels(algoLabels)
    xtickangle(45)
    yticks(1:5)
    yticklabels(paramLabels)
    colormap(magma)
    z = colorbar;
    ylabel(z, 'F1/Param.', 'FontSize', figureDetails.fontSize, 'FontWeight', 'bold')
    set(gca, 'FontSize', figureDetails.fontSize)

    subplot(1, 2, 2)
    imagesc(flip(intercepts))
    title('Intercepts', 'FontSize', figureDetails.fontSize, 'FontWeight', 'bold')
    xticks(1:input.nAlgos)
    xticklabels(algoLabels)
    xtickangle(45)
    yticks(1:5)
    %yticklabels(paramLabels)
    yticklabels({'', '', '', '', ''})
    colormap(magma)
    z = colorbar;
    ylabel(z, '\Delta F1/Param.', 'FontSize', figureDetails.fontSize, 'FontWeight', 'bold')
    set(gca, 'FontSize', figureDetails.fontSize)

    print(sprintf('%s/figs/SupFig163-Resource-%i-%i-%i-%i-%i_%i', ...
        HOME_DIR2, ...
        input.gDate, ...
        input.gRun, ...
        input.nDatasets, ...
        workingOnServer), ...
        '-dpng')

elseif doLinearFit == 2
    x = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3];
    slopes = zeros(input.nParams, input.nAlgos);
    intercepts = zeros(input.nParams, input.nAlgos);
    deltas = zeros(input.nParams, input.nAlgos);
    allP = zeros(input.nParams, input.nAlgos);
    allAdjRSq = zeros(input.nParams, input.nAlgos);
    allRMSE = zeros(input.nParams, input.nAlgos);
    allMSE = zeros(input.nParams, input.nAlgos);

    fig162 = figure(162);
    clf
    set(fig162, 'Position', [1, 1, 900, 1200])
    for algoi = 1: input.nAlgos
        fprintf('[INFO] Algo: %s ...\n', algoLabels{algoi})
        for parami = 1:input.nParams
            count = count + 1;
            %fprintf('[INFO] Param: %s ...\n', paramLabels{parami})
            datasets = (parami-1)*10 + (1:10);

            %Param. Mod. Values
            if parami == 1
                myXTickLabels = {num2str(noiseVal1), num2str(noiseVal2), num2str(noiseVal3)};
            elseif parami == 2
                myXTickLabels = {num2str(ewVal1), num2str(ewVal2), num2str(ewVal3)};
            elseif parami == 3
                myXTickLabels = {num2str(impVal1), num2str(impVal2), num2str(impVal3)};
            elseif parami == 4
                myXTickLabels = {num2str(htrVal1), num2str(htrVal2), num2str(htrVal3)};
            elseif parami == 5
                myXTickLabels = {num2str(backVal1), num2str(backVal2), num2str(backVal3)};
            end
            y = squeeze(joinBundledResults3_mod(datasets, :, :)); %Selects datasets based on the paramter being modulated
            %myText = sprintf('p = %.4f, Adj. r^2 = %.4f', p, mdl.Rsquared.Adjusted);
            %disp(myText)
            %yR = reshape(squeeze(mean(squeeze(y(:, algoi, :)), 1)), [], 1); %reshaped
            yR = reshape(squeeze(y(:, algoi, :)), [], 1); %reshaped
            %yR = squeeze(y(:, algoi, :));
            yR0 = zeros(size(yR));

            mdl = fitlm(x, yR, 'linear');

            %[p, F, r] = coefTest(mdl, [1, 0]);
            [p, F, r] = coefTest(mdl);

            allP(parami, algoi) = p;
            allAdjRSq(parami, algoi) = mdl.Rsquared.Adjusted;

            slopes(parami, algoi) = mdl.Coefficients.Estimate(2);
            intercepts(parami, algoi) = mdl.Coefficients.Estimate(1);

            allRMSE(parami, algoi) = mdl.RMSE;
            allMSE(parami, algoi) = mdl.MSE;

            subplot(input.nAlgos, input.nParams, count)
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
            %ylim([-0.1, 1.1])
            %             title(sprintf('%s', algoLabels{algoi}), ...
            %                 'FontSize', figureDetails.fontSize-1, ...
            %                 'FontWeight', 'bold')
            title('')
            if parami == 1
                ylabel({sprintf('[%s]', algoLabels{algoi}); 'Adj. F1'}, ...
                    'FontSize', figureDetails.fontSize-1, ...
                    'FontWeight', 'bold')
            else
                ylabel('')
            end
            xticklabels(myXTickLabels)
            if algoi == 10
                xlabel(sprintf('%s', paramLabels{parami}), ...
                    'FontSize', figureDetails.fontSize-3, ...
                    'FontWeight', 'bold')
            else
                xlabel('')
            end
            %             if algoi == 10 && parami == 1
            %                 ylabel('Adj. F1 Score', ...
            %                     'FontSize', figureDetails.fontSize-3, ...
            %                     'FontWeight', 'bold')
            %             else
            %                 ylabel('')
            %             end

            fprintf('Algo: %s, Param: %s, y = %.4fx + %.4f, r^{2}=%.4f, Adj. r^{2}=%.4f, p=%.4f, h=%i, alpha=%.4f\n', ...
                algoLabels{algoi}, ...
                paramLabels{parami}, ...
                mdl.Coefficients.Estimate(2), ...
                mdl.Coefficients.Estimate(1), ...
                mdl.Rsquared.Ordinary, ...
                mdl.Rsquared.Adjusted, ...
                p, ...
                p<alpha, ...
                alpha);

            if parami == 5 && algoi == 10
                lgd = legend(linFitLabels, ...
                    'Location', 'northeast');
                lgd.FontSize = figureDetails.fontSize-6;
            else
                legend('hide')
            end
            set(gca, 'FontSize', figureDetails.fontSize-3)
        end
    end
    print(sprintf('%s/figs/SupFig162-synthPhysDependenceLinearFits-%i-%i-%i-%i-%i_%i', ...
        HOME_DIR2, ...
        input.gDate, ...
        input.gRun, ...
        input.nDatasets, ...
        workingOnServer), ...
        '-dpng')

    sigP = zeros(size(allP));
    sigP(allP >= alpha) = 0;
    sigP(allP < alpha) = 1;
end

% x1 = joinBundledResults3_base(datasets, :, :); %"Paired"
% x2 = joinBundledResults3_base; %"Unpaired"
%
% myX1 = x1(:, :) - x1(:, :);
%
% N = size(joinBundledResults3_base, 1);
% K = length(datasets);
% myX2 = x2(randsample(N, K), :) - x2(randsample(N, K), :);
%
% myX3 = x2(randsample(N, N), :) - x2(randsample(N, N), :);
%
% if parami == 1
%     myY = squeeze(y(:, :, 3)) - squeeze(y(:, :, 1))/(noiseVal3-noiseVal1);
% elseif parami == 2
%     myY = squeeze(y(:, :, 3)) - squeeze(y(:, :, 1))/(ewVal3-ewVal1);
% elseif parami == 3
%     myY = squeeze(y(:, :, 3)) - squeeze(y(:, :, 1))/(impVal3-impVal1);
% elseif parami == 4
%     myY = squeeze(y(:, :, 3)) - squeeze(y(:, :, 1))/(htrVal3-htrVal1);
% elseif parami == 5
%     myY = squeeze(y(:, :, 3)) - squeeze(y(:, :, 1))/(backVal3-backVal1);
% end
%
% %myYR = reshape(myY, [], 1);
% %myXR = reshape(myX, [], 1);
% if doANOVA
%     disp('[INFO] One-Way ANOVA ...')
%     [p, ~, stats] = anova1(myY);
%
%     for typei = 1:nTests
%         %typei = 2
%         %subplot(input.nParams, nTests, parami*typei)
%         if typei == 1
%             ctype = 'tukey-kramer'
%         elseif typei == 2
%             ctype = 'bonferroni'
%         elseif typei == 3
%             ctype = 'dunn-sidak'
%         elseif typei == 4
%             ctype = 'lsd'
%         elseif typei == 5
%             ctype = 'scheffe'
%         end
%         [c, ~, ~, ~] = multcompare(stats, ...
%             'CType', ctype, ...
%             'Alpha', alpha)
%         title(sprintf('Param: %s', paramLabels{parami}))
%         yticklabels(algoLabels)
%         set(gca, 'ydir', 'reverse')
%         allC(parami, typei, :, :) = c;
%     end
%     disp('[INFO] ... done.')
% end
%
% if doTTest
%     disp('[INFO] Two-Sample T-Test ...')
%     % Two-sample T-Test
%     myX0 = myX3;
%     varY = var(myY);
%     varX = var(myX0);
%
%     % Rule of Thumb with variance
%     if max(varY, varX)/min(varY, varX) < 4
%         disp('[INFO] ... w/ equal variances ...')
%         [h, p, ci, stats] = ttest2(myY, myX0, ...
%             'Alpha', alpha, ...
%             'Vartype', 'equal')
%     else
%         disp('[INFO] ... w/ unequal variances ...')
%         [h, p, ci, stats] = ttest2(myY, myX0, ...
%             'Alpha', alpha, ...
%             'Vartype', 'unequal')
%     end
%     disp('[INFO] ... done.')
% end
disp('... done.')
% end