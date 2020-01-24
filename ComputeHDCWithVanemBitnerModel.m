

% Define a joint probability distribution (taken from Vanem and Bitner-
% Gregersen, 2012):
PM.name = 'Vanem and Bitner-Gregersen (2012), doi: 10.1016/j.apor.2012.05.006';
PM.modelType = 'CMA';
PM.distributions = {'weibull'; 'lognormal'};
PM.isConditionals = {[0 0 0]; [1 1]};
PM.coeffs = {{2.776 1.471 0.8888}; 
                             { @(x1)0.1000 + 1.489 * x1^0.1901;
                               @(x1)0.0400 + 0.1748 * exp(-0.2243*x1)}
                            };
PM(1).labels = {'Significant wave height, h_s (m)';
                                'Zero-up-crossing period, t_z (s)'};
PM.gridCenterPoints = {0:0.05:24; 0:0.05:18};

% Define the exceedance probabilities.
alphas = [0.5 0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001]';

figShapeIs1p471 = figure();
subplot(3, 2, 1);
hold on

C1 = cell(size(alphas));
C2 = cell(size(alphas));
maxHs = nan(size(alphas));
hsAlpha = nan(size(alphas));
pMarginal = nan(size(alphas));

for i = 1 : length(alphas)
    % Calculate the highest density contour.
    [fm, x1Hdc, x2Hdc] = computeHdc(PM, alphas(i), PM.gridCenterPoints, 0);
    
    C1{i} = x1Hdc{1};
    C2{i} = x2Hdc{1};
    maxHs(i) = max(x1Hdc{1});
    pMarginal(i) = 1 - cdf('weibull', maxHs(i) - 0.8888, 2.776, 1.471);
    hsAlpha(i) = icdf('weibull', 1 - alphas(i), 2.776, 1.471) + 0.8888;
    plot(C2{i}, C1{i});
end
ylim([0  22]);
xlim([0 18]);
xlabel(PM(1).labels{2})
ylabel(PM(1).labels{1})
legendCell = cellstr(num2str(alphas, '%1.0e'));
lgd = legend(legendCell, 'location', 'northwest');
lgd.Title.String = '\alpha-values';
legend box off

subplot(3, 2, 2);
hold on
plot(alphas, maxHs, '--k');
plot(alphas, maxHs, 'ok');
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas(2:end)));
set(gca, 'xscale', 'log')
ylim([0  22]);
xlabel('\alpha (-)');
ylabel('Max. h_s-value along the contour, h_{s,HDC} (m)');

subplot(3, 2, 3);
hold on
plot(alphas, pMarginal, '--k');
plot(alphas, pMarginal, 'ok');
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas(2:end)));
set(gca, 'xscale', 'log')
set(gca, 'YDir', 'reverse');
set(gca, 'yscale', 'log');
xlabel('\alpha (-)');
ylabel('Marginal exceedance prob. for h_{s,HDC}, p_M (-)');

subplot(3, 2, 4)
hold on
plot(alphas, hsAlpha, '--k');
plot(alphas, hsAlpha, 'ok');
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas(2:end)));
set(gca, 'xscale', 'log')
ylim([0  22]);
xlabel('\alpha (-)');
ylabel('h_{s,M} (-)');

subplot(3, 2, 5)
hold on
plot(alphas, maxHs ./ hsAlpha, '--k');
plot(alphas, maxHs ./ hsAlpha, 'ok');
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas(2:end)));
set(gca, 'xscale', 'log')
xlabel('\alpha (-)');
ylabel('h_{s,HDC} / h_{s,M}(-)');

subplot(3, 2, 6)
hold on
plot(alphas, alphas ./ pMarginal, '--k');
plot(alphas, alphas ./ pMarginal, 'ok');
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas(2:end)));
set(gca, 'xscale', 'log')
xlabel('\alpha (-)');
ylabel('\alpha_{HDC} / \alpha_{M} (-)');

betas =  [1.1, 1.471, 2 3 4];
betaMarkers = {'-rs', '-bo', '-k^', '-m+', '-g*', '-cd'};
C1 = cell(length(betas), length(alphas));
C2 = cell(length(betas), length(alphas));
maxHs = nan(length(betas), length(alphas));
hsAlpha = nan(length(betas), length(alphas));
pMarginal = nan(length(betas), length(alphas));

for i = 1:length(betas)
    % Define a joint probability distribution, the model is based on the model 
    % from Vanem and Bitner-Gregersen, 2012, however, the shape parameter 
    % is changed.
    PM.name = 'altered Vanem and Bitner-Gregersen (2012), doi: 10.1016/j.apor.2012.05.006';
    PM.modelType = 'CMA';
    PM.distributions = {'weibull'; 'lognormal'};
    PM.isConditionals = {[0 0 0]; [1 1]};
    PM.coeffs = {{2.776 betas(i) 0.8888}; 
                                 { @(x1)0.1000 + 1.489 * x1^0.1901;
                                   @(x1)0.0400 + 0.1748 * exp(-0.2243*x1)}
                                };
    PM(1).labels = {'Significant wave height, h_s (m)';
                                    'Zero-up-crossing period, t_z (s)'};
    PM.gridCenterPoints = {0:0.05:40; 0:0.05:30};
    for j = 1 : length(alphas)
        % Calculate the highest density contour.
        [fm, x1Hdc, x2Hdc] = computeHdc(PM, alphas(j), PM.gridCenterPoints, 0);
        C1{i, j} = x1Hdc{1};
        C2{i, j} = x2Hdc{1};
        maxHs(i,j) = max(x1Hdc{1});
        pMarginal(i, j) = 1 - cdf('weibull', maxHs(i, j) - 0.8888, 2.776, betas(i));
        hsAlpha(i, j) = icdf('weibull', 1 - alphas(j), 2.776, betas(i)) + 0.8888;
    end
end

figDifferentShapeValues = figure();
subplot(1, 3, 1)
hold on
j = 5; %10^-6
for i = 1:length(betas)
    betaMarker = betaMarkers{i};
    plot(C2{i, j}, C1{i, j}, betaMarker(1:2));
end
xlabel(PM(1).labels{2})
ylabel(PM(1).labels{1})
legendCell = cellstr(num2str(betas', 'beta = %2.3f'));
legend(legendCell, 'location', 'northwest', 'orientation', 'vertical')
legend box off


subplot(1, 3, 2)
hold on
for i = 1:length(betas)
    plot(alphas', maxHs(i, :) ./ hsAlpha(i, :), betaMarkers{i});
end
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas(2:end)));
set(gca, 'xscale', 'log')
xlabel('\alpha (-)');
ylabel('h_{s,HDC} / h_{s,M} (-)');
legendCell = cellstr(num2str(betas', 'beta = %2.3f'));
legend(legendCell, 'location', 'northeast', 'orientation', 'vertical')
legend box off

subplot(1, 3, 3)
hold on
for i = 1:length(betas)
    plot(alphas', alphas' ./ pMarginal(i, :), betaMarkers{i});
end
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas(2:end)));
set(gca, 'xscale', 'log')
ylim([1 15])
xlabel('\alpha (-)');
ylabel('\alpha_{HDC} / \alpha_{M} (-)');
legendCell = cellstr(num2str(betas', 'beta = %2.3f'));
lgd = legend(legendCell, 'location', 'northwest', 'orientation', 'vertical');
lgd.NumColumns = 2;
legend box off
