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
PM(1).labels = {'Significant wave height (m)';
                                'Zero-upcrossing period (s)'};
PM.gridCenterPoints = {0:0.1:24; 0:0.1:18};

% Define the exceedance probabilities.
alphas = [0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001]';

figContours = figure();
subplot(1, 2, 1);
hold on

C1 = cell(size(alphas));
C2 = cell(size(alphas));
maxHs = nan(size(alphas));

for i = 1 : length(alphas)
    % Calculate the highest density contour.
    [fm, x1Hdc, x2Hdc] = computeHdc(PM, alphas(i), PM.gridCenterPoints, 0);
    
    C1{i} = x1Hdc{1};
    C2{i} = x2Hdc{1};
    maxHs(i) = max(x1Hdc{1});
    plot(C2{i}, C1{i});
end
ylim([0  22]);
xlim([0 18]);
xlabel(PM(1).labels{2})
ylabel(PM(1).labels{1})
legend(num2str(alphas, '%10.0e'), 'location', 'northwest');
legend box off

subplot(1, 2, 2);
hold on
plot(alphas, maxHs, '--k');
plot(alphas, maxHs, 'ok');
set(gca, 'XDir', 'reverse');
set(gca, 'xtick', flip(alphas));
set(gca, 'xscale', 'log')
ylim([0  22]);
xlabel('\alpha (-)');
ylabel('Maximum value h_s-value along the contour (m)');
