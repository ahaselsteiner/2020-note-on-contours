% Example sea state contour.

% This is the Vanem and Bitner-Gregersen model adjusted for Tp instead of
% Tz. For a JONSWAP spectrum with mu = 3.3 tp = 1.2796 * tp; .
PM.name = 'Vanem and Bitner-Gregersen for Tp (2012)';
PM.modelType = 'CMA';
PM.distributions = {'weibull'; 'lognormal'};
PM.isConditionals = {[0 0 0]; [1 1]};
PM.coeffs = {{2.776 1.471 0.8888}; 
                             { @(x1)0.1000 + 1.489 * x1^0.1901;
                               @(x1)1.2796 * (0.0400 + 0.1748 * exp(-0.2243*x1))}
                            };
PM.labels = {'Significant wave height (m)';
                             'Spectral peak period (s)'};
PM.gridCenterPoints = {0:0.05:20; 0:0.05:18};


nYears = 50;
stateDuration = 6;
alpha = 1 / (nYears * 365.25 * 24 / stateDuration);
[hsIFORM, tpIFORM] = computeIformContour(PM, alpha, 360); 
[fm, hsHDC, tpHDC] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);

figContour = figure('position', [100 100 350 300]);
hold on
plot([tpIFORM tpIFORM(1)], [hsIFORM hsIFORM(1)], '-k');
plot(tpHDC{1}, hsHDC{1}, '--k');
xlabel(PM.labels{2})
ylabel(PM.labels{1})

tp = [0 : 0.5 : 20];
hs = computeHsToReachResponse(tp, 25);
plot(tp, hs, '-r');
xlim([0 18]);
ylim([0 20]);

legend({'50-year IFORM contour', '50-year HD contour', 'Failure surface'}, 'location', 'northwest');
legend box off
box off


function r = ross2019Response(hs, tp)
    alpha = 2;
    beta = 0.007;
    tp0 = 7;
 
  
    r = alpha * hs ./ (1 + beta .* (tp - tp0).^2);
end

function hs = computeHsToReachResponse(tp, r)
    alpha = 2;
    beta = 0.007;
    tp0 = 7;
    
    % r = alpha * hs ./ (1 + beta .* (tp - tp0).^2);
    hs = r * (1 + beta .* (tp - tp0).^2) ./ alpha;
end