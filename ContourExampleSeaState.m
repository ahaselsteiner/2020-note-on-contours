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

% Compute the contours.
[hsIFORM, tpIFORM] = computeIformContour(PM, alpha, 360); 
[h, t] = WblLogN_simulate(1 / alpha * 200);
[hsAE, tpAE] = direct_sampling_contour(h, t, alpha, 5);
[hsISORM, tpISORM] = computeIsormContour(PM, alpha, 360); 
[fm, hsHD, tpHD] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);



figFourContours = figure('position', [100 100 350 300]);
hold on
h = zeros(8, 1);
% IFORM
h(1) = plot([tpIFORM tpIFORM(1)], [hsIFORM hsIFORM(1)], '-b');
responseIFORM = ross2020Response(hsIFORM, tpIFORM);
[maxResponseIFORM, iMaxIFORM] = max(responseIFORM);
h(2) = plot(tpIFORM(iMaxIFORM), hsIFORM(iMaxIFORM), 'xb', 'linewidth', 2);

% AE
h(3) = plot([tpAE; tpAE(1)], [hsAE; hsAE(1)], '-k');
responseAE = ross2020Response(hsAE, tpAE);
[maxResponseAE, iMaxAE] = max(responseAE);
h(4) = plot(tpAE(iMaxAE), hsAE(iMaxAE), 'xk', 'linewidth', 2);

% ISORM
h(5) = plot([tpISORM tpISORM(1)], [hsISORM hsISORM(1)], '--b');
responseISORM = ross2020Response(hsISORM, tpISORM);
[maxResponseISORM, iMaxISORM] = max(responseISORM);
h(6) = plot(tpISORM(iMaxISORM), hsISORM(iMaxISORM), 'xb', 'linewidth', 2);

% HD
h(7) = plot(tpHD{1}, hsHD{1}, '--k');
responseHD = ross2020Response(hsHD{1}, tpHD{1});
[maxResponseHD, iMaxHD] = max(responseHD);
h(8) = plot(tpHD{1}(iMaxHD), hsHD{1}(iMaxHD), 'xk', 'linewidth', 2);

xlabel(PM.labels{2})
ylabel(PM.labels{1})
legend(h([1,3,5,7,8]), {'50-year IFORM contour', '50-year AE contour', ...
    '50-year ISORM contour', '50-year HD contour', 'Maximum response'}, ...
    'location', 'northwest');
legend box off
box off

contourName = {'IFORM'; 'Angular exceedance'; 'ISORM'; 'Highest density'};
maxResponse = [maxResponseIFORM; maxResponseAE; maxResponseISORM; maxResponseHD];
maxResponseHs = [hsIFORM(iMaxIFORM); hsAE(iMaxAE); hsISORM(iMaxISORM); hsHD{1}(iMaxHD)];
maxResponseTp = [tpIFORM(iMaxIFORM); tpAE(iMaxAE); tpISORM(iMaxISORM); tpHD{1}(iMaxHD)];
T = table(contourName, maxResponse, maxResponseHs, maxResponseTp)


figContourWithFailureSurface = figure('position', [100 100 350 300]);
hold on
%plot([tpIFORM tpIFORM(1)], [hsIFORM hsIFORM(1)], '-b');
h(1) = plot([tpAE; tpAE(1)], [hsAE; hsAE(1)], '-k');
%plot([tpISORM tpISORM(1)], [hsISORM hsISORM(1)], '--b');
h(2) = plot(tpHD{1}, hsHD{1}, '--k');
xlabel(PM.labels{2})
ylabel(PM.labels{1})
tp = [0 : 0.5 : 20];
hs = computeHsToReachResponse(tp, 25);
h(3) = plot(tp, hs, '-r');
h(4) = plot(tpAE(iMaxAE), hsAE(iMaxAE), 'xk', 'linewidth', 2);
h(5) = plot(tpHD{1}(iMaxHD), hsHD{1}(iMaxHD), 'xk', 'linewidth', 2);

xlim([0 18]);
ylim([0 20]);
lCell = {'50-year AE contour', '50-year HD contour', ...
    'Failure surface, r_{SDOF}=25', 'Maximum response'};
legend(h(1:4), lCell, 'location', 'northoutside', 'numcolumns', 2);
legend box on
box off



figResponseSurfaces = figure('position', [100 100 350 300]);
hold on
xlabel(PM.labels{2})
ylabel(PM.labels{1})
responseValues = [5 : 5 : 35];
legendCell = cell(length(responseValues), 1)
for i = 1 : length(responseValues)
    tp = [0 : 0.5 : 20];
    hs = computeHsToReachResponse(tp, responseValues(i));
    plot(tp, hs);
    legendCell{i} = num2str(responseValues(i),'r_{SDOF}=%-d');
end
legend(legendCell, 'orientation', 'vertical', 'location', 'northoutside', 'NumColumns', 4);
xlim([0 18]);
ylim([0 20]);
box off
legend box on

function r = ross2020Response(hs, tp)
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