% Example sea state contour.

% Vanem and Bitner-Gregersen (DOI: 10.1016/j.apor.2012.05.006)
PM.name = 'Vanem and Bitner-Gregerse (2012), DOI: 10.1016/j.apor.2012.05.006';
PM.modelType = 'CMA';
PM.distributions = {'weibull'; 'lognormal'};
PM.isConditionals = {[0 0 0]; [1 1]};
PM.coeffs = {{2.776 1.471 0.8888}; 
                             { @(x1)0.1000 + 1.489 * x1^0.1901;
                               @(x1)0.0400 + 0.1748 * exp(-0.2243*x1)}
                            };
PM.labels = {'Significant wave height (m)';
                             'Zero-up-crossing period (s)'};
PM.gridCenterPoints = {0:0.05:20; 0:0.05:18};

% We will use tz instead of tp and assume that tp = 1.2796 * tz;
tztpCoeff = 1.2796;

nYears = 50;
stateDuration = 6;
alpha = 1 / (nYears * 365.25 * 24 / stateDuration);

% Compute the contours.
[hsIFORM, tzIFORM] = computeIformContour(PM, alpha, 360); 
tpIFORM = tztpCoeff * tzIFORM;
[h, t] = WblLogN_simulate(1 / alpha * 200);
[hsAE, tzAE] = direct_sampling_contour(h, t, alpha, 5);
tpAE = tztpCoeff * tzAE;
[hsISORM, tzISORM] = computeIsormContour(PM, alpha, 360); 
tpISORM = tztpCoeff * tzISORM;
[fm, hsHD, tzHD] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);
hsHD = hsHD{1};
tpHD = tztpCoeff * tzHD{1};

% Compute the true response using Monte Carlo simulation.
disp('Now starting MC simulations for all sea state approach.');
tic
[allSeaStateResponseMC, rStd] = allSeaStateApproachE1MonteCarlo(alpha, 1 / alpha * 100, 20);
tic
% Compute the true response using numerical integration
disp('Now starting numerical intergration for all sea state approach.');
tic
allSeaStateResponse = allSeaStateApproachE1NumericIntegration(alpha);
tic

% Compute the failure surface.
tpFailureSurface = [0 : 0.5 : 20];
hsFailureSurface = hsToReachRoss2020Response(tpFailureSurface, 19);

figFourContours = figure('position', [100 100 500 450]);
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
h(7) = plot(tpHD, hsHD, '--k');
responseHD = ross2020Response(hsHD, tpHD);
[maxResponseHD, iMaxHD] = max(responseHD);
h(8) = plot(tpFailureSurface, hsFailureSurface, '-r');
h(9) = plot(tpHD(iMaxHD), hsHD(iMaxHD), 'xk', 'linewidth', 2);

%xlabel(PM.labels{2})
xlabel('Spectral peak period (s)');
ylabel(PM.labels{1})
xlim([0 20]);
ylim([0 20])
legend(h([1,3,5,7,8,9]), {'50-year IFORM contour', '50-year AE contour', ...
    '50-year ISORM contour', '50-year HD contour', ...
    'Failure surface, r_{SDOF}=19', 'Maximum response'}, ...
    'location', 'northoutside', 'numcolumns', 3);
box off



contourName = {'IFORM'; 'Angular exceedance'; 'ISORM'; 'Highest density'; 'All sea state'};
maxResponse = [maxResponseIFORM; maxResponseAE; maxResponseISORM; maxResponseHD; allSeaStateResponse];
maxResponseHs = [hsIFORM(iMaxIFORM); hsAE(iMaxAE); hsISORM(iMaxISORM); hsHD(iMaxHD); NaN];
maxResponseTp = [tpIFORM(iMaxIFORM); tpAE(iMaxAE); tpISORM(iMaxISORM); tpHD(iMaxHD); NaN];
T = table(contourName, maxResponse, maxResponseHs, maxResponseTp)


figContourWithFailureSurface = figure('position', [100 100 500 450]);
hold on
%plot([tpIFORM tpIFORM(1)], [hsIFORM hsIFORM(1)], '-b');
h(1) = plot([tpAE; tpAE(1)], [hsAE; hsAE(1)], '-k');
%plot([tpISORM tpISORM(1)], [hsISORM hsISORM(1)], '--b');
h(2) = plot(tpHD, hsHD, '--k');
%xlabel(PM.labels{2})
xlabel('Spectral peak period (s)');
ylabel(PM.labels{1})
h(3) = plot(tpFailureSurface, hsFailureSurface, '-r');
h(4) = plot(tpAE(iMaxAE), hsAE(iMaxAE), 'xk', 'linewidth', 2);
h(5) = plot(tpHD(iMaxHD), hsHD(iMaxHD), 'xk', 'linewidth', 2);

xlim([0 20]);
ylim([0 20]);
lCell = {'50-year AE contour', '50-year HD contour', ...
    'Failure surface, r_{SDOF}=19', 'Maximum response'};
legend(h(1:4), lCell, 'location', 'northoutside', 'numcolumns', 2);
legend box on
box off


figResponseSurfaces = figure('position', [100 100 500 450]);
hold on
%xlabel(PM.labels{2})
xlabel('Spectral peak period (s)');
ylabel(PM.labels{1})
responseValues = [5 : 5 : 35];
legendCell = cell(length(responseValues), 1)
h = zeros(length(responseValues), 1);
for i = 1 : length(responseValues)
    tp = [0 : 0.5 : 20];
    hs = hsToReachRoss2020Response(tp, responseValues(i));
    h(i) = plot(tp, hs);
    legendCell{i} = num2str(responseValues(i),'r_{SDOF}=%-d');
end
plot([7 7], [0 20], '--k');
legend(h, legendCell, 'orientation', 'vertical', 'location', 'northoutside', 'NumColumns', 4);
xlim([0 18]);
ylim([0 20]);
box off
legend box on
