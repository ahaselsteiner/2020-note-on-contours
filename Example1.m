% Script for Example 1 presented in "Marginal and Total Exceedance 
% Probabilities for Environmental Contours" by Ed Mackay and Andreas
% F. Haselsteiner.

addpath('example1-subfunctions')
addpath('compute-hdc')

% Define a joint distribution for the sea state. The used model was 
% proposed by Vanem and Bitner-Gregersen, DOI: 10.1016/j.apor.2012.05.006 .
PM.name = 'Vanem and Bitner-Gregerse (2012), DOI: 10.1016/j.apor.2012.05.006';
PM.modelType = 'CMA';
PM.distributions = {'weibull'; 'lognormal'};
PM.isConditionals = {[0 0 0]; [1 1]};
PM.coeffs = {{2.776 1.471 0.8888}; 
                             { @(x1)0.1000 + 1.489 * x1.^0.1901;
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

% Compute the contours: IFORM, DS, ISORM and HD.
[hsIFORM, tzIFORM] = computeIformContour(PM, alpha, 360); 
tpIFORM = tztpCoeff * tzIFORM;
% Compute a line for each 5 degree, 1 degree would lead to numerical
% issues at the region where the contour's curvature is high.
[hsDS, tzDS] = computeDsContour(PM, 1 / alpha * 200, alpha, 5);
tpDS = tztpCoeff * tzDS;
[hsISORM, tzISORM] = computeIsormContour(PM, alpha, 360); 
tpISORM = tztpCoeff * tzISORM;
[fm, hsHD, tzHD] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);
hsHD = hsHD{1};
tpHD = tztpCoeff * tzHD{1};


% Compute the true response using numerical integration.
disp('Now starting numerical intergration for all sea state approach.');
tic
allSeaStatesResp = allSeaStateApproachE1NumericIntegration(alpha);
tic

% Compute the failure surface.
tpFailureSurface = [0 : 0.5 : 25];
hsFailureSurface = hsToReachRoss2020Response(tpFailureSurface, 15);

% Plot the contours and the failure surface.
figFourContours = figure('position', [100 100 500 450]);
hold on
h = zeros(8, 1);
h(1) = plot([tpIFORM tpIFORM(1)], [hsIFORM hsIFORM(1)], '-b'); % IFORM
responseIFORM = ross2020Response(hsIFORM, tpIFORM);
[maxRespIFORM, iMaxIFORM] = max(responseIFORM);
h(2) = plot(tpIFORM(iMaxIFORM), hsIFORM(iMaxIFORM), '+b', 'linewidth', 1);
h(3) = plot([tpDS; tpDS(1)], [hsDS; hsDS(1)], '-k'); % DS
responseDS = ross2020Response(hsDS, tpDS);
[maxRespDS, iMaxDS] = max(responseDS);
h(4) = plot(tpDS(iMaxDS), hsDS(iMaxDS), '+k', 'linewidth', 1);
h(5) = plot([tpISORM tpISORM(1)], [hsISORM hsISORM(1)], '--b'); % ISORM
responseISORM = ross2020Response(hsISORM, tpISORM);
[maxRespISORM, iMaxISORM] = max(responseISORM);
h(6) = plot(tpISORM(iMaxISORM), hsISORM(iMaxISORM), '+b', 'linewidth', 1);
h(7) = plot(tpHD, hsHD, '--k'); % HD
responseHD = ross2020Response(hsHD, tpHD);
[maxRespHD, iMaxHD] = max(responseHD);
h(8) = plot(tpFailureSurface, hsFailureSurface, '-r', 'linewidth', 2); % Failure surface
h(9) = plot(tpHD(iMaxHD), hsHD(iMaxHD), '+k', 'linewidth', 1);

xlabel('Spectral peak period (s)');
ylabel(PM.labels{1})
xlim([0 25]);
ylim([0 20])
legend(h([1,3,5,7,8,9]), {'IFORM', 'Direct sampling', ...
    'ISORM', 'Highest density', ...
    'Failure surface', 'Max.response'}, ...
    'location', 'northwest');
legend boxoff
box off

% Create a table with the mail results.
contourName = {'IFORM'; 'Direct sampling'; 'ISORM'; 'Highest density'; 'All sea states'};
maxResponse = [maxRespIFORM; maxRespDS; maxRespISORM; maxRespHD; allSeaStatesResp];
maxResponseHs = [hsIFORM(iMaxIFORM); hsDS(iMaxDS); hsISORM(iMaxISORM); hsHD(iMaxHD); NaN];
maxResponseTp = [tpIFORM(iMaxIFORM); tpDS(iMaxDS); tpISORM(iMaxISORM); tpHD(iMaxHD); NaN];
T = table(contourName, maxResponse, maxResponseHs, maxResponseTp)
