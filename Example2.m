% Script for Example 2 presented in "Marginal and Total Exceedance 
% Probabilities for Environmental Contours" by Ed Mackay and Andreas
% F. Haselsteiner.

addpath('compute-hdc')
addpath('example2-subfunctions')

% Define the contour's exceedance probability.
ALPHA = 1 / (1 * 365.25 * 24/3);

% We are using a complex distribution with many parameter to describe 
% directional significant wave height. Thus, we will use numerical 
% integration instead of using the analytical expressions when we calculate
% contours.
% Define a grid at which the PDF will be evaluated:
m = 20;
dx = 0.02; %0.02 is sufficient for ALPHA = 1e-4
hx = -m : dx : m;
hy = -m : dx : m;
dx = hx(2) - hx(1);
dy = hy(2) - hy(1);
[HX, HY] = meshgrid(hx, hy);

% Calculate density and check integral.
fxy = directionalDistributionDensityCartesian(HX, HY);
fxy(HX == 0 & HY == 0) = 0;
err = sum(fxy(:)) * dx * dy - 1; % Error in numerical integration.
%disp(['Error in numerical integration: ' num2str(err)]);

% Simulate random points for DS contour
N = 1e6;
[theta, Hs] = directionalDistributionSimulate(N);
xr = Hs .* cos(theta);
yr = Hs .* sin(theta);

% Calculate contours.
[xcont_IF, ycont_IF] = count2iform(hx, hy, fxy, ALPHA, 0.25);
[xcont_IS, ycont_IS] = count2isorm(hx, hy, fxy, ALPHA, 0.25);
[xcont_HD, ycont_HD] = hdr2D(hx, hy, fxy, ALPHA);
[xcont_DS, ycont_DS] = dsContour2D(xr, yr, ALPHA, 2);

% Calculate the maximum responses along the contours.
[maxResponseIF, i] = max(directionalResponse(xcont_IF, ycont_IF));
xMaxResponseIF = xcont_IF(i);
yMaxResponseIF = ycont_IF(i);
[maxResponseIS, i]= max(directionalResponse(xcont_IS, ycont_IS));
xMaxResponseIS = xcont_IS(i);
yMaxResponseIS = ycont_IS(i);
[maxResponseHD, i]= max(directionalResponse(xcont_HD{1}, ycont_HD{1}));
xMaxResponseHD = xcont_HD{1}(i);
yMaxResponseHD = ycont_HD{1}(i);
[maxResponseDS, i] = max(directionalResponse(xcont_DS, ycont_DS));
xMaxResponseDS = xcont_DS(i);
yMaxResponseDS= ycont_DS(i);

% Transform to Hs-theta.
Hs_IF=sqrt(xcont_IF.^2 + ycont_IF.^2);
Hs_IS=sqrt(xcont_IS.^2 + ycont_IS.^2);
Hs_DS=sqrt(xcont_DS.^2 + ycont_DS.^2);
Hs_HD=sqrt(xcont_HD{1}.^2 + ycont_HD{1}.^2);
theta_IF=atan2(ycont_IF,xcont_IF);
theta_IS=atan2(ycont_IS,xcont_IS);
theta_DS=atan2(ycont_DS,xcont_DS);
theta_HD=atan2(ycont_HD{1},xcont_HD{1});
[theta_IF,ind]=sort(theta_IF);
Hs_IF=Hs_IF(ind);
[theta_IS,ind]=sort(theta_IS);
Hs_IS=Hs_IS(ind);
[theta_DS,ind]=sort(theta_DS);
Hs_DS=Hs_DS(ind);
[theta_HD,ind]=sort(theta_HD);
Hs_HD=Hs_HD(ind);

HsMaxResponseIF = sqrt(xMaxResponseIF^2 + yMaxResponseIF^2);
HsMaxResponseIS = sqrt(xMaxResponseIS^2 + yMaxResponseIS^2);
HsMaxResponseDS = sqrt(xMaxResponseDS^2 + yMaxResponseDS^2);
HsMaxResponseHD = sqrt(xMaxResponseHD^2 + yMaxResponseHD^2);
thetaMaxRespIF = atan2(yMaxResponseIF, xMaxResponseIF);
thetaMaxRespIS = atan2(yMaxResponseIS, xMaxResponseIS);
thetaMaxRespDS = atan2(yMaxResponseDS, xMaxResponseDS);
thetaMaxRespHD = atan2(yMaxResponseHD, xMaxResponseHD);

% Compute true response.
allSeaStateResponse = allSeaStateApproachE2(ALPHA);

figure('position',[0, 0, 600, 600])
h = zeros(9, 1);
h(1) = polarplot(theta_IF,Hs_IF, 'b-','linewidth', 1); hold on
h(2) = polarplot(theta_IS,Hs_IS, 'k-','linewidth', 1);
h(3) = polarplot(theta_DS,Hs_DS, 'k--','linewidth', 1);
h(4) = polarplot(theta_HD,Hs_HD, 'b--','linewidth', 1);

[rcx, rcy] = computeResponseSurfaceDirectional(allSeaStateResponse);
hs = sqrt(rcx.^2 + rcy.^2);
waveAngle = atan2d(rcy, rcx);
h(5) = polarplot(deg2rad([waveAngle; waveAngle(1)]), [hs; hs(1)], '-r', 'linewidth', 2);

h(6) = polarplot(thetaMaxRespIF, HsMaxResponseIF, 'xb');
h(7) = polarplot(thetaMaxRespIS, HsMaxResponseIS, 'xk');
h(8) = polarplot(thetaMaxRespDS, HsMaxResponseDS, 'xk');
h(9) = polarplot(thetaMaxRespHD, HsMaxResponseHD, 'xb');

title('$H_s$ [m] vs. $\theta$ [deg]', 'interpreter', 'latex')
ax = gca;
set(ax,'TickLabelInterpreter', 'latex')
set(ax,'RTick',0 : 2 : 12)
L = legend(h([1,2,3,4,5,6]), ...
    'IFORM', 'ISORM', 'DS', 'HD', 'Failure surface', 'Max response');
set(L, 'Box', 'off')
set(L,'Position', [0.8 0.07 0.18 0.13])

% Compute probability of failures if structures were designed such that
% their capacity were exactly the maximum response of the contour.
pfIF = (1 - longTermResponseCdfE2(maxResponseIF)) / ALPHA;
pfDS = (1 - longTermResponseCdfE2(maxResponseDS)) / ALPHA;
pfIS = (1 - longTermResponseCdfE2(maxResponseIS)) / ALPHA;
pfHD = (1 - longTermResponseCdfE2(maxResponseHD)) / ALPHA;

contourName =   {'IFORM';         'DS contour';    'ISORM';         'HD contour';    'All sea state'};
maxResponse =   [maxResponseIF;   maxResponseDS;   maxResponseIS;   maxResponseHD;   allSeaStateResponse];
maxResponseHs = [HsMaxResponseIF; HsMaxResponseDS; HsMaxResponseIS; HsMaxResponseHD; NaN];
maxRespTheta =  [thetaMaxRespIF;  thetaMaxRespDS;  thetaMaxRespIS;  thetaMaxRespHD;  NaN];
probOfFailure = [pfIF;            pfDS;            pfIS;            pfHD;            1];

maxRespTheta = maxRespTheta / pi * 180;
maxRespTheta(maxRespTheta < 0) = 360 + maxRespTheta(maxRespTheta < 0);

T = table(contourName, maxResponse, maxResponseHs, maxRespTheta, probOfFailure)
