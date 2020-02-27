clear
%plotsettings
%figurefolder='C:\Users\em604\OneDrive - University of Exeter\Extremes work\Papers\7 - Contours\figures';
%load('FigureFormat_ACP_v3.mat');
formatEPS.FontMode='auto';
formatEPS.LockAxesTicks='on';
formatEPS.LineMode='auto';
formatEPS.Bounds='loose';

% range and resolution for density
m=20;
dx=0.02; %0.02 is sufficient for alpha = 1e-4
x=-m:dx:m;
y=-m:dx:m;
dx=x(2)-x(1);
dy=y(2)-y(1);
[X,Y]=meshgrid(x,y);

% calculate density and check integral
fxy=directional_distribution_density_cartesian(X,Y);
fxy(X==0&Y==0)=0;
sum(fxy(:))*dx*dy-1 % error in numerical integration

% simulate random points for DS contour
N=1e6;
[theta,Hs]=directional_distribution_simulate(N);
xr=Hs.*cos(theta);
yr=Hs.*sin(theta);

% calculate contours
alpha=1e-4;
P=1-alpha;
[xcont_IF,ycont_IF]=count2iform(x,y,fxy,P);
[xcont_IS,ycont_IS]=count2isorm(x,y,fxy,P);
[xcont_HD,ycont_HD]=hdr2D(x,y,fxy,alpha);
[xcont_DS, ycont_DS]=direct_sampling_contour(xr,yr,P,2);

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

% transform to Hs-theta
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
thetaMaxResponseIF = atan2(yMaxResponseIF, xMaxResponseIF);
thetaMaxResponseIS = atan2(yMaxResponseIS, xMaxResponseIS);
thetaMaxResponseDS = atan2(yMaxResponseDS, xMaxResponseDS);
thetaMaxResponseHD = atan2(yMaxResponseHD, xMaxResponseHD);

% Compute true response.
allSeaStateResponse = allSeaStateApproachE2(alpha);

% figure
% hold on
% % scatter(xr,yr)
% % contour(x,y,log10(fxy),-6:0)
% plot(xcont_IF,ycont_IF,'r')
% plot(xcont_IS,ycont_IS,'g')
% plot(xcont_DS,ycont_DS,'k')
% plot(xcont_HD{1},ycont_HD{1},'b')

figure('position',[-800,300,400,400])
hold on; grid on
plot(theta_IF,Hs_IF,'r')
plot(theta_DS,Hs_DS,'k--')
plot(theta_IS,Hs_IS,'b','linewidth',2)
plot(theta_HD,Hs_HD,'g--','linewidth',2)
L=legend('IFORM','DS','ISORM','HD');
set(L,'Box','off')

figure('position',[-800,300,400,400])
h = zeros(9, 1);
h(1) = polarplot(theta_IF,Hs_IF,'b-','linewidth',1); hold on
h(2) = polarplot(theta_IS,Hs_IS,'k-','linewidth',1);
h(3) = polarplot(theta_DS,Hs_DS,'k--','linewidth',1);
h(4) = polarplot(theta_HD,Hs_HD,'b--','linewidth',1);

[rcx, rcy] = computeResponseSurfaceDirectional(allSeaStateResponse);
hs = sqrt(rcx.^2 + rcy.^2);
waveAngle = atan2d(rcy, rcx);
h(5) = polarplot(deg2rad([waveAngle; waveAngle(1)]), [hs; hs(1)], '-r', 'linewidth', 2);

h(6) = polarplot(thetaMaxResponseIF, HsMaxResponseIF, '+k');
h(7) = polarplot(thetaMaxResponseIS, HsMaxResponseIS, '+k');
h(8) = polarplot(thetaMaxResponseDS, HsMaxResponseDS, '+k');
h(9) = polarplot(thetaMaxResponseHD, HsMaxResponseHD, '+k');

title('$H_s$ [m] vs. $\theta$ [deg]')
ax=gca;
set(ax,'TickLabelInterpreter','latex')
set(ax,'RTick',0:2:12)
L=legend(h([1,2,3,4,5,6]),'IFORM','ISORM','DS','HD', 'Failure surface', 'Max response');
set( L, 'Box','off')
set(L,'Position',[0.8 0.07 0.18 0.13])
%hgexport(gcf, [figurefolder '\DirectionCont.eps'], formatEPS);


pfIF = (1 - longTermResponseCdfE2(maxResponseIF)) / alpha;
pfDS = (1 - longTermResponseCdfE2(maxResponseDS)) / alpha;
pfIS = (1 - longTermResponseCdfE2(maxResponseIS)) / alpha;
pfHD = (1 - longTermResponseCdfE2(maxResponseHD)) / alpha;
contourName = {'IFORM'; 'DS contour'; 'ISORM'; 'HD contour'; 'All sea state'};
maxResponse =      [maxResponseIF; maxResponseDS; maxResponseIS; maxResponseHD; allSeaStateResponse];
maxResponseHs =    [HsMaxResponseIF; HsMaxResponseDS; HsMaxResponseIS; HsMaxResponseHD; NaN];
maxResponseTheta = [thetaMaxResponseIF; thetaMaxResponseDS; thetaMaxResponseIS; thetaMaxResponseHD; NaN] / pi * 180;
probOfFailure = [pfIF; pfDS; pfIS; pfHD; 1];
T = table(contourName, maxResponse, maxResponseHs, maxResponseTheta, probOfFailure)
