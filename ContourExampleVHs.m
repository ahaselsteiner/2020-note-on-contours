% Example wind-wave contour.

PM.name = 'V-Hs model from X for dataset D (FINO1)';
PM.modelType = 'CMA';
PM.distributions = {'exponentiated-weibull'; 'exponentiated-weibull'};
PM.isConditionals = {[0 0 0 ]; [1 1 1]};
PM.coeffs = {
    {10.0 2.42 0.761}; 
    { 
    @(x1) (0.394 + 0.0178 * x1^1.88) / (2.0445^(1 / (0.582 + 1.90 / (1 + exp(-0.248 * (x1 - 8.49))))));
    @(x1) 0.582 + 1.90 / (1 + exp(-0.248 * (x1 - 8.49)));
    @(x1) 5}
    };
PM.labels = {'Wind speed (m/s)';
    'Significant wave height (m)'};
PM.gridCenterPoints = {0:0.05:50; 0:0.05:30};
%PM = getProbabilisticModel(9);


nYears = 50;
stateDuration = 1;
alpha = 1 / (nYears * 365.25 * 24 / stateDuration);
[vIFORM, hsIFORM] = computeIformContour(PM, alpha, 360); 
[fm, vHDC, hsHDC] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);

figContour = figure('position', [100 100 350 300]);
hold on

plot([vIFORM vIFORM(1)], [hsIFORM hsIFORM(1)], '-k');
plot(vHDC{1}, hsHDC{1}, '--k');
xlabel(PM(1).labels{1})
ylabel(PM(1).labels{2})

v = [0:1:60];
hs = computeVToReachHorizontalForce(v, 1, 5.5*10^(-2), 3.8*10^(-4));
plot(v, hs, '-r');
xlim([0 60]);
ylim([0 21]);

legend({'50-year IFORM contour', '50-year HD contour', 'Failure surface'});
legend box off
box off


function r = totalHorziontalForce(v, hs, a, b)
    r = a*hs + b*v.^2;
end

function hs = computeVToReachHorizontalForce(v, r, a, b)
    % r = a*hs + b*v.^2;
    hs = (r - b * v.^2) / a;
end