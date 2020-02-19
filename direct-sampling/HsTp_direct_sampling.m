% Hs Weibull parameters.
alpha = 2.776;
beta = 1.471;
gamma = 0.8888;

% Tp conditional parameters.
tztpCoeff = 1.2796;
a1 = 0.100;
a2 = 1.489;
a3 = 0.190;
b1 = 0.040;
b2 = 0.175;
b3 = -0.224;

% Number of points and non-exceedance probabilities for contours.
N = 1e5;
P = 1 - [1e-1 1e-2 1e-3];

% Simulate data from joint distribution.
disp('Simulating data...')
tic
[h, t] = WblLogN_simulate(N, alpha, beta, gamma, ...
    tztpCoeff, a1, a2, a3, b1, b2, b3);
toc

% Calculate AE contour.
disp('Calculating contour...')
tic
[tcont, hcont] = direct_sampling_contour(t, h, P, 5);
toc

% Calculate empirical joint occurrence.
disp('Counting joint occurrence...')
tic
dh = 0.2;
dt = 0.2;
tedges = floor(min(t)) : dt : ceil(max(t));
hedges = floor(min(h)) : dt : ceil(max(h));
tc = tedges(1 : end - 1) + dt / 2;
hc = hedges(1 : end - 1) + dh / 2;
count = histcounts2(h, t, hedges, tedges);
count(count==0) = NaN;
toc

% Plot contour.
figure
hold on
nanimage(tc, hc, log10(count))
plot([tcont; tcont(1,:)], [hcont; hcont(1,:)], 'r-')
xlabel('Spectral peak period (s)');
ylabel('Significant wave height (m)');