function fxy = directionalDistributionDensityCartesian(hx, hy)

hs = sqrt(hx.^2 + hy.^2);
theta = atan2(hy, hx);

% Parameters.
alpha0 = 3.75;
beta0 = 1.91;
a_alpha = [0.69	-0.42	-0.32	-0.53	-0.18	0.14	0.06	0.06];
a_beta = [0.24	-0.08	-0.01	-0.11	-0.0004	0.06	0.06	0.004];
b_alpha = [-0.28	-1.64	-0.4	0.19	0.22	0.14	0.04	-0.03];
b_beta = [-0.13	-0.17	-0.03	0.03	0.003	0.02	0.02	-0.01];
w = [0.21 0.79];
mu = [2.1 5.54];
k = [0.74 13.11];
gamma = 0.5;

% Directional PDF.
ftheta = 0;
for i = 1 : 2
    ftheta = ftheta + w(i) * exp(k(i) * cos(theta - mu(i))) / ...
        (2 * pi * besseli(0, k(i)));
end

% Weibull parameters.
alpha = alpha0;
beta = beta0;
for j = 1 : length(a_alpha)
    alpha = alpha + a_alpha(j) * cos(j * theta) + b_alpha(j) * sin(j * theta);
    beta = beta + a_beta(j) * cos(j * theta) + b_beta(j) * sin(j * theta);
end

% Conditional weibull PDF.
alpha = alpha / 2;
fHs = wblpdf(hs - gamma, alpha, beta);

% Joint PDF.
f = fHs .* ftheta;
fxy = f ./ hs;
