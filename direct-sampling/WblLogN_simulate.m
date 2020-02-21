function [hs, tz] = WblLogN_simulate(N, alpha, beta, gamma, ...
    a1, a2, a3, b1, b2, b3)
    if ~exist('alpha', 'var')
        % Hs Weibull parameters.
        alpha = 2.776;
        beta = 1.471;
        gamma = 0.8888;

        % Tp conditional parameters.
        a1 = 0.100;
        a2 = 1.489;
        a3 = 0.190;
        b1 = 0.040;
        b2 = 0.175;
        b3 = -0.224;
    end

    hs = wblrnd(alpha, beta, N, 1) + gamma;

    mu = a1 + a2 * hs.^a3;
    sigma = b1 + b2 * exp(b3 * hs);

    tz = lognrnd(mu, sigma, N, 1);
end