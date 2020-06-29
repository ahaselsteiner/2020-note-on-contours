function p = longTermResponseCdfE1TwoPeaks(x)
    p = integral2(@(hs, tz)rCdfTimesSeaStateDensity(x, hs, tz), 0, 40, 0, 40);
end

function r  = rCdfTimesSeaStateDensity(x, hs, tz)
    % Hs Weibull parameters.
    alpha = 2.776;
    beta = 1.471;
    gamma = 0.8888;

    % Tp conditional parameters.
    a1 = 0.100;
    a2 = 1.489;
    a3 = 0.1901;
    b1 = 0.040;
    b2 = 0.1748;
    b3 = -0.2243;

    %Hs = TranslatedWeibull(alpha, beta, gamma);
    %fHs = Hs.pdf(hs);
    fHs = wblpdf(hs - gamma, alpha, beta);

    mu = a1 + a2 * hs.^a3;
    sigma = b1 + b2 * exp(b3 * hs);
    fTp = lognpdf(tz, mu, sigma);

    fjoint = fHs .* fTp;

    tztpCoeff = 1.2796;
    tp =  tztpCoeff * tz;
    r  = responseTwoPeaksCdf(x, hs, tp) .* fjoint;
end