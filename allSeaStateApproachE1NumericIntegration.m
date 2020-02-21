function r = allSeaStateApproachE1NumericIntegration(alpha)
    %TR = 50;
    %TSS = 6;
    %alpha = 1 / (TR * 365.25 * 24 / TSS);
    
    r = fzero(@(x)responseCdf(x) - (1 - alpha), 20);
end

function p = responseCdf(x)
    p = integral2(@(hs, tz)rCdfTimesSeaStateDensity(x, hs, tz), 0, 30, 0, 40);
end

function r  = rCdfTimesSeaStateDensity(x, hs, tz)
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

    Hs = TranslatedWeibull(alpha, beta, gamma);
    fHs = Hs.pdf(hs);

    mu = a1 + a2 * hs.^a3;
    sigma = b1 + b2 * exp(b3 * hs);
    fTp = lognpdf(tz, mu, sigma);

    fjoint = fHs .* fTp;

    tztpCoeff = 1.2796;
    tp =  tztpCoeff * tz;
    r  = ross2020ResponseCdf(x, hs, tp) .* fjoint;
    %r = fjoint;
end