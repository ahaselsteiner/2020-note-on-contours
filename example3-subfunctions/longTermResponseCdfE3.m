function p = longTermResponseCdfE3(x, phi, a, b)
    minInt = -30;
    maxInt = 30;
    p = integral2(@(hx, hy)rCdfTimesSeaStateDensity(x, hx, hy, phi, a, b), ...
        minInt, maxInt, minInt, maxInt);
end

function r  = rCdfTimesSeaStateDensity(x, hx, hy, phi, a, b)
    fjoint = directionalDistributionDensityCartesian(hx, hy);
    r  = directionalResponseCdf(x, hx, hy, phi, a, b) .* fjoint;
end