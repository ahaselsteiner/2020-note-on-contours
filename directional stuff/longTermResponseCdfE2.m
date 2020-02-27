function p = longTermResponseCdfE2(x)
    minInt = -30;
    maxInt = 30;
    p = integral2(@(hx, hy)rCdfTimesSeaStateDensity(x, hx, hy), ...
        minInt, maxInt, minInt, maxInt);
end

function r  = rCdfTimesSeaStateDensity(x, hx, hy)
    fjoint = directional_distribution_density_cartesian(hx, hy);
    r  = directionalResponseCdf(x, hx, hy) .* fjoint;
    %r = fjoint;
end