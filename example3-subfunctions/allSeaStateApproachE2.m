function r = allSeaStateApproachE2(alpha, phi, a, b)
    r = fzero(@(x)longTermResponseCdfE2(x, phi, a, b) - (1 - alpha), 100);
end
