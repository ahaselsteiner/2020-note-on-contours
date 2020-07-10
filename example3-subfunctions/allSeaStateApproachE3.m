function r = allSeaStateApproachE3(alpha, phi, a, b)
    r = fzero(@(x)longTermResponseCdfE3(x, phi, a, b) - (1 - alpha), 100);
end
