function r = allSeaStateApproachE2(alpha)
    r = fzero(@(x)longTermResponseCdfE2(x) - (1 - alpha), 100);
end
