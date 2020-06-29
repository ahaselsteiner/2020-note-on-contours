function r = allSeaStateApproachE1Ross2020(alpha)
    r = fzero(@(x)longTermResponseCdfE1Ross2020(x) - (1 - alpha), 20);
end
