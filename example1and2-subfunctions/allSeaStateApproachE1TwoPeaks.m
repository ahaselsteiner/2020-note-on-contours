function r = allSeaStateApproachE1TwoPeaks(alpha)
    r = fzero(@(x)longTermResponseCdfE1TwoPeaks(x) - (1 - alpha), 20);
end
