function r = responseTwoPeaks(hs, tp)
    % RESPONSETWOPEAKS is a deterministic response function with two peaks
    % that represents an extension to the deterministic response 
    % function that is is given in the publication by Ross et al. (2020), 
    % DOI: 10.1016/j.oceaneng.2019.106194 at page 10.
    
    a1 = 4;
    a2 = 1.1;
    b1 = 0.1;
    b2 = 0.05;
    te1 = 25;
    te2 = 12.5;

   
    r = hs .* (a1 ./ (1 + b1 .* (tp - te1).^2) + a2 ./ (1 + b2 .* (tp - te2).^2));
end
