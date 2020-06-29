function hs = hsToReachRoss2020Response(tp, r)
    % HSTOREACHROSS2020RESPONSE returns the hs value that is required to
    % reach a given response, r. It uses the deterministic response
    % function that is given in the publication by Ross et al. (2020), 
    % DOI: 10.1016/j.oceaneng.2019.106194 at page 10.
    
    a1 = 4;
    a2 = 1.1;
    b1 = 0.1;
    b2 = 0.05;
    te1 = 25;
    te2 = 12.5;

    % The response function reads
    % r = hs .* (a1 ./ (1 + b1 .* (tp - te1).^2) + a2 ./ (1 + b2 .* (tp - te2).^2));
    % such that hs is:    
    hs = r ./ (a1 ./ (1 + b1 .* (tp - te1).^2) + a2 ./ (1 + b2 .* (tp - te2).^2));
end