function hs = hsToReachRoss2020Response(tp, r)
    % HSTOREACHROSS2020RESPONSE returns the hs value that is required to
    % reach a given response, r. It uses the deterministic response
    % function that is given in the publication by Ross et al. (2020), 
    % DOI: 10.1016/j.oceaneng.2019.106194 at page 10.
    
    alpha = 2;
    beta = 0.007;
    tp0 = 30;
    
    % The response function reads r = alpha * hs ./ (1 + beta .* (tp - tp0).^2)
    % such that hs is:    
    hs = r * (1 + beta .* (tp - tp0).^2) ./ alpha;
end