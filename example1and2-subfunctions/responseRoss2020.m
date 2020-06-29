function r = ross2020Response(hs, tp)
    % ROSS2020RESPONSE is the deterministic response function
    % that is given in the publication by Ross et al. (2020), 
    % DOI: 10.1016/j.oceaneng.2019.106194 at page 10.

    alpha = 2;
    beta = 0.007;
    tp0 = 30;
   
    r = alpha * hs ./ (1 + beta .* (tp - tp0).^2);
end