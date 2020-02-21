function hs = hsToReachRoss2020Response(tp, r)
    alpha = 2;
    beta = 0.007;
    tp0 = 7;
    
    % r = alpha * hs ./ (1 + beta .* (tp - tp0).^2);
    hs = r * (1 + beta .* (tp - tp0).^2) ./ alpha;
end