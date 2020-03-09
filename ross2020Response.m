function r = ross2020Response(hs, tp)
    alpha = 2;
    beta = 0.007;
    tp0 = 30;
   
    r = alpha * hs ./ (1 + beta .* (tp - tp0).^2);
end