function p = ross2020ResponseCdf(x, hs, tp)
    alpha = 2;
    beta = 0.007;
    tp0 = 7;
   
    r = alpha * hs ./ (1 + beta .* (tp - tp0).^2);
    
    p = double(r <= x);
end