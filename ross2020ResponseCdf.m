function p = ross2020ResponseCdf(x, hs, tp)
    r = ross2020Response(hs, tp);
    p = double(r <= x);
end