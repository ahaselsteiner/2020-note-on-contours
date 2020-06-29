function p = directionalResponseCdf(x, hx, hy, phi, a, b)  
    r = directionalResponse(hx, hy, phi, a, b);
    p = double(r <= x);
end