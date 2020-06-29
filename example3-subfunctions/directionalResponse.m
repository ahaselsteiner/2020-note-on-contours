function r = directionalResponse(hx, hy, phi, a, b)
    hxPrime = hx * cos(phi) + hy * sin(phi);
    hyPrime = -1 * hx * sin(phi) + hy * cos(phi);
    r = sqrt(a * hxPrime.^2 + b * hyPrime.^2);
end