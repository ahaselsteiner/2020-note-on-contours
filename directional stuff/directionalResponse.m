function r = directionalResponse(hx, hy)
    theta = -310 / 180 * pi;
    a = 1.3;
    b = 5;
    hxPrime = hx * cos(theta) - hy * sin(theta);
    hyPrime = hx * sin(theta) + hy * cos(theta);
    r = sqrt(a * hxPrime.^2 + b * hyPrime.^2);
end