function [rcx, rcy] = computeResponseSurfaceDirectional(r, phi, a, b)    
    phi
    a
    b
    hxPrime = -15 : 0.01 : 15;
    hxPrime = hxPrime';
    hyPrime1 = nan(length(hxPrime), 1);
    hyPrime2 = nan(length(hxPrime), 1);
    for i = 1 : length(hxPrime)

        x = hxPrime(i);
        y = sqrt((r^2 - a*x^2) / b);
        if isreal(y)
            hyPrime1(i) = real(y);
            hyPrime2(i) = -1 * hyPrime1(i);
        end
    end

    % Rotate the ellipse.
    rcxPrime = [hxPrime; flip(hxPrime)];
    rcyPrime = [hyPrime1; flip(hyPrime2)];
    rcx = rcxPrime * cos(phi) - rcyPrime * sin (phi);
    rcy = rcxPrime * sin(phi) + rcyPrime * cos(phi);


    rcx = rcx(~isnan(rcy));
    rcy = rcy(~isnan(rcy));
end
