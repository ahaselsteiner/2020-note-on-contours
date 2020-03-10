function [rcx, rcy] = computeResponseSurfaceDirectional(r)
    hxPrime = -15 : 0.01 : 15;
    hxPrime = hxPrime';
    hyPrime1 = nan(length(hxPrime), 1);
    hyPrime2 = nan(length(hxPrime), 1);
    for i = 1 : length(hxPrime)
        a = 1.3;
        b = 5;
        x = hxPrime(i);
        y = sqrt((r^2 - a*x^2) / b);
        if isreal(y)
            hyPrime1(i) = real(y);
            hyPrime2(i) = -1 * hyPrime1(i);
        end
    end

    % Rotate the ellipse.
    theta = 310 / 180 * pi;
    rcxPrime = [hxPrime; flip(hxPrime)];
    rcyPrime = [hyPrime1; flip(hyPrime2)];
    rcx = rcxPrime * cos(theta) - rcyPrime * sin (theta);
    rcy = rcxPrime * sin(theta) + rcyPrime * cos(theta);


    rcx = rcx(~isnan(rcy));
    rcy = rcy(~isnan(rcy));
end
