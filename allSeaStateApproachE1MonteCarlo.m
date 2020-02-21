function [rMean, rStd] = allSeaStateApproachE1MonteCarlo(alpha, N, B)
    responses = nan(B, 1); 
    for i = 1 : B
        if N > 1e6
            fprintf('   bootstrap %u/%u\n', i, B)
        end
        responses(i) = allSeaStateApproachResponseE1(alpha, N);
    end
    rMean = mean(responses);
    rStd = std(responses);
end
function r  = allSeaStateApproachResponseE1(alpha, N)
    [hs, tz] = WblLogN_simulate(N);
    tztpCoeff = 1.2796;
    tp =  tztpCoeff * tz;
    responses = ross2020Response(hs, tp);
    r = quantile(responses, 1 - alpha);
end