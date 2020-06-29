function p = responseTwoPeaksCdf(x, hs, tp)
    % RESPONSETWOPEAKSCDF returns the short term CDF of the response.
    % Because the response is deterministic the short term CDF returns 
    % either 0 or 1.
    
    r = responseTwoPeaks(hs, tp);
    p = double(r <= x); 
end