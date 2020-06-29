function p = responseRoss2020Cdf(x, hs, tp)
    % RESPONSEROSS2020CDF returns the short term CDF of the response.
    % Because the response is deterministic the short term CDF returns 
    % either 0 or 1.
    
    r = responseRoss2020(hs, tp);
    p = double(r <= x); 
end