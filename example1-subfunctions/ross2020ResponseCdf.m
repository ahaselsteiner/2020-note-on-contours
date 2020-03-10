function p = ross2020ResponseCdf(x, hs, tp)
    % ROSS2020RESPONSECDF returns the short term CDF of the response.
    % Because the response is deterministic the short term CDF returns 
    % either 0 or 1.
    
    r = ross2020Response(hs, tp);
    p = double(r <= x); 
end