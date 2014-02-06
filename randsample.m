function [ x ] = randsample(v,probs) 
% Samples a single element from the vector v with probability
% as specified in probs.

    probs = probs / sum(probs); % Normalise
    c = cumsum(probs);
    r = rand(1);
    
    i = 1;
    while (r >= c(i)) 
        i = i+1;
    end
    
    x = v(i);
end