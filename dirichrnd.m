function [ x ] = dirichrnd(a) 
    x = a;
    for (i = 1:length(a))
        x(i) = gamrnd1(a(i),1);
    end
    x = x / sum(x);
end