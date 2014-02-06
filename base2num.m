function [ base ] = base2num(b)

    A = zeros(1,84);
    A(65) = 1;
    A(67) = 2;
    A(71) = 3;
    A(84) = 4;

    base = A(b + 0);
end