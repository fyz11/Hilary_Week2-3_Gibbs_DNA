function prob = felsenstein_likelihood(seq,start_positions,k,p)

% Mar. 2013
% Joe Herman <herman@stats.ox.ac.uk>
% Adam Novak <novak@stats.ox.ac.uk>
% with thanks to Andreas Harris <andreas.harris@exeter.ox.ac.uk>
        
% Tree represented as
%            1
%     0.1 /     \0.3
%      2          3
% 0.4/  \0.2  0.2/  \0.14
%   4    5      6    7
% 
% Tree = { [0.1 0.3] [0.4 0.2] [0.2 0.14] };
% Children = { [2 3] [4 5] [6 7] };

% i.e. the Tree cell array has left and right branch lengths for each
% internal node

% seq will be a cell array of sequences passed in

% Here's tree.tree converted into our format
Tree = { [0.512 0.712] [0.297 0.597] [0.365 0.464] [0.464 0.564] [0.143 0.343] };
Children = { [6 2] [3 4] [7 5] [10 11] [8 9] };

nseqs = length(seq);
column = zeros(nseqs,1);
for (i = 1:nseqs) 
    % select the relevant bits of the sequences as specified by
    % start_positions
    column(i) = seq{i}(start_positions(i)+k-1);
end

% Create the rate matrix according to the Felsenstein 1981 (F81) model
Q = [p; p; p; p];
Q = Q - diag(diag(Q));
Q = Q - diag(sum(Q'));
[V,D] = eig(Q); 
D = diag(D);
Vinv = inv(V);

prior = p; % stationary dist

% Recursive function to compute the likelihood of observing the
% characters below a particular internal node i of the tree.
% Returns a vector of length 4, containing the likelihoods if
% there is an A,C,G or T at node i.
function [ x ] = fels(i)

    if (i >= nseqs)
        x = zeros(4,1);
        x(column(i-nseqs+1)) = 1;
    else      
        x = ((V * diag(exp(Tree{i}(1) * D)) * Vinv) * fels(Children{i}(1))) .* ...
        ((V * diag(exp(Tree{i}(2) * D)) * Vinv) * fels(Children{i}(2)));
    end
end

% Compute the posterior probability of the characters in this column
% as the dot product of the prior probability for each character at
% the root, and the likelihood of the observed characters at the
% tips given each possible root character (ACGT).
prob = prod(prior * fels(1));

end