function [ Z, S, mu, min_ent_M, min_ent_s, max_lr_M,max_lr_s, posterior_mean_M, information ]  = find_motifs(sequence_file,K, n_iterations,burn_in,a, mu_start, mu_unknown, beta)
% This code will run the Gibbs sampler motif detection algorithm of
% Lawrence et al. (1993) on a set of sequences inputted as a FASTA file. 
%
% Joe Herman, Feb. 2013 (herman@stats.ox.ac.uk)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES ON ARGUMENTS
%
%% sequence_file: a FASTA-formatted file containing the input sequences
%% K:             the length of the motif 
%% n_iterations:  number of iterations for which Gibbs sampler
%                 should be run
%% burn_in:       number of iterations to allow for the burn-in
%                 phase, while the MCMC is converging.
%% a:             a constant multiplier for the uniform prior on the
%                 motif
%% mu_start:      the starting value of mu
%% mu_unknown:    0 if mu is fixed to mu_start, 1 if mu is unknown and
%                 is to be estimated by MCMC.
%% beta:          a length-two vector, containing prior parameters
%                 for mu (used only if mu_unknown == 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = ones(1,4) * a;
% The prior distribution for each column of the motif

[head, seqs] = fastaread1(sequence_file);
% read in a set of sequences from a FASTA formatted file

nseqs = length(seqs);

for (i = 1:nseqs) 
    seqs{i} = base2num(seqs{i});
    % Turns a sequence of ACGT... into a numerical sequence of 1234...
end

background = compute_background(seqs);
% Compute the background distribution

z = ones(nseqs,1); 
% Vector of zeros and ones to indicate which sequences contain the motif.
% This is to be used when we believe that some of the sequences do not
% contain the motif, otherwise all entries will stay as 1.
Z = ones(nseqs,n_iterations-1); 
% This matrix (with the capitalised name) stores the value of
% z for each iteration of the simulation, for analysis purposes. 

mu = mu_start;
% The probability of each sequence containing a motif
Mu = ones(n_iterations,1);
Mu(1) = mu;
% The value of mu from each iteration.

s = ones(nseqs,1); 
% The start position of the motif in each sequence
S = ones(nseqs,n_iterations); 
% A matrix that stores the s from each iteration.

for (i = 1:nseqs) % For each sequence
    L_i = length(seqs{i});
    s(i) = randi([1,L_i-K+1],1);
    S(:,i) = s(i);
end

min_ent = Inf; % Keep a record of the minimum PWM entropy 
min_ent_M = zeros(4,K); % And the corresponding PWM
min_ent_s = ones(nseqs,1); 
background_entropy = -sum(background .* log(background));
entropy = zeros(n_iterations,1);
% The above variables keep track of the entropy of the PWM, and store the
% PWM and starting positions for the minimum entropy PWM.

max_lr = 0; % The maximum observed likelihood ratio
max_lr_M = zeros(4,K);
max_lr_s = ones(nseqs,1); 
lr = ones(n_iterations,1); % Record the likelihood ratio 
% The above variables are used for the purposes of finding the PWM that
% maximises the likelihood ratio between the motif model and the random
% background model, as described in the accompanying explanatory document.

posterior_mean_M = zeros(4,K);
% The posterior mean PWM is updated once the Gibbs sampler has converged
% i.e. after the burn_in is over.

background_M = repmat(background,1,K);
% This is a background PWM derived from the inputted background
% frequencies. 

for (iter = 1:n_iterations) 
    if (mod(iter,10)==0) % Every ten iterations
        iter
        % Print out where we're up to in the simulation
    end
    if (iter > 1) 
        Z(:,iter-1) = z;
        % Store the previous value of the presence/absence vector
        S(:,iter-1) = s;
        % Store the previous value of the starting positins

        if (mu_unknown)
            mu(iter) = sample_mu(z,beta);
        else
            mu(iter) = mu(iter-1);
        end
    end
    for (i = 1:nseqs) % For each sequence
    
        L_i = length(seqs{i}); 
        
        M = sample_M(seqs,K,z,s,alpha); 
        % Compute the current PWM for the motif from all sequences.
        % M will be a 4 x K array giving the probability of
        % each base at each location along the motif.
    
        prob = zeros(L_i-K+1,1);
        background_prob = zeros(L_i-K+1,1);

        for (j = 1:(L_i-K+1)) % For each potential starting position
            
            prob(j) = likelihood(seqs{i},j,M,K);
            % likelihood of the motif starting here under the model 
            % defined by M
            background_prob(j) = likelihood(seqs{i},j,background_M,K);
            % likelihood of motif starting here under the random
            % background model
        end
        
        % Now sample a new start position from the full conditional
        
        [likelihood_ratio, s(i)] = sample_s(prob,background_prob,mu(iter));
        % NB if a zero is sampled, it corresponds to the motif not being
        % present in this particular sequence.
        
        if (s(i) == 0) 
            % Then there is no motif in this sequence
            z(i) = 0;
        else
            z(i) = 1;
            % Form the product of the likelihood ratios for all of the
            % sequences
            lr(iter) = lr(iter) * likelihood_ratio;
        end
       
    end
    
    if (iter > burn_in)
        % i.e. if the chain has converged to the stationary distribution
        posterior_mean_M = posterior_mean_M + M/(n_iterations-burn_in);
    end
    % Now update the min entropy and max likelihood ratio PWMs etc.
    entropy(iter) = mean(-sum(M .* log(M))); % Mean entropy per site
    if (iter > 1 && entropy(iter) < min_ent)
        min_ent = entropy(iter);
        min_ent_M = M;
        min_ent_s = s;
    end
    if (iter > 1 && lr(iter) > max_lr)
        max_lr = lr(iter);
        max_lr_M = M;
        max_lr_s = s;
    end
    % Store the last values:
    Z(:,iter) = z;
    S(:,iter) = s;
end

% Compute the average information per site (a vector of length
% n_iterations)
information = background_entropy - entropy;

end

function [ background ] = compute_background(seqs) 
% This function returns a length-four column vector
% containing the background model for the sequences.
% The code defined below just uses a uniform distribution,
% but a better approach would be to compute the relative
% frequencies of each base in the sequences, and to use
% these as the background.
      
    background = zeros(4,1);
    for i = 1:length(seqs)
        countA = sum(sum( seqs{i}==1 ) );
        countC = sum(sum( seqs{i}==2 ) );
        countG = sum(sum( seqs{i}==3 ) );
        countT = sum(sum( seqs{i}==4 ) );
        background = background + [ countA; countC; countG; countT ];
    end
    background = background / sum(background );
    
end

function [ M ] = sample_M(seqs,K,z,s,alpha)
% Samples a new PWM from the full conditional
    M = zeros(4,K);
    % for task 1, z is fixed
    
    for slideAloneMotif = 1:K
        
        %count the frequencies
        counts = zeros(1,4);
        length_seqs = length(seqs);
        for seqs_index = 1:length_seqs
            EachSeqRow = seqs{ seqs_index };
            s_i = s(seqs_index);
            
            if (s_i ~= 0)
                motifCandidate_InEachSeq = EachSeqRow ( s_i + slideAloneMotif - 1 );
                %motifCandidate_InEachSeq;
                if motifCandidate_InEachSeq == 1
                    counts = counts + [1 0 0 0];
                end
                if motifCandidate_InEachSeq == 2
                    counts = counts + [0 1 0 0];
                end
                if motifCandidate_InEachSeq == 3
                    counts = counts + [0 0 1 0];
                end
                if motifCandidate_InEachSeq == 4
                    counts = counts + [0 0 0 1];
                end
            end
        end        
        %use the frequecies to update dir parameters

        alpha_new = alpha + counts;
        M(:,slideAloneMotif) = dirichrnd( alpha_new );

    end
        
end

function [ likelihood_ratio, s_i ] =  sample_s(prob,background_prob,mu)
    % 'prob' and 'background_prob' are vectors containing the
    % probability of the motif starting at each possible site
    % from 1 to L_i - K + 1 in sequence i, under the motif model M,
    % and the background model G, respectively. 

    % This function returns a start position sampled according to 
    % its posterior probability, along with the corresponding
    % likelihood ratio.
    % If this function returns s_i = 0, it means that sequence i
    % does not contain a copy of the motif.
    vector_loglikelyhood_ratio =  log(prob) - log(background_prob);
    vector_likelihood_ratio = exp( vector_loglikelyhood_ratio );
    
    % set z to 0 for a particular probability (i.e. just return s_i = 0;
    length_vector = length(vector_likelihood_ratio); 
    z_set_zero_probability = (length_vector*(1-mu))/(length_vector*(1-mu) + mu*sum(vector_likelihood_ratio));
    
    if (rand() < z_set_zero_probability)
        s_i = 0;
        likelihood_ratio =1;
    else
        winner = gendist(vector_likelihood_ratio',1,1);
        s_i=winner;
        likelihood_ratio = vector_likelihood_ratio (winner);
    end   
end



function [ mu ] = sample_mu(z,beta)
% Samples mu from its full conditional, given the current values
% for z, and the prior parameters, beta.

% z is a vector of counts;
% mu is a value;
    z_counts = [sum(z), length(z) - sum(z)];
    beta = beta + z_counts;
    mu = betarnd(beta(1),beta(2),1,1); % samples one value from the beta distribution 

end

function [ p ] = likelihood(sequence,s_i,M,K)
% Computes the likelihood of a subsequence of length K,
% beginning at s_i according to the model specified by M
sub_sequence=sequence(s_i:s_i+K-1);
pm=0;
for i=1:K
    pm=pm+log(M(sub_sequence(i),i));
end

p=exp(pm);

end

    

