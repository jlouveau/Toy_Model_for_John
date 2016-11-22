function [ mutant ] = energy_affecting_CDR_mutation( wildtype, k, sigma, mu)
%   If the wildtype cell saw the first antigen in the previous round of
%   selection (or was selected to become founder), wildtype(nb_Ag+2) = 1.
%   Then Ev1 and Ec are increased by values randomly chosen from the skewed
%   gevd distribution. 
%   Ev2 is redrawn from uniform distribution between 0 and 1.

mutant = wildtype;

if wildtype(length(mutant)) == 1
    mutant(1) = mutant(1) + gevrnd(k, sigma, mu, 1, 1);
    mutant(2) = rand;
    mutant(3) = mutant(3) + gevrnd(k, sigma, mu, 1, 1);
else
    mutant(1) = rand;
    mutant(2) = mutant(2) + gevrnd(k, sigma, mu, 1, 1);    
    mutant(3) = mutant(3) + gevrnd(k, sigma, mu, 1, 1);
end
end