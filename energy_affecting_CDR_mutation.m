function [ mutant ] = energy_affecting_CDR_mutation( mutant, kappa, sigma, mu, nb_Ag)
%   If the wildtype cell saw the first antigen in the previous round of
%   selection (or was selected to become founder), wildtype(nb_Ag+2) = 1.
%   Then Ev1 and Ec are increased by values randomly chosen from the skewed
%   gevd distribution. 
%   Ev2 is redrawn from uniform distribution between 0 and 1.
%disp('energy changing_mutations');
a = gevrnd(kappa, sigma, mu, 1, 1);
 mutant(nb_Ag + 1) = mutant(3) + a;
%  if a > 0
%      disp(['energy changing_mutations is beneficial' num2str(a)]);
%  end
if mutant(nb_Ag +2) == 1
    mutant(1) = mutant(1) + gevrnd(kappa, sigma, mu, 1, 1);
    mutant(2) = rand;   
else
    mutant(1) = rand;
    mutant(2) = mutant(2) + gevrnd(kappa, sigma, mu, 1, 1);    
   
end

end