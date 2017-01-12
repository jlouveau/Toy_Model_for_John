function [ mutant ] = energy_affecting_CDR_mutation(mutant, lognormal_params, gevd_params, distribution, nb_Ags)
%energy_affecting_CDR_mutation( mutant, kappa, sigma, mu, nb_Ag)
%   If the wildtype cell saw the first antigen in the previous round of
%   selection (or was selected to become founder), wildtype(nb_Ag+2) = 1.
%   Then Ev1 and Ec are increased by values randomly chosen from the skewed
%   gevd distribution. 
%   Ev2 is redrawn from uniform distribution between 0 and 1.
%disp('energy changing_mutations');

if distribution == 0
    %%gevd
    size_step = gevd_params.step;
    a = gevrnd(gevd_params.kappa, gevd_params.sigma, gevd_params.mu, 1, 1);
    b = gevrnd(gevd_params.kappa, gevd_params.sigma, gevd_params.mu, 1, 1);

    mutant(nb_Ags + 1) = mutant(nb_Ags + 1) + size_step*a;

    if mutant(nb_Ags +2) == 1
        mutant(1) = mutant(1) + size_step*b;
        mutant(2) = rand;   
    else
        mutant(1) = rand;
        mutant(2) = mutant(2) + size_step*b;    

    end
else
    %%lognormal
    a = lognormal_params.factor * (lognormal_params.offset - exp( normrnd(lognormal_params.mu, lognormal_params.sigma)));
    b = lognormal_params.factor * (lognormal_params.offset - exp( normrnd(lognormal_params.mu, lognormal_params.sigma)));
    
    mutant(nb_Ags + 1) = mutant(nb_Ags + 1) + a;
    
    if mutant(nb_Ags +2) == 1
        mutant(1) = mutant(1) + b;
        mutant(2) = rand;   
    else
        mutant(1) = rand;
        mutant(2) = mutant(2) + b;    

    end
end

end