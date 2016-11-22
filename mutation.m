function [ mutant ] = mutation( wildtype, activation_energy, threshold_energy, p_CDR, p_FR_lethal, overlap)
%   "mutation" takes a cell and mutates it. 
%   The type of mutation is determined and its effect on the energies for 
%   the two Ags is reported in the modified mutant matrix.

p_CDR_lethal = 0.3;
p_CDR_silent = 0.5;
rand_CDR = rand;
rand_type = rand;
k = -0.7;
sigma = 1.2;
mu = -1.5;

% size of mutant is (1, nb_Ag +2)
mutant = wildtype;

if rand_CDR <= p_CDR
    %% Mutation in the Complementarity-Determining Region
    if rand_type < p_CDR_lethal
        mutant = [];
    else if rand_type > p_CDR_lethal + p_CDR_silent
            mutant = energy_affecting_CDR_mutation(wildtype, k, sigma, mu);    
        end
    end
else
    %% Mutation in the Framework Region
    if rand_type < p_FR_lethal
        mutant = [];
    else
        mutant = energy_affecting_FR_mutation(wildtype,  overlap, activation_energy, threshold_energy);    
    end
end

