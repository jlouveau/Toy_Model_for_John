function [ mutant ] = mutation( wildtype, activation_energy, threshold_energy, p_CDR, p_FR_lethal, overlap, p_CDR_lethal, p_CDR_silent,kappa, sigma, mu )
%   "mutation" takes a cell and mutates it. 
%   The type of mutation is determined and its effect on the energies for 
%   the two Ags is reported in the modified mutant matrix.

% size of mutant is (1, nb_Ag +2)
mutant = wildtype;

rand_CDR = rand;
rand_type = rand;

if rand_CDR <= p_CDR
    %% Mutation in the Complementarity-Determining Region
    if rand_type < p_CDR_lethal
        mutant = [];
    else if rand_type > p_CDR_lethal + p_CDR_silent
            mutant = energy_affecting_CDR_mutation(wildtype, kappa, sigma, mu);    
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

