function [ mutant ] = mutation(wildtype, params)
%mutation( wildtype, activation_energy, threshold_energy, p_CDR, p_FR_lethal, overlap, p_CDR_lethal, p_CDR_silent,kappa, sigma, mu, nb_Ag)
%   "mutation" takes a cell and mutates it. 
%   The type of mutation is determined and its effect on the energies for 
%   the two Ags is reported in the modified mutant matrix.

% size of mutant is (1, nb_Ag +3)
% mutant = wildtype;
% 
% mutant = wildtype;

rand_CDR = rand;
rand_type = rand;

mutant = wildtype;

if rand_CDR <= params.variable_params.p_CDR
    %% Mutation in the Complementarity-Determining Region
    if rand_type <= params.experimental_params.p_CDR_lethal
        %disp('lethal');
        mutant = [];
    elseif rand_type > params.experimental_params.p_CDR_lethal + params.experimental_params.p_CDR_silent
        mutant = energy_affecting_CDR_mutation(mutant, params.algorithm_constants.lognormal_params, params.algorithm_constants.gevd_params, params.algorithm_constants.distribution, params.variable_params.nb_Ags);    
        mutant(params.variable_params.nb_Ags + 3) = mutant(params.variable_params.nb_Ags + 3) +1; %we're not counting silent mutations
    end
else
    %% Mutation in the Framework Region
    if rand_type < p_FR_lethal
        mutant = [];
    else
        mutant = energy_affecting_FR_mutation(mutant, params.variable_params.overlap, params.algorithm_constants.energy_params);    
    end
end

%%adds mutation
 
end

