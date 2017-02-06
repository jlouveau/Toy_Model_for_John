function [ mutant ] = energy_affecting_FR_mutation(wildtype, flexibility_params, distribution_FR, nb_Ags)
%energy_affecting_FR_mutation(wildtype,  overlap, activation_energy, threshold_energy)
%   "energy_affecting_FR_mutation" modifies the flexibility from the
%   wildtype to mutant.
%   The weighted energy for Ag_i: Ewi = (1-overlap)*[Evi + flex * Evj]+ overlap*(1-flex)*Ec.
%   delta_Ewi = (1-overlap)*delta_flex*Evi - overlap*delta_flex*Ec

mutant = wildtype; 
 

if distribution_FR == 0
    delta_flex = normrnd(flexibility_params.mu, flexibility_params.sigma);
else
    delta_flex = flexibility_params.min + rand*(flexibility_params.max - flexibility_params.min);
end

mutant(nb_Ags + 2) = mutant(nb_Ags + 2) + delta_flex;

if mutant(nb_Ags + 2) > flexibility_params.max
    mutant(nb_Ags + 2) = flexibility_params.max;
elseif mutant(nb_Ags + 2) < flexibility_params.min
    mutant(nb_Ags + 2) = flexibility_params.min;   
end

% case where Ec is dominant: increasing flexibility will lower both
% weighted energies, meaning that the cell will likely not survive
% selection. If Ev dominates, increasing flexibility is beneficial.

% conserved_part = params.variable_params.overlap*wildtype(3);
% variable1 = (1 - params.variable_params.overlap)*wildtype(1);
% variable2 = (1 - params.variable_params.overlap)*wildtype(2);
% 
% if conserved_part > variable1 && conserved_part > variable2
%     
% mutant = wildtype;                                                                                                                                                                                                                                                                                              
%     Ew1 = (1 - overlap)*wildtype(1) + overlap*wildtype(3);
%     Ew2 = (1 - overlap)*wildtype(2) + overlap*wildtype(3);
%     if Ew1 < Ew2
%         mutant = flexibility(mutant, Ew1, 1, Ew2, 2, energy_params.activation_energy, energy_params.threshold_energy);
%     else
%         mutant = flexibility(mutant, Ew2, 2, Ew1, 1, energy_params.activation_energy, energy_params.threshold_energy);
%     end
% end