function [ mutant ] = energy_affecting_FR_mutation(wildtype,  overlap, activation_energy, threshold_energy)
%   "energy_affecting_FR_mutation" increases the flexibility from the
%   wildtype to mutant.
%   The weighted energy for Ag_i: Ewi = (1-overlap)*Evi + overlap*Ec.
%   If Ewi < activation_energy, then modify Ec and Evi with the hypothesis
%   that they are as likely to be modified by a FR mutation.

mutant = wildtype;                                                                                                                                                                                                                                                                                              
    Ew1 = (1 - overlap)*wildtype(1) + overlap*wildtype(3);
    Ew2 = (1 - overlap)*wildtype(2) + overlap*wildtype(3);
    if Ew1 < Ew2
        mutant = flexibility(mutant, Ew1, 1, Ew2, 2, activation_energy, threshold_energy);
    else
        mutant = flexibility(mutant, Ew2, 2, Ew1, 1, activation_energy, threshold_energy);
    end
end