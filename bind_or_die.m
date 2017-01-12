function [survival, B_cell, weighted_energy] = bind_or_die(B_cell, params)
%   "bind_or_die: tests the B-cell against ONE of the antigens present:
%   Ag_number is the number of antigens.
%   B-cell is a vector of length nb_Ag+1, giving the energies (between 0 and 1) of the
%   B-cell for the variable part of each antigen, followed by the energy
%   for the conserved part of the antigens.
%   B-cell: [ Ev1   Ev2    Ec   Ag_index].
%   survival is 0 if B-cell dies, otherwise is 1.
%   conc is antigen concentration, it is a parameter affecting the
%   auxiliary variable 'factor' multiplicatively.
%   overlap is the parameter between 0 and 1 which determines the size of
%   the region shared (or conserved) by the two Ags. 
%   activation_energy is the value between 0 and 1 which gives threshold
%   binding energy.
%   energy_scale scales up activation_energy, it affects the auxiliary
%   variable 'factor' exponentially.

Ag_index = randi(params.variable_params.nb_Ags); %randomly choose antigen
B_cell(params.variable_params.nb_Ags +2) = Ag_index;

Ev = B_cell(Ag_index); %find energy for variable part of randomly chosen antigen
Ec = B_cell(params.variable_params.nb_Ags + 1); %find energy for conserved part
weighted_energy = (1 - params.variable_params.overlap)*Ev + params.variable_params.overlap*Ec; %find overall affinity for randomly chosen antigen

Factor = params.algorithm_constants.AM_constants.conc * exp(params.algorithm_constants.energy_params.energy_scale*(weighted_energy - params.algorithm_constants.energy_params.activation_energy));
langmuir = Factor/(1+Factor); %ratio of bound Ags over total number of B cell receptors (c.f. Shenshen)
    
r = rand;
if r >= langmuir  %B-cell dies!
    survival = 0;
else %B-cell binds!
    survival = 1;
end
end