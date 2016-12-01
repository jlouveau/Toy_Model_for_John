function [ new_exit_cells, B_cells_trial ] = GC_cycle( B_cells_trial, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, overlap, nb_Ag, energy_scale, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu)
%   GC_cycle will modify the B_cells_trial matrix (remove those that have lethal
%   mutation etc...) and it will create the output of new_exit_cells.
%   GC reaction is in three steps:
% 1) Dark zone: 2 divisions with potential mutation occur.
%   the function "division_and_mutation" is called twice and modifies the
%   size and content of B_cells_trial
% 2) Light zone: the function "select" discards the worst performers and
%   all cells which don't meet the activation_energy.
% 3) Recycling: removes cells from the GC and adds them to new_exit_cells.

%%DARK ZONE: mutation
%% for each B cell, determine if there is mutation, whether it's in the CDR or FR and the type of mutation. Then change the affinities (or delete) accordingly.

% 1st division + SHM
daughters1 = division_and_mutation(B_cells_trial, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, overlap, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu );

%2nd division + SHM
B_cells_trial = division_and_mutation(daughters1, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, overlap, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu );

%%LIGHT ZONE: selection
%% B cells that have affinity at least higher than a threshold and are in the top portion remain. 

B_cells_trial = select(B_cells_trial, nb_Ag, energy_scale, overlap, conc, activation_energy, t_cell_selection);

%%RECYCLE
%% randomly pick exit_cells from the selected b_cells.
[ B_cells_trial, new_exit_cells ] = recycling(B_cells_trial, p_recycle, t_cell_selection);

end

