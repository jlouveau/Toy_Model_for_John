function [new_exit_cells, B_cells_trial] = GC_cycle( B_cells_trial, cycle_number, params)
%GC_cycle( B_cells_trial, cycle_number, nb_cycle_max, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, overlap, nb_Ag, energy_scale, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu, nb_max_B_cells)
%   GC_cycle will modify the B_cells_trial matrix (remove those that have lethal
%   mutation etc...) and it will create the output of new_exit_cells.
%   GC reaction is in three steps:
% 1) Dark zone: 2 divisions with potential mutation occur.
%   the function "division_and_mutation" is called twice and modifies the
%   size and content of B_cells_trial
% 2) Light zone: the function "select" discards the worst performers and
%   all cells which don't meet the activation_energy.
% 3) Recycling: removes cells from the GC and adds them to new_exit_cells.

% UPDATE B_cells_trial only if its final size is below < nb_max_B_cells

%%DARK ZONE: mutation
%% for each B cell, determine if there is mutation, whether it's in the CDR or FR and the type of mutation. Then change the affinities (or delete) accordingly.

% 1st division + SHM
daughters1 = division_and_mutation(B_cells_trial, params);
%disp('size after first division_mutation '); disp(size(B_cells_trial));

%2nd division + SHM
%daughters2 = daughters1;
daughters2 = division_and_mutation(daughters1, params );
%disp(['size after second division_mutation ' num2str(size(B_cells_trial,1))]);

%%LIGHT ZONE: selection
%% B cells that have affinity at least higher than a threshold and are in the top portion remain. 
daughters2 = select(daughters2, params);
%disp(['size after selection ' num2str(size(B_cells_trial,1))]);

%%RECYCLE
%% randomly pick exit_cells from the selected b_cells.
[ daughters2, new_exit_cells ] = recycling(daughters2, params.algorithm_constants.AM_constants.p_recycle);
%[ B_cells_trial ] = recycling(B_cells_trial, p_recycle);
%disp(['size after recycling ' num2str(size(B_cells_trial,1))]);

if size(B_cells_trial,1) <= params.algorithm_constants.AM_constants.nb_max_B_cells && cycle_number <= params.algorithm_constants.AM_constants.nb_cycle_max
    B_cells_trial = daughters2;
end

