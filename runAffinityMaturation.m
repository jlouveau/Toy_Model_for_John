function [B_cells, number_recycled_b_cells, number_exit_cells, final_cycles, success ] = runAffinityMaturation(B_cells, number_recycled_b_cells, number_exit_cells, initial_cycle_number, params)
%runAffinityMaturation(B_cells, number_recycled_b_cells, number_exit_cells, nb_trial_max, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, initial_cycle_number, overlap, nb_max_B_cells, nb_cycle_max, nb_Ag, energy_scale, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu)
%   "runAffinityMaturation" runs "runTrial" for each trial until the max
%   number of trials is reached.
%   The B_cells, number_recycled_b_cells and number_exit_cells 
%   are updated.

% B_cells_trial: size = (n_max_Bcells, nb_Ag +5). Starts with the
% replicated founders.
% number_recycled_b_cells_trial starts as zeros(1, n_cycle_max);
% number_exit_cells_trial starts as zeros(1, n_cycle_max);
% initial_cycle_number = 2;

trial_number = 1;
success = 0;
trial_max = params.algorithm_constants.AM_constants.nb_trial_max;
final_cycles = zeros(trial_max,1);

while trial_number <= trial_max
    
    disp(['TRIAL NUMBER ' num2str(trial_number)]);
    
    B_cells_trial = squeeze(B_cells(trial_number,:,:));
    number_recycled_b_cells_trial = number_recycled_b_cells(trial_number,:);
    number_exit_cells_trial = number_exit_cells(trial_number,:);
    
    [B_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, final_cycle, success ] = runTrial(success, B_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, initial_cycle_number, params );

    final_cycles(trial_number) = final_cycle;
    
    gusnumber = min(size(B_cells_trial,1), params.algorithm_constants.AM_constants.nb_max_B_cells);
    
    for j = 1:gusnumber %nb_B_cells I CHANGED THIS!!!
        for k = 1:size(B_cells_trial,2) %energies & mutations
            B_cells(trial_number, j,k) = B_cells_trial(j,k);
        end
    end
    
    for i = 1:final_cycle 
        number_recycled_b_cells(trial_number,i) = number_recycled_b_cells_trial(i);
        number_exit_cells(trial_number,i) = number_exit_cells_trial(i);
    end

    trial_number = trial_number +1;
    
end
end