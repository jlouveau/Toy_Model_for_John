function [B_cells, number_recycled_b_cells, number_exit_cells, final_cycles, success ] = runAffinityMaturation(B_cells, number_recycled_b_cells, number_exit_cells, nb_trial_max, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, initial_cycle_number, overlap, nb_max_B_cells, nb_cycle_max, nb_Ag, energy_scale, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu)
%   "runAffinityMaturation" runs "runTrial" for each trial until the max
%   number of trials is reached.
%   The B_cells, exit_cells, number_recycled_b_cells and number_exit_cells 
%   are updated.

% B_cells_trial: size = (n_max_Bcells, nb_Ag +3). Starts with the
% replicated founders.
% exit_cells_trial starts as zeros(n_cycle_max, floor(n_max_Bcells/4), nb_Ag +3);
% number_recycled_b_cells_trial starts as zeros(1, n_cycle_max);
% number_exit_cells_trial starts as zeros(1, n_cycle_max);
% initial_cycle_number = 2;

trial_number = 1;
success = 0;
final_cycles = zeros(nb_trial_max,1);

while trial_number <= nb_trial_max
    
    disp(['TRIAL NUMBER ' num2str(trial_number)]);
    
    B_cells_trial = squeeze(B_cells(trial_number,:,:));
    number_recycled_b_cells_trial = number_recycled_b_cells(trial_number,:);
    %exit_cells_trial = squeeze(exit_cells(trial_number, :, :,:));
    number_exit_cells_trial = number_exit_cells(trial_number,:);
    
    %[B_cells_trial, number_recycled_b_cells_trial, exit_cells_trial, number_exit_cells_trial, final_cycle ] = runTrial( B_cells_trial, exit_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, overlap, nb_max_B_cells, nb_cycle_max, initial_cycle_number, nb_Ag, energy_scale, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu );
    [B_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, final_cycle, success ] = runTrial(success, B_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, overlap, nb_max_B_cells, nb_cycle_max, initial_cycle_number, nb_Ag, energy_scale, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu );

    final_cycles(trial_number) = final_cycle;
    
    gusnumber = min(size(B_cells_trial,1), nb_max_B_cells);
    
    for j = 1:gusnumber %nb_B_cells I CHANGED THIS!!!
        for k = 1:size(B_cells_trial,2) %energies & mutations
            B_cells(trial_number, j,k) = B_cells_trial(j,k);
        end
    end
    
    for i = 1:final_cycle 
        number_recycled_b_cells(trial_number,i) = number_recycled_b_cells_trial(i);
        number_exit_cells(trial_number,i) = number_exit_cells_trial(i);
        
%         for j = 1:size(exit_cells_trial,2) %nb_B_cells
%             for k = 1:size(B_cells_trial,3) %energies
%                 exit_cells(trial_number,i,j,k) = exit_cells_trial(i,j,k);
%             end
%         end
    end

    
%     if final_cycle <= nb_cycle_max && number_recycled_b_cells_trial(final_cycle) >= nb_max_B_cells
%         for i = final_cycle : nb_cycle_max
%             number_recycled_b_cells(trial_number,i) = nb_max_B_cells;
%         end
%     end
    
    trial_number = trial_number +1;
    
end
end