function [B_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, final_cycle, success ] = runTrial(success, B_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, initial_cycle_number, params )
%runTrial( success, B_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, overlap, nb_max_B_cells, nb_cycle_max, initial_cycle_number, nb_Ag, energy_scale, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu  )  
%   The function "runTrial" runs "GC_cycle" for each cycle until the end of
%   the GC reaction. 
%   The GC_B_cells at the end of the GC reaction are stored
%   in B_cells_trial.
%   For each cycle, the cells that exit the GC are stored
%   in exit_cells_trial.
%   The number of exit cells and number of GC B cells are also recorded for
%   each cycle.
%   Final_cycle gives the end of the GC reaction. 
%   The GC reaction is terminated if the GC is extinguished or if all the
%   antigens have been consumed. 

cycle_number = initial_cycle_number;

%%GC REACTION:
%% the founders seed the GC and undergo AM. GC_cycle modifies b_cells and adds the new plasma cells. In this toy model, exit_cells designed the cells that exit the GC at the end of a cycle (memory + plasma cells).
% B_cells_trial = zeros(n_max_Bcells, nb_Ag + 5); a B_cell is [Ev1 Ev2 Ec
% Ag_index nb_mutations]

while 1
    cycle_number = cycle_number +1;
    %disp(['cycle_number ' num2str(cycle_number)]);
    [new_exit_cells, B_cells_trial] = GC_cycle( B_cells_trial, cycle_number, params);
    number_recycled_b_cells_trial(cycle_number) = size(B_cells_trial,1); %n GC cells
    number_exit_cells_trial(cycle_number) = size(new_exit_cells, 1);%n exit cells
    
    if number_recycled_b_cells_trial(cycle_number) > params.algorithm_constants.AM_constants.nb_max_B_cells
        success = success + 1;
        break
    end
    if cycle_number > params.algorithm_constants.AM_constants.nb_cycle_max || isempty(B_cells_trial) 
        %disp(['size when break ' num2str(size(B_cells_trial,1))]);
        break
    end

end

final_cycle = cycle_number - 1;
     
end
