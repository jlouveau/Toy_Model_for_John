function [B_cells_trial] = select(B_cells_trial, params)
%select(B_cells_trial, nb_Ag, energy_scale, overlap, conc, activation_energy, t_cell_selection)
%   The function "select" takes a collection of B-cells, binds them (or not) 
%   with antigens by calling the function "bind_or_die", applies T-cell 
%   selection, and discards any with low binding energy.
%   B_cells is a matrix whose rows are instances of B-cell (c.f. bind_or_die)
%   T-cell factor is the proportion of B-cell-Ag complexes which receive
%   T-cell help (why should it be a ratio, rather than an absolute?)

nb_Ag = params.variable_params.nb_Ags;

%B_cells_trial = zeros(n_max_Bcells, nb_Ag + 5);
%[Ev1 Ev2 Ec flex Ag_index nb_CDR_energyaffecting_mutations nb_FR_mutations]
nb_B_cells = size(B_cells_trial);
nb_B_cells = nb_B_cells(1);

selected_B_cells = zeros(nb_B_cells, nb_Ag + 6); % selected_B_cell: [ [B_Cell] Ew]

%% Ag binding
index = 0;
for i = 1 : nb_B_cells
    [survival, B_cell, weighted_energy] = bind_or_die(B_cells_trial(i,:), params);
    if survival == 1
        index = index + 1;
        selected_B_cells(index,:) = [B_cell weighted_energy];
    end
end
selected_B_cells = selected_B_cells(1:index, :);
%index is the number of B cells that survive binding

%B_cells_trial = selected_B_cells(:, 1:nb_Ag+3);
%disp(['n_bound ' num2str(size(B_cells_trial,1))]);

%% T cell help
selected_B_cells = sortrows(selected_B_cells, -(nb_Ag+6)); %sort by weighted_energy in descending order
cutoff = floor(index * params.algorithm_constants.AM_constants.t_cell_selection);
selected_B_cells = selected_B_cells(1:cutoff, :);
%B_cells_trial = selected_B_cells(:, 1:nb_Ag + 5); % top ordered B cells
%cutoff is the number of B cells that survive T cell selection

%disp(['n_helped ' num2str(size(B_cells_trial,1))]);

nb_selected = 1;
%nb_selected is the number of B cells that meet the energy activation
%threshold

%% Activation threshold
while (nb_selected <= cutoff && selected_B_cells(nb_selected, params.variable_params.nb_Ags + 5) >= params.algorithm_constants.energy_params.activation_energy)
    nb_selected = nb_selected + 1;
end
%B_cells_trial = selected_B_cells(1:nb_selected - 1, 1:nb_Ag + 5);
B_cells_trial = selected_B_cells(:, 1:nb_Ag + 5);
end
