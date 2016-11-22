function [ B_cells_trial, new_exit_cells ] = recycling(B_cells_trial, p_recycle, t_cell_selection)
% "recycling" chooses cells that exit the GC, adds them to new_exit_cells and
% remove those from the B_cells_trial matrix.
% the size of B_cells_trial is reduced by n_exit.

n_selected = floor(t_cell_selection*size(B_cells_trial,2));
n_exit = floor((1 - p_recycle)*n_selected);
new_exit_cells = zeros(1, n_exit, size(B_cells_trial,3));

for k = 1:n_exit
    ind = randi(size(B_cells_trial,2)); %randomly pick a cell in B_cells_trial
    for l = 1:size(B_cells_trial,3) %energies
        new_exit_cells(:,k,l) = B_cells_trial(:, ind, l); %add the chosen b cells to the exit b cells
    end
    B_cells_trial(:, ind, :) = []; % remove that b_cell from the list of GC b cells.
end
end