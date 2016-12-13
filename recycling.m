function [ B_cells_trial, new_exit_cells ] = recycling(B_cells_trial, p_recycle)
% "recycling" chooses cells that exit the GC, adds them to new_exit_cells and
% remove those from the B_cells_trial matrix.
% the size of B_cells_trial is reduced by n_exit.

% exit_cells_trial = zeros(n_cycle_max, floor(n_max_Bcells/4), n_Ag + 2);
% B_cells_trial = zeros(n_max_Bcells, nb_Ag + 2);

n_selected = floor(size(B_cells_trial,1));
n_exit = floor((1 - p_recycle)*n_selected);
new_exit_cells = zeros(n_exit, size(B_cells_trial,2));
%disp(['n_selected ' num2str(n_selected)]);

for k = 1:n_exit
    ind = randi(size(B_cells_trial,1)); %randomly pick a cell in B_cells_trial
    for l = 1:size(B_cells_trial,2) %energies
        new_exit_cells(k,l) = B_cells_trial(ind, l); %add the chosen b cells to the exit b cells
    end
    B_cells_trial(ind, :) = []; % remove that b_cell from the list of GC b cells.
end
%disp(['n_recycled ' num2str(size(B_cells_trial,1))]);
end