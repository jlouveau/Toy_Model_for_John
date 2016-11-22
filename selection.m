function [ b_cells_trial ] = selection( b_cells_trial, conc, e_act, t_cell_selection) 
% each b_cell sees one antigen per cycle, chosen randomly. 
% the probability to internalize the Ag is of a langmuir form. Among these
% successful b_cells, only the best ones are selected by T cells.
% Amongst these B cells that are selected, only those that have affinity for
% the specific antigen above a_act are actually allowed to survive because
% the division rate depends on the affinity. 
% b_cells_trial = zeros(1, n_Ag, n_max_Bcells);

%disp(['number of b cells after selection by threshold ' num2str(length(b_cells))]);

conc_Ag = double(conc/size(b_cells_trial,2));
dominant_Ag = zeros(1,size(b_cells_trial,3));

%% Internalization of antigen: case where each b cell sees one antigen per cycle
%first step 
%is to pick the antigen that the b cell sees for each b cell it
%will be stored in dominant_Ag
%second step 
%is to stochastically choose which b cells survive using the
%langmuir probability with the affinity for the dominant_Ag

for n = 1:size(b_cells_trial,3)
    rand_Ag = rand;
    proba = rand;
    if rand_Ag < 0.5
        index_Ag = 1;
    else
        index_Ag = 2;
    end
    dominant_Ag(1,n) = index_Ag;
    factor = conc_Ag*exp(b_cells_trial(1,index_Ag,n) - e_act);
    langmuir = factor/(1+factor);
    
    if proba >= langmuir
        %b_cell didn't internalize the Ag, so it dies
        b_cells_trial(1,:,n) = [NaN NaN];
    end 
end

%remove the b_cells that failed to internalize the Ag
b_cells_trial = b_cells_trial(:, :, isnan(b_cells_trial(1,1,:)) < 1);  


%% T cell selection: 
%first step is
%to sort the b_cells by increasing value of the affinity for the antigen
%seen. 
%second step is 
%to select the top performing b cells.
%third step is 
%to remove from this selection the b cells with affinity
%below a_act since they won't be allowed to divide.

% Notice that all three steps are done on the internalized_b_cells matrix
% which contains the affinities for the dominant_Ag only

internalized_b_cells = zeros(1,1,size(b_cells_trial,3));

for n = 1:size(b_cells_trial,3)
    internalized_b_cells(1,:,n) = b_cells_trial(1, dominant_Ag(1,n),n);
end

%step 1
[sorted_internalized_b_cells, indexes] = sort(internalized_b_cells, 3, 'descend'); %%%%%%%%%%%%%%%%%%how it works???

%step 2
n_selected = floor(t_cell_selection*size(internalized_b_cells,3));
selected_internalized_b_cells = zeros(1,1,n_selected);
sorted_b_cells = zeros(1,size(b_cells_trial,2), n_selected);

for j = 1:n_selected
    selected_internalized_b_cells(1,1,j) = sorted_internalized_b_cells(1,1,j);
    sorted_b_cells(:,:,j) = b_cells_trial(1,:, indexes(1,1,j)); 
end

%step 3
for j = 1:n_selected   
    if selected_internalized_b_cells(1,1,j) < e_act
        sorted_b_cells(1,:,j) = [NaN NaN];
    end
end   
b_cells_trial = sorted_b_cells;
b_cells_trial = b_cells_trial(:, :, isnan(b_cells_trial(1,1,:)) < 1); 

end