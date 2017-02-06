function [ daughters ] = division_and_mutation(B_cells_trial, params)
%division_and_mutation( B_cells_trial, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, overlap, p_CDR_lethal, p_CDR_silent, kappa, sigma, mu, nb_Ag  )
%   "division_and_mutation" divides the cells and then mutates cells with
%   probability p_mut. When a cell is undergoes mutation, the function
%   "mutation" is called.


%% Division
% size of daughters = (2*nb_max_B_cells, nb_Ag +5)
daughters = cat(1,B_cells_trial, B_cells_trial);

%% Somatic Hyper-Mutation 
for n = 1:size(daughters,1)
    rand_mut = rand;   
    if rand_mut < params.experimental_params.p_mut %mutation occurs
        %daughters(n, nb_Ag + 5) = daughters(n, nb_Ag + 5) + 1;
        mutant = mutation(daughters(n,:), params );
%        daughters(n,:) = mutant(1,:);
        if ~isempty(mutant)
            daughters(n,:) = mutant(:);
        else            
            daughters(n,:) = NaN(1,size(daughters,2));
        end
        
    end
end

%% Removes the cells that got lethal mutations (empty mutant)
daughters = daughters(isnan(daughters(:,1)) < 1, :);

end
