function [B_cells] = replication(founder_B_cells, params)
%replication(founder_B_cells, rep, nb_trial_max, nb_max_B_cells, nb_Ag)
% "replication" takes the founder cells and replicates them without
% mutation.

B_cells = zeros(params.algorithm_constants.AM_constants.nb_trial_max, params.algorithm_constants.AM_constants.nb_max_B_cells, params.variable_params.nb_Ags + 5);
nb_founders = size(founder_B_cells,1);

for f = 1 : nb_founders
    f_start = (f-1)*2^params.experimental_params.rep + 1;   
    for b = f_start:f_start + 2^params.experimental_params.rep - 1
        for n = 1 : params.algorithm_constants.AM_constants.nb_trial_max
            B_cells(n,b,:) = founder_B_cells(f,:);
        end
    end
end
end