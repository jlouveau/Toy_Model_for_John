function [B_cells] = replication(founder_B_cells, rep, nb_trial_max, nb_max_B_cells, nb_Ag)
% "replication" takes the founder cells and replicates them without
% mutation.

B_cells = zeros(nb_trial_max, nb_max_B_cells, nb_Ag + 2);
nb_founders = size(founder_B_cells,1);

for f = 1 : nb_founders
    f_start = (f-1)*2^rep + 1;   
    for b = f_start:f_start + 2^rep - 1
        for n = 1 : nb_trial_max
            B_cells(n,b,:) = founder_B_cells(f,:);
        end
    end
end
end