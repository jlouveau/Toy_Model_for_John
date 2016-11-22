function [ founder_B_cells ] = create_founders(indices, nb_Ag, activation_energy )
%   "create_founders" determines the energies for each founder cell given
%   the indices vector.
%   indices gives, for each founder cell, the Ag_index for which the 
%   founder cell meets the activation energy.
%   [Ev1 Ev2 Ec Ag_index].
%   If founder meets activation_energy for Ag_1 then:
%(1) weighted_energy_1 = (1-overlap)*Ev1 + overlap*Ec = activation_energy;
%(2) weighted_energy_2 = (1-overlap)*Ev2 + overlap*Ec <= activation_energy;
%   Under the hypothesis that shared and variable regions are equally
%   likely to bind to the Ag, Ev1 = Ec.
%   So Ev1 = Ec = activation_energy 
%   and Ev2 should be drawn from [0 activation_energy]

nb_founders = length(indices);
nb_founders = nb_founders(1);

founder_B_cells = zeros(nb_founders, nb_Ag + 2);

for f = 1 : nb_founders
    rand_f = rand;
    %in case of 2 Ags if <0.5 then meets activation for Ag_1
        if indices(f) < 0.5
            %Ag_index = 1;
            founder_B_cells(f, nb_Ag + 2) = 1; 
            founder_B_cells(f,1) = activation_energy; 
            founder_B_cells(f,2) = rand_f*activation_energy; 
            founder_B_cells(f,nb_Ag + 1) = activation_energy;
        else
            %Ag_index = 2;
            founder_B_cells(f, nb_Ag + 2) = 2; 
            founder_B_cells(f,1) = rand_f*activation_energy; 
            founder_B_cells(f,2) = activation_energy; 
            founder_B_cells(f,nb_Ag + 1) = activation_energy;
        end
end   

end

