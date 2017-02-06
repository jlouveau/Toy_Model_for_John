function [ founder_B_cells ] = create_founders(indices, nb_Ag, energy_params, distrib_flex, flex_params, init )
%   "create_founders" determines the energies for each founder cell given
%   the indices vector.
%   indices gives, for each founder cell, the Ag_index for which the 
%   founder cell meets the activation energy.
% [Ev1 Ev2 Ec flex Ag_index nb_CDR_energyaffecting_mutations nb_FR_mutations].

%   If founder meets activation_energy for Ag_1 then:
%(1) weighted_energy_1 = (1-overlap)*(Ev1 + flex*E0) + overlap*(1-flex)*Ec = activation_energy;
%(2) weighted_energy_2 = (1-overlap)*(Ev2 + flex*E0) + overlap*(1-flex)*Ec <= activation_energy;
%   Under the hypothesis that shared and variable regions are equally
%   likely to bind to the Ag, Ev1 = Ec.
%   So Ev1 = Ec = activation_energy 
%   and Ev2 should be drawn from [0 activation_energy]

nb_founders = length(indices);
nb_founders = nb_founders(1);

if nb_Ag ==2
    
    founder_B_cells = zeros(nb_founders, nb_Ag + 5);
    
    for f = 1 : nb_founders
        rand_f = rand;
        %in case of 2 Ags if <0.5 then meets activation for Ag_1
        if indices(f) < 0.5
            %Ag_index = 1;
            founder_B_cells(f, nb_Ag + 3) = 1;
            founder_B_cells(f,1) = energy_params.activation_energy;
            founder_B_cells(f,2) = rand_f*energy_params.activation_energy;            
        else
            %Ag_index = 2;
            founder_B_cells(f, nb_Ag + 3) = 2;
            founder_B_cells(f,1) = rand_f*energy_params.activation_energy;
            founder_B_cells(f,2) = energy_params.activation_energy;
        end
        
        if init.level_to_threshold_energy > 0
            founder_B_cells(f,nb_Ag + 1) = energy_params.threshold_energy + init.epsilon;
        else
            founder_B_cells(f,nb_Ag + 1) = energy_params.activation_energy + rand*(energy_params.threshold_energy-energy_params.activation_energy) - init.epsilon;
        end

        if distrib_flex == 0 %(normal)
             founder_B_cells(f, nb_Ag + 2) = init.initial_flex_center + normrnd(flex_params.mu, flex_params.sigma);
        else
            founder_B_cells(f, nb_Ag + 2) = flex_params.min + rand*(flex_params.max - flex_params.min);
        end
        
    end
else
    disp('not enough antigens');
    exit
end
end

