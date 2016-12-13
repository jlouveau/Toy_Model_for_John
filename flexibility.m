function [ mutant ] = flexibility (mutant, low_Ew, low_index, high_Ew, high_index, activation_energy, threshold_energy)
% Problem with this function: if the delta has the effect of making the
% other Ew cross the limit. While loop?


delta_act = activation_energy - low_Ew;
delta_thr = threshold_energy - high_Ew;

if delta_act > 0 % low_Ew < activation_energy
    mutant(low_index) = mutant(low_index) + delta_act;
    mutant(1) = mutant(1) + delta_act;
end
if delta_thr < 0 % high_Ew > threshold_energy
    mutant(high_index) = mutant(high_index) + delta_thr;
    mutant(1) = mutant(1) + delta_thr;
end    
end

