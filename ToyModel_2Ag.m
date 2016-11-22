clear all;
close all;
clc;

%%parameters
nb_Ag = 2;
nb_founders = 3;
rep = 9;
nb_max_B_cells = nb_founders*2^rep;
nb_cycle_max = 250;
nb_trial_max = 50;
activation_energy = 0.1;
threshold_energy = 0.7;
energy_scale = 10;
conc = 1;
overlap = 0.3;

p_mut = 0.2; %per division
p_CDR = 0.3;
p_FR_lethal = 0.8;
p_recycle = 0.7;
t_cell_selection = 0.35;

%% Creating matrices for founder B cells, B cells, exit cells and record the 
% number of GC cells that are recycled and those who exit the GC for each 
% trial and each cycle.
% founder_B_cells is a matrix whose rows are instances of B_cell [Ev1 Ev2 Ec
% Ag_index] where Ag_index is the index of the Ag for which the founder
% cell meets the activation_energy.
%
% B_cells is the matrix of GC B cells = [trial B_cell_index energy] 

indices = rand(nb_founders,1); %gives the Ag for which each founder B cell meets the activation_energy 

exit_cells = zeros(nb_trial_max, nb_cycle_max, floor(nb_max_B_cells/4),  nb_Ag + 2);
number_recycled_b_cells = zeros(nb_trial_max, nb_cycle_max);
number_exit_cells = zeros(nb_trial_max, nb_cycle_max);

tic;

%%INITIALIZATION + proliferation: 
%% each founder needs to meet a_act for one Ag
%cycle 1: creates founders
cycle_number = 1;
% founder_B_cells = zeros(nb_founders, nb_Ag + 2);
founder_B_cells = create_founders(indices, nb_Ag, activation_energy );
 
number_recycled_b_cells(:,cycle_number) = nb_founders;
number_exit_cells(:,cycle_number) = 0;

%cycle 2: replication without mutation
cycle_number = cycle_number + 1;

% B_cells = zeros(nb_trial_max, nb_max_B_cells, nb_Ag + 2);
B_cells = replication(founder_B_cells, rep, nb_trial_max, nb_max_B_cells, nb_Ag);

number_recycled_b_cells(:,cycle_number) = size(B_cells,2);
number_exit_cells(:,cycle_number) = 0;

%% Stochastic process: reproduce the GC reaction many times.
initial_cycle_number = 2 ;
[B_cells, exit_cells, number_recycled_b_cells, number_exit_cells ] = runAffinityMaturation(B_cells, exit_cells, number_recycled_b_cells, number_exit_cells, nb_trial_max, conc, activation_energy, threshold_energy, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, cycle_number, overlap, nb_max_B_cells, nb_cycle_max, nb_Ag, energy_scale);

toc; 
%% Analyze trials
[ pop_time, total_exit_cells, neutralized, breadth, energy, affinity ] = analysis( number_recycled_b_cells, number_exit_cells, exit_cells, nb_trial_max, activation_energy, nb_cycle_max, p_mut, p_recycle, t_cell_selection, conc, p_CDR);


