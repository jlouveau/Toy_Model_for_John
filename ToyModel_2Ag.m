%% 1 or 2 antigens only!!!!!
%%
clear all;
close all;
clc;

%% 3 CDRs of about 10 residues each on the heavy chain. 
% The FR represents 60-85% of the chain.
% Shenshen used 46 CDR residues, 18 conserved (35%).

experimental_params = struct('p_mut', 0.2, 'rep', 9, 'p_CDR_lethal', 0.3, 'p_CDR_silent', 0.5);
variable_params = struct('nb_Ags', 2, 'overlap', 1, 'nb_founders', 3, 'p_CDR', 1, 'p_FR_lethal', 0.8);
initial_nb_Bcells = variable_params.nb_founders * 2^experimental_params.rep;

AM_constants = struct('nb_cycle_max',250, 'nb_trial_max', 10, 'nb_max_B_cells', initial_nb_Bcells, 'p_recycle', 0.7, 't_cell_selection', 0.7, 'conc', 1.2);
lognormal_params = struct('offset', 3.0, 'sigma', 0.5, 'mu', 1.9, 'factor', 1, 'activation_energy', 0.1, 'energy_scale', 0.05);
gevd_params = struct('kappa', -0.7, 'sigma', 1.2 , 'mu', -1.5, 'step', 5, 'activation_energy', 0.1, 'energy_scale', 0.05);

%% choose lognormal or gevd
distribution = 0; %0 for gevd, 1 for lognormal

if distribution == 0
    energy_params = struct('activation_energy', gevd_params.activation_energy, 'threshold_energy', 0.7, 'energy_scale', gevd_params.energy_scale);
else
    energy_params = struct('activation_energy', lognormal_params.activation_energy, 'threshold_energy', 0.7, 'energy_scale', lognormal_params.energy_scale);   
end

%%
algorithm_constants = struct('AM_constants', AM_constants, 'energy_params', energy_params, 'gevd_params', gevd_params, 'lognormal_params', lognormal_params, 'distribution', distribution);
params = struct('experimental_params', experimental_params, 'variable_params', variable_params, 'algorithm_constants', algorithm_constants);

%% Creating matrices for founder B cells, B cells, exit cells and record the 
% number of GC cells that are recycled and those who exit the GC for each 
% trial and each cycle.
% founder_B_cells is a matrix whose rows are instances of B_cell [Ev1 Ev2 Ec
% Ag_index nb_mutations] where Ag_index is the index of the Ag for which the founder
% cell meets the activation_energy.
%
% B_cells is the matrix of GC B cells = [trial B_cell_index energies nb_mutations] 

indices = rand(variable_params.nb_founders,1); %gives the Ag for which each founder B cell meets the activation_energy 

number_recycled_b_cells = zeros(AM_constants.nb_trial_max, AM_constants.nb_cycle_max);
number_exit_cells = zeros(AM_constants.nb_trial_max, AM_constants.nb_cycle_max);

tic;

%%INITIALIZATION + proliferation: 
%% each founder needs to meet a_act for one Ag
%cycle 1: creates founders
cycle_number = 1;
% founder_B_cells = zeros(nb_founders, nb_Ag + 3);
founder_B_cells = create_founders(indices, variable_params.nb_Ags, energy_params.activation_energy);
 
number_recycled_b_cells(:,cycle_number) = variable_params.nb_founders;
number_exit_cells(:,cycle_number) = 0;

%cycle 2: replication without mutation
cycle_number = cycle_number + 1;

% B_cells = zeros(nb_trial_max, nb_max_B_cells, nb_Ag + 3);
B_cells = replication(founder_B_cells, params);

number_recycled_b_cells(:,cycle_number) = size(B_cells,2);   
number_exit_cells(:,cycle_number) = 0;      

%% Stochastic process: reproduce the GC reaction many times.
initial_cycle_number = 2;
[B_cells, number_recycled_b_cells, number_exit_cells, final_cycles, success ] = runAffinityMaturation(B_cells, number_recycled_b_cells, number_exit_cells, initial_cycle_number, params);
toc; 
%% Analyze trials
[muts] = analysis(success, B_cells, number_recycled_b_cells, final_cycles, params);



