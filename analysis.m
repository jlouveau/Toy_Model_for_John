function [muts_CDR muts_FR] = analysis(success, B_cells, number_recycled_b_cells, final_cycles, params)
%analysis( success, B_cells, overlap, number_recycled_b_cells, nb_trial_max, nb_max_B_cells, p_mut, p_recycle, t_cell_selection, conc, p_CDR, final_cycles, nb_Ag, nb_cycle_max)
%% GC population over time
figure();
% x = linspace(2,200,1000);
% curve1 = nb_max_B_cells * (4 * (1-p_mut)^2 * p_recycle * t_cell_selection * conc /(1+conc) ).^(x-2);
% plot(x, curve1, ':');

for i = 1:params.algorithm_constants.AM_constants.nb_trial_max
    for j = 1:final_cycles(i)
        hold on; plot(number_recycled_b_cells(i,1:j)); 
    end
end
title({['Population of GC b cells for 2 Ags with ' num2str(params.variable_params.overlap*100) '% shared residues']; [' conc = ' num2str(params.algorithm_constants.AM_constants.conc) ' proba CDR = ' num2str(params.variable_params.p_CDR)]}, 'Fontweight', 'bold');
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',10);

% %% Final GC population
% % of interest for debugging with 3 cycles max
% final_pop = zeros(nb_max_B_cells,1);
% for i = 1:nb_trial_max
%     final_pop(number_recycled_b_cells(i, final_cycles(i)),1) = final_pop(number_recycled_b_cells(i, final_cycles(i)),1)+1;
% end
% figure(); bar(final_pop, 0.8);
% title('Number of GC B cells at the end ', 'Fontweight', 'bold');
% set(gca,'FontSize',10);

%% GC success rate
success = success / params.algorithm_constants.AM_constants.nb_trial_max;
disp(['ratio of GCs that are successful ' num2str(success)]);

%% final cycle
figure(); histogram(final_cycles, 'Normalization', 'probability', 'BinWidth', 5);
title('Number of cycles until GC reaction stops ', 'Fontweight', 'bold');
set(gca,'FontSize',10);

%% mutations
gusmuts_CDR = zeros(30,1);
gusmuts_FR = zeros(30,1);
x = 0:29;
total_final_cells_successful_GCs = 0;
for i = 1:params.algorithm_constants.AM_constants.nb_trial_max
    if number_recycled_b_cells(i, final_cycles(i)) > 0
        for j = 1:number_recycled_b_cells(i, final_cycles(i))
            k = B_cells(i,j,params.variable_params.nb_Ags+4);
            l = B_cells(i,j,params.variable_params.nb_Ags+5);
            gusmuts_CDR(k+1,1) = gusmuts_CDR(k+1,1)+1;
            gusmuts_FR(l+1,1)  = gusmuts_FR(l+1,1)+1;
            total_final_cells_successful_GCs = total_final_cells_successful_GCs +1;
        end
    end
end
muts_CDR = gusmuts_CDR/total_final_cells_successful_GCs;
muts_FR = gusmuts_FR/total_final_cells_successful_GCs;

figure(); bar(muts_CDR, 0.5);
title('Number of energy affecting CDR mutations for the final B cells ', 'Fontweight', 'bold');
set(gca,'FontSize',10);

figure(); bar(x,muts_FR);
title('Number of FR mutations for the final B cells ', 'Fontweight', 'bold');
set(gca,'FontSize',10);

%% energies at the end of the GC reaction for successful GCs
e_conserved = zeros(total_final_cells_successful_GCs,1);
%e_variable = zeros(nb_trial_max,nb_Ag);
index = 1;
for i = 1:params.algorithm_constants.AM_constants.nb_trial_max
    if number_recycled_b_cells(i, final_cycles(i)) > 0
        for j = 1:number_recycled_b_cells(i, final_cycles(i))
            e_conserved(index) = B_cells(i,j,params.variable_params.nb_Ags+1);
            index = index +1;
        end
    end
end
figure(); histogram(e_conserved, 'Normalization', 'probability', 'BinWidth', 0.01);
title('Conserved energies at the end of the GC reaction for successful GCs ', 'Fontweight', 'bold');
set(gca,'FontSize',10);

% %% affinities at the end of the GC reaction for successful GCs
aff_conserved = exp(e_conserved);

figure(); histogram(aff_conserved, 'Normalization', 'probability', 'BinWidth', 0.01);
title('Conserved affinities at the end of the GC reaction for successful GCs ', 'Fontweight', 'bold');
set(gca,'FontSize',10);

% %% flexibilities at the end of the GC reaction for successful GCs
% flexs = zeros(total_final_cells_successful_GCs,1);
% index = 1;
% for i = 1:params.algorithm_constants.AM_constants.nb_trial_max
%     if number_recycled_b_cells(i, final_cycles(i)) > 0
%         for j = 1:number_recycled_b_cells(i, final_cycles(i))
%             flexs(index) = B_cells(i,j,params.variable_params.nb_Ags+1);
%             index = index +1;
%         end
%     end
% end
% figure(); histogram(flexs, 'Normalization', 'probability', 'BinWidth', 0.01);
% title('Flexibilities at the end of the GC reaction for successful GCs ', 'Fontweight', 'bold');
% set(gca,'FontSize',10);
end
